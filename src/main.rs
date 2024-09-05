use gloo::timers::callback::Interval;
use rand::prelude::*;
use wasm_bindgen::JsValue;
use yew::prelude::*;

fn main() {
    let document = gloo::utils::document();
    let target_element = document.get_element_by_id("onthis").unwrap();
    yew::Renderer::<App>::with_root(target_element).render();
}

const W: usize = 200;
const H: usize = 150;
const P_FOOD: f64 = 0.3;

const MAXHP: usize = 15;
const RECOVER: usize = 5;

const CHAN_NUM: usize = 50;

const MIX_NUM: usize = 5;
const MUT_NUM: usize = 5;
const MUT_CPY: usize = 2;
const MUT_MUCH: usize = 100;

const TRAIN_STEPMAX: usize = 1000;
const TRAIN_NUM: usize = 50;

#[derive(Debug, Clone, PartialEq, Properties)]

struct Field {
    field: [[bool; W]; H],
}

pub fn log<T: AsRef<str>>(str: T) {
    web_sys::console::log_1(&JsValue::from_str(str.as_ref()))
}

#[derive(Debug, Clone, PartialEq)]
struct Chan {
    pos: (usize, usize),
    hp: usize,
    act: Gene,
    mem: u8,
    ori: Ori,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Ori {
    U,
    D,
    R,
    L,
}

const MOVE_NUM: usize = 2_usize.pow((3 * 3) as u32 - 1) * (2_usize.pow(8));

#[derive(Debug, Clone, PartialEq)]
struct Gene {
    act: [Action; MOVE_NUM],
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Action {
    Forword,
    Back,
    TurnRight,
    TurnLeft,
}

impl Gene {
    fn get_action<F>(&self, surr: F, mem: u8) -> &Action
    where
        F: Fn((isize, isize)) -> bool,
    {
        let mut index: u16 = mem as u16;
        for i in -1..=1 {
            for j in -1..=1 {
                if surr((i, j)) {
                    let s = 3 * i + j;
                    let s = if s < 0 { -s + 3 } else { s - 1 };
                    index |= 1 << (8 + s);
                }
            }
        }
        &self.act[index as usize]
    }

    fn random(rng: &mut ThreadRng) -> Self {
        let mut act = [Action::Forword; MOVE_NUM];
        for i in 0..MOVE_NUM {
            act[i] = match rng.gen_range(0..4) {
                0 => Action::Forword,
                1 => Action::Back,
                2 => Action::TurnLeft,
                3 => Action::TurnRight,
                _ => unreachable!(),
            };
        }
        Self { act }
    }

    fn mutate(&mut self, num: usize) {
        let mut rng = thread_rng();
        for _ in 0..num {
            let i = rng.gen_range(0..MOVE_NUM);
            let o = match rng.gen_range(0..4) {
                0 => Action::Forword,
                1 => Action::Back,
                2 => Action::TurnLeft,
                3 => Action::TurnRight,
                _ => unreachable!(),
            };
            self.act[i] = o;
        }
    }

    fn mix(&self, other: &Self) -> Self {
        let mut act = [Action::Forword; MOVE_NUM];
        let mut rng = thread_rng();
        for i in 0..MOVE_NUM {
            act[i] = if rng.gen_bool(0.5) {
                self.act[i]
            } else {
                other.act[i]
            }
        }
        Self { act }
    }
}

fn goto(pos: (usize, usize), diff: (isize, isize)) -> Option<(usize, usize)> {
    let i = {
        let i = pos.0 as isize + diff.0;
        if i < 0 || H as isize <= i {
            return None;
        }
        i as usize
    };
    let j = {
        let j = pos.1 as isize + diff.1;
        if j < 0 || W as isize <= j {
            return None;
        }
        j as usize
    };
    Some((i, j))
}

fn rotate_from_chan(diff: (isize, isize), ori: Ori) -> (isize, isize) {
    match ori {
        Ori::U => diff,
        Ori::D => (-diff.0, -diff.1),
        Ori::R => (-diff.1, diff.0),
        Ori::L => (diff.1, -diff.0),
    }
}

#[derive(Debug, Clone, PartialEq, Properties)]
struct OneGame {
    remain: usize,
    field: Field,
    chans: Vec<Result<Chan, (usize, Gene)>>,
}

impl OneGame {
    fn new() -> Self {
        let mut rng = thread_rng();
        let mut field = [[false; W]; H];
        for i in 0..H {
            for j in 0..W {
                if rng.gen_bool(P_FOOD) {
                    field[i][j] = true;
                }
            }
        }
        let chans = (0..CHAN_NUM)
            .map(|_| {
                let pos = (rng.gen_range(0..H), rng.gen_range(0..W));
                Ok(Chan {
                    pos,
                    hp: MAXHP,
                    act: Gene::random(&mut rng),
                    mem: 0,
                    ori: Ori::U,
                })
            })
            .collect::<Vec<_>>();

        Self {
            remain: CHAN_NUM,
            field: Field { field },
            chans,
        }
    }
    fn from_gene(genes: Vec<Gene>) -> Self {
        assert_eq!(genes.len(), CHAN_NUM);
        let mut rng = thread_rng();
        let mut field = [[false; W]; H];
        for i in 0..H {
            for j in 0..W {
                if rng.gen_bool(P_FOOD) {
                    field[i][j] = true;
                }
            }
        }
        let mut chans = vec![];
        for gene in genes {
            let pos = (rng.gen_range(0..H), rng.gen_range(0..W));
            let chan = Chan {
                pos,
                hp: MAXHP,
                act: gene,
                mem: 0,
                ori: Ori::U,
            };
            chans.push(Ok(chan));
        }
        Self {
            remain: CHAN_NUM,
            field: Field { field },
            chans,
        }
    }

    fn is_end(&self) -> bool {
        self.remain == 0
    }

    fn end(&mut self) -> Self {
        let mut genes = nice_gene(&self);
        // top n の mix
        let mix_genes = {
            let mut v = vec![];
            for i in 0..MIX_NUM {
                for j in 0..MIX_NUM {
                    if i != j {
                        v.push(genes[i].mix(&genes[j]));
                    }
                }
            }
            v
        };
        // top n の mutate
        let mut_genes = {
            let mut v = vec![];
            for i in 0..MUT_NUM {
                for _ in 0..MUT_CPY {
                    let mut g = genes[i].clone();
                    g.mutate(MUT_MUCH);
                    v.push(g);
                }
            }
            v
        };

        for _ in 0..mix_genes.len() {
            genes.pop();
        }

        for _ in 0..mut_genes.len() {
            genes.pop();
        }

        genes.extend(mix_genes);
        genes.extend(mut_genes);
        OneGame::from_gene(genes)
    }
}

fn nice_gene(
    OneGame {
        remain,
        field,
        chans,
    }: &OneGame,
) -> Vec<Gene> {
    let mut gene_order: Vec<_> = chans
        .iter()
        .map(|result| match result {
            Ok(chan) => ((0, chan.hp), chan.act.clone()),
            Err((order, gene)) => ((1, *order), gene.clone()),
        })
        .collect();
    gene_order.sort_by_key(|(o, g)| *o);
    gene_order.into_iter().map(|(o, g)| g).collect()
}

fn train() -> OneGame {
    log("start");
    let mut game = OneGame::new();
    'a: for i in 0..TRAIN_NUM {
        log(format!("train {i}"));
        for _ in 0..TRAIN_STEPMAX {
            game.step();
            if game.is_end() {
                game = game.end();
                continue 'a;
            }
        }
        break 'a;
    }
    game.end()
}

impl OneGame {
    fn step(&mut self) {
        let OneGame {
            remain,
            field,
            chans,
        } = self;
        for chan_opt in chans {
            let Ok(chan) = chan_opt else {
                continue;
            };
            let surr = |diff: (isize, isize)| -> bool {
                let diff = rotate_from_chan(diff, chan.ori);
                let Some((i, j)) = goto(chan.pos, diff) else {
                    return false;
                };
                field.field[i][j]
            };
            let act = chan.act.get_action(surr, chan.mem);
            match act {
                Action::Forword => {
                    let diff = rotate_from_chan((1, 0), chan.ori);
                    if let Some(pos) = goto(chan.pos, diff) {
                        chan.pos = pos;
                    }
                }
                Action::Back => {
                    let diff = rotate_from_chan((1, 0), chan.ori);
                    if let Some(pos) = goto(chan.pos, diff) {
                        chan.pos = pos;
                    }
                }
                Action::TurnRight => {
                    chan.ori = match chan.ori {
                        Ori::U => Ori::R,
                        Ori::D => Ori::L,
                        Ori::R => Ori::D,
                        Ori::L => Ori::U,
                    };
                }
                Action::TurnLeft => {
                    chan.ori = match chan.ori {
                        Ori::U => Ori::L,
                        Ori::D => Ori::R,
                        Ori::R => Ori::U,
                        Ori::L => Ori::D,
                    };
                }
            }
            if field.field[chan.pos.0][chan.pos.1] {
                chan.hp = std::cmp::max(chan.hp + RECOVER, MAXHP);
                field.field[chan.pos.0][chan.pos.1] = false;
            }
            if chan.hp == 0 {
                *chan_opt = Err((*remain, chan.act.clone()));
                *remain -= 1;
            } else {
                chan.hp -= 1;
            }
        }
    }
}

const LEN: usize = 5;

fn onegame_view(
    OneGame {
        remain: _,
        field,
        chans,
    }: &OneGame,
) -> Html {
    html! {
        <>
    <svg version="1.1"
     width={(W * LEN).to_string()} height={(H * LEN).to_string()}
     xmlns="http://www.w3.org/2000/svg">
        {field_view(field)}
        {for chans.iter().filter_map(|chan|{
            if let Ok(chan) = chan {
                Some(chan_view(chan))
            } else {
                None
            }
        })}
        </svg>
        </>
    }
}

fn rect(pos: (usize, usize), diff: (usize, usize), fill: String) -> Html {
    html! {
        <rect y={pos.0.to_string()} x={pos.1.to_string()} height={diff.0.to_string()} width={diff.1.to_string()} fill={fill} />
    }
}

fn field_view(Field { field }: &Field) -> Html {
    let mut v: Vec<Html> = vec![];
    for i in 0..H {
        for j in 0..W {
            if field[i][j] {
                v.push(rect((i * LEN, j * LEN), (LEN, LEN), "blue".to_string()))
            }
        }
    }
    v.into_iter().collect()
}

fn circle(pos: (usize, usize), r: usize, fill: String) -> Html {
    html! {<circle cy={pos.0.to_string()} cx={pos.1.to_string()} r={r.to_string()} fill={fill} />}
}

fn chan_view(
    Chan {
        pos,
        hp,
        act,
        mem,
        ori,
    }: &Chan,
) -> Html {
    circle(
        (pos.0 * LEN + LEN / 2, pos.1 * LEN + LEN / 2),
        LEN / 2,
        "green".to_string(),
    )
}

#[derive(Debug)]
struct App {
    game: OneGame,
    #[allow(dead_code)]
    interval: Interval,
}

enum Msg {
    Tick,
    End,
}

impl Component for App {
    type Message = Msg;
    type Properties = ();
    fn create(ctx: &Context<Self>) -> Self {
        let callback = ctx.link().callback(|_| Msg::Tick);
        let interval = Interval::new(10, move || callback.emit(()));
        Self {
            game: train(),
            interval,
        }
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        if self.game.is_end() {
            ctx.link().send_message(Msg::End);
        }
        onegame_view(&self.game)
    }
    fn update(&mut self, _ctx: &Context<Self>, msg: Self::Message) -> bool {
        match msg {
            Msg::Tick => {
                self.game.step();
                true
            }
            Msg::End => {
                self.game = self.game.end();
                true
            }
        }
    }
}

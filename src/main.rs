use gloo::timers::callback::Interval;
use rand::prelude::*;
use wasm_bindgen::JsValue;
use yew::prelude::*;

// なんか一気にlogが描画されるせいでちゃんと動いてるのかわかりにくい。
fn main() {
    let mut rng = thread_rng();
    let genes = train(&mut rng);
    let document = gloo::utils::document();
    let target_element = document.get_element_by_id("onthis").unwrap();
    yew::Renderer::<App>::with_root_and_props(target_element, StartAppGame { genes }).render();
}

const W: usize = 200;
const H: usize = 150;
const P_FOOD: f64 = 0.1;

const MAXHP: usize = 15;
const RECOVER: usize = 5;

const CHAN_NUM: usize = 50;

const MIX_NUM: usize = 5;
const MUT_NUM: usize = 5;
const MUT_CPY: usize = 2;
const MUT_MUCH: usize = 100;

const TRAIN_STEPMAX: usize = 1000;
const TRAIN_NUM: usize = 20;
const TRAIN_B: usize = 20;

// 表示するときのマスの大きさ
const LEN: usize = 5;

pub fn log<T: AsRef<str>>(str: T) {
    web_sys::console::log_1(&JsValue::from_str(str.as_ref()))
}

#[derive(Debug, Clone, PartialEq, Properties)]
struct Field {
    field: [[bool; W]; H],
}

impl Field {
    fn random(rng: &mut ThreadRng) -> Self {
        let mut field = [[false; W]; H];
        for row in field.iter_mut() {
            for item in row {
                if rng.gen_bool(P_FOOD) {
                    *item = true;
                }
            }
        }
        Self { field }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Chan {
    pos: (usize, usize),
    hp: usize,
    gene: Gene,
    mem: u8,
    ori: Ori,
}

impl Chan {
    fn spawn(gene: Gene, rng: &mut ThreadRng) -> Self {
        let pos = (rng.gen_range(0..H), rng.gen_range(0..W));
        Chan {
            pos,
            hp: MAXHP,
            gene,
            mem: 0,
            ori: Ori::U,
        }
    }
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
    gene: [Action; MOVE_NUM],
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
        &self.gene[index as usize]
    }

    fn random(rng: &mut ThreadRng) -> Self {
        let mut gene = [Action::Forword; MOVE_NUM];
        for gene_i in gene.iter_mut() {
            *gene_i = match rng.gen_range(0..4) {
                0 => Action::Forword,
                1 => Action::Back,
                2 => Action::TurnLeft,
                3 => Action::TurnRight,
                _ => unreachable!(),
            };
        }
        Self { gene }
    }

    fn mutate(&mut self, num: usize, rng: &mut ThreadRng) {
        for _ in 0..num {
            let i = rng.gen_range(0..MOVE_NUM);
            let o = match rng.gen_range(0..4) {
                0 => Action::Forword,
                1 => Action::Back,
                2 => Action::TurnLeft,
                3 => Action::TurnRight,
                _ => unreachable!(),
            };
            self.gene[i] = o;
        }
    }

    fn mix(&mut self, other: &Self, rng: &mut ThreadRng) {
        for i in 0..MOVE_NUM {
            if rng.gen_bool(0.5) {
                self.gene[i] = other.gene[i];
            }
        }
    }
}

fn goto(pos: (usize, usize), diff: (isize, isize)) -> (usize, usize) {
    let i = (pos.0 as isize + diff.0).rem_euclid(H as isize);
    let j = (pos.1 as isize + diff.1).rem_euclid(W as isize);
    (i as usize, j as usize)
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
    fn random(rng: &mut ThreadRng) -> Self {
        let field = Field::random(rng);
        // log("A");
        let chans = (0..CHAN_NUM)
            .map(|_| Ok(Chan::spawn(Gene::random(rng), rng)))
            .collect::<Vec<_>>();
        // log("B");
        Self {
            remain: CHAN_NUM,
            field,
            chans,
        }
    }
    fn from_gene(genes: Vec<Gene>, rng: &mut ThreadRng) -> Self {
        assert_eq!(genes.len(), CHAN_NUM);
        let field = Field::random(rng);
        let chans = genes
            .into_iter()
            .map(|gene| Ok(Chan::spawn(gene, rng)))
            .collect();
        Self {
            remain: CHAN_NUM,
            field,
            chans,
        }
    }

    fn is_end(&self) -> bool {
        self.remain == 0
    }

    fn get_genes_ordered(&self) -> Vec<Gene> {
        let OneGame {
            remain: _,
            field: _,
            chans,
        } = self;

        let mut gene_order: Vec<_> = chans
            .iter()
            .map(|result| match result {
                Ok(chan) => ((0, chan.hp), chan.gene.clone()),
                Err((order, gene)) => ((1, *order), gene.clone()),
            })
            .collect();
        gene_order.sort_by_key(|(o, _)| *o);
        gene_order.into_iter().map(|(_, g)| g).collect()
    }
}

fn next_gene(genes: Vec<Gene>, rng: &mut ThreadRng) -> Vec<Gene> {
    let mut new_genes = Vec::with_capacity(genes.len());

    // top n の mix
    for i in 0..MIX_NUM {
        for j in 0..MIX_NUM {
            if i != j {
                let mut new_gene = genes[i].clone();
                new_gene.mix(&genes[j], rng);
                new_genes.push(new_gene);
            }
        }
    }
    // top n の mutate
    for i in 0..MUT_NUM {
        for _ in 0..MUT_CPY {
            let mut new_gene = genes[i].clone();
            new_gene.mutate(MUT_MUCH, rng);
            new_genes.push(new_gene);
        }
    }

    let l = genes.len() - new_genes.len();
    new_genes.extend(genes[0..l].to_owned());
    new_genes
}

fn train(rng: &mut ThreadRng) -> Vec<Gene> {
    // log("start");
    let mut game = OneGame::random(rng);
    'a: for i in 0..TRAIN_NUM {
        // log(format!("train {i}"));
        for _ in 0..TRAIN_STEPMAX {
            game.step();
            if game.is_end() {
                game = OneGame::from_gene(next_gene(game.get_genes_ordered(), rng), rng);
                continue 'a;
            }
        }
        break 'a;
    }
    game.get_genes_ordered()
}

fn traint_from_gene(genes: Vec<Gene>, num: usize, rng: &mut ThreadRng) -> Vec<Gene> {
    let mut game = OneGame::from_gene(genes, rng);
    'a: for i in 0..num {
        // log(format!("train {i}"));
        for _ in 0..TRAIN_STEPMAX {
            game.step();
            if game.is_end() {
                game = OneGame::from_gene(next_gene(game.get_genes_ordered(), rng), rng);
                continue 'a;
            }
        }
        break 'a;
    }
    game.get_genes_ordered()
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
                let (i, j) = goto(chan.pos, diff);
                field.field[i][j]
            };
            let act = chan.gene.get_action(surr, chan.mem);
            match act {
                Action::Forword => {
                    let diff = rotate_from_chan((1, 0), chan.ori);
                    chan.pos = goto(chan.pos, diff);
                }
                Action::Back => {
                    let diff = rotate_from_chan((-1, 0), chan.ori);
                    chan.pos = goto(chan.pos, diff);
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
                *chan_opt = Err((*remain, chan.gene.clone()));
                *remain -= 1;
            } else {
                chan.hp -= 1;
            }
        }
    }
}

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
        hp: _,
        gene: _,
        mem: _,
        ori: _,
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

#[derive(Debug, Clone, PartialEq, Properties)]
struct StartAppGame {
    genes: Vec<Gene>,
}

impl Component for App {
    type Message = Msg;
    type Properties = StartAppGame;
    fn create(ctx: &Context<Self>) -> Self {
        let callback = ctx.link().callback(|_| Msg::Tick);
        let interval = Interval::new(10, move || callback.emit(()));
        let StartAppGame { genes } = ctx.props();
        let mut rng = thread_rng();
        let game = OneGame::from_gene(genes.to_vec(), &mut rng);
        Self {
            game: game.clone(),
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
                let mut rng = thread_rng();
                let genes = self.game.get_genes_ordered();
                let new_genes = traint_from_gene(genes, TRAIN_B, &mut rng);
                self.game = OneGame::from_gene(new_genes, &mut rng);
                true
            }
        }
    }
}

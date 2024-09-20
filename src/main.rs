use anyhow::bail;
use gloo::file::callbacks::FileReader;
use gloo::storage::Storage;
use gloo::timers::callback::Interval;
use rand::prelude::*;
use serde::Serialize;
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::ops::{Add, AddAssign, Deref, DerefMut, Index, IndexMut};
use std::rc::Rc;
use wasm_bindgen::JsValue;
use yew::prelude::*;

// なんか一気にlogが描画されるせいでちゃんと動いてるのかわかりにくい。
fn main() {
    let document = gloo::utils::document();
    let target_element = document.get_element_by_id("onthis").unwrap();
    yew::Renderer::<App>::with_root_and_props(target_element, ()).render();
}

const W: usize = 200;
const H: usize = 150;
const P_FOOD: f64 = 0.1;

const MAXHP: usize = 15;
const RECOVER: usize = 5;

const CHAN_NUM: usize = 50;

// MIX_NUM * (MIX - 1) + MUT_NUM * MUT_CPY <= CHAN_NUM
const MIX_NUM: usize = 5;
const MUT_NUM: usize = 5;
const MUT_CPY: usize = 2;
const MUT_MUCH: usize = 100;

const TRAIN_STEPMAX: usize = 1000;
const TRAIN_NUM: usize = 10;

const LOCAL_STORAGE_GENE_KEY: &str = "GENE";
const LOCAL_STORAGE_NUM_KEY: &str = "NUM";

// 表示するときのマスの大きさ
const LEN: usize = 5;

pub fn log<T: AsRef<str>>(str: T) {
    web_sys::console::log_1(&JsValue::from_str(str.as_ref()))
}

#[derive(Debug, Clone, PartialEq)]
struct Field {
    field: [[bool; W]; H],
}

impl Field {
    fn random<T>(rng: &mut T) -> Self
    where
        T: BorrowMut<ThreadRng>,
    {
        let mut field = [[false; W]; H];
        for row in field.iter_mut() {
            for item in row {
                if rng.borrow_mut().gen_bool(P_FOOD) {
                    *item = true;
                }
            }
        }
        Self { field }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Pos(usize, usize);

impl Pos {
    fn generate_random<T>(rng: &mut T) -> Self
    where
        T: BorrowMut<ThreadRng>,
    {
        Pos(
            rng.borrow_mut().gen_range(0..H),
            rng.borrow_mut().gen_range(0..W),
        )
    }
}

impl Add<(isize, isize)> for Pos {
    type Output = Self;
    fn add(self, rhs: (isize, isize)) -> Self::Output {
        Pos(
            (self.0 as isize + rhs.0).rem_euclid(H as isize) as usize,
            (self.1 as isize + rhs.1).rem_euclid(W as isize) as usize,
        )
    }
}

impl AddAssign<(isize, isize)> for Pos {
    fn add_assign(&mut self, rhs: (isize, isize)) {
        *self = *self + rhs;
    }
}

impl Index<Pos> for Field {
    type Output = bool;
    fn index(&self, index: Pos) -> &Self::Output {
        &self.field[index.0][index.1]
    }
}

impl IndexMut<Pos> for Field {
    fn index_mut(&mut self, index: Pos) -> &mut Self::Output {
        &mut self.field[index.0][index.1]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Mem(bool, bool, bool, bool, bool, bool);

impl From<Mem> for u8 {
    fn from(value: Mem) -> Self {
        let f = |b: bool, i: u8| -> u8 {
            if b {
                1 << i
            } else {
                0
            }
        };
        f(value.0, 0)
            + f(value.1, 1)
            + f(value.2, 2)
            + f(value.3, 3)
            + f(value.4, 4)
            + f(value.5, 5)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Move {
    Forword,
    Back,
    TurnRight,
    TurnLeft,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Action(Mem, Move);

impl Action {
    fn generate_random<T>(rng: &mut T) -> Self
    where
        T: BorrowMut<ThreadRng>,
    {
        let u: u8 = rng.borrow_mut().gen();
        u.into()
    }
}

impl From<Action> for u8 {
    fn from(Action(value0, value1): Action) -> Self {
        let v0: u8 = value0.into();
        let v1 = match value1 {
            Move::Forword => 0,
            Move::Back => 1,
            Move::TurnRight => 2,
            Move::TurnLeft => 3,
        } * (1 << 6);
        v0 + v1
    }
}

impl From<u8> for Action {
    fn from(value: u8) -> Self {
        let f = |i: u8| -> bool { value & (1 << i) != 0 };
        let move_action = match (f(6), f(7)) {
            (false, false) => Move::Forword,
            (true, false) => Move::Back,
            (false, true) => Move::TurnRight,
            (true, true) => Move::TurnLeft,
        };
        let mem = Mem(f(0), f(1), f(2), f(3), f(4), f(5));
        Action(mem, move_action)
    }
}

const MOVE_NUM: usize = 2_usize.pow((3 * 3) as u32 - 1) * (2_usize.pow(6));

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
struct Gene(#[serde(with = "serde_bytes")] Vec<u8>);

impl Gene {
    fn get_action<F>(&self, surr: F, mem: Mem) -> Action
    where
        F: Fn((isize, isize)) -> bool,
    {
        // x\y -1  0  1
        // -1   7  6  5
        //  0   4  x  0
        //  1   1  2  3
        let mut surr_index: u8 = 0;
        for i in -1..=1 {
            for j in -1..=1 {
                if (i != 0 && j != 0) && surr((i, j)) {
                    let s = 3 * i + j;
                    let s = if s < 0 {
                        (-s + 3) as u32
                    } else {
                        (s - 1) as u32
                    };
                    surr_index |= 1_u8.checked_shl(s).unwrap();
                }
            }
        }
        let index: u16 = surr_index as u16 + u8::from(mem) as u16;
        let u = self.0[index as usize];
        u.into()
    }

    fn generate_random<T>(rng: &mut T) -> Self
    where
        T: BorrowMut<ThreadRng>,
    {
        let gene = (0..MOVE_NUM).map(|_| rng.borrow_mut().gen()).collect();
        Self(gene)
    }

    fn mutate<T>(&mut self, rng: &mut T)
    where
        T: BorrowMut<ThreadRng>,
    {
        for _ in 0..MUT_MUCH {
            let i = rng.borrow_mut().gen_range(0..MOVE_NUM);
            self.0[i] = rng.borrow_mut().gen();
        }
    }

    fn mix<T>(&mut self, other: &Self, rng: &mut T)
    where
        T: BorrowMut<ThreadRng>,
    {
        for i in 0..MOVE_NUM {
            if rng.borrow_mut().gen_bool(0.5) {
                self.0[i] = other.0[i];
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Chan {
    pos: Pos,
    hp: usize,
    gene: Gene,
    mem: Mem,
    ori: Ori,
}

impl Chan {
    fn spawn<T>(gene: Gene, rng: &mut T) -> Self
    where
        T: BorrowMut<ThreadRng>,
    {
        Chan {
            pos: Pos::generate_random(rng),
            hp: MAXHP,
            gene,
            mem: Mem::default(),
            ori: Ori::U,
        }
    }

    fn fix_act(&mut self, Action(mem, move_action): Action) {
        match move_action {
            Move::Forword => {
                let diff = rotate_from_chan((-1, 0), self.ori);
                self.pos += diff;
            }
            Move::Back => {
                let diff = rotate_from_chan((1, 0), self.ori);
                self.pos += diff;
            }
            Move::TurnRight => {
                self.ori = match self.ori {
                    Ori::U => Ori::R,
                    Ori::D => Ori::L,
                    Ori::R => Ori::D,
                    Ori::L => Ori::U,
                };
            }
            Move::TurnLeft => {
                self.ori = match self.ori {
                    Ori::U => Ori::L,
                    Ori::D => Ori::R,
                    Ori::R => Ori::U,
                    Ori::L => Ori::D,
                };
            }
        }
        self.mem = mem;
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Ori {
    U,
    D,
    R,
    L,
}

fn rotate_from_chan(diff: (isize, isize), ori: Ori) -> (isize, isize) {
    match ori {
        Ori::U => diff,
        Ori::D => (-diff.0, -diff.1),
        Ori::R => (diff.1, -diff.0),
        Ori::L => (-diff.1, diff.0),
    }
}

#[derive(Debug, Clone, PartialEq, Properties)]
struct OneGame {
    remain: usize,
    field: Field,
    chans: Vec<Result<Chan, (usize, Gene)>>,
}

impl OneGame {
    fn random<T>(rng: &mut T) -> Self
    where
        T: BorrowMut<ThreadRng>,
    {
        let field = Field::random(rng);
        // log("A");
        let chans = (0..CHAN_NUM)
            .map(|_| Ok(Chan::spawn(Gene::generate_random(rng), rng)))
            .collect::<Vec<_>>();
        // log("B");
        Self {
            remain: CHAN_NUM,
            field,
            chans,
        }
    }

    fn from_gene<T>(genes: Vec<Gene>, rng: &mut T) -> Self
    where
        T: BorrowMut<ThreadRng>,
    {
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

fn next_gene<T>(genes: Vec<Gene>, rng: &mut T) -> Vec<Gene>
where
    T: BorrowMut<ThreadRng>,
{
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
            new_gene.mutate(rng);
            new_genes.push(new_gene);
        }
    }

    let l = genes.len() - new_genes.len();
    new_genes.extend(genes[0..l].to_owned());
    new_genes
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
                field[chan.pos + diff]
            };
            let act = chan.gene.get_action(surr, chan.mem);
            chan.fix_act(act);
            if field[chan.pos] {
                chan.hp = std::cmp::max(chan.hp + RECOVER, MAXHP);
                field[chan.pos] = false;
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

fn train<T>(rng: &mut T) -> Vec<Gene>
where
    T: BorrowMut<ThreadRng>,
{
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

fn train_from_gene<T>(genes: Vec<Gene>, num: usize, rng: &mut T) -> Vec<Gene>
where
    T: BorrowMut<ThreadRng>,
{
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

struct StartScene {}

#[derive(Debug, Clone, PartialEq, Properties)]
struct StartSceneProps {
    on_choose: Callback<(usize, Vec<Gene>)>,
}

enum StartSceneMsg {
    Random,
    Read,
}

impl Component for StartScene {
    type Message = StartSceneMsg;
    type Properties = StartSceneProps;
    fn create(ctx: &Context<Self>) -> Self {
        Self {}
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        let on_ramdom = ctx.link().callback(|_: MouseEvent| StartSceneMsg::Random);
        let on_read_localstorage = ctx.link().callback(|_: MouseEvent| StartSceneMsg::Read);
        html! {
            <>
            <button onclick={on_ramdom}> {"random"} </button>
            <button onclick={on_read_localstorage}> {"read from local storage"} </button>
            </>
        }
    }
    fn update(&mut self, ctx: &Context<Self>, msg: Self::Message) -> bool {
        let StartSceneProps { on_choose } = ctx.props();
        let mut rng = thread_rng();
        let res = read();
        if let Err(ref err) = &res {
            log(format!("{err}"))
        }
        let n = match (msg, res) {
            (StartSceneMsg::Read, Ok((u, gene))) => {
                let mut v = Vec::with_capacity(CHAN_NUM);
                for _ in 0..CHAN_NUM {
                    let mut gene: Gene = gene.clone();
                    gene.mutate(&mut rng);
                    v.push(gene)
                }
                (u, v)
            }
            _ => (
                0,
                (0..CHAN_NUM)
                    .map(|_| Gene::generate_random(&mut rng))
                    .collect(),
            ),
        };
        on_choose.emit(n);
        true
    }
}

struct GameWatchScene {
    game: OneGame,
    interval: Interval,
}

#[derive(Debug, Clone, PartialEq, Properties)]
struct GameWatchProps {
    genes: Vec<Gene>,
    on_end: Callback<()>,
}

#[derive(Debug, Clone, PartialEq)]
enum GameWatchMsg {
    Tick,
}

impl Component for GameWatchScene {
    type Message = GameWatchMsg;
    type Properties = GameWatchProps;
    fn create(ctx: &Context<Self>) -> Self {
        let GameWatchProps { genes, on_end: _ } = ctx.props();
        let mut rng = thread_rng();
        let callback = ctx.link().callback(|_| GameWatchMsg::Tick);
        let interval = Interval::new(10, move || callback.emit(()));
        let game = OneGame::from_gene(genes.to_vec(), &mut rng);
        Self { game, interval }
    }
    fn view(&self, _ctx: &Context<Self>) -> Html {
        html! {
            <>
            {onegame_view(&self.game)}
            </>
        }
    }
    fn update(&mut self, ctx: &Context<Self>, msg: Self::Message) -> bool {
        match msg {
            GameWatchMsg::Tick => {
                self.game.step();
                if self.game.is_end() {
                    let GameWatchProps { genes: _, on_end } = ctx.props();
                    on_end.emit(());
                }
                true
            }
        }
    }
}

fn save((u, gene): (usize, Gene)) -> Result<(), gloo::storage::errors::StorageError> {
    gloo::storage::LocalStorage::set(LOCAL_STORAGE_NUM_KEY, u)?;
    gloo::storage::LocalStorage::set(LOCAL_STORAGE_GENE_KEY, gene)?;
    Ok(())
}

fn read() -> Result<(usize, Gene), gloo::storage::errors::StorageError> {
    let u = gloo::storage::LocalStorage::get(LOCAL_STORAGE_NUM_KEY)?;
    let gene = gloo::storage::LocalStorage::get(LOCAL_STORAGE_GENE_KEY)?;
    Ok((u, gene))
}
struct TrainScene {
    // start_genes: serde_json::Value,
    genes: Vec<Gene>,
    num: usize,
    interval: Interval,
}

#[derive(Debug, Clone, PartialEq, Properties)]
struct TrainProps {
    start_genes: Vec<Gene>,
    on_train_end: Callback<Vec<Gene>>,
    train_num: usize,
}

enum TrainMsg {
    TrainStart,
}

impl Component for TrainScene {
    type Message = TrainMsg;
    type Properties = TrainProps;
    fn create(ctx: &Context<Self>) -> Self {
        let callback = ctx.link().callback(|_| TrainMsg::TrainStart);
        let interval = Interval::new(10, move || callback.emit(()));
        let TrainProps {
            start_genes,
            on_train_end,
            train_num,
        } = ctx.props();
        // let json = serde_json::to_value(start_genes).unwrap();
        Self {
            // start_genes: json,
            genes: start_genes.clone(),
            num: 0,
            interval,
        }
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        let TrainProps {
            start_genes,
            on_train_end,
            train_num,
        } = ctx.props();
        html! {
            <>
            {self.num} {"/"} {train_num}
            // <JsonFileSaveView json_value={self.start_genes.clone()} />
            </>
        }
    }
    fn update(&mut self, ctx: &Context<Self>, msg: Self::Message) -> bool {
        match msg {
            TrainMsg::TrainStart => {
                let TrainProps {
                    start_genes,
                    on_train_end,
                    train_num,
                } = ctx.props();
                let mut rng = thread_rng();
                self.genes = train_from_gene(self.genes.clone(), 1, &mut rng);
                self.num += 1;
                if self.num == *train_num {
                    on_train_end.emit(self.genes.clone())
                }
                true
            }
        }
    }
}

#[derive(Debug)]
enum Scene {
    Start,
    Game,
    Train,
}

#[derive(Debug)]
struct App {
    genes: Vec<Gene>,
    num: usize,
    scene: Scene,
}

enum Msg {
    Start((usize, Vec<Gene>)),
    GameEnd,
    TrainEnd(Vec<Gene>),
}

impl Component for App {
    type Message = Msg;
    type Properties = ();
    fn create(ctx: &Context<Self>) -> Self {
        let mut rng = thread_rng();
        let genes = (0..CHAN_NUM)
            .map(|_| Gene::generate_random(&mut rng))
            .collect();
        Self {
            genes,
            num: 0,
            scene: Scene::Start,
        }
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        match self.scene {
            Scene::Start => {
                let on_choose = ctx.link().callback(Msg::Start);
                html! {
                    <StartScene on_choose={on_choose}/>
                }
            }
            Scene::Game => {
                let on_end = ctx.link().callback(|_| Msg::GameEnd);
                html! {
                    <GameWatchScene
                        genes={self.genes.clone()}
                        on_end={on_end}
                    />
                }
            }
            Scene::Train => {
                let on_train_end = ctx.link().callback(Msg::TrainEnd);
                html! {
                    <TrainScene
                        start_genes={self.genes.clone()}
                        on_train_end={on_train_end}
                        train_num={TRAIN_NUM}
                    />
                }
            }
        }
    }
    fn update(&mut self, _ctx: &Context<Self>, msg: Self::Message) -> bool {
        match msg {
            Msg::Start((num, genes)) => {
                self.num = num;
                self.genes = genes;
                self.scene = Scene::Game;
            }
            Msg::GameEnd => {
                if let Err(err) = save((self.num, self.genes[0].clone())) {
                    log(format!("{err}"))
                };
                self.scene = Scene::Train;
            }
            Msg::TrainEnd(genes) => {
                self.num += 1;
                if let Err(err) = save((self.num, self.genes[0].clone())) {
                    log(format!("{err}"))
                };
                log(format!("{}", self.num));
                self.genes = genes;
                self.scene = Scene::Game;
            }
        }
        true
    }
}

// #[derive(Debug, Clone, PartialEq, Properties)]
// pub struct JsonFileSaveProps {
//     pub json_value: serde_json::Value,
// }

// #[function_component(JsonFileSaveView)]
// pub fn json_file_save_view(JsonFileSaveProps { json_value }: &JsonFileSaveProps) -> Html {
//     let head_string = "data:text/json;charset=utf-8,";
//     let data = json_value.to_string();
//     html! {
//         <a href={format!("{}{}", head_string, data)} download="data.json"> {"save as json"}</a>
//     }
// }

// pub enum JsonFileReadMsg {
//     Read(DragEvent),
//     LoadEnd(Result<String, anyhow::Error>),
// }

// #[derive(Debug, Clone, PartialEq, Properties)]
// pub struct JsonFileReadProps {
//     pub on_drop_json: Callback<serde_json::Value>,
// }

// #[derive(Debug, Default)]
// pub struct JsonFileReadView {
//     reader: Option<FileReader>,
// }

// impl Component for JsonFileReadView {
//     type Message = JsonFileReadMsg;
//     type Properties = JsonFileReadProps;
//     fn create(_ctx: &Context<Self>) -> Self {
//         Self::default()
//     }
//     fn view(&self, ctx: &Context<Self>) -> Html {
//         html! {
//             <>
//             <div id="drop-container"
//                 ondrop={ctx.link().callback(|event: DragEvent|{
//                     event.prevent_default();
//                     JsonFileReadMsg::Read(event)
//                 })}
//                 ondragover={Callback::from(|event: DragEvent| {
//                     event.prevent_default();
//                 })}
//                 ondragenter={Callback::from(|event: DragEvent| {
//                     event.prevent_default();
//                 })}
//             > <p> {"drop here"} </p> </div>
//             </>
//         }
//     }
//     fn update(&mut self, ctx: &Context<Self>, msg: Self::Message) -> bool {
//         match msg {
//             JsonFileReadMsg::Read(dragevent) => {
//                 let read = move |event: DragEvent| -> Result<FileReader, anyhow::Error> {
//                     let Some(data_transfer) = event.data_transfer() else {
//                         bail!("data transfer fail")
//                     };
//                     let Some(files) = data_transfer.files() else {
//                         bail!("files fail")
//                     };
//                     let Some(file) = files.get(0) else {
//                         bail!("file fail")
//                     };
//                     let file: gloo::file::File = file.into();
//                     let link = ctx.link().clone();
//                     let task = gloo::file::callbacks::read_as_text(&file, move |res| {
//                         link.send_message(JsonFileReadMsg::LoadEnd(res.map_err(|e| e.into())))
//                     });
//                     Ok(task)
//                 };
//                 match read(dragevent) {
//                     Ok(task) => {
//                         self.reader = Some(task);
//                     }
//                     Err(err) => {
//                         log(format!("{err:?}"));
//                     }
//                 }
//                 true
//             }
//             JsonFileReadMsg::LoadEnd(res) => {
//                 match res {
//                     Ok(string) => match serde_json::from_str(&string) {
//                         Ok(val) => {
//                             ctx.props().on_drop_json.emit(val);
//                         }
//                         Err(err) => {
//                             log(format!("{err:?}"));
//                         }
//                     },
//                     Err(err) => {
//                         log(format!("{err:?}"));
//                     }
//                 }
//                 true
//             }
//         }
//     }
// }

#[cfg(test)]
mod tests {
    use std::u8;

    use super::*;
    #[test]
    fn test1() {
        let pos0 = Pos(0, 0);
        let pos1 = Pos(H - 1, W - 1);
        assert_eq!(pos0 + (-1, -1), pos1);
        assert_eq!(pos1 + (1, 1), pos0);

        // x\y -1, 0, 1
        // -1
        //  0      ↑
        //  1
        let diff = (1, 0);
        assert_eq!(rotate_from_chan(diff, Ori::U), (1, 0));
        assert_eq!(rotate_from_chan(diff, Ori::D), (-1, 0));
        assert_eq!(rotate_from_chan(diff, Ori::L), (0, 1));
        assert_eq!(rotate_from_chan(diff, Ori::R), (0, -1));
    }
    #[test]
    fn memact() {
        for i in 0..u8::MAX {
            let act: Action = i.into();
            let u: u8 = act.into();
            assert_eq!(i, u);
        }
    }
    #[test]
    fn train_test() {
        let a: u8 = u8::MAX;
        let b: u8 = u8::MAX;
        let c: u16 = (a as u16) + ((b as u16) << 8);
        let mut rng = thread_rng();
        train(&mut rng);
    }
}

use std::time::Instant;

pub mod cli;
use cli::*;

pub mod count;
pub mod consts;

fn main() {

    println!("[TOOL NAME] Intrahost Variant\n");

    let start = Instant::now();

    let args = cli::parse_args();
    match args.mode {
        Mode::Count(count_args) => count::count(count_args),
    }

    println!("Hello, world!");

    let end = Instant::now();
    eprintln!("\n\nDone in {}s", end.duration_since(start).as_secs());

}

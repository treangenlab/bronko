use std::time::Instant;

pub mod cli;
use cli::*;

use crate::consts::BRONKO_VERSION;

pub mod consts;
pub mod call;
pub mod io;
pub mod lcb;
pub mod build;
pub mod util;

fn main() {
    println!("bronko v{}", BRONKO_VERSION);
    println!("Developed by Ryan Doughty (Rice University)");
    println!("Correspondence: rdd4@rice.edu, treangen@rice.edu\n");

    let start = Instant::now();

    let args = cli::parse_args();
    match args.mode {
        Mode::Call(call_args) => call::call(call_args),
        Mode::Build(build_args) => build::build(build_args),
    }

    let end = Instant::now();
    eprintln!("\nbronko v{} finished in {}s", BRONKO_VERSION, end.duration_since(start).as_secs_f32());
}

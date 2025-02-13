use crate::consts::*;
use clap::{Args, Subcommand, Parser};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode{
    //end-to-end --> runs everything
    Count(CountArgs)
}

#[derive(Args, Default)]
pub struct CountArgs{
    #[clap(short, long="input", help_heading = "INPUT", help="Input fastq(.gz) file")]
    pub input: String,

    #[clap(short, long="kmer-size", default_value_t = DEFAULT_KMER_SIZE, help_heading="KMER", help="Kmer size")]
    pub kmer: u32,

    //Debug mode
    #[clap(long="debug", help="Debug output")]
    pub debug: bool,

    //Verbose mode (prints most checkpoints)
    #[clap(long="verbose", help="Verbose output (warning: very verbose)")]
    pub verbose: bool,
}


pub fn parse_args() -> Cli {
    Cli::parse()
}


use crate::consts::*;
use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Cli {
    #[clap(subcommand)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    //end-to-end --> runs everything
    Build(BuildArgs), //Build --> parses a genome(s) and builds the data strucutre that we want to map to later
    Query(QueryArgs), //Query --> takes the genomes and sequencing data and does the intrahost variation analysis
}


#[derive(Args, Default)]
pub struct BuildArgs {
    // the fasta file to be used as a reference
    #[clap(
        short,
        long = "genome",
        help_heading = "INPUT",
        help = "Genome fasta(.gz) file"
    )]
    pub genome: String,

    // the kmer size
    #[clap(short, long="kmer-size", default_value_t = DEFAULT_KMER_SIZE, help_heading="KMER", help="Kmer size")]
    pub kmer: usize,

    //Number of threads
    #[clap(short, long="threads", default_value_t=1, help="Number of threads")]
    pub threads: usize,

    //Debug mode
    #[clap(long = "debug", help = "Debug output")]
    pub debug: bool,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

#[derive(Args)]
pub struct QueryArgs {

    // REFERENCE INPUT
    #[clap(num_args=1.., short='g', long = "genomes", help_heading = "REFERENCE INPUT", help = "Genome fasta(.gz) files to use as references")]
    pub genomes: Vec<String>,

    // todo: add database version (pre built)

    // SEQUENCING DATA INPUT
    #[clap(num_args=1.., short = 'r', long = "reads", help_heading = "READS INPUT", help = "Input single-end reads (fastq/gzip)")]
    pub reads: Vec<String>,

    #[clap(num_args=1..,short='1', long="first-pairs", help_heading = "READS INPUT", help = "First pairs for raw paired-end reads (fastq/gzip)")]
    pub first_pair: Vec<String>,

    #[clap(num_args=1.., short='2', long="second-pairs", help_heading="READS INPUT", help = "Second pairs for raw paired-end reads (fastq/gzip)")]
    pub second_pair: Vec<String>,

    //GENERAL ALGORITHM PARAMETERS
    //kmer size
    #[clap(short, long="kmer-size", default_value_t = DEFAULT_KMER_SIZE, help_heading="ALGORITHM", help="Kmer size used for analysis")]
    pub kmer: usize,

    //minimum kmers in initial filtering
    #[clap(long="min-kmers", default_value_t = MIN_KMER_COUNT,  help_heading="ALGORITHM", help="Minimum times a kmer must occur in sequencing data to be used")]
    pub min_kmers: usize,

    //VARIANT CALLING PARAMETERS
    #[clap(long="min-af", default_value_t = DEFAULT_MIN_AF, help_heading="VARIANT CALLING PARAMETERS", help="Minimum minor allele frequency to be reported")]
    pub min_af: f64,

    #[clap(long="no-end-filter", default_value_t = DEFAULT_NO_FILTER_ENDS, help_heading="VARIANT CALLING PARAMETERS", help="Do not filter variants from beginning and end k bases of each segment")]
    pub no_end_filter: bool,

    #[clap(long="no-strand-filter", default_value_t = DEFAULT_NO_STRAND_FILTER, help_heading="VARIANT CALLING PARAMETERS", help="Do not filter variants that are present on one strand but not the other")]
    pub no_strand_filter: bool,

    #[clap(long="n-per-strand", default_value_t = DEFAULT_N_KMERS_PER_STRAND, help_heading="VARIANT CALLING PARAMETERS", help="Minimum number of unique kmers to observe on each strand to call a variant at any site")]
    pub n_per_strand: usize,

    //OUTPUT PARAMETERS
    //todo add output locations, output formats
    #[clap(long="output", help_heading="OUTPUT", help="Folder to output all resulting files")]
    pub output: String,

    // OTHER PARAMETERS
    //Number of threads
    #[clap(short, long="threads", default_value_t=1, help="Number of threads")]
    pub threads: usize,

    //Debug mode
    #[clap(long = "debug", help = "Debug output")]
    pub debug: bool,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

pub fn parse_args() -> Cli {
    Cli::parse()
}

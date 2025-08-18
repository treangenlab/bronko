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
    pub first_pairs: Vec<String>,

    #[clap(num_args=1.., short='2', long="second-pairs", help_heading="READS INPUT", help = "Second pairs for raw paired-end reads (fastq/gzip)")]
    pub second_pairs: Vec<String>,

    //GENERAL ALGORITHM PARAMETERS
    //kmer size
    #[clap(short, long="kmer-size", default_value_t = DEFAULT_KMER_SIZE, help_heading="ALGORITHM", help="Kmer size used for analysis")]
    pub kmer: usize,

    //minimum kmers in initial filtering
    #[clap(long="min-kmers", default_value_t = MIN_KMER_COUNT,  help_heading="ALGORITHM", help="Minimum times a kmer must occur in sequencing data to be used")]
    pub min_kmers: usize,

    //using all the buckets for a kmer
    #[clap(long="use-full-kmer", default_value_t = DEFAULT_USE_FULL_KMER, help_heading="ALGORITHM", help="Use the entire kmer length for variant positions rather than having [--n-fixed] bases on each end")]
    pub use_full_kmer: bool,

    //number of buckets to ignore on ends of kmers 
    #[clap(long="n-fixed", default_value_t = DEFAULT_N_FIXED, help_heading="ALGORITHM", help="Number of fied positions at the end of each kmer")]
    pub n_fixed: usize,

    //VARIANT CALLING PARAMETERS
    //minimum allele frequency to be reported
    #[clap(long="min-af", default_value_t = DEFAULT_MIN_AF, help_heading="VARIANT CALLING PARAMETERS", help="Minimum minor allele frequency to be reported")]
    pub min_af: f64,

    //do not filter variants that are in the first and last k bases in any given reference segment
    #[clap(long="no-end-filter", default_value_t = DEFAULT_NO_FILTER_ENDS, help_heading="VARIANT CALLING PARAMETERS", help="Do not filter variants from beginning and end k bases of each segment")]
    pub no_end_filter: bool,

    //do not do strand filtering 
    #[clap(long="no-strand-filter", default_value_t = DEFAULT_NO_STRAND_FILTER, help_heading="VARIANT CALLING PARAMETERS", help="Do not filter variants that are present on one strand but not the other")]
    pub no_strand_filter: bool,

    //the number of kmers per strand that are required to call a variant
    #[clap(long="n-per-strand", default_value_t = DEFAULT_N_KMERS_PER_STRAND, help_heading="VARIANT CALLING PARAMETERS", help="Min number of unique kmers to observe to call a variant at any site (needed on both strands if strand filter active)")]
    pub n_per_strand: usize,

    #[clap(long="min_strand_difference", default_value_t = DEFAULT_PERCENT_STRAND_DIFF, help_heading="VARIANT CALLING PARAMETERS", help="Minimum percent one strand's depth must be of the total depth in order to bypass strand filtering")]
    pub min_strand_diff: f64,

    //min depth to call a variant
    #[clap(long="min-depth", default_value_t = DEFAULT_MIN_DEPTH, help_heading="VARIANT CALLING PARAMETERS", help="Minimum depth at an allele to call a minor variant (default=100*min_kmers)")]
    pub min_depth: usize,

    //OUTPUT PARAMETERS
    //todo add output locations, output formats
    #[clap(long="output", help_heading="OUTPUT", help="Folder to output all resulting files")]
    pub output: String,

    //EXTRA OUTPUT FORMATS
    #[clap(long="pileup", default_value_t = DEFAULT_TSV_PILEUP, help_heading="EXTRA OUTPUT FILES", help="Also output a tsv of the approximate pileup for each sample and reference")]
    pub output_pileup: bool,

    #[clap(long="alignment", default_value_t = DEFAULT_ALIGNMENT, help_heading="EXTRA OUTPUT FILES", help="Output an multifasta containing the alignment of all samples to the reference and themselves (only major variant positions)")]
    pub output_alignment: bool,

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

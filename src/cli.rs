use crate::consts::*;
use clap::{Args, Parser, Subcommand};
use clap::builder::styling::{AnsiColor, Effects, Styles};


fn custom_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Blue.on_default() | Effects::BOLD)
        .usage(AnsiColor::White.on_default() | Effects::BOLD)
        .literal(AnsiColor::White.on_default() | Effects::BOLD)
        .placeholder(AnsiColor::White.on_default())
}

#[derive(Parser)]
#[command(author, version, about, long_about = None, styles=custom_styles())]
#[command(propagate_version = true)]
pub struct Cli {
    #[clap(subcommand)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    Build(BuildArgs), //Build --> parses a genome(s) and builds the data strucutre that we want to map to later
    Call(CallArgs), //Call --> takes the genomes and sequencing data and does the viral variation analysis
}

#[derive(Args, Default)]
#[clap(about="Create an bronko index of existing viral references for a given species", arg_required_else_help = true)]
pub struct BuildArgs {

    //REFERENCE INPUT
    // the fasta file to be used as a reference
    #[clap(num_args=1.., short='g', long = "genomes", help_heading = "REFERENCE INPUT", help = "Genome files to be built into index (fasta/gzip)")]
    pub genomes: Vec<String>,

    #[clap(long="file-input", help_heading="REFERENCE INPUT", help = "Path to .txt file containing paths to each file, one line per file, each containing a single genome")]
    pub file_input: Option<String>,

    // INDEXING PARAMETERS
    // the kmer size
    #[clap(short, long="kmer-size", default_value_t = DEFAULT_KMER_SIZE, help_heading="KMER", help="Kmer size")]
    pub kmer: usize,

    //OUTPUT 
    #[clap(short='o', long="output", default_value = DEFAULT_INDEX_OUTPUT, help_heading="OUTPUT", help="Name of index file (.bkdb will be added)")]
    pub output: String,

    // OTHER PARAMETERS
    //Number of threads
    #[clap(short, long="threads", default_value_t=4, help="Number of threads")]
    pub threads: usize,

    //Debug mode
    #[clap(long = "debug", help = "Debug output")]
    pub debug: bool,

    //Verbose mode (prints most checkpoints)
    #[clap(long = "verbose", help = "Verbose output (warning: very verbose)")]
    pub verbose: bool,
}

#[derive(Args)]
#[clap(about="Perform rapid viral variant calling of viral sequencing data", arg_required_else_help = true)]
pub struct CallArgs {

    // REFERENCE INPUT
    #[clap(num_args=1.., short='g', long = "genomes", help_heading = "REFERENCE INPUT", help = "Genome fasta(.gz) files to use as references (bronko build will be called)")]
    pub genomes: Option<Vec<String>>,

    #[clap(short='d', long= "db", help_heading="REFERENCE INPUT", help = "Use a prebuilt bronko db (.bkdb) of genomes of interest")]
    pub db: Option<String>, 

    // SEQUENCING DATA INPUT
    #[clap(num_args=1.., short = 'r', long = "reads", help_heading = "READS INPUT", help = "Input single-end reads (fastq/gzip)")]
    pub reads: Vec<String>,

    #[clap(num_args=1..,short='1', long="first-pairs", help_heading = "READS INPUT", help = "First pairs for raw paired-end reads (fastq/gzip)")]
    pub first_pairs: Vec<String>,

    #[clap(num_args=1.., short='2', long="second-pairs", help_heading="READS INPUT", help = "Second pairs for raw paired-end reads (fastq/gzip)")]
    pub second_pairs: Vec<String>,

    #[clap(long="file-input", help_heading="READS INPUT", help = "Path to .txt file containing paths to each fastq file, one line per file (paired end reads should be tab delimited on same line)")]
    pub file_input: Option<String>,

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
    #[clap(long="n-fixed", default_value_t = DEFAULT_N_FIXED, help_heading="ALGORITHM", help="Number of fixed positions at the end of each kmer that cannot contribute to pileup")]
    pub n_fixed: usize,

    //specify pattern for position skipping; 1 to include, 0 to skip
    #[clap(long="specify-pattern", help_heading="ALGORITHM", help="Explicit 1/0 keep-mask matching kmer size exactly, e.g. k=15, pattern=\"110000000000001\" -> keeps only the first two and last position. Overrides --n-fixed/--use-full-kmer if set.")]
    pub specify_pattern: Option<String>,

    //VARIANT CALLING PARAMETERS
    //minimum allele frequency to be reported
    #[clap(long="min-af", default_value_t = DEFAULT_MIN_AF, help_heading="VARIANT CALLING PARAMETERS", help="Minimum minor allele frequency to be reported")]
    pub min_af: f64,

    //do not filter variants that are in the first and last k bases in any given reference segment
    #[clap(long="no-end-filter", default_value_t = DEFAULT_NO_FILTER_ENDS, help_heading="VARIANT CALLING PARAMETERS", help="Do not filter variants from beginning and end k bases of each segment")]
    pub no_end_filter: bool,

    //do not do strand filtering 
    #[clap(long="no-strand-filter", default_value_t = DEFAULT_NO_STRAND_FILTER, help_heading="VARIANT CALLING PARAMETERS", help="Do not utilize SOR test to filter variants that are present on one strand but not the other")]
    pub no_strand_filter: bool,

    //do not filter variants with extreme strand balance issues
    #[clap(long="no-strand-balance-filter", default_value_t = DEFAULT_NO_STRAND_BALANCE_FILTER, help_heading="VARIANT CALLING PARAMETERS", help="Allow variants with extreme strand disbalance pass without SOR check (will not matter if --no-strand-filter present)")]
    pub no_strand_balance_filter: bool,

    //strand balance ratio for strand to be ignored
    #[clap(long="balance-ratio", default_value_t = DEFAULT_STRAND_BALANCE_RATIO, help_heading="VARIANT CALLING PARAMETERS", help="Percent of total depth that one strand must be under to be considered unbalanced (must be [0.0-1.0])")]
    pub strand_balance_ratio: f64,

    //the number of kmers per strand that are required to call a variant
    #[clap(long="n-per-strand", default_value_t = DEFAULT_N_KMERS_PER_STRAND, help_heading="VARIANT CALLING PARAMETERS", help="Min number of unique kmers to observe to call a variant at any site (needed on both strands if strand filter active)")]
    pub n_per_strand: usize,

    #[clap(long="strand_odds", default_value_t = DEFAULT_MAX_STRAND_ODDS, help_heading="VARIANT CALLING PARAMETERS", help="Maximum strand odds ratio for a variant to pass strand filtering")]
    pub strand_odds_max: f64,

    //min depth to call a variant
    #[clap(long="min-depth", default_value_t = DEFAULT_MIN_DEPTH, help_heading="VARIANT CALLING PARAMETERS", help="Minimum total depth at an allele to call a minor variant (default=100*min_kmers)")]
    pub min_depth: usize,

    //minimum reads/kmers supporting a minor variant
    #[clap(long="min-variant-depth", default_value_t = MIN_KMER_COUNT, help_heading = "VARIANT CALLING PARAMETERS", help="Minimum depth of a minor variant to be called present (default=min_kmers)")]
    pub min_variant_depth: usize,

    //noise multiplier
    #[clap(long="noise-multiplier", default_value_t = DEFAULT_NOISE_MULTIPLIER, help_heading = "VARIANT CALLING PARAMETERS", help="How much greater (1x, 1.5x, etc) the minor allele frequency of a variant must be above estimated baseline noise in that region (must be > 1.0x). Note that for variants under 1%, multiplier will be increased exponentially up to +0.5 more")]
    pub variant_multiplier: f64,

    //INDEL PARAMETERS
    //disable indel reconstruction entirely (consensus will not contain indels)
    #[clap(long="no-indels", default_value_t = DEFAULT_NO_INDELS, help_heading="INDEL DETECTION", help="Disable indel detection and reconstruction (the consensus sequence will not contain indels)")]
    pub no_indels: bool,

    //max ratio of low/high total depth between neighbors to flag a breakpoint
    #[clap(long="indel-max-ratio", default_value_t = DEFAULT_INDEL_MAX_RATIO, help_heading="INDEL DETECTION", help="Flag an indel breakpoint when min/max total depth between two neighboring positions falls below this ratio (a sharp coverage cliff)")]
    pub indel_max_ratio: f64,

    //min total depth at a position to consider a breakpoint
    #[clap(long="indel-min-depth", default_value_t = DEFAULT_INDEL_MIN_DEPTH, help_heading="INDEL DETECTION", help="Minimum total depth to identify breakpoints for indel detection")]
    pub indel_min_depth: usize,

    //max total depth after a coverage drop to consider a breakpoint
    #[clap(long="indel-max-drop-depth", default_value_t = DEFAULT_INDEL_MAX_DROP_DEPTH, help_heading="INDEL DETECTION", help="Maximum total depth to be considered an indel breakpoint (the low side of a coverage cliff must be under this)")]
    pub indel_max_drop_depth: usize,

    #[clap(long="indel-width", default_value_t = DEFAULT_INDEL_WIDTH, help_heading="INDEL DETECTION", help="Maximum distance between a coverage drop and the following rise to be paired as an indel for reconstruction")]
    pub indel_width: usize,


    //OUTPUT PARAMETERS
    //todo add output locations, output formats
    #[clap(short, long="output", help_heading="OUTPUT", default_value = DEFAULT_OUT_FOLDER, help="Folder to output all resulting files")]
    pub output: String,

    //EXTRA OUTPUT FORMATS
    #[clap(long="pileup", default_value_t = DEFAULT_TSV_PILEUP, help_heading="EXTRA OUTPUT FILES", help="Also output a tsv of the approximate pileup for each sample and reference")]
    pub output_pileup: bool,

    #[clap(long="alignment", default_value_t = DEFAULT_ALIGNMENT, help_heading="EXTRA OUTPUT FILES", help="Output an multifasta containing the alignment of all samples to the reference and themselves (only major variant positions)")]
    pub output_alignment: bool,

    #[clap(long="keep-kmer-info", default_value_t = DEFAULT_KEEP_KMER_INFO, help_heading="EXTRA OUTPUT FILES", help="Keep kmer count information and temporary files (deleted by default)")]
    pub keep_kmer_counts: bool, 

    #[clap(long="consensus", default_value_t = DEFAULT_CONSENSUS, help_heading="EXTRA OUTPUT FILES", help="Return consensus sequences for each virus (.fasta format)")]
    pub output_consensus: bool,

    // OTHER PARAMETERS
    //Number of threads
    #[clap(short, long="threads", default_value_t=4, help="Number of threads")]
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

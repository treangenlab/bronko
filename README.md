# bronko - ultra rapid and precise viral variant calling

## Introduction
bronko is a viral variant caller that can rapidly detect most major and minor variants in any given sequencing dataset. bronko bypasses read mapping and variant calling and directly outputs a VCF from reads and a reference genome.  

bronko also allows users to run multiple samples of the same species and build an alignment that can be inputted directly into phylogenetics software

### Why use bronko over existing viral variant callers?
1. **Simplicity** -- bronko bypasses most steps of typical variant calling (indexing, read mapping, sam manipulation, variant calling, alignment, etc) and packages everything into one command. Given just (sets of) reads and a reference, we will return vcf files, pileups, and can even give the option to return an alignment that can be used for downstream phylogenetics
2. **Consistently ultrafast** -- bronko is around an order of magnitude faster than existing tools for read mapping (bowtie2, ie) and then variant calling (lofreq, ivar), particularly when sequencing depth is super high/the number of reads increases. Additionally, although multithreaded, bronko does not need large amounts of threads to run quickly, allowing lightweight deployment
3. **Comparable accuracy and sensitivity** -- On our benchmarks so far, bronko achieves reasonable consistency with both ivar and lofreq on variant calling, even sometimes outperforming them on recall

### Why not to use bronko
1. If you are interested in identifying indels in viral genomes, currently we do not report on indel presence, but we are working on this problem. 

### How does bronko work
byronko bypasses readmapping by directly mapping kmers with small edit distance to a pileup representing the forward and reverse strands. It then uses the depth information and number of kmers mapping to each position/base to perform variant calling. 

## Some notes before running
Please perform quality control on your samples before running through this tool. In particular, please remove any primer sequences that were used. It is also helpful to have reasonable base quality thresholds (>25 or >30) as bronko does not take into account base quality information. 

## Requirements
The only non-rust requirement is KMC3 (https://github.com/refresh-bio/KMC), so please follow the instructions on their github to download the software on your system (either through conda or by downloading the source code directly)

## How to run
Once fully released I will put everything on conda, but at the moment, you will have to pull down the source code by cloning this repository, then can use this command from the repo to test if things are working. It should download dependencies (aside from KMC3) from cargo directly:

```
cargo run -- query --help
```

If this works, then you should be able to run this command to do a basic run for single SE reads:

```
cargo run -- query -g /path/to/reference.fasta -r /path/to/single_end_reads.fastq -o /path/to/output
```

or a single set of paired end reads: 
```
cargo run -- query -g /path/to/reference.fasta -1 /path/to/reads_R1.fastq -2 /path/to/reads_R1.fastq -o /path/to/output
```

or both at the same time:
```
cargo run -- query -g /path/to/reference.fasta -r /path/to/single_end_reads.fastq -1 /path/to/reads_R1.fastq -2 /path/to/reads_R1.fastq -o /path/to/output
```

You can also run multiple files in the same directory at the same time by using * for any of the reads (keep single end reads in the -r and paired in the -1/-2). If the number of paired end reads do not match then an error will be thrown. 


## Issues
bronko is still in early versions and we welcome any feedback. If you run into any issues, please email Ryan Doughty at rdd4@rice.edu or raise an issue on github

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/treangenlab/bronko/rust.yml) ![GitHub License](https://img.shields.io/github/license/treangenlab/bronko) ![Conda Downloads](https://img.shields.io/conda/d/bioconda/bronko)

# bronko - ultra-rapid mapping free viral variant calling

## Introduction
**bronko** is a viral variant caller that can rapidly detect most major and minor variants in any given sequencing dataset. bronko bypasses read mapping and variant calling and directly outputs a VCF from reads and a reference genome.  

bronko also allows users to run multiple samples of the same species and build an alignment that can be inputted directly into phylogenetics software

### Why use bronko over existing viral variant callers?
1. **Simplicity** -- bronko bypasses most steps of typical variant calling (indexing, read mapping, sam manipulation, variant calling, alignment, etc) and packages everything into 1-2 commands
2. **Consistently ultrafast** -- bronko is between 10-1000x faster than existing pipelines for read mapping (bowtie2, ie) and then variant calling (lofreq, ivar, ie), depending on the sequencing depth / number of reads as well as the number of threads available
3. **Comparable accuracy and sensitivity** -- On our benchmarks so far, bronko achieves reasonable consistency with both ivar and lofreq on variant calling for SNVs, even sometimes outperforming on recall in some cases

### Why not to use bronko
1. If you are interested in identifying novel indels in viral genomes, currently we do not report on indel presence, but we are working on this problem. For historical sequencing data, we attempt to bypass this problem by letting you use multiple reference genomes and then selecting the one with highest identity to each sample.
2. If you are interested in super (<0.5% MAF) frequency mutations. Generally we have found that performance degrades below that threshold, but we are working on addressing this. 

### How does bronko work
byonko bypasses readmapping by directly mapping kmers with small edit distance to a pileup representing the forward and reverse strands. It then uses the depth information and number of kmers mapping to each position/base to perform variant calling, similarly to existing tools. The lack of a formal pileup makes the process much more efficient, and thus we are able to achieve similar results in a fraction of the time.  

## Some notes before running
Please perform quality control on your samples before running through this tool. In particular, please remove any primer sequences that were used. It is also helpful to have reasonable base quality thresholds (>25 or >30) as bronko is unable to take this into account due to the transformation into kmer space. 

## Requirements
The only non-rust requirement is KMC3 (https://github.com/refresh-bio/KMC), so please follow the instructions on their github to download the software on your system (either through conda or by downloading the source code directly). If you are downloading from conda (recommended), you do not need to worry about this. 

## Installation Instructions

#### Conda Installation (Recommended)
The best and recommended way to use bronko is through mamba/conda and the bioconda channel. If you do not have these, please follow this link to install them through miniforge: https://github.com/conda-forge/miniforge 

Once installed, you should be able to download bronko using the following command
```
conda create -n bronko bronko
conda activate bronko 
```

You should be able to test if installed successfully by running
```
bronko --help
```

## Basic Usage
For full details, please see the wiki.

### Building bronko database: 
bronko enables users to build indexes containing many strains of a single species. This can be used to run against a heterogenous dataset and the best reference will be chosen for each sample. If you have a single reference that you know you want to use, that is okay as well. 

Basic usage is as follows:
```
bronko build -g /path/to/reference_genomes/*.fa -o bronko_db
```

This will output a `bronko_db.bkdb` index that you can use for downstream variant calling. The only parameter to change how this function is run is the kmer size k. 

We are working on building a number of prebuilt databases for viruses of interest, so once these are available we will release them here. We also provide additional recommendations on how to build a bronko database to achieve the best results in the wiki. 

### Variant calling 

To perform mapping / reference selection / variant calling, use the `bronko call` command:

If using a prebuilt database, use the -d/--db command
```
bronko call -d /path/to/bronko_db.bkdb -r /path/to/single_end_reads.fastq -o /path/to/output
```

or similarly, on paired end reads: 
```
bronko call -d /path/to/bronko_db.bkdb -1 /path/to/reads/*_R1.fastq -2 /path/to/reads/*_R2.fastq -o /path/to/output
```

or both single end and paired end reads:
```
bronko call -d /path/to/bronko_db.bkdb -r /path/to/single_end_reads.fastq -1 /path/to/reads/*_R1.fastq -2 /path/to/reads/*_R2.fastq -o /path/to/output
```
Note: Keep single end reads in the -r and paired in the -1/-2. If the number of paired end reads do not match then an error will be thrown. 

If you do not want to use a pre-built database, you can replace `-d` with `-g` and provide genomes in fasta format. This will still build a database, but not save it. If you are using a single genome, this may be sufficient.

#### Other Parameters
There are many parameters that can be used to play with the stringency of variant calling. See the wiki here for a description.

## Issues
bronko is still in early versions and we welcome any feedback, suggestions, or feature requests. If you run into any issues, please email Ryan Doughty at rdd4@rice.edu or raise an issue on github

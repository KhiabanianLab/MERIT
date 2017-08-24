
# MERIT

MERIT - *Mutation Error Rate Identification Toolkit* - is a comprehensive pipeline designed for in-depth quantification of erroneous HTS calls developed in [Khiabanian Lab](www.khiabanian-lab.org) by Mohammad Hadigol. 

## Introduction

MERIT uses SAMtools to identify all positions with alternate alleles from the aligned, indexed sequencing reads. It then probes SAMtools mpileup data  to extract the reference and alternate alleles' depths, Phred quality, and position-in-read information for all variants, even when they are present in only a single read amongst tens of thousands. Finally, it obtains the genomic context of the variants from the reference genome, including the nucleotides immediately at their 5' and 3', and estimates error rates for 96 possible single nucleotide substitutions as well as four single-base and 16 double-base insertions/deletions (indels). An optional annotation step is also available. 

## Workflow

The MERIT pipeline is organized into 6 main steps:

- Step 1: Run samtools mpileup to detect variants and generate Pileup and BCF/VCF 
- Step 2: Optional annotation via snpEff and SnpSift 
- Step 3: Convert VCf to human readable file via bcftools query and generate SampleName.var 
- Step 4: Extract Phred scores and more from Pileup and store in SampleName.var.qual 
- Step 5: Adding context to variants and generate SampleName.var.ctx 
- Step 6: Estimate context-specific error rates 

As input, MERIT takes a coordinate-sorted indexed BAM file (See BAM directory for example input BAM) and it produces context-specific error rates as output (See RES directory). 

## Dependencies

The following programs must be in your `PATH`:

- Python 2.7.11
- [Samtools and Bcftools 1.3.1](http://www.htslib.org) 

Optional (required if running annotation step 2):

- Java 1.8
- [SnpEff and SnpSift 4.1](http://snpeff.sourceforge.net) 

The following Python package is required:

- [numpy](http://www.numpy.org)

## Installation

Just clone this repository and run MERIT as described below. 


## Reference file

MERIT needs an indexed reference fasta file (whatever you mapped your BAM files to) -- e.g., hg19.fa and hg19.fa.fai.


## Usage example

Running MERIT on merge reads of example BAM file; UltraII_TP53.2.ME:

```
python MERIT_1.0.py  -N UltraII_TP53.2.ME
```

Try:

```
python MERIT_1.0.py  --help
```

for detailed description of input parameters and their default values. 


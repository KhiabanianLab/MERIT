
# MERIT

MERIT - *Mutation Error Rate Identification Toolkit* - is a comprehensive pipeline designed for in-depth quantification of erroneous HTS calls developed in [Khiabanian Lab](http://khiabanian-lab.org/) by Mohammad Hadigol. 

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

## Description of MERIT's results 

For each sample, MERIT reports the following in /RES/SampleName directory: 

1. SampleName.var.ctx: this tab separated file includes all the mismatches (variants) in all sequencing reads relative to the reference genome. Columns are chromosome (chr), position of the variant (pos), reference sequence (ref), context-specific (nucleotides immediately at 5' and 3') reference sequence (ref_ctx), alternative allele (alt), context-specific (nucleotides immediately at 5' and 3') alternative allele (alt_ctx), total coverage in R1 reads (total_depth_frw), total coverage in R2 reads (total_depth_rev), variant allele coverage in R1 reads (alt_depth_frw), variant allele coverage in R2 reads (alt_depth_rev), alternative allele frequency as freq=100*(alt_depth_frw + alt_depth_rev)/(total_depth_frw + total_depth_rev), mode of the position-in-read of the all the variant alleles at that position in R1 reads (pos_inread_alt_frw), mode of the position-in-read of the all the variant alleles at that position in R2 reads (pos_inread_alt_rev), if indel, type of the inserted/deleted bases (indel_bases), number of repeated indel bases (for homopolymeric inddels > 1) (No. of hp repeats), annotation (effect). 

2. In /RES/SampleName/err_count, the following comma separated files include error counts (context-specific total depth of reference and alternative allele):

  - SNV_err_count.csv: 96 possible context-specific base substitution error counts 
  - DEL_1nt_err_count.csv: 4 single-base deletion error counts 
  - INS_1nt_err_count.csv: 4 single-base insertion error counts 
  - DEL_1nt_hp_err_count.csv: single-base homopolymer deletion error counts 
  - INS_1nt_hp_err_count.csv: single-base homopolymer insertion error counts 
  
if option `--two_nt_indel` is enabled, the following four additional files for double-base indels are also generated:

   - DEL_2nt_err_count.csv: 4 double-base deletion error counts
   - INS_2nt_err_count.csv: 4 double-base insertion error counts
   - DEL_2nt_hp_err_count.csv: double-base homopolymer deletion error counts
   - INS_2nt_hp_err_count.csv: double-base homopolymer insertion error counts

3. And finally in directory /RES/SampleName/err_count, the corresponding context specific error rates estimated using Binomial model are presented. 

We note that the number of homopolymeric repeats considered in error rate estimation is controlled via input option `--max_hp_no`. In estimating these error rates, `--cut_off_freq` and `--min_depth` filters are applied to SampleName.var.ctx when counting and estimating error rates. 

For detailed description of input parameters and their default values try: 

```
python MERIT_1.0.py  --help
  -h, --help            show this help message and exit
  -B BAM_DIR, --BAM_dir BAM_DIR
                        Directory of input BAM file
  -R RES_DIR, --RES_dir RES_DIR
                        Directory of results
  -I POSITIONS, --positions POSITIONS
                        BED file containing a list of regions or sites where
                        pileup or BCF should be generated
  -H REF, --REF REF     Reference genome
  -N NAME, --Name NAME  Sample name
  --max_hp_no MAX_HP_NO
                        Maximum number of homopolymer repeats to be considered
                        in indel error rate analysis
  --cut_off_freq CUT_OFF_FREQ
                        Cut off frequency (percentage) to identify errors
  --min_depth MIN_DEPTH
                        Minimum total depth of variant to be considered in
                        error rate analysis
  --qual-offset QUAL_OFFSET
                        Quality offset
  --ann ANN             Comma-separated list of annotating VCFs with which to
                        provide additional annotation in SnpSift. Absolute
                        path should be provided.
  --steps STEPS         steps to run
  --memory MEMORY       the memory for the (SnpEff) Java virtual machine in
                        gigabytes
  -d MAX_DEPTH, --max_depth MAX_DEPTH
                        Max depth option -d for samtools mpileup
  -q MIN_MQ, --min_MQ MIN_MQ
                        Minimum mapping quality for an alignment to be used
  -Q MIN_BQ, --min_BQ MIN_BQ
                        Minimum base quality for a base to be used
  -e EXT_PROB, --ext_prob EXT_PROB
                        Phred-scaled gap extension sequencing error
                        probability
  -F GAP_FRAC, --gap_frac GAP_FRAC
                        Minimum fraction of gapped reads
  --tandem_qual TANDEM_QUAL
                        Coefficient for modeling homopolymer errors
  -L MAX_IDEPTH, --max_idepth MAX_IDEPTH
                        Skip INDEL calling if the average per-input-file depth
                        is above INT
  -o OPEN_PROB, --open_prob OPEN_PROB
                        Phred-scaled gap open sequencing error probability
  --verbose, -v         Print commands used in each step (default: off)
  --two_nt_indel        Consider two nt indels in error rate analysis
                        (default: off)
  --noclean             do not delete temporary intermediate files (default:
                        off)

```

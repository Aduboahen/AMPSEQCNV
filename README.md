# AMPSEQCNV
Determine Copy Number Variation from AMPSEQ data

## Overview
AMPSEQCNV is a tool designed to analyze amplicon sequencing data for copy number variation (CNV) detection. It performs GC bias correction on amplicon coverage data and calculates CNV by comparing the corrected coverage against a reference sample. 
<!-- The tool is particularly useful for researchers working with targeted amplicon sequencing data, enabling them to identify CNVs in specific genomic regions of interest. -->

## Features
- **GC Bias Correction**: Corrects for GC content bias in amplicon coverage data.
- **Per-Amplicon Coverage Summary**: Generates a summary of coverage for each amplicon.
- **CNV Calculation**: Computes CNVs by comparing corrected coverage data against a reference sample.
- **Visualization**: Optionally plots CNV results for easier interpretation.

## Requirements
- Required Python packages listed in `requirements.txt`
- Sorted BAM files of amplicon sequencing data
- Reference genome in FASTA format (indexed by faidx)
- BED file defining target amplicon regions

## How it works
1. **GC Bias Correction**: The `gc_bias_correction.py` script takes sorted BAM files, a reference genome in FASTA format, and a BED file of target amplicon regions. It computes GC content for each amplicon using a specified window size and applies LOESS smoothing to correct for GC bias.
2. **CNV Calculation**: The `cnv_calculator.py` script reads the corrected coverage data from the TSV file. It requires a reference sample ID for normalisation and calculates CNVs by finding log2 of the ratio of the corrected coverage against the reference sample.

```
   Raw GC values & coverage
   ┌───────────────┐
   │ GC%  Coverage │
   ├──────┬────────┤
   │ 0.35 │   22   │
   │ 0.40 │   28   │
   │ 0.45 │   35   │
   │ ...  │  ...   │
   └──────┴────────┘
           │
           ▼
     LOESS smoothing
           │
           ▼
   Expected coverage curve
   (GC → coverage)
           │
           ▼
  Interpolation
           │
           ▼
 Region grouping ─────────────┐
                              ▼
           ┌─────────────────────────────┐
           │ For each region:            │
           │  corrected = Σ(obs) / Σ(exp)│
           │                             │
           └─────────────────────────────┘
                              │
                              ▼
                 Corrected coverage values
         { region1: 0.98, region2: 1.05, ... }
```

## Installation

```bash
git clone	[AMPSEQ](https://github.com/Aduboahen/AMPSEQCNV.git)

cd AMPSEQCNV

mamba create -n ampseqcnv

mamba install -r requirements.txt
```

## Usage

### GC bias correction and per-amplicon coverage summary

python ./gc_bias_correction.py -h

usage: gc_bias_correction.py [-h] [--bams BAMS] [--fasta FASTA] [--bed BED]
														 [--window WINDOW] [--frac FRAC]
														 [--summary SUMMARY]


options: \
  &emsp; -h --help &ensp; show this help message and exit \
  &emsp; --bams BAMS &ensp; Path to sorted BAM files. \
  &emsp; --fasta FASTA &ensp; Reference genome in FASTA format (indexed by faidx). \
  &emsp; --bed BED &ensp; BED file of target amplicon regions. \
  &emsp; --window WINDOW &ensp; Fragment window size to compute GC (300). \
  &emsp; --frac FRAC &ensp; LOESS smoothing parameter (0.3). \
  &emsp; --summary SUMMARY &ensp; Output TSV path for per-amplicon summary (amplicon_coverage.tsv)


### CNV calculation from corrected coverage data

python cnv_calculator.py -h

usage: cnv_calculator.py [-h] [--data-path DATA_PATH]
                         [--reference-sample REFERENCE_SAMPLE]
                         [--output-path OUTPUT_PATH] [--chromosome CHROMOSOME]
												 [--plot-results]

options: \
  &emsp; -h, --help &ensp; show this help message and exit \
  &emsp; --data-path &ensp; DATA_PATH Path to the TSV file containing the corrected coverage data (default: amplicon_covearge.tsv) \
  &emsp; --reference-sample &ensp; REFERENCE_SAMPLE Sample ID of the reference sample to use for normalization \
  &emsp; --output-path &ensp; OUTPUT_PATH  Path to save the CNV results (default: cnv_results.tsv) \
  &emsp; --plot-results &ensp; Whether to plot the CNV results
  &emsp; --chromosome &ensp; CHROMOSOME Specific chromosome to plot (optional, default: None)
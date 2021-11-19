# LV_calling

A downstream analysis pipeline focusing on large gene modification (large variants) after CRISPR/Cas9 editing. The pipeline takes the output from longread_umi[[1]](#1) pipeline using SMRTseq data and generates a series of files and figures for large gene modification analysis.

Updated date: 2021-11-19

NEXT plan:
- Provide example data
- Package the script to conda

**Table of contents**
- [Schematics](#schematics)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Usage](#usage)
- [References](#references)

## Schematics
![Schematics](Schematics.png)

## Installation 

Note that our pipeline is used for downstream analysis after longread_umi pipeline. If you need to learn how the initial steps work, please refer to https://github.com/SorenKarst/longread_umi

1. Requirements & Dependencies
`python` >= 3.6.0
`minimap2` >= 2.17 [[2]](#2)

2. Install python packages
```bash
conda install numpy scikit-learn pandas scipy
conda install -c conda-forge matplotlib
```
These packages are applied in distribution generation.

3. Download the scripts
```bash
git clone https://github.com/baolab-rice/LV_calling.git
```
or
Download zip file from our github page.

4. The SMRTseq_data_processing .py in the scrips folder is the main script to be used. 

## Quick-start

## Usage
```bash
usage: SMRTseq_data_processing.py [-h] -d DIRECTORY -o OUTPUT [-st {1,2}]
                                  [-os {1,2}] [-r READS] [-m] [-g REFERENCE]
                                  [-ld] [-ld_ps LARGE_DELETION_PARAMETERS]
                                  [-ls] [-ld_c]
                                  [-ld_cps LARGE_DELETION_CLUSTERING_PARAMETERS]

longumi_read downstraem analysis: Version 4.0 this script is for extracting data
from longumi_read pipeline output. You may use 'seqtk seq -a <fastq> > <fasta>
to convert file if raw reads need to be processed.'

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Directory of outout folder with longumi_read pipeline,
                        the output folder should contain a raconx2 subfolder.
  -o OUTPUT, --output OUTPUT
                        Output directory. Output the completly organized file.
  -st {1,2}, --stats {1,2}
                        Output the stats filem (default=2): -st 1: UMI_ID
                        Read_IDs Consensus_read_sequence; -st 2: UMI_ID
                        Read_count Consensus_read_sequence
  -os {1,2}, --output_style {1,2}
                        If also involve raw reads, can also produce a file
                        contaning all reads. -os 1: UMI_ID Read_ID
                        Read_sequence Centroid_ID. -os 2: UMI Consensus_seq
                        Reads_seqs
  -r READS, --reads READS
                        Input PicBio ccs fasta file. (converted with seqtk)
  -m, --mapping         Mapping all filetered read to reference genome using
                        minimap2. For the large deletion analysis option,
                        could ONLY use minimap2.
  -g REFERENCE, --reference REFERENCE
                        Alignment refence genome.
  -ld, --large_deletion
                        Large deletion calling, devide LDs as small INDELs or
                        unmodified, 50bp-200bp, and >200bp.
  -ld_ps LARGE_DELETION_PARAMETERS, --large_deletion_parameters LARGE_DELETION_PARAMETERS
                        Large deletion analysis parameters: FORMAT:
                        cut_site_pos+c+tolerance_bp+t Default: 2778c10t
  -ls, --large_insertion
                        Large insertion (>=50bp) calling.
  -ld_c, --large_deletion_clustering
                        Large deletion clustering based on deletion size and
                        deletion start site.
  -ld_cps LARGE_DELETION_CLUSTERING_PARAMETERS, --large_deletion_clustering_parameters LARGE_DELETION_CLUSTERING_PARAMETERS
                        Large deletion clustering parameters: FORMAT: deletion
                        _size_tolenrance+d+deletion_position_tolerance+l
                        Default: 10d10l
  -a, --all             Process all reads.

```

### Process
[Gathering umibins list...]\
[Reading from umibins...]\
[Polishing umis...]\
[Reading from raw fasta file...]\
[Mapping reads with reads IDs...]\
[Writing into files...]\
[Generating fasta file for consensus seqs...]\
[Alignment using minimap2...]\
[Large deletions calling...]\
[Generating distribution figure for large deletions (>200bp)...]\
[Clustering large deletions (>200bp)...]\
[Generating distribution figure for clustered large deletions (>200bp)...]\
[Large insertions calling...]\
Program finished.

## References
<a id="1">[1]</a> 
SM Karst, RM Ziels, RH Kirkegaard, EA Sørensen, D. McDonald, Q Zhu, R Knight, & M Albertsen. (2020). Enabling high-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. bioRxiv, 6459039. https://github.com/SorenKarst/longread_umi
<a id="1">[2]</a> 
Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191 https://github.com/lh3/minimap2

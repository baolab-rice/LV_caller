## SMRTseq_data_processing.py

###Schematics
![Schematics](Schematics.png)

Version 4.0 Update 20210812
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

blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 database.2bit query.fa output.psl\

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
Program finished.\


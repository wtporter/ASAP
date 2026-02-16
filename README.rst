.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN

Amplicon Sequencing Analysis Pipeline (ASAP)
========================================

OVERVIEW:
------------------------------
ASAP is a start-to-finish Nextflow pipeline that supports analysis of amplicon sequencing reads. The pipeline ingests demultiplexed sequencing reads (Illumina, PacBio, or ONT), performs read-level quality control, reference-based alignment, and produces outputs including quality-control metrics, SNPs, iSNVs, and consensus sequences. These outputs are either viewable through built-in workflows, or raw data can be easily accessed for further analysis and the creation of project-specific outputs. ASAP provides a flexible, open-source alternative usable by individuals with limited bioinformatics experience, while offering advanced customization and data access for experienced users. Previously, ASAP has been used for pathogen and AMR detection and for viral whole-genome assembly, providing a comprehensive platform for amplicon sequencing analysis workflows.

INSTALL:
------------------------------
Not Updated...

USAGE:
========================================
*nextflow run \*/ASAP/nextflow/main.nf*
------------------------------

*Input/Output Options:*
------------------------------
* --read_dir [string, null]: Directory containing input FASTQ files.
* --json [string, null]: Path to a JSON file with sample information.
* --outdir [string] [default: ASAP_Results]: Output directory.
* --file_name [string] [default: ASAP_Output]: Meaningful name for reports.
* --technology [string] [default: illumina]: Technology to guide workflow (QC and alignment). Accepted: illumina, ont, ont.v14, pacbio.

*Quality Control & Alignment:*
------------------------------
* --adapter_fasta [string, null] [default: ../asap/illumina_adapters_all.fasta]: FASTA file containing adapter sequences for trimming.
* --fastp_extra_args [string]: Extra args for fastp or fastplong.
* --aligner [string] [default: bowtie2]: Read alignment tool. Accepted: bowtie2, bwa, minimap2.
* --aligner_extra_args [string]: Extra args for chosen aligner.

*Primer Masking & Trimming:*
------------------------------
Note: Primer files must have 6 columns (reference, start, end, name, score, strand) and no headers.

* --primer_file [string, null]: Path to primer file for removal via masking or iVar trimming.
* --mask_primers [boolean]: Mask primers using custom ASAP script. Default: true if bed file is specified.
* --wiggle [integer] [default: 9]: (Masking Only) Wiggle distance for primer masking.
* --mask_bam [boolean] [default: true]: (Masking Only) Whether to mask primers in the BAM file.
* --primer_only [boolean] [default: false]: (Masking Only) Only perform primer masking and exit.

*Variant Calling & Consensus:*
------------------------------
* --proportion [number] [default: 0.1]: Min proportion required to call a variant/iSNV (ivar -t).
* --consensus_proportion [number] [default: 0.8]: Min proportion to call a base, else 'N' (ivar -c).
* --depth [integer] [default: 100]: Minimum read depth required to consider a position covered.
* --min_base_qual [integer] [default: 5]: Min Phred quality score for variant or consensus calling.

*BAM Processor & JSON Options:*
------------------------------
* --asap_snps [boolean] [default: true]: Enable/Disable ASAP SNP and Consensus profiling.
* --breadth [number] [default: 0.8]: Min breadth of coverage to consider an amplicon present.
* --mutation_depth [integer] [default: 5]: Min reads required to call a mutation.
* --whole_genome [boolean] [default: false]: If true, optimized for WGS (omits consensus/depth arrays in JSON).
* --out_file [string]: Name of the output HTML file.

*ASAP Tools (Downstream R-Processing):*
------------------------------
* --asaptools_processing [boolean] [default: true]: Enable downstream R-script processing.  
* --asaptools_cov_table [boolean] [default: true]: Enable Coverage Table generation.  
* --asaptools_snp_table [boolean] [default: true]: Enable SNP Table generation.  
* --asaptools_genbank_location [string]: Path to the reference GenBank file for SNP annotation.  
* --asaptools_max_sample_snp_count [integer] [default: 50]: Filter samples with >X SNPs in the table.  

*Help Options:*  
------------------------------
* --help [boolean, string]: Show the help message for all top level parameters. Pass a specific parameter (e.g., --help technology) for full details.
* --helpFull [boolean]: Show the help message for all non-hidden parameters.
* --showHidden [boolean]: Show all hidden parameters in the help message (use with --help or --helpFull).


LICENSE:
--------

Copyright |copy| The Translational Genomics Research Institute See the
included "LICENSE" document.

CONTACT:
--------
Tanner Porter (tporter@tgen.org)
| TGen North
| 3051 W Shamrell Blvd Ste 106
| Flagstaff, AZ 86001-9435

Darrin Lemmer (dlemmer@tgen.org)
| TGen North
| 3051 W Shamrell Blvd Ste 106
| Flagstaff, AZ 86001-9435

REFERENCES:
-----------
Not Updated.....

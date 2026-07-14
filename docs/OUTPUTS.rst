..
   Outputs, Data & Performance — detailed reference for ASAP.
   For an overview and quick start, see the main README: ../README.rst

This page documents the full ASAP output directory layout, the ``.Rdata`` frames
provided for custom analysis, and performance/tuning considerations. For an overview
and quick start, see the `main README <../README.rst>`_.

.. contents:: On this page
   :local:
   :depth: 1

Output Structure
================

.. code-block:: none

   <outdir>/
   ├── pipeline_info/
   │   ├── sample_read_type_summary.tsv        # Sample IDs and read type (PE / SE)
   │   ├── dag_<name>_<timestamp>.png          # Pipeline DAG
   │   ├── report_<name>_<timestamp>.html      # Per-process execution report
   │   ├── trace_<name>_<timestamp>.txt        # Per-task resource usage
   │   └── timeline_<name>_<timestamp>.html    # Interactive job timeline
   │
   ├── reference/
   │   ├── reference.fasta                     # Extracted amplicon reference sequences
   │   └── bwa_index/ or bt2_index/            # Aligner index files
   │
   ├── sample_info/<sample>/
   │   ├── fastqc_initial/                     # Pre-trim FastQC reports
   │   ├── fastp/ or fastp_long/               # Trimmed reads + QC HTML / JSON
   │   ├── fastqc_post/                        # Post-trim FastQC reports
   │   ├── bwa/ or bowtie2/ or minimap2/       # Sorted, indexed BAM + flagstat
   │   ├── xml/
   │   │   └── <sample>.xml                    # Per-sample ASAP results
   │   ├── rdata/
   │   │   ├── <sample>_XML_Data.Rdata         # Per-sample R data object
   │   │   └── <sample>_Summary.csv            # Per-sample summary table
   │   ├── mask_primers/                       # Primer-masked BAM + masking log/TSV
   │   ├── identity_filter/                    # Identity-filtered BAM + log
   │   ├── smor/ or smor_correction/           # SMOR-processed BAM + log
   │   ├── ivar_trim/                          # iVAR-trimmed BAM
   │   ├── ivar_variants/                      # iVAR variant TSV
   │   └── ivar_consensus/                     # iVAR consensus FASTA
   │
   └── sample_reports/
       ├── <name>_analysis.xml                 # Combined XML (all samples)
       ├── <name>_report.html                  # Combined HTML report
       ├── multiqc/
       │   └── <name>_multiqc.html             # Aggregated QC report
       ├── rdata/
       │   ├── <name>_ASAP_Data.Rdata          # Combined R data (all samples)
       │   └── SNP_Amino_Acid_Table.Rdata      # Amino acid changes (GenBank required)
       ├── general_reports/
       │   ├── <name>_Summary.csv              # Combined summary table
       │   └── <name>_coverage_table.xlsx      # Coverage depth table
       ├── plots/
       │   ├── *.html                          # QC figures (interactive)
       │   └── *.jpg                           # QC figures (static)
       ├── snp_reports/
       │   ├── *.csv                           # SNP / iSNV table + amino acid changes
       │   └── *.xlsx                          # SNP table Excel (if --asaptools_snp_table_xls)
       └── fasta/
           └── *.fasta                         # Consensus FASTA exports


Custom Analysis with Rdata Frames
=================================

Every table and figure ASAP Tools produces is built from a small set of tidy data frames
saved as ``.Rdata`` objects under ``sample_reports/rdata/``. These are provided so you can
load the same data into an interactive R session and build **custom plots or analyses**
without re-parsing XML. Load an object with ``load()`` (it restores the named frames into
your environment):

.. code-block:: r

   load("sample_reports/rdata/MyRun_ASAP_Data.Rdata")   # -> final_asap, final_snps, final_array
   load("sample_reports/rdata/SNP_Amino_Acid_Table.Rdata")  # -> Amino_Acids, Gene_SNPS (GenBank runs)

   library(tidyverse)
   # Example: depth-of-coverage trace for one sample/amplicon
   final_array %>%
     filter(name == "Sample01", assay_name == "gyrA") %>%
     ggplot(aes(position, depth)) + geom_line()

``<name>_ASAP_Data.Rdata`` — combined across all samples
--------------------------------------------------------

**final_asap** — one row per sample × assay (amplicon); the master summary frame. Key
columns: ``run``, ``name`` (sample ID), ``name_short`` (shortened plot label),
``assay_name``, ``assay_gene``, ``assay_type``, ``amplicon_number``; read-fate counts
(``total_reads``, ``trimmed_reads``, ``mapped_reads``, ``unassigned_reads``,
``unmapped_reads``, ``amplicon_reads``, ``aligned_reads``, ``primer_reads``,
``no_primer_reads``, ``identity_input``, ``identity_discarded``, ``smor_*``); coverage
(``breadth``, ``avg_depth``); and the ``consensus_seq`` plus comma-delimited per-position
arrays (``depths``, ``proportions``, ``n_reads``, ``quality_discards``). The per-position
arrays are pre-expanded for you in ``final_array`` — prefer that for plotting.

**final_snps** — one row per detected SNP / iSNV per sample × assay. Key columns:
``run``, ``name``, ``name_short``, ``assay_name``, ``amplicon_number``, ``snp_name``,
``snp_position``, ``snp_reference``, ``snp_depth``, ``snp_proportion`` (allele
frequency, %), ``location_depth`` (total depth at the position), and the codon-merge /
linked-SNP annotation columns (``codon_merge_*``, ``linked_snp_*``) when
``--codon_correction`` / ``--discover_roi`` are enabled.

**final_array** — long/tidy per-position table (one row per sample × assay × position):
``run``, ``name``, ``name_short``, ``assay_name``, ``position``, ``depth``,
``proportion``, ``n_reads``, ``quality_discards``. This is the frame behind the coverage,
percent-N, and depth plots and is the easiest starting point for custom per-position
figures.

``SNP_Amino_Acid_Table.Rdata`` — GenBank runs only
--------------------------------------------------

Written only when a GenBank reference (or ``--asaptools_genbank_location``) supplies CDS
annotations.

**Amino_Acids** — one row per coding SNP: ``assay_name``, ``SNP`` (genome-level variant,
e.g. ``A1401G``), ``Product`` (gene/product name), and ``AA`` (amino acid change).

**Gene_SNPS** — the genome-to-gene coordinate mapping for each SNP (gene-relative
position), joined to ``assay_name`` — useful for annotating variants against gene
coordinates rather than genome coordinates.

.. note::

   The per-sample ``sample_info/<sample>/rdata/<sample>_XML_Data.Rdata`` files hold the
   same three frames before cohort merging, under the names ``ASAP``, ``SNPS``, and
   ``array_info`` — handy for inspecting a single sample.


Performance Considerations
==========================

**SLURM resource scaling:** All processes use exponential retry. Memory and wall-time
double on each retry attempt (``memory = { N.GB * 2**(task.attempt-1) }``), preventing
transient resource limits from aborting runs.

**R post-processing memory:** ``PROCESS_COMBINE_RDATA`` (default 100 GB) and
``PROCESS_SNPS_TO_AMINOACIDS`` (default 30 GB) are the most memory-intensive steps,
as they load all per-sample Rdata objects simultaneously. For large runs (100+ samples),
ensure sufficient memory is available on the target SLURM partition.

**Shared environment cache (``--env_dir``):** By default Conda/Singularity environments
are cached under Nextflow's working directory. Set ``--env_dir <path>`` to redirect the
caches to ``<path>/conda`` and ``<path>/singularity`` — point multiple runs (or all users
on a cluster) at one shared location to build each environment once and reuse it, rather
than rebuilding per run. Unset leaves Nextflow's default behavior unchanged.

**SMOR and identity filtering:** These steps add computational overhead but substantially
improve variant call quality. For routine screening, they can be omitted; for
high-resolution iSNV detection or resistance profiling, they are strongly recommended.

**Primer masking:** Runs on a single CPU and completes quickly relative to alignment.
The BAM re-sort step (``--mask_bam true``) adds a small overhead but ensures downstream
tools receive a correctly sorted BAM.

**Read-based phasing (``--discover_roi``):** This is the most compute-intensive option
inside ASAP BAM processing. It first makes an extra pass over the BAM to build a
per-fragment allele table, then runs a **pairwise (O(n²)) co-occurrence search** over
every qualifying SNP pair per amplicon. Cost therefore grows with read depth,
quadratically with the number of candidate SNPs, and with **read/fragment length**:
only SNP pairs a single fragment can physically span are evaluated, so short Illumina
reads skip most distant pairs, whereas long reads (ONT / PacBio) and large paired-end
inserts bring far more pairs within reach and make the search substantially denser.
Noisy or highly diverse amplicons (and low variant-frequency thresholds, which admit
many candidates) are the expensive cases. The dominant lever is
``--discover_roi_min_snp_perc``: it pre-filters the SNP set before the O(n²) search, so
raising it sharply reduces runtime on noisy data. By contrast,
codon-aware linkage (``--codon_correction``) is bounded to same-codon pairs and is
comparatively cheap. Both are off by default and add no cost unless enabled.


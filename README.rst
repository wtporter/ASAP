.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN

.. code-block:: none

    █████╗  ███████╗  █████╗  ██████╗
   ██╔══██╗ ██╔════╝ ██╔══██╗ ██╔══██╗
   ███████║ ███████╗ ███████║ ██████╔╝
   ██╔══██║ ╚════██║ ██╔══██║ ██╔═══╝
   ██║  ██║ ███████║ ██║  ██║ ██║
   ╚═╝  ╚═╝ ╚══════╝ ╚═╝  ╚═╝ ╚═╝
   ──────────────────────────────────────────────────────────────────────────────────────────
   Amplicon Sequencing Analysis Pipeline (ASAP)
   ──────────────────────────────────────────────────────────────────────────────────────────

**Pre-release** | Authors: Darrin Lemmer, W. Tanner Porter, *et al.*,
Pathogen & Microbiome Division, Translational Genomics Research Institute (TGen) |
License: |copy| TGen North (non-commercial)

**ASAP** is a start-to-finish Nextflow pipeline for targeted amplicon sequencing analysis.
The pipeline ingests demultiplexed sequencing reads from multiple platforms (Illumina,
Oxford Nanopore, PacBio), performs read-level quality control, reference-based alignment,
and produces quality-control metrics, SNPs, iSNVs, and consensus sequences. Outputs are
viewable through built-in reporting workflows or accessible as structured data for
project-specific downstream analysis. ASAP provides a flexible, open-source platform
suitable for users with limited bioinformatics experience, while offering advanced
customization for expert users. Previously, ASAP has been applied to pathogen and
antimicrobial resistance (AMR) detection and viral whole-genome assembly.

**Intended Uses:**

1. Pathogen detection and AMR profiling from targeted amplicon panels
2. Viral whole-genome assembly from tiled amplicon sequencing (e.g., SARS-CoV-2, RSV)
3. Multi-reference analysis for species differentiation or cross-panel quality control
4. High-resolution iSNV detection in mixed infections or heteroresistant populations

----

.. contents:: Table of Contents
   :depth: 2
   :backlinks: none

----

Quick Start
===========

**1. Install** — clone the repo and create the Nextflow environment (one time):

.. code-block:: bash

   git clone https://github.com/wtporter/ASAP.git
   cd ASAP
   conda env create -f ASAP_nextflow_env.yml   # Nextflow 25.10.4 + nf-test 0.9.5
   conda activate ASAP_nextflow_env
   cd nextflow
   nextflow -version                            # verify the environment

Singularity containers and all per-step Conda environments are resolved automatically on
first run — no further setup required. The **nf-schema** plugin is downloaded on first run.
See Requirements_ for supported platforms and pinned tool versions, and
`docs/INSTALLATION.rst <docs/INSTALLATION.rst>`_ for detailed setup (HPC/SLURM, shared
environment caches, and troubleshooting).

**2. Run** — a minimal Illumina + GenBank analysis on SLURM:

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./reads/ \
     --reference_input "./refs/*.gb" \
     --outdir          ASAP_Results \
     --file_name       MyRun \
     -profile slurm

**3. Explore parameters** — print the full, schema-driven parameter list:

.. code-block:: bash

   nextflow run main.nf --help

More runnable recipes are in Examples_. A step-by-step walkthrough of every stage and
every parameter is in `docs/PIPELINE_STEPS.rst <docs/PIPELINE_STEPS.rst>`_.

----

Overview
========

Targeted sequencing has become an essential tool across clinical research, ecology, and
infectious disease surveillance. By amplifying targeted genomic regions, this approach
provides high coverage depth and increased sensitivity in complex samples at a fraction of
the cost of non-targeted approaches. Target-specific primers and probes can be designed to
amplify highly specific genomic regions, enabling high taxonomic resolution, resistance
marker detection, virulence factor identification, or tiled coverage across a gene or
full genome.

ASAP is designed for **target-specific amplicon analysis** — aligning reads against one or
more references and using the resulting alignment files to synthesize final datasets
including amplicon presence/absence, coverage statistics, consensus sequences,
single-nucleotide polymorphism (SNP) tables, and intra-host single-nucleotide variants
(iSNVs).

**Highlights:**

- Supports Illumina short reads (paired-end and single-end), Oxford Nanopore, and PacBio
- Multi-reference analysis: multiple amplicons or species analyzed in a single run
- Primer masking, percent-identity filtering, and SMOR error correction for high-resolution variant calling
- Supports SNP analysis including automatic SNP to amino acid conversions, presence-true absence detection, codon-aware SNP linkage, and read-based phasing
- Various automatic plots and reports for quick reference and visualization of data.
- Stored Rdata files for custom data visualizations.
- Optional iVAR integration for primer trimming, variant calling, and consensus generation
- HPC-ready via SLURM with exponential retry and configurable resource scaling


----

Pipeline Summary
================

.. code-block:: none

   Reference Input (FASTA / GenBank / Excel / JSON)
             │
             ▼
   PREPARE_ASAP_JSON ──► GENERATE_REFERENCE_FASTA ──────────────┐
                                                                │
   Read Files (FASTQ) ──────────────────────────────────────────┤
             │                                                  │
             ▼                                                  │
   [Optional] FastQC (initial QC)                               │
             │                                                  │
             ▼                                                  │
   ┌─────────────────────────┐                                  │
   │ Illumina:  fastp        │                                  │
   │ ONT/PacBio: fastplong   │                                  │
   └─────────────────────────┘                                  │
             │                                                  │
   [Optional] FastQC (post-trim QC)                             │
             │                                                  │
             ▼                                                  ▼
   ┌──────────────────────────────────────────────────────────────┐
   │                      Alignment                               │
   │   Illumina:     Bowtie2 (default) or BWA-MEM                 │
   │   ONT/PacBio:   minimap2 (auto-selected)                     │
   └──────────────────────────────────────────────────────────────┘
             │
             ▼
   [Optional] Primer Masking        ← maskPrimers.py + BED (or primer CSV → GENERATE_PRIMER_BED)
             │
   [Optional] Percent-Identity Filtering  ← identityFilter.py
             │
   [Optional] SMOR Masking / SMOR Correction ← generateSMORbam.py
             │
             ├─────────────────────────────────────────────────────┐
             ▼                                                     ▼
   ┌──────────────────────────┐                    ┌───────────────────────────┐
   │  ASAP BAM Processing     │                    │  iVAR [optional]          │
   │  (ASAPBamProcessor.py)    │                    │  ─ Primer trimming        │
   │  → per-sample XML        │                    │  ─ Variant calling (.tsv) │
   └──────────────────────────┘                    │  ─ Consensus FASTA        │
             │                                     └───────────────────────────┘
             ▼
   [Optional] ASAP Tools R Post-Processing
     ├── Per-sample XML → Rdata  (parallel)
     ├── Combine all Rdata       (gather)
     ├── Coverage Table          (Excel)
     ├── QC Figures              (HTML + JPG)
     ├── Consensus FASTA export
     ├── SNP → Amino Acid translation  (if GenBank provided)
     └── SNP / iSNV Table        (CSV ± Excel)
             │
   [Optional] OUTPUT_COMBINER + FORMAT_OUTPUT
             │  (combined XML → HTML report)
             ▼
   [Optional] MultiQC (aggregated QC report)   ← skipped by default


----

Reference Input Formats
=======================

ASAP is designed to analyze multiple amplicon targets simultaneously in a single run.
This supports:

- **Species co-analysis** — e.g., RSV-A and RSV-B from the same sequencing run
- **Multi-gene AMR panels** — multiple resistance genes across *M. tuberculosis*
- **Cross-contamination QC** — include neighboring-species references to detect and quantify cross-contamination between samples

Multi-reference analysis can be enabled by:

- Providing multiple sequences in a multi-FASTA or JSON reference file
- Specifying multiple GenBank files (``"./refs/*.gb"``)
- Using the ASAP Excel template with multiple amplicon rows
- Subsetting output annotation to specific positions via ``--asaptools_positions_of_interest``

ASAP accepts four reference formats via ``--reference_input``:

+---------+----------------------------------------------------------------+----------------+-----------------+-------------------------------------------+
| Format  | Extension(s)                                                   | Files accepted | AA annotation   | Notes                                     |
+=========+================================================================+================+=================+===========================================+
| GenBank | ``.gb``, ``.gbk``, ``.gbf``, ``.gbb``, ``.gbff``, ``.genbank`` | 1 or more      | Yes (CDS-based) | **Recommended** — full features, standard |
+---------+----------------------------------------------------------------+----------------+-----------------+-------------------------------------------+
| JSON    | ``.json``                                                      | 1              | Yes (full)      | Native format; custom significance rules  |
+---------+----------------------------------------------------------------+----------------+-----------------+-------------------------------------------+
| FASTA   | ``.fasta``, ``.fa``                                            | 1              | No              | Presence/absence only; quick runs         |
+---------+----------------------------------------------------------------+----------------+-----------------+-------------------------------------------+
| Excel   | ``.xlsx``, ``.xls``                                            | 1              | No              | Spreadsheet-based panel entry             |
+---------+----------------------------------------------------------------+----------------+-----------------+-------------------------------------------+

All non-JSON formats are converted to an internal JSON assay description by
``prepareJSONInput_nextflow.py`` before processing.

----

Common Tasks
============

The most-used ASAP features and the parameters that switch them on. Each is optional and
off (or at its default) unless set, and they combine freely — see Examples_ for full
commands and `docs/PIPELINE_STEPS.rst <docs/PIPELINE_STEPS.rst>`_ for what each does and
every tuning parameter.

.. list-table::
   :header-rows: 1
   :widths: 44 56

   * - Goal
     - Key parameter(s)
   * - Remove primer bases from aligned reads
     - ``--primer_file <bed|csv>`` ``--mask_primers``
   * - Keep only primer-overlapping reads
     - ``--primer_only``
   * - Exclude off-target / near-neighbor reads
     - ``--identity 0.97``
   * - Resolve low-frequency variants (Illumina overlap)
     - ``--smor_correction`` (or ``--smor`` for full-overlap assays)
   * - Tune SNP / iSNV calling
     - ``--proportion`` ``--mutation_depth`` ``--depth``
   * - Call variants and consensus with iVar
     - ``--ivar``
   * - Annotate amino-acid changes
     - use a GenBank reference (or ``--asaptools_genbank_location``)
   * - Annotate positions of interest
     - ``--asaptools_positions_of_interest <csv>``
   * - Link same-codon SNPs / phase SNPs onto reads
     - ``--codon_correction`` / ``--discover_roi``
   * - Combine per-sample results into one HTML report
     - ``--combine_output``
   * - Analyze long reads (ONT / PacBio)
     - ``--technology ont`` (or ``pacbio``)
   * - Skip R post-processing (XML + HTML only)
     - ``--asaptools_processing false``

----

Examples
========

These commands assume the ``ASAP_nextflow_env`` environment is active and you are in the
``nextflow`` directory (see `Quick Start`_ for setup).

**Illumina paired-end reads with a GenBank reference (SLURM) — minimal:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./reads/ \
     --reference_input "./refs/*.gb" \
     --outdir          ASAP_Results \
     --file_name       MyRun \
     -profile slurm

**JSON reference passed directly (no conversion):**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./reads/ \
     --reference_input assay.json \
     --outdir          Results \
     --file_name       MyRun \
     -profile slurm

**TB amplicon panel — primer masking, identity filtering, SMOR correction, full asaptools:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir          ./reads/ \
     --reference_input   ./refs/H37Rv.gb \
     --outdir            TB_Results \
     --file_name         TB_Run \
     --aligner           bwa \
     --primer_file       ./primers/TB_primers.bed \
     --mask_primers      true \
     --identity          0.97 \
     --smor_correction   true \
     --fastp_extra_args  "-l 100" \
     --proportion        0.01 \
     --asaptools_positions_of_interest ./genes/H37Rv_genes.csv \
     -profile slurm

**SARS-CoV-2 tiled amplicon — iVAR variant calling and consensus:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./sc2_reads/ \
     --reference_input ./refs/SC2_Reference.json \
     --outdir          SC2_Results \
     --file_name       SC2_Run \
     --aligner         bwa \
     --primer_file     ./primers/SC2_primers.bed \
     --mask_primers    true \
     --ivar            true \
     --combine_output  true \
     -profile slurm

**RSV multi-GenBank reference (RSV-A + RSV-B simultaneously):**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./rsv_reads/ \
     --reference_input "./refs/RSV*.gb" \
     --outdir          RSV_Results \
     --file_name       RSV_Run \
     --aligner         bwa \
     --combine_output  true \
     -profile slurm

**ONT long reads with minimap2:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir        ./ont_reads/ \
     --reference_input ./refs/reference.gb \
     --outdir          ONT_Results \
     --file_name       ONT_Run \
     --technology      ont \
     --combine_output  true \
     -profile slurm

**ASAP Tools disabled — XML and HTML report only:**

.. code-block:: bash

   nextflow run main.nf \
     --read_dir              ./reads/ \
     --reference_input       ./refs/reference.json \
     --outdir                Results \
     --file_name             Run \
     --asaptools_processing  false \
     --combine_output        true \
     -profile slurm

**Resume a failed or interrupted run:**

.. code-block:: bash

   nextflow run main.nf [params] -resume

----

Outputs at a Glance
===================

A completed run writes everything under ``<outdir>/``:

.. code-block:: none

   <outdir>/
   ├── pipeline_info/     # run reports: DAG, execution/timeline/trace, sample read-type summary
   ├── reference/         # extracted amplicon reference FASTA + aligner index
   ├── sample_info/       # per sample: trimmed reads, QC, BAMs, per-sample XML + Rdata
   └── sample_reports/    # cohort level: combined XML/HTML report, tables, plots, FASTA, SNP reports

See `docs/OUTPUTS.rst <docs/OUTPUTS.rst>`_ for the complete directory tree, the ``.Rdata``
frame schemas for custom analysis, and performance/tuning notes.

----

Requirements
============

- **Nextflow** 25.10.4 — provided via ``ASAP_nextflow_env.yml``
- **nf-schema** plugin 2.5.1 — loaded automatically via ``nextflow.config`` on first run
- **nf-test** 0.9.5 — provided via ``ASAP_nextflow_env.yml``; required only to run the test suite
- **Singularity / Apptainer** (for containerized alignment and QC tools)
- **Conda / Mamba** (environments are built automatically from module YMLs — no manual setup required)
- A reference file in FASTA, GenBank, Excel (.xlsx), or JSON format

**Execution profiles:**

+---------------+---------------------+---------------------+-----------------------------------+
| Profile       | Executor            | Containers          | Best For                          |
+===============+=====================+=====================+===================================+
| ``slurm``     | SLURM (child jobs)  | Singularity + Conda | Production HPC runs               |
+---------------+---------------------+---------------------+-----------------------------------+
| ``conda``     | Local (current node)| Conda               | Interactive ``srun`` or laptop    |
+---------------+---------------------+---------------------+-----------------------------------+

**Key tool versions:**

+-----------+---------+---------------------------------------------+
| Tool      | Version | Purpose                                     |
+===========+=========+=============================================+
| fastp     | 0.23.4  | Illumina read QC and adapter trimming       |
+-----------+---------+---------------------------------------------+
| fastplong | 0.4.1   | ONT / PacBio read QC                        |
+-----------+---------+---------------------------------------------+
| Bowtie2   | 2.x     | Short-read alignment                        |
+-----------+---------+---------------------------------------------+
| BWA-MEM   | 0.7.x   | Short-read alignment (alternative)          |
+-----------+---------+---------------------------------------------+
| minimap2  | 2.x     | Long-read alignment                         |
+-----------+---------+---------------------------------------------+
| SAMtools  | 1.x     | BAM manipulation and indexing               |
+-----------+---------+---------------------------------------------+
| iVAR      | 1.4.4   | Primer trimming, variant calling, consensus |
+-----------+---------+---------------------------------------------+
| FastQC    | 0.12.1  | Per-sample read QC                          |
+-----------+---------+---------------------------------------------+
| MultiQC   | latest  | Aggregated QC reporting                     |
+-----------+---------+---------------------------------------------+

----

Documentation
=============

- `docs/INSTALLATION.rst <docs/INSTALLATION.rst>`_ — full setup guide: prerequisites,
  execution profiles, shared environment caches, and troubleshooting.
- `docs/PIPELINE_STEPS.rst <docs/PIPELINE_STEPS.rst>`_ — detailed per-step reference and
  the full parameter tables for every stage.
- `docs/OUTPUTS.rst <docs/OUTPUTS.rst>`_ — complete output layout, ``.Rdata`` frames for
  custom R analysis, and performance considerations.
- `docs/TESTING.rst <docs/TESTING.rst>`_ — running the nf-test end-to-end and R unit-test
  suites, test tags, and per-test descriptions.
- `DEVELOPER.md <DEVELOPER.md>`_ — architecture, Python/R module internals, the output
  schema, and how to extend the pipeline.
- ``nextflow run main.nf --help`` — the complete, schema-driven parameter list.

----

License
=======

Copyright |copy| The Translational Genomics Research Institute (TGen).
See the included ``LICENSE`` document.

Available for academic and research use under a license from TGen that is free for
non-commercial use. Distributed on an "AS IS" basis without warranties or conditions
of any kind, either express or implied.

----

Contact
=======

| TGen North
| 3051 W Shamrell Blvd Ste 106
| Flagstaff, AZ 86001-9435

| Darrin Lemmer — dlemmer@tgen.org
| W. Tanner Porter — tporter@tgen.org

Issues and feature requests:
https://github.com/TGenNorth/ASAP/issues

----

References
==========

Key tools ASAP builds on (in-text citations appear in
`docs/PIPELINE_STEPS.rst <docs/PIPELINE_STEPS.rst>`_):

1. Grubaugh ND, Gangavarapu K, Quick J, *et al.* An amplicon-based sequencing framework
   for accurately measuring intrahost virus diversity using PrimalSeq and iVar.
   *Genome Biology*. 2019;20(1):8.
2. Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. *Nature Methods*.
   2012;9:357–359.
3. Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler Aligner.
   *Bioinformatics*. 2009;25(14):1754–1760.
4. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*.
   2018;34(18):3094–3100.
5. Di Tommaso P, Chatzou M, Floden EW, *et al.* Nextflow enables reproducible
   computational workflows. *Nature Biotechnology*. 2017;35(4):316–319.
6. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor.
   *Bioinformatics*. 2018;34(17):i884–i890.
7. Li H, Handsaker B, Wysoker A, *et al.* The Sequence Alignment/Map Format and SAMtools.
   *Bioinformatics*. 2009;25(16):2078–2079.
8. Andrews S. FastQC: A quality control tool for high throughput sequence data. 2010.
   http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
9. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for
   multiple tools and samples in a single report. *Bioinformatics*. 2016;32(19):3047–3048.
10. fastplong and the SMOR method — citations forthcoming.

----

AI Development Assistance
==========================

Portions of this pipeline — including module logic, subworkflow design, parameter
schema, test suite, and documentation — were developed with assistance from
**Claude Sonnet 4.6** (Anthropic). AI assistance was used as a coding and design
collaborator under active human direction and review. All scientific decisions,
parameter choices, and pipeline architecture reflect the work of the authors.

Pipeline Steps — Detailed Reference
===================================

Full per-step reference for the ASAP pipeline: what each step does, when to use the
optional steps, and every parameter each step exposes. For a high-level overview and
quick start, see the `main README <../README.rst>`_. For the complete, schema-driven
parameter list, run ``nextflow run main.nf --help``.

.. contents:: On this page
   :local:
   :depth: 1

Step 1 — Quality Control
-------------------------

**Illumina reads** are processed by ``fastp`` [6]_, which performs adapter trimming,
quality filtering, and per-sample HTML/JSON reports.

**ONT and PacBio reads** are processed by ``fastplong`` [CITATION]_, a long-read
variant of fastp.

``FastQC`` [8]_ can run on reads before and after trimming for per-sample QC
assessment, but is **skipped by default** (``--skip_fastqc true``). Pass
``--skip_fastqc false`` to enable both the initial and post-trim FastQC runs
(see `Step 9`_).

+------------------------+--------------+----------------------------------------------------------+
| Parameter              | Default      | Description                                              |
+========================+==============+==========================================================+
| ``--technology``       | ``illumina`` | Platform: ``illumina``, ``ont``, ``pacbio``              |
+------------------------+--------------+----------------------------------------------------------+
| ``--adapter_fasta``    | bundled      | Adapter FASTA; bundled Illumina adapters used by default |
+------------------------+--------------+----------------------------------------------------------+
| ``--fastp_extra_args`` | ``""``       | Additional fastp flags (e.g. ``-l 100`` for min length)  |
+------------------------+--------------+----------------------------------------------------------+
| ``--skip_fastqc``      | ``true``     | Skip both FastQC runs; set ``false`` to enable FastQC    |
+------------------------+--------------+----------------------------------------------------------+

Step 2 — Alignment
-------------------

Trimmed reads are aligned to ``reference.fasta`` (extracted from the assay JSON):

- **Illumina:** ``bowtie2`` (default) or ``bwa mem``
- **ONT / PacBio:** ``minimap2`` (selected automatically by ``--technology``)

All aligners produce a coordinate-sorted, indexed BAM published to
``sample_info/<sample>/bwa/``, ``bowtie2/``, or ``minimap2/`` respectively.
BWA and Bowtie2 also emit ``flagstat`` files for MultiQC.

+-----------------------------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Parameter                         | Default     | Description                                                                                                                                                                                                                                                                                           |
+===================================+=============+=======================================================================================================================================================================================================================================================================================================+
| ``--aligner``                     | ``bowtie2`` | Aligner: ``bowtie2``, ``bwa``,``minimap2``                                                                                                                                                                                                                                                            |
+-----------------------------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``--aligner_extra_args``          | ``""``      | Additional arguments passed to the aligner                                                                                                                                                                                                                                                            |
+-----------------------------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``--filter_secondary_alignments`` | ``true``    | Drop secondary/supplementary alignment records (``samtools view -F 0x900``) immediately after alignment, so downstream read counts (``mapped_reads``, ``amplicon_reads``, ``aligned_reads``, depth/breadth) stay consistent. Set to ``false`` to retain all alignment records emitted by the aligner. |
+-----------------------------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Step 3 — Primer Masking *(optional)*
--------------------------------------

*When to use:* amplicon assays where primers are included in the sequenced region.
Masking prevents primer-derived base calls from inflating or distorting SNP frequencies
or consensus sequences.

Primer-derived base calls are masked in aligned reads. A BED-format TSV file specifying
primer coordinates is required — but ``--primer_file`` also accepts a **primer CSV**
(see `Generating a primer BED from a CSV`_), in which case ASAP locates the primers
against the reference and builds the BED for you.

For each read, if the alignment start (R1) or end (R2) falls within ``--wiggle`` bases
of a primer boundary, the primer region is masked: base quality scores are set to 0
and (by default) bases are replaced with ``N``. A per-amplicon log file tallies reads
with and without detected primer sequences, confirming correct masking.

Masking also reports the **number of masked reads per primer per reference**: a
per-sample table ``<sample>_masked_reads_per_primer.tsv`` (under
``sample_info/<sample>/mask_primers/``) and a combined
``<name>_masked_reads_per_primer.tsv`` across all samples (under
``sample_reports/general_reports/``), each with columns ``sample_id, ref_name,
primer_name, direction, masked_reads``. Primers that masked zero reads are included.

+--------------------------+-----------+-------------------------------------------------------+
| Parameter                | Default   | Description                                           |
+==========================+===========+=======================================================+
| ``--primer_file``        | ``null``  | Path to a primer BED **or** a primer CSV (required to |
|                          |           | enable masking). A ``.csv`` is auto-detected and      |
|                          |           | converted to a BED (see below).                       |
+--------------------------+-----------+-------------------------------------------------------+
| ``--primer_max_mismatch``| ``2``     | CSV→BED only: max mismatches allowed when locating    |
|                          |           | each primer in the reference (ambiguity codes count   |
|                          |           | as mismatches)                                        |
+--------------------------+-----------+-------------------------------------------------------+
| ``--mask_primers``       | ``null``  | Enable primer masking (auto-enabled when              |
|                          |           | ``--primer_file`` is provided; set ``false`` to force |
|                          |           | disable)                                              |
+--------------------------+-----------+-------------------------------------------------------+
| ``--wiggle``             | ``9``     | Bases outside primer boundary to include in mask      |
+--------------------------+-----------+-------------------------------------------------------+
| ``--mask_bam``           | ``true``  | Replace masked bases with ``N`` in BAM sequence field |
+--------------------------+-----------+-------------------------------------------------------+
| ``--primer_only``        | ``false`` | Retain only primer-overlapping reads; discard all     |
|                          |           | others after masking                                  |
+--------------------------+-----------+-------------------------------------------------------+

.. _Generating a primer BED from a CSV:

**Generating a primer BED from a CSV**

If you only have primer sequences (not coordinates), pass a CSV to ``--primer_file``
with a ``.csv`` extension and three columns: ``primer_name``, ``direction``
(``F``/``R``), and ``sequence``. Before masking, ASAP runs ``GENERATE_PRIMER_BED``,
which searches the pipeline reference (both strands, allowing
``--primer_max_mismatch`` mismatches/indels) for each primer and writes a
pipeline-ready 6-column BED. Three files are published to ``<outdir>/primer_bed/``:

- ``<name>_primers.bed`` — the generated 6-column primer BED.
- ``<name>_primer_search_results.csv`` — the full search table (found and not-found
  primers, both orientations, with per-match mismatch details for QC).
- ``<name>_primer_match_summary.csv`` — one row per primer with the **number of
  matches per reference sequence** (one column per reference) and a
  ``total_matches`` column; primers that matched nothing show an all-zero row.

The same generated BED is reused by primer masking, iVar trimming, and the SNP table.
A ready-made BED (any non-``.csv`` file) is used directly, unchanged.

Step 4 — Percent-Identity Filtering *(optional)*
-------------------------------------------------

*When to use:* when near-neighbor organisms co-amplify with the target and off-target
reads must be excluded before variant calling.

Read-level filtering removes reads that do not meet a minimum percent identity to their
aligned reference amplicon. Identity is computed via local Smith-Waterman alignment.
This is especially valuable for:

- **Near-neighbor co-infections** — e.g., distinguishing *M. tuberculosis* from
  non-tuberculous mycobacteria when primers amplify both.
- **Resistance gene specificity** — ensuring only reads from the precise target gene
  are used for resistance calling, preventing false calls from paralogs or related genes.

+--------------------+----------+-------------------------------------------------------------+
| Parameter          | Default  | Description                                                 |
+====================+==========+=============================================================+
| ``--identity``     | ``null`` | Minimum fractional identity threshold (e.g. ``0.97`` = 97%) |
+--------------------+----------+-------------------------------------------------------------+
| ``--filter_pairs`` | ``true`` | If either mate of a read pair fails the identity check on a |
|                    |          | checked reference, discard both mates. Set to ``false`` to  |
|                    |          | filter each mate independently.                             |
+--------------------+----------+-------------------------------------------------------------+

Step 5 — SMOR Processing *(optional)*
--------------------------------------

*When to use:* high-resolution applications where sequencing noise would otherwise
obscure true low-frequency variants —

- High-resolution AMR profiling in heteroresistant *M. tuberculosis* infections
- Detection of low-frequency viral variants in mixed infections
- Distinguishing true iSNVs from sequencing artifacts at low allele frequencies

Two complementary approaches leverage paired-end read overlap for Illumina error
correction, without requiring additional wet-lab steps.

**SMOR Masking (``--smor true``)**

Designed for assays where both reads in a pair are expected to fully overlap.
Non-overlapping read regions are masked, and discordant overlapping positions between
R1 and R2 are also masked. This approach yields an approximately 8-fold decrease in
overall error rate relative to unmasked data [10]_, enabling resolution of minor
allele frequencies that would otherwise be indistinguishable from sequencing noise.

**SMOR Correction (``--smor_correction true``)**

Does not require full overlap. Within overlapping regions, discordant positions are
resolved by selecting the base with the higher Phred quality score (threshold: Q10
default). When reads agree, quality scores are combined to produce higher-confidence
calls. This approach is advantageous when reads partially overlap and quality degrades
toward the ends of R1 or R2.

+---------------------------------------------+-----------+--------------------------------------------------------+
| Parameter                                   | Default   | Description                                            |
+=============================================+===========+========================================================+
| ``--smor``                                  | ``false`` | SMOR masking (full-overlap assays)                     |
+---------------------------------------------+-----------+--------------------------------------------------------+
| ``--smor_correction``                       | ``false`` | SMOR correction (partial-overlap assays)               |
+---------------------------------------------+-----------+--------------------------------------------------------+
| ``--smor_correction_qual_diff_threshold``   | ``10``    | Phred quality difference required between R1/R2        |
|                                             |           | bases at a mismatch for the higher-quality base        |
|                                             |           | to be selected during SMOR correction; otherwise       |
|                                             |           | the position is masked with ``--fill_character``       |
+---------------------------------------------+-----------+--------------------------------------------------------+
| ``--smor_correction_agreement_method``      | ``sum``   | How agreeing R1/R2 base qualities are combined         |
|                                             |           | during SMOR correction: ``sum`` = min(q1+q2, 60),      |
|                                             |           | ``max`` = max(q1, q2)                                  |
+---------------------------------------------+-----------+--------------------------------------------------------+

Step 6 — ASAP BAM Processing
------------------------------

The core analysis step. ``ASAPBamProcessor.py`` reads the assay JSON and the aligned
(optionally masked/filtered/SMOR'd) BAM to produce a per-sample XML containing:

- Aligned read counts per amplicon
- Per-position depth, breadth of coverage, and consensus sequence
- SNPs / iSNVs detected above the proportion and depth thresholds, with per-base distributions
- Per-SNP QC data: a strand/mate distribution (per-base R1-vs-R2 read counts, plus an
  SE count for single-end reads) and per-base quality metrics (mean, median, min, max),
  supporting the strand-bias and base-quality plots in `Step 8`_
- Significance calls based on rules defined in the assay JSON (if provided)
- Optional codon-aware SNP linkage (``<codon_merge>``) and read-based phasing of
  linked SNPs (``<linked_snps>``)

+----------------------------------+----------+------------------------------------------------------------+
| Parameter                        | Default  | Description                                                |
+==================================+==========+============================================================+
| ``--depth``                      | ``100``  | Minimum read depth to consider a position covered          |
+----------------------------------+----------+------------------------------------------------------------+
| ``--breadth``                    | ``0.8``  | Minimum breadth of coverage to call an amplicon present    |
+----------------------------------+----------+------------------------------------------------------------+
| ``--proportion``                 | ``0.1``  | Minimum allele frequency to call a SNP / iSNV              |
+----------------------------------+----------+------------------------------------------------------------+
| ``--mutation_depth``             | ``5``    | Minimum read count to call a SNP / iSNV                    |
+----------------------------------+----------+------------------------------------------------------------+
| ``--min_base_qual``              | ``5``    | Minimum Phred base quality score                           |
+----------------------------------+----------+------------------------------------------------------------+
| ``--consensus_proportion``       | ``0.8``  | Minimum frequency to call a consensus base (else ``N``)    |
+----------------------------------+----------+------------------------------------------------------------+
| ``--fill_character``             | ``N``    | Character written at masked / gap positions (used by SMOR  |
|                                  |          | masking and bam_processor)                                 |
+----------------------------------+----------+------------------------------------------------------------+
| ``--fill_gaps``                  | ``n``    | Character written at zero-coverage positions in consensus  |
|                                  |          | sequence                                                   |
+----------------------------------+----------+------------------------------------------------------------+
| ``--mark_deletions``             | ``_``    | Character written at deletion positions in consensus       |
+----------------------------------+----------+------------------------------------------------------------+
| ``--suppress_per_base``          | ``false``| Suppress all per-position arrays (see --prune_per_base)    |
+----------------------------------+----------+------------------------------------------------------------+
| ``--prune_per_base``             | ``false``| Retain per-position arrays only where depth >= --depth     |
+----------------------------------+----------+------------------------------------------------------------+
| ``--asap_snps``                  | ``true`` | Enable ASAP BAM processing (set ``false`` to skip)         |
+----------------------------------+----------+------------------------------------------------------------+
| ``--combine_output``             | ``true`` | Combine per-sample XMLs and generate HTML report           |
+----------------------------------+----------+------------------------------------------------------------+
| ``--stylesheet``                 | bundled  | XSLT stylesheet for HTML report generation                 |
+----------------------------------+----------+------------------------------------------------------------+
| ``--codon_correction``           | ``false``| Annotate same-codon SNP pairs with read-level allele-      |
|                                  |          | linkage info (``<codon_merge>``)                           |
+----------------------------------+----------+------------------------------------------------------------+
| ``--codon_correction_error``     | ``0.05`` | Frequency tolerance (0–1) for 'complete' vs. 'partial'     |
|                                  |          | codon-merge linkage                                        |
+----------------------------------+----------+------------------------------------------------------------+
| ``--codon_correction_min_reads`` | ``10``   | Minimum spanning reads to confirm codon-level SNP linkage  |
+----------------------------------+----------+------------------------------------------------------------+
| ``--discover_roi``               | ``false``| Read-based phasing: annotate each SNP with                 |
|                                  |          | other variants phased onto the same reads                  |
|                                  |          | (``<linked_snps>``)                                        |
+----------------------------------+----------+------------------------------------------------------------+
| ``--discover_roi_min_perc``      | ``0.1``  | Minimum co-occurrence proportion (0–1) to report a linked  |
|                                  |          | SNP                                                        |
+----------------------------------+----------+------------------------------------------------------------+
| ``--discover_roi_min_reads``     | ``10``   | Minimum co-occurring read count to report a linked SNP     |
+----------------------------------+----------+------------------------------------------------------------+
| ``--discover_roi_min_snp_perc``  | ``0.05`` | Minimum variant frequency (0–1) for a SNP to be considered |
|                                  |          | in discover-roi linkage analysis                           |
+----------------------------------+----------+------------------------------------------------------------+

**SNP Call Quality Control**

Every SNP is emitted with two per-base QC annotations that help distinguish true variants
from sequencing or library artifacts:

- **Base quality** (``<base_quality>``) — for each observed base, the mean, median, min,
  and max Phred quality of the supporting reads. A variant supported only by low-quality
  bases (especially relative to the reference base at the same position) is a likely
  artifact rather than a real iSNV.
- **Strand / mate distribution** (``<base_strand_distribution>``) — for each observed
  base, the number of supporting reads split by mate: ``R1`` vs ``R2`` (with an ``SE``
  count for single-end data). Depending on amplification techniques, a genuine variant is 
  typically seen on both mates in proportion to coverage; a call that appears almost 
  exclusively on R1 or R2 indicates **strand bias** and is treated with suspicion.

These annotations feed the base-quality comparison and strand-bias plots produced by ASAP
Tools (see `Step 8`_), so calls can be reviewed visually across all samples and amplicons.

**Codon-Aware SNP Linkage & Read-Based Phasing**

Two optional, GenBank-aware analyses annotate each SNP with how its variant allele
co-occurs with other variants on the same sequencing read/fragment:

*Codon-aware linkage* (``--codon_correction``) — for SNP pairs that fall within the
same codon (per GenBank CDS annotations), each SNP is annotated with a
``<codon_merge>`` element describing how the pair's alleles co-occur on the same
fragments. Each observed combination of the two positions' bases is reported as a
``<combo>`` and classified as:

- ``reference`` — neither variant (both positions match the reference)
- ``variant`` — both variants together, i.e. the actual combined codon change
- ``discordant`` — only one of the two variants present (e.g. sequencing noise)

The pair is also flagged as ``complete`` (both SNPs' frequencies agree closely
enough that they represent the same underlying change) or ``partial``
(frequencies diverge, suggesting independent or partially-linked events) via the 
``--codon_correction_error`` parameter. This flows through to ASAP Tools' 
SNP / Amino Acid table (see `Step 8`_).

*Read-based phasing* (``--discover_roi``) — more general than codon-aware linkage,
this phases *any* SNP against the reads that span it, annotating it with a
``<linked_snps>`` list of other SNPs (regardless of codon membership) whose variant
alleles fall on the same reads/fragments. Because phasing is established directly from
individual reads rather than inferred, it is useful for resolving candidate haplotype
structure or regions of interest spanning multiple variants. This process calculates the
prevalence of linked SNPs relative to the number of spanning reads, providing a realistic
prevalence values.

Step 7 — iVAR Processing *(optional)*
---------------------------------------

iVAR [1]_ provides an alternative or complementary variant calling and consensus
generation workflow. Enable the full iVAR workflow with ``--ivar true``, or enable
individual steps independently.

+---------------------------------+-----------+-----------------------------------------------------+
| Parameter                       | Default   | Description                                         |
+=================================+===========+=====================================================+
| ``--ivar``                      | ``false`` | Enable all iVAR steps (trim + variants + consensus) |
+---------------------------------+-----------+-----------------------------------------------------+
| ``--ivar_trim``                 | ``false`` | iVAR primer trimming only (requires primer BED)     |
+---------------------------------+-----------+-----------------------------------------------------+
| ``--ivar_variants``             | ``false`` | iVAR variant calling only                           |
+---------------------------------+-----------+-----------------------------------------------------+
| ``--ivar_consensus``            | ``false`` | iVAR consensus calling only                         |
+---------------------------------+-----------+-----------------------------------------------------+
| ``--ivar_trim_extra_args``      | ``""``    | Additional ``ivar trim`` arguments                  |
+---------------------------------+-----------+-----------------------------------------------------+
| ``--ivar_variants_extra_args``  | ``""``    | Additional ``ivar variants`` arguments              |
+---------------------------------+-----------+-----------------------------------------------------+
| ``--ivar_consensus_extra_args`` | ``""``    | Additional ``ivar consensus`` arguments             |
+---------------------------------+-----------+-----------------------------------------------------+

.. _Step 8:

Step 8 — ASAP Tools R Post-Processing *(optional)*
----------------------------------------------------

A suite of R scripts transforms per-sample XML outputs into tabular summaries,
figures, and FASTA files. Processing follows a fan-out / gather pattern:

.. code-block:: none

   PROCESS_XML_R         (per sample, parallel)  →  sample_info/<id>/rdata/
         │
         ▼
   PROCESS_COMBINE_RDATA (gather all)             →  sample_reports/rdata/ (Rdata)
         │                                           sample_reports/general_reports/ (CSV)
         ├── PROCESS_GENERATE_COV_TABLE           →  sample_reports/general_reports/ (Excel)
         ├── PROCESS_QC_PLOTS                     →  sample_reports/plots/ (JPG, +HTML if enabled)
         ├── PROCESS_FASTP_PANEL                  →  sample_reports/plots/ (JPG, +HTML if enabled)
         ├── PROCESS_SNP_PLOTS                    →  sample_reports/plots/ (JPG, +HTML if enabled)
         ├── PROCESS_GENERATE_FASTA               →  sample_reports/fasta/ (FASTA)
         ├── PROCESS_SNPS_TO_AMINOACIDS           →  sample_reports/rdata/ (Rdata)
         │                                           sample_reports/snp_reports/ (CSV)
         └── PROCESS_GENERATE_SNP_TABLE           →  sample_reports/snp_reports/ (CSV ± Excel)

+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| Parameter                             | Default   | Description                                                                                                                |
+=======================================+===========+============================================================================================================================+
| ``--asaptools_processing``            | ``true``  | Enable R post-processing                                                                                                   |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_cov_table``             | ``true``  | Generate coverage depth table                                                                                              |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_qc_plots``              | ``true``  | Generate QC figures (depth of coverage, % N bases, breadth heatmap, alignment summary, read funnel) and the FastP QC panel |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_snp_plots``             | ``true``  | Generate SNP figures (position prevalence, proportion density, strand bias, base quality, GenBank-driven genome track)     |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_generate_fasta``        | ``true``  | Export consensus FASTA files                                                                                               |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_snp_table``             | ``true``  | Generate SNP / iSNV table                                                                                                  |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_snp_table_xls``         | ``false`` | Also export SNP table as Excel                                                                                             |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_interactive_plots``     | ``false`` | Also export the self-contained interactive HTML widgets for the QC and SNP figures. JPGs are always produced; the HTML     |
|                                       |           | widgets are expensive to render, so they are off by default.                                                               |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_positions_of_interest`` | ``null``  | CSV of genomic positions to annotate in outputs                                                                            |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_genbank_location``      | ``null``  | GenBank file for amino acid annotation                                                                                     |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_snp_proportion``        | ``null``  | Override allele frequency threshold for SNP table                                                                          |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_max_sample_snp_count``  | ``10000`` | Max SNPs per sample before flagging as noisy                                                                               |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_samples_to_remove``     | ``null``  | Sample IDs to exclude from combined outputs                                                                                |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+
| ``--asaptools_breadth_threshold``     | ``null``  | Minimum breadth for FASTA export (uses ``--breadth`` if unset)                                                             |
+---------------------------------------+-----------+----------------------------------------------------------------------------------------------------------------------------+

**Outputs produced**

The R suite provides ready-to-use outputs for sample QC and downstream analysis without
requiring users to write custom analysis code:

**Coverage Table** — per-sample, per-amplicon statistics including breadth of coverage,
average depth, aligned reads, and depth at positions of interest. Depth thresholds flag
positions without adequate coverage, distinguishing true absence of SNPs from
insufficient data.

**QC Figures** — static JPG figures (and, when ``--asaptools_interactive_plots true``,
self-contained interactive HTML versions) produced by three independent, individually
toggleable steps, all published to ``sample_reports/plots/``:

- *QC Plots* (``--asaptools_qc_plots``) — reference depth of coverage, percent masked
  bases ("N"s), a per-sample/per-amplicon breadth-of-coverage heatmap, an alignment
  summary (read fate through trimming and alignment), and a read-fate funnel per amplicon.
- *SNP Plots* (``--asaptools_snp_plots``) — SNP position prevalence and proportion
  density across the panel, strand-bias assessment, base-quality comparison (called SNP
  vs. reference base), and — when a GenBank reference is provided and the genome is
  under 2 Mb — a combined gene-track / SNP-density / per-sample heatmap figure.
- *FastP Panel* (bundled with ``--asaptools_qc_plots``) — pre/post-trim read counts,
  filtering breakdown, Q20/Q30 rates, GC content, duplication rate, insert size, and
  mean base quality per sequencing cycle, from the fastp/fastplong JSON reports.

**SNP / iSNV Table** — a comprehensive table of detected variants across all samples and
amplicons. User-defined frequency and depth thresholds separate true iSNVs from noise.
When a GenBank reference is provided, coding-region SNPs are annotated with the
resulting amino acid change, enabling immediate identification of resistance-conferring
mutations. When ``--codon_correction`` finds a same-codon SNP pair with ``complete``
linkage, the table reports a single combined row (e.g. ``T5118A|T5119A``) with one
amino acid annotation for the joint codon change instead of two separate
single-position rows; ``partial`` pairs show both the individual rows and the
combined row.

**Consensus FASTA Export** — consensus sequences for each sample and amplicon, exported
at a user-defined breadth-of-coverage threshold. Suitable for downstream phylogenetic
analysis or genome assembly. If no sample/assay combination meets the breadth threshold,
a placeholder ``*_WARNING_no_sequences_passed_breadth_filter.fasta`` file is written
instead of silently producing no output.

.. _Step 9:

Step 9 — MultiQC
-----------------

MultiQC [9]_ aggregates reports from FastQC (initial and post-trim),
fastp/fastplong JSON, and SAMtools flagstat into a single interactive HTML report,
published to ``<outdir>/sample_reports/multiqc/``.

**FastQC and MultiQC are skipped by default** (``--skip_fastqc`` and
``--skip_multiqc`` both default to ``true``). To enable them, pass
``--skip_fastqc false`` and/or ``--skip_multiqc false``. When FastQC is skipped but
MultiQC runs, the MultiQC report is built from fastp JSON and flagstat only.

References
==========

.. [1] Grubaugh ND, Gangavarapu K, Quick J, *et al.* An amplicon-based sequencing
       framework for accurately measuring intrahost virus diversity using PrimalSeq
       and iVar. *Genome Biology*. 2019;20(1):8.

.. [6] Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ
       preprocessor. *Bioinformatics*. 2018;34(17):i884–i890.

.. [8] Andrews S. FastQC: A quality control tool for high throughput sequence data.
       2010. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

.. [9] Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results
       for multiple tools and samples in a single report.
       *Bioinformatics*. 2016;32(19):3047–3048.

.. [10] [SMOR citation — TGen internal or forthcoming publication]

.. [CITATION] fastplong — citation pending.

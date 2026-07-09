# ASAP Developer Reference

This document is the technical companion to [README.rst](README.rst), aimed at
contributors who need to understand pipeline internals, module responsibilities,
and the per-sample output schema. README.rst covers installation and running the
pipeline; this document covers how it works and how to extend it.

---

## 1. Overview & Audience

ASAP (Amplicon Sequencing Analysis Pipeline) is a Nextflow workflow for
processing targeted amplicon sequencing data. The processing chain is:

```
FASTQ  →  QC (fastp/fastplong + FastQC)  →  Alignment  →  optional filtering
       →  ASAP BAM processing (per-sample XML)
       →  optional R post-processing (tables, figures, FASTA, AA annotations)
       →  HTML report
```

The pipeline is pre-release and under active development at TGen North. It is
built around three main sub-systems:

| Sub-system | Directory | Language | Role |
|---|---|---|---|
| Core Python | `asap/` | Python 3 | Preprocessing, BAM analysis, output |
| Nextflow workflow | `nextflow/` | Nextflow DSL2 | Orchestration, conda/SLURM management |
| R downstream tools | `asap_tools/` | R | Tables, figures, FASTA export, AA annotation |

---

## 2. Repository Layout

```
ASAP/
├── asap/                        # Python package (core analysis + pre/post-processing scripts)
│   ├── __init__.py              # Package version (1.9.0) and author metadata
│   ├── assayInfo.py             # Assay data model, JSON de/serialization, FASTA generation
│   ├── prepareJSONInput_nextflow.py  # Convert FASTA/GenBank/Excel to ASAP JSON
│   ├── maskPrimers.py           # Primer masking (soft-clip BAM bases in primer regions)
│   ├── identityFilter.py        # Percent-identity filtering of alignments
│   ├── generateSMORbam.py       # SMOR consensus BAM generation (overlap merging)
│   ├── generateSMORbam_correction.py  # SMOR with quality-score-driven mismatch correction
│   ├── newBamProcessor.py       # Core BAM analysis: pileup → SNP calls → per-sample XML
│   ├── outputCombiner.py        # Merge per-sample XMLs into analysis-level XML
│   ├── formatOutput.py          # XSLT transform XML → HTML report
│   └── outputData.py            # Minimal XML utility functions (shared helpers)
│
├── asap_tools/                  # R downstream processing pipeline
│   ├── process_xml.R            # Per-sample entry: XML → Rdata
│   ├── process_combine_rdata.R  # Gather: merge all per-sample Rdata → combined Rdata
│   ├── process_asaptools_cov_table.R     # Coverage summary Excel report
│   ├── process_asaptools_generate_figures.R  # Interactive/static QC plots
│   ├── process_asaptools_fasta_export.R  # Consensus FASTA export per assay
│   ├── process_asaptools_snps_amino_acids.R  # SNP → gene position → AA translation
│   ├── process_asaptools_snp_table.R    # Final SNP table (CSV ± Excel)
│   ├── process_primers_to_bed.R # Primer CSV → primer search → 6-col BED (self-contained)
│   └── asap_tools_functions/    # Shared R functions (sourced by scripts above)
│       ├── _read.ASAP.individual.R      # Parse per-amplicon metrics from XML
│       ├── _read.ASAP.snps.individual.R # Parse SNP-level data including codon_merge XML
│       ├── _expand_codon_merges.R       # Expand codon_merge annotations into combined SNP rows
│       ├── _genome.snp.to.gene.snp.R    # Map genomic SNP positions to gene coordinates
│       ├── _snps.to.amino.R             # Translate codon changes to amino acid changes
│       ├── _ASAP.get.depth.R            # Unpack comma-delimited depth arrays
│       ├── _ASAP.get.nreads.R           # Unpack N-read count arrays
│       ├── _ASAP.get.proportions.R      # Unpack allele proportion arrays
│       └── _ASAP.get.quality.discards.R # Unpack quality-discard count arrays
│
├── nextflow/                    # Nextflow DSL2 pipeline
│   ├── main.nf                  # Main workflow orchestration
│   ├── nextflow.config          # Profile definitions (slurm, conda, singularity)
│   ├── nextflow_schema.json     # nf-schema parameter schema (validation + docs)
│   ├── bin/                     # Symlinks/copies of all scripts above (Nextflow PATH)
│   ├── modules/
│   │   ├── asap/main.nf         # Nextflow processes for all ASAP Python scripts
│   │   ├── asap_tools/main.nf   # Nextflow processes for all R asaptools scripts
│   │   ├── bowtie2/             # Bowtie2 alignment
│   │   ├── bwa/                 # BWA-MEM alignment
│   │   ├── fastp/               # Illumina read trimming
│   │   ├── fastplong/           # ONT/PacBio read trimming
│   │   ├── fastqc/              # Per-sample read QC
│   │   ├── ivar/                # iVar variant calling + primer trimming + consensus
│   │   ├── minimap2/            # Long-read / universal alignment
│   │   └── multiqc/             # Aggregated QC report
│   └── tests/                   # nf-test end-to-end test suite
│       └── ASAP_EtE.nf.test     # Tagged test definitions (sc2, tb, ont, ivar, …)
│
├── output_transforms/           # XSLT stylesheets for HTML report generation
├── assay_data/                  # Example assay JSONs, Excel templates, reference FASTAs
│   └── references/              # Companion reference FASTA files for bundled assays
├── tests/                       # Python pytest unit test suite
├── asap_tutorials/              # Tutorial data and example runs
├── web_resources/               # Static assets for HTML reports (logos, icons)
└── ASAP_nextflow_env.yml        # Conda YAML for the Nextflow runner environment
```

---

## 3. Pipeline Architecture

### 3.1 Main Workflow Orchestration (`nextflow/main.nf`)

See the ASCII diagram in README.rst (Pipeline Summary section) for a visual
overview. The main.nf workflow is organized into five phases:

**Reference setup**: The input reference is auto-detected by file extension
(`.gb`/`.gbk`/`.genbank` → GenBank; `.fasta`/`.fa` → FASTA; `.xlsx`/`.xls`
→ Excel; `.json` → direct JSON passthrough). GenBank and FASTA/Excel inputs
flow through `PREPARE_ASAP_JSON` to produce `assay_input.json`, which is
then used by `GENERATE_REFERENCE_FASTA` to extract a FASTA alignment target.
If `params.codon_correction` is enabled, a GenBank reference is mandatory;
the workflow errors immediately if only a FASTA/Excel/JSON reference is
provided.

**QC and trimming**: Read technology (`params.technology = illumina | ont |
ont.v14 | pacbio`) controls the trimming branch. Illumina reads go through
`RUN_FASTP` (adapter file from `params.adapter_fasta` or the bundled
`illumina_adapters_all.fasta`); long reads go through `FASTPLONG`. FastQC
runs before and after trimming. Paired-end vs. single-end detection is
automatic.

**Alignment and filtering chain**: Alignment uses Bowtie2 (Illumina default),
BWA (`params.aligner = bwa`), or minimap2 (auto-selected for ONT/PacBio, or
forced with `params.aligner = minimap2`). After alignment, a conditional
cascade applies optional filters in sequence:

1. `MASK_PRIMERS` — if `params.mask_primers` or `params.primer_file` is set.
   When `params.primer_file` is a `.csv`, `GENERATE_PRIMER_BED` runs first to
   build the BED (see §3.3); the resulting value channel
   (`effective_primer_bed_ch`) is shared by masking, `IVAR_TRIM`, and
   `PROCESS_GENERATE_SNP_TABLE`.
2. `IDENTITY_FILTER` — if `params.identity != null`
3. `SMOR` or `SMOR_CORRECTION` — if `params.smor` or `params.smor_correction`
   (`SMOR_CORRECTION` takes precedence over `SMOR`)

A snapshot of the pre-filter (original) BAM is passed alongside the filtered
BAM to `PROCESS_BAM` so that per-amplicon `aligned_reads` counts reflect
total mapped reads, not post-filter counts.

**ASAP / iVAR fork**: The filtered BAM is consumed by two optional parallel
paths: the ASAP path (`params.asap_snps = true`, default) runs
`PROCESS_BAM`; the iVAR path (`params.ivar`, `params.ivar_trim`,
`params.ivar_variants`, `params.ivar_consensus`) runs the iVar tool suite.

**asaptools fan-out / gather**: When both `params.asap_snps` and
`params.asaptools_processing` are true, per-sample XMLs feed into
`PROCESS_XML_R` (runs in parallel, one job per sample). Outputs are gathered
by `PROCESS_COMBINE_RDATA` into a single combined Rdata, which then drives
up to five independent downstream processes:

| Process | Param gate |
|---|---|
| `PROCESS_GENERATE_COV_TABLE` | `params.asaptools_cov_table` |
| `PROCESS_QC_PLOTS` | `params.asaptools_qc_plots` |
| `PROCESS_GENERATE_FASTA` | `params.asaptools_generate_fasta` |
| `PROCESS_SNPS_TO_AMINOACIDS` | `params.asaptools_snp_table` + GenBank file present |
| `PROCESS_GENERATE_SNP_TABLE` | `params.asaptools_snp_table` |

**Output consolidation**: If `params.combine_output` is true,
`OUTPUT_COMBINER` merges all per-sample XMLs into a single
`{file_name}_analysis.xml`, which `FORMAT_OUTPUT` transforms into an HTML
report via the XSLT stylesheet at `params.stylesheet` (default:
`default_stylesheet/ASAP_fulldetails_web.xsl`). `MULTIQC` aggregates FastQC,
trimmer JSON, and alignment flagstats into a QC report.

### 3.2 Process Table — `nextflow/modules/asap/main.nf`

| Process | Script | Key inputs | Key outputs | Governing params |
|---|---|---|---|---|
| `PREPARE_ASAP_JSON` | `prepareJSONInput_nextflow.py` | FASTA/GenBank/Excel reference(s) | `assay_input.json` | format auto-detected |
| `GENERATE_REFERENCE_FASTA` | `assayInfo.py` | `assay_input.json` | `reference.fasta` | always |
| `MASK_PRIMERS` | `maskPrimers.py` | BAM+BAI, primer BED | masked BAM+BAI, `primer_masking_stats.tsv`, `{id}_masked_reads_per_primer.tsv` (per-primer/ref counts; also gathered into `sample_reports/general_reports/{name}_masked_reads_per_primer.tsv` via `collectFile` in `main.nf`) | `params.primer_file`, `params.wiggle`, `params.mask_bam`, `params.primer_only` |
| `IDENTITY_FILTER` | `identityFilter.py` | BAM+BAI | filtered BAM+BAI, `identity_filter_stats.tsv` | `params.identity`, `params.filter_pairs` |
| `SMOR` | `generateSMORbam.py` | BAM+BAI | SMOR BAM+BAI, `smor_stats.tsv` | `params.smor`, `params.fill_character` |
| `SMOR_CORRECTION` | `generateSMORbam_correction.py` | BAM+BAI | SMOR BAM+BAI, `smor_stats.tsv` | `params.smor_correction`, `params.qual_diff_threshold` |
| `PROCESS_BAM` | `newBamProcessor.py` | filtered BAM, original BAM, fastp JSON, stats TSVs, assay JSON, optional GenBank | `{sample_id}.xml` | `params.depth/breadth/proportion/mutation_depth/min_base_qual/consensus_proportion/fill_gaps/mark_deletions/whole_genome/codon_correction*/discover_roi*` |
| `OUTPUT_COMBINER` | `outputCombiner.py` | collected XMLs | `{file_name}_analysis.xml` | `params.combine_output` |
| `FORMAT_OUTPUT` | `formatOutput.py` | analysis XML, XSLT stylesheet | HTML report | `params.stylesheet`, `params.out_file` |

### 3.3 Process Table — `nextflow/modules/asap_tools/main.nf`

| Process | Script | Key inputs | Key outputs | Governing params |
|---|---|---|---|---|
| `GENERATE_PRIMER_BED` | `process_primers_to_bed.R` | primer CSV (`primer_name,direction,sequence`), reference FASTA | `{name}_primers.bed`, `{name}_primer_search_results.csv`, `{name}_primer_match_summary.csv` (wide: matches per primer per reference) (→ `<outdir>/primer_bed/`) | `params.primer_file` ends in `.csv`, `params.primer_max_mismatch`, `task.cpus` |
| `PROCESS_XML_R` | `process_xml.R` | per-sample XML, min proportion | `{id}_XML_Data.Rdata`, `{id}_Summary.csv` | `params.proportion` |
| `PROCESS_COMBINE_RDATA` | `process_combine_rdata.R` | all Rdata files (gathered), optional POI CSV | `{name}_ASAP_Data.Rdata`, `{name}_Summary.csv` | `params.asaptools_positions_of_interest`, `params.file_name` |
| `PROCESS_GENERATE_COV_TABLE` | `process_asaptools_cov_table.R` | combined Rdata | `*_Coverage_Report.xlsx` | `params.asaptools_cov_table`, `params.depth` |
| `PROCESS_QC_PLOTS` | `process_asaptools_generate_figures.R` | combined Rdata | `*.html`, `*.jpg` (3 × 2) | `params.asaptools_qc_plots` |
| `PROCESS_GENERATE_FASTA` | `process_asaptools_fasta_export.R` | combined Rdata | `*_{assay}.fasta` (one per assay) | `params.asaptools_generate_fasta`, `params.asaptools_breadth_threshold` |
| `PROCESS_SNPS_TO_AMINOACIDS` | `process_asaptools_snps_amino_acids.R` | combined Rdata, GenBank file(s) | `SNP_Amino_Acid_Table.Rdata` | `params.asaptools_snp_table` + GenBank present |
| `PROCESS_GENERATE_SNP_TABLE` | `process_asaptools_snp_table.R` | combined Rdata, AA Rdata, optional GenBank/BED/POI | `*_SNP_Table_*.csv` (4 files), optional `*.xlsx` | `params.asaptools_snp_table`, `params.asaptools_snp_proportion`, `params.asaptools_max_sample_snp_count`, `params.depth`, `params.asaptools_snp_table_xls` |

---

## 4. Python Modules (`asap/`)

### `__init__.py`

Package metadata: `__version__ = "1.9.0"`, author Darrin Lemmer at TGen North.

### `assayInfo.py`

**Role**: Defines the core assay data model and provides JSON serialization/
deserialization. Used as a library by `newBamProcessor.py` and as a CLI that
converts an assay JSON to FASTA.

**Library API** (consumed by `newBamProcessor.py` and `prepareJSONInput_nextflow.py`):

- `parseJSON(file_path)` → list of `Assay` objects. A custom JSON hook
  (`_json_decode`) reconstructs typed objects by pattern-matching dict keys.
- `parseOperation(file_path)` → list of `Operation` objects (AND/OR/NOT logic
  trees used for cross-amplicon significance evaluation).
- `writeJSON(assay_data, filename)` → serializes to JSON.
- `generateReference(assay_list)` → yields FASTA records.

**Key classes**: `Assay` (name, type, target, percid), `Target` (function,
gene\_name, start/end positions, reverse\_comp, amplicons), `Amplicon` (sequence,
variant\_name, SNPs, significance, percid), `SNP` (position, reference, variant,
name, significance), `Significance` (message, resistance string),
`Operation`/`ITEM` (logic tree for `evaluateOperation` in `newBamProcessor.py`).

Valid `assay_type` values: `"presence/absence"`, `"SNP"`, `"gene variant"`,
`"ROI"`, `"mixed"`. Valid `target_function` values: `"species ID"`,
`"strain ID"`, `"resistance type"`, `"virulence factor"`.

**CLI**: `python -m asap.assayInfo <json_file>` prints assay sequences as FASTA.

### `prepareJSONInput_nextflow.py`

**Role**: Converts assay definitions from FASTA, GenBank, or Excel format into
the ASAP JSON schema expected by `newBamProcessor.py`.

**CLI**: `prepareJSONInput_nextflow.py -o <out.json> (-f <fasta> | -g <gb> [<gb>…] | -x <excel> [-w worksheet])`

**Core logic**: For FASTA input, each sequence ID becomes an assay name; for
GenBank, the filename is the assay name; for Excel, a row-based parser constructs
assays with multiple amplicons and SNPs per target. Delegates to
`assayInfo.writeJSON()` for output.

### `maskPrimers.py`

**Role**: Removes primer sequences from BAM alignments by zeroing Phred quality
scores and optionally converting bases to `'N'` for the primer-overlapping region.

**CLI**: `maskPrimers.py -b <bam> -p <primer_bed> [-o <out.bam>] [--wiggle INT] [--mask-bam | --no-mask-bam] [--primer-only | --no-primer-only]`

**Core logic**: Loads primer coordinates (format: CHROM, start, end, name,
unused, strand `+`/`-`) into a NumPy structured array. For each read,
checks whether its start falls within a forward primer window or its end falls
within a reverse primer window (with `--wiggle` tolerance). Matched primer
regions are masked by zeroing the quality array; with `--mask-bam`, the
corresponding sequence bases are also replaced with `'N'`. Uses
`get_aligned_pairs()` to accurately map reference to query coordinates.

**Outputs**: `primer_masking.tsv` (per-read log), `primer_masking_stats.tsv`
(per-reference aggregate — consumed by `newBamProcessor.py` via `--primer-stats`),
and `primer_masking_primer_stats.tsv` (per-reference **per-primer** masked-read
counts: `ref_name, primer_name, direction, masked_reads`; includes primers that
masked zero reads). The pipeline prepends a `sample_id` column to the last file for
its per-sample and combined masked-reads-per-primer reports.

### `identityFilter.py`

**Role**: Filters reads from a BAM by a minimum percent-identity threshold,
with optional pair-aware filtering (both mates dropped if one fails).

**CLI**: `identityFilter.py -b <bam> -i <threshold> [-r <ref_name>…] [-o <out.bam>] [--filter-pairs | --no-filter-pairs]`

**Core logic**: Two-pass: Pass 1 collects query names of reads failing the
identity check on specified references (computed as exact base matches /
aligned length via `get_aligned_pairs(with_seq=True)`, which works from the
MD tag). Pass 2 writes the output, marking failed reads as unmapped or (with
`--filter-pairs`) dropping both mates of a failing pair.

**Outputs**: `identity_filter_stats.tsv` (per-reference counts — consumed by
`newBamProcessor.py` via `--identity-stats`).

### `generateSMORbam.py`

**Role**: Generates SMOR (Stacked/Merged Overlapping Read) consensus alignments
by merging overlapping paired-end reads into synthetic single-fragment records.

**CLI**: `generateSMORbam.py -b <bam> [-o <out.bam>] [-c <fill_char>] [-q <min_base_qual>] [-w]`

**Core logic**: Groups name-sorted reads by query name. For each overlapping
pair, `_get_consensus()` iterates both reads' aligned pairs in lock-step: where
bases agree the base is used; where they disagree, `--fill-character` is written
and the lower quality is retained. The result is a new synthetic read spanning
the full fragment, with a dynamically built CIGAR list.

**Outputs**: `smor_stats.tsv` (per-reference: input reads, pairs dropped,
consensus reads, singleton reads — consumed by `newBamProcessor.py` via
`--smor-stats`).

### `generateSMORbam_correction.py`

**Role**: Variant of SMOR generation that selects between mismatching bases
based on Phred quality difference rather than masking all mismatches.

**CLI**: `generateSMORbam_correction.py -b <bam> [-o <out.bam>] [-c <fill_char>] [-q <qual_diff_threshold>] [-l <logfile>]`

**Core logic**: Differs from `generateSMORbam.py` at the disagreement step:
if the quality difference between mismatching bases exceeds
`--qual-diff-threshold`, the higher-quality base is selected instead of the
fill character; otherwise the fill character is used. Quality scores at
agreeing positions are summed (capped at 60). CIGAR operations are explicitly
tracked to prevent coordinate shifts from deletions.

### `outputCombiner.py`

**Role**: Merges per-sample XML files from a directory into a single
run-level `<analysis>` XML document.

**CLI**: `outputCombiner.py -n <run_name> -x <xml_dir> [-o <out.xml>]`

**Core logic**: Creates `<analysis run_name="...">` root element, iterates
alphabetically sorted XMLs in `<xml_dir>`, parses each, and appends the root
element as a child. Pretty-prints via minidom.

### `formatOutput.py`

**Role**: Applies an XSLT stylesheet to the combined analysis XML to produce
an HTML report (or other text output).

**CLI**: `formatOutput.py -s <stylesheet.xsl> -x <analysis.xml> [-o <out_file>] [-d <out_dir>] [-t]`

**Core logic**: Uses lxml to load XML and XSLT; registers a custom XPath
`distinct-values()` function for stylesheet use. Applies the transformation
and writes the result as HTML (or plain text with `-t`).

### `outputData.py`

**Role**: Minimal XML helper module — `createXMLNode`, `writeOutput`,
`parse`, `_write_parameters`, `_parse_parameters`. Provides basic
ElementTree wrappers for pipeline-internal metadata records. Not used in the
main analysis path.

### `newBamProcessor.py` (deep reference)

`newBamProcessor.py` (1,634 lines, current version) is the core analysis
engine. Its structure:

#### Pileup → SNP calling (lines 53–354)

- `_process_pileup` (line 131): Main per-amplicon analysis loop. Iterates
  the pysam pileup, counts bases per position (with quality filtering and
  deletion handling), computes consensus sequence, and builds an SNP list.
  For each position in the assay's SNP dictionary (`snp_dict`, from
  `_create_snp_dict`), it records depth, base distribution, and proportion
  regardless of whether a call was made. De-novo SNPs (positions with depth
  > threshold and a variant > mutdepth at proportion > threshold, not in
  the assay definition) are auto-named `{ref}{position}{variant}`.
- `_create_snp_dict` (line 313): Builds a position-keyed dict from the
  assay's `Amplicon.SNPs`; handles wildcard variants (`any`) by expanding
  to all three non-reference bases.
- `_add_snp_node` (line 331): Serializes a SNP dict to an XML `<snp>`
  element with `<snp_call>`, optional `<significance>`, and
  `<base_distribution>` children. Returns the created node so callers can
  append further children (`<codon_merge>`, `<linked_snps>`).

#### GenBank CDS parsing (lines 356–506)

- `CdsFeature` (line 360): namedtuple — `name`, `strand`, `codon_boundaries`
  (list of `(amp_start, amp_end)` 0-based amplicon-local tuples, one per codon).
- `_load_genbank_records` (line 365): LRU-cached skbio GenBank reader; avoids
  re-parsing the same file for every amplicon.
- `_parse_genbank_cds` (line 374): Extracts CDS features from one or more
  GenBank files that overlap a given amplicon sequence. Locates the amplicon
  by exact string match (with reverse-complement fallback); falls back to
  `local_pairwise_align_nucleotide` for small GenBank records (<100 kb) where
  exact match fails. Handles compound CDS locations (joins), `codon_start`
  offsets, and both forward and reverse strand features. Returns a list of
  `CdsFeature` objects with amplicon-local codon boundary coords. Returns `[]`
  (with a `logging.warning`) if the amplicon cannot be located or no CDS
  features overlap.

#### Read-level allele linkage primitives (lines 507–627)

These three functions provide the shared machinery used by both
`_apply_codon_correction` and `_apply_discover_roi`:

- `_build_fragment_allele_table` (line 507): Single BAM pass. Records the
  base each DNA fragment (read pair merged by query name) carries at each
  requested 0-based amplicon-local position. Tracks three data structures:
  `pos_table = {pos: {qname: base}}`, `reach = {pos: (min_frag_start,
  max_frag_end)}` (the reference span any fragment touching `pos` can reach,
  used to skip SNP pairs no fragment could bridge), and `masked = {pos:
  {qname, …}}` (fragments with `'N'` at `pos`, e.g. from primer masking).
  `masked` and `pos_table` are kept strictly disjoint: if either mate provides
  a real base call at a position, the fragment is treated as called, not masked.
- `_tally_allele_linkage` (line 589): Given a `pos_table` and a set of
  positions, counts how often each allele combination co-occurs across all
  those positions on the same fragment. Only fragments present at ALL requested
  positions are counted. Returns a `Counter` mapping
  `frozenset{(pos, base), …}` → fragment count. Iterates the shallowest
  position's table to minimize work.
- `_amp_to_translated` / `_translated_to_amp` (lines 622, 629): Convert
  between 0-based amplicon-local positions and gene-relative translated
  (1-based) positions, accounting for negative offsets (amplicons that start
  before the gene boundary).
- `_snp_variant_freq` (line 636): Returns `(freq, count)` of the variant
  allele for a SNP dict entry.

#### Codon-aware linkage (lines 629–767, 931–948)

- `_apply_codon_correction` (line 646): Annotates pairs of SNPs that fall in
  the same codon (per `CdsFeature.codon_boundaries`) with read-level
  allele-linkage information. For each codon containing **exactly two**
  polymorphic positions (codons with ≥3 SNPs are skipped), it:
  1. Calls `_tally_allele_linkage` over the two amplicon-local positions.
  2. Skips if the dominant combo has fewer than `min_reads` fragments.
  3. Classifies `linkage` as `"complete"` if `|freq_A − freq_B|
     ≤ error_threshold`, else `"partial"`.
  4. Labels each observed allele combination as `"reference"` (both reference
     alleles), `"variant"` (both called variant alleles), or `"discordant"`
     (anything else).
  5. Appends a `codon_merges` list entry to **both** SNP dicts in place.
  Codon boundaries are deduplicated across all `CdsFeature` objects before
  processing, so overlapping CDS features (e.g. SARS-CoV-2 ORF1ab/ORF1a) do
  not produce duplicate annotations.

- `_add_codon_merges_node` (line 957): Serializes the `codon_merges` list to
  `<codon_merge>` XML children of a `<snp>` element; each
  `<codon_merge>` gets one `<combo>` child per observed allele combination.

#### Discover-ROI (lines 769–928)

- `_apply_discover_roi` (line 795): For each SNP in the list (filtered to
  those with variant frequency ≥ `min_snp_perc`), scans all other SNPs on
  the same amplicon and measures read-level co-occurrence. Skips SNP pairs
  that are already annotated via `codon_merges` (those are handled with the
  richer `<codon_merge>` annotation). For each qualifying pair, splits the
  anchor SNP's variant-supporting reads into four buckets: `linked`
  (co-carries the other SNP's variant), `standalone` (spans both positions
  but carries a different allele), `masked` (`'N'` at the other SNP's
  position), `non-overlapping` (read doesn't reach the other position). Three
  distinct percentages are reported: `linked_pct` = linked / all-variant-reads,
  `percentage_linked` = linked / comparable-reads (linked + standalone), and
  `spanning_depth` = all fragments spanning both positions regardless of allele.
  A pair is only reported if `linked_count ≥ min_reads` and
  `linked_pct / 100 ≥ min_perc`.

- `_add_linked_snps_node` (line 935): Serializes the `linked_snps` list to a
  `<linked_snps>` XML child containing one `<linked_snp>` element per entry
  (attributes: `name`, `ref_pos`, `spanning_depth`, `linked_depth`,
  `comparable_depth`, `percentage_linked`, `snp_percentage_linked`), plus a
  `<full_distribution>` child text node with the four-bucket breakdown.

#### Output assembly (lines 1130–1616)

`main()` (line 1156) is the CLI entry point. Key steps:

1. Opens BAM via pysam; extracts sample name from the `RG` read group header
   (falls back to BAM filename stem). Computes `mapped_reads`,
   `unmapped_reads`, and `unassigned_reads` (lines ~782–798) from a single,
   consistent BAM snapshot — `--original-bam` (pre-ASAP-filter, straight off
   the aligner) when provided, else the `-b` input BAM. Do not compute these
   three from different BAMs (e.g. one pre-filter, one post-identity-filter)
   — `identityFilter.py` re-flags failing reads as unmapped rather than
   deleting them (see §4 `identityFilter.py`), so mixing snapshots
   double-attributes those reads. `mapped_reads` = primary-mapped count
   (`samtools flagstat` "primary mapped", excludes secondary/supplementary).
   `unmapped_reads` = pysam `.unmapped`, the *total* unmapped-read count.
   `unassigned_reads` = pysam `.nocoordinate`, the subset of unmapped reads
   with no alignment coordinate at all (both mates failed to align) —
   already included inside `unmapped_reads`, not a separate pool. See §6.3
   for how the R-side Alignment Summary plot accounts for this overlap.
2. For each amplicon: runs pileup → calls `_process_pileup`; if
   `codon_correction` or `discover_roi` is enabled and there are ≥2 SNPs with
   depth > 0, does a single BAM pass via `_build_fragment_allele_table`
   (shared by both features); applies `_apply_codon_correction` and/or
   `_apply_discover_roi`; serializes SNPs with `_add_snp_node` (then appends
   `<codon_merge>` and/or `<linked_snps>` children as needed).
3. Evaluates `Operation` objects across the assembled XML (cross-amplicon
   significance logic) and writes `<operation>` elements.
4. Calls `_write_output` (line 1533) to serialize as pretty-printed XML or JSON.

`cast_json_output_types` (line 1561) is the JSON `object_hook` that converts
string-valued XML attributes to typed Python values after xmltodict parsing.
See Section 5 for the full type table.

---

## 5. Per-Sample Output Schema (XML / JSON)

`newBamProcessor.py` produces one XML file per sample. When
`--output-format json` is used, the XML is converted via xmltodict +
`cast_json_output_types`.

### 5.1 XML Element Hierarchy

`mapped_reads`/`unmapped_reads`/`unassigned_reads` are all computed from one
BAM snapshot (see §3.3 step 1). `unassigned_reads` (no alignment coordinate
at all) is a *subset* of `unmapped_reads` (total unmapped count) — not an
additional, disjoint category. To get a mutually-exclusive "truly unmapped
but not already counted as unassigned" figure, compute
`unmapped_reads − unassigned_reads`.

```
<sample name="…" mapped_reads="…" unmapped_reads="…" unassigned_reads="…"
        depth_filter="…" breadth_filter="…" proportion_filter="…"
        mutation_depth_filter="…" json_file="…" bam_file="…"
        [total_reads="…"] [trimmed_reads="…"] [SMOR="True"]>

  <assay name="…" type="…" function="…" gene="…" start="…" end="…">

    <amplicon reads="…" [variant="…"] [aligned_reads="…"]
              [primer_reads="…"] [no_primer_reads="…"]
              [identity_input="…"] [identity_discarded="…"]
              [smor_input="…"] [smor_pairs_dropped="…"]
              [smor_consensus_reads="…"] [smor_singleton_reads="…"]>

      <significance [flag="no coverage|low coverage|insufficient breadth of coverage"]
                    [resistance="…"]>…text…</significance>       <!-- optional -->

      <snp name="…" position="…" depth="…" reference="…">
        <snp_call count="…" percent="…">VARIANT_BASE</snp_call>
        <significance [flag="…"] [resistance="…"] [level="low|high"]>…</significance>  <!-- optional -->
        <base_distribution A="…" T="…" C="…" G="…" [_="…"] […]/>  <!-- optional -->

        <!-- present if --codon-correction matched this SNP to a same-codon partner -->
        <codon_merge linked_snp="OTHER_SNP_NAME" spanning_depth="…"
                     linkage="complete|partial"
                     snp_percentage_linked="…" total_percentage_depth_linked="…">
          <combo bases="X|Y" count="…" percent="…" type="reference|variant|discordant"/>
          …  <!-- one <combo> per observed allele combination -->
        </codon_merge>
        …  <!-- one <codon_merge> per codon containing this SNP (usually 1) -->

        <!-- present if --discover-roi found co-occurring variants -->
        <linked_snps>
          <linked_snp name="…" ref_pos="…" spanning_depth="…"
                      linked_depth="…" comparable_depth="…"
                      percentage_linked="…" snp_percentage_linked="…">
            <full_distribution>linked …% (N=…), standalone …% (N=…), …</full_distribution>
          </linked_snp>
          …
        </linked_snps>
      </snp>
      …  <!-- one <snp> per SNP/iSNV per amplicon -->

      <breadth>…</breadth>
      <average_depth>…</average_depth>
      <consensus_sequence>…</consensus_sequence>         <!-- omitted if --whole-genome -->
      <gapfilled_consensus_sequence>…</gapfilled_consensus_sequence>
      <depths>pos1,pos2,…</depths>                       <!-- comma-separated -->
      <proportions>pos1,pos2,…</proportions>
      <quality_discards>pos1,pos2,…</quality_discards>
      <n_reads>pos1,pos2,…</n_reads>
      <ref_positions>pos1,pos2,…</ref_positions>
    </amplicon>
    …
  </assay>
  …

  <operation flag="…" message="…"/>  <!-- one per triggered cross-amplicon logic rule -->
</sample>
```

### 5.2 JSON Type Casting (`cast_json_output_types`)

`cast_json_output_types` is called as an `object_hook` for every JSON object
(dict) produced by xmltodict. Attribute names are prefixed with `@`, child
text with `#text`. The function applies type conversions and list-wrapping for
consistency.

| Element | Key | Cast |
|---|---|---|
| `<sample>` | `@breadth_filter` | `float` |
| `<sample>` | `@depth_filter` | `int` |
| `<sample>` | `@mapped_reads` | `int` |
| `<sample>` | `@proportion_filter` | `float` |
| `<sample>` | `@unassigned_reads` | `int` |
| `<sample>` | `@unmapped_reads` | `int` |
| `<sample>` | `assay` | always a `list` (wraps single-assay result) |
| `<assay>` | `amplicon` | always a `list` |
| `<amplicon>` | `@reads` | `int` |
| `<amplicon>` | `breadth` | `float` |
| `<amplicon>` | `average_depth` | `float` |
| `<amplicon>` | `depths` | `list[int]` (split on `,`) |
| `<amplicon>` | `proportions` | `list[float]` (split on `,`) |
| `<amplicon>` | `snp` | always a `list` |
| `<snp>` | `@depth` | `int` |
| `<snp>` | `@position` | `int` |
| `<snp>` | `base_distribution` | `dict[str, int]` (values cast to `int`) |
| `<codon_merge>` | `@spanning_depth` | `int` (generic — matches any element) |
| `<codon_merge>` | `codon_merge` | always a `list` |
| `<combo>` | `@count` | `int` (generic — matches any element with `@count`) |
| `<combo>` | `@percent` | `float` (generic — matches any element with `@percent`) |
| `<snp_call>` | `@count` | `int` (same generic rule) |
| `<snp_call>` | `@percent` | `float` (same generic rule) |
| `<significance>` | `@resistance` | `list[str]` (split on `,`) |
| `<significance>` | (if string) | wrapped in `{'#text': value}` |
| any | `@changes` | `int` |
| any | `mutation` | always a `list` |

**Known gaps** (current behavior, not bugs to fix here):

- `linked_snp` (the `<linked_snps>` container) is **not** wrapped in a list.
  If a SNP has exactly one `<linked_snp>` child, xmltodict returns it as a
  dict rather than a one-element list in JSON output.
- `<linked_snp>` numeric attributes (`@linked_depth`, `@comparable_depth`,
  `@percentage_linked`, `@snp_percentage_linked`) are **not** cast to numbers.
  They remain strings in JSON. Note that `@spanning_depth` **is** cast to
  `int` (it matches the generic `@spanning_depth` check that was added for
  `<codon_merge>`) — creating a subtle inconsistency between attributes on the
  same element.

### 5.3 Worked Example — `<codon_merge>`

From a SARS-CoV-2 run with `--codon-correction` enabled: positions T5118A and
T5119A both fall in the same ORF1ab codon. The `<snp>` element for T5118A
includes:

```xml
<snp name="T5118A" position="5118" depth="1727" reference="T">
  <snp_call count="97" percent="5.619...">A</snp_call>
  <base_distribution T="1630" A="97"/>
  <codon_merge linked_snp="T5119A"
               spanning_depth="1718"
               linkage="complete"
               snp_percentage_linked="5.6"
               total_percentage_depth_linked="5.6">
    <combo bases="T|T" count="1621" percent="94.4" type="reference"/>
    <combo bases="A|A" count="97"   percent="5.6"  type="variant"/>
  </codon_merge>
</snp>
```

Interpretation:

- `spanning_depth=1718`: 1718 fragments had confident base calls at **both**
  positions 5118 and 5119.
- `linkage="complete"`: `|freq(T5118A) − freq(T5119A)|` ≤ `error_threshold`
  (0.05 default) — both SNPs change together at the same frequency.
- `combo type="reference"` (T|T, 94.4%): the vast majority of fragments carry
  both reference alleles — wild-type codon.
- `combo type="variant"` (A|A, 5.6%): fragments carrying both variant alleles —
  the actual codon change.
- `snp_percentage_linked=5.6`: variant_combo_count / spanning_depth × 100.
- `total_percentage_depth_linked=5.6`: variant_combo_count / this SNP's own
  depth (1727) × 100 ≈ same here because depth ≈ spanning_depth.
- No `"discordant"` combo: every fragment either carries both variants or
  neither — the mutation is perfectly co-occurring.

T5119A carries a symmetric `<codon_merge>` entry pointing back to T5118A with
the same `spanning_depth`, `linkage`, `snp_percentage_linked`, and `combos`.

In JSON (after `cast_json_output_types`):

```json
{
  "@name": "T5118A", "@position": 5118, "@depth": 1727, "@reference": "T",
  "snp_call": {"@count": 97, "@percent": 5.619, "#text": "A"},
  "codon_merge": [{
    "@linked_snp": "T5119A",
    "@spanning_depth": 1718,
    "@linkage": "complete",
    "@snp_percentage_linked": "5.6",
    "@total_percentage_depth_linked": "5.6",
    "combo": [
      {"@bases": "T|T", "@count": 1621, "@percent": 94.4, "@type": "reference"},
      {"@bases": "A|A", "@count": 97,   "@percent": 5.6,  "@type": "variant"}
    ]
  }]
}
```

Note `codon_merge` is a list; `@spanning_depth` is `int`; `@count`/`@percent`
on each `combo` are `int`/`float`.

### 5.4 Worked Example — `<linked_snps>`

From a synthetic test case illustrating the `--discover-roi` four-bucket
breakdown. Anchor SNP s_a is A→T at position 101 (10 variant reads out of 15
total depth); candidate SNP s_b is A→C at position 111 (14 total depth):

- 6 reads: T@101 + C@111 → **linked**
- 3 reads: T@101 + A@111 → **standalone** (confident call, different allele)
- 1 read (too short): T@101, never reaches 111 → **non-overlapping**
- 5 reads: A@101 + A@111 → reference, not s_a's variant (span both positions)

Result: `spanning_depth = 14` (6+3+5 fragments with confident calls at both
positions); `comparable_count = 9` (6+3 reads carrying s_a's variant and
reaching 111); `percentage_linked = 66.7` (6/9×100); `snp_percentage_linked =
60.0` (6/10×100, denominator is all of s_a's variant reads including the
short one).

```xml
<snp name="s_a" position="101" depth="15" reference="A">
  <snp_call count="10" percent="66.67">T</snp_call>
  <linked_snps>
    <linked_snp name="s_b" ref_pos="A111"
                spanning_depth="14" linked_depth="6"
                comparable_depth="9" percentage_linked="66.7"
                snp_percentage_linked="60.0">
      <full_distribution>
        linked 60.0% (N=6), standalone 30.0% (N=3),
        masked 0.0% (N=0), non-overlapping 10.0% (N=1)
      </full_distribution>
    </linked_snp>
  </linked_snps>
</snp>
```

In JSON: `@spanning_depth` is cast to `int` (generic rule); `@linked_depth`,
`@comparable_depth`, `@percentage_linked`, `@snp_percentage_linked` remain
**strings** (not yet cast — known gap, see §5.2).

---

## 6. R `asaptools` Pipeline (`asap_tools/`)

### 6.1 Data Flow

```
process_xml.R (per sample, parallel)
  sources: _read.ASAP.individual.R
           _read.ASAP.snps.individual.R
           _ASAP.get.{depth,nreads,proportions,quality_discards}.R
  outputs: <id>_XML_Data.Rdata   (ASAP, SNPS, array_info)
           <id>_Summary.csv

process_combine_rdata.R (gather all samples)
  outputs: <name>_ASAP_Data.Rdata   (final_asap, final_snps, final_array)
           <name>_Summary.csv
      │
      ├─→ process_asaptools_cov_table.R
      │     outputs: <prefix>_Coverage_Report.xlsx
      │
      ├─→ process_asaptools_generate_figures.R
      │     outputs: <prefix>_{coverage_depth,n_reads_prop,SNP_prop}.{jpg,html}
      │
      ├─→ process_asaptools_fasta_export.R
      │     outputs: <prefix>_<assay>.fasta (one per assay)
      │
      ├─→ process_asaptools_snps_amino_acids.R
      │     sources: _genome.snp.to.gene.snp.R
      │              _snps.to.amino.R
      │              _expand_codon_merges.R
      │     outputs: SNP_Amino_Acid_Table.Rdata  (Amino_Acids, Gene_SNPS)
      │
      └─→ process_asaptools_snp_table.R
            sources: _expand_codon_merges.R
            inputs:  SNP_Amino_Acid_Table.Rdata (from above)
            outputs: <prefix>_SNP_Table_{Included,All}_Samples.csv
                     <prefix>_SNP_Linelist_{Included,All}_Samples.csv
                     <prefix>_SNP_Table_Final.xlsx  (if --xls)
```

**Note — `process_primers_to_bed.R` is separate from this flow.** It is a
pre-alignment utility (Nextflow process `GENERATE_PRIMER_BED`), not part of the
combined-Rdata post-processing chain above. Unlike the other scripts it is
**self-contained** (it inlines `find.primers` / `find.diff.in.seq` / a 6-column
`create.asap.bed.file` rather than `source()`-ing helpers), so it runs identically
standalone or from `bin/`. It reads a primer CSV (`primer_name, direction,
sequence`), searches the reference with `Biostrings::matchPattern` (both strands,
`--primer_max_mismatch` mismatches/indels), and writes a headerless 6-column BED
(0-based start; strand `F→+`, `R→-`), a full search-results CSV, and a wide
`*_primer_match_summary.csv` (one row per primer; one match-count column per reference
sequence plus `total_matches`). The primer loop
is parallelized over primers via `foreach %dopar%` with forked workers
(`registerDoParallel(task.cpus)`); the reference is shared copy-on-write. Dependencies
(`tidyverse`, `Biostrings`, `foreach`, `data.table`, `doParallel`) are all already in
`modules/asap_tools/r_env.yml`.

CLI: `process_primers_to_bed.R <primer_csv> <reference_fasta> <prefix> [max_mismatch=2] [cores=1]`

### 6.2 Per-Sample XML Parsing

**`_read.ASAP.individual.R`** parses sample-level and amplicon-level attributes
into a tall data frame (one row per amplicon). Columns include: sample metadata
(`run`, `name`, `total_reads`, `trimmed_reads`, `mapped_reads`, …), assay
metadata (`assay_name`, `assay_type`, `assay_function`, `assay_gene`), amplicon
metrics (`breadth`, `avg_depth`, `amplicon_reads`, `aligned_reads`, …), and the
serialized position arrays as comma-delimited strings (`depths`, `proportions`,
`quality_discards`, `n_reads`, `consensus_seq`).

**`_read.ASAP.snps.individual.R`** parses SNP-level data into one row per SNP.
The critical codon-merge transformation flattens the XML hierarchy: for each
`<snp>` element, all `<codon_merge>` children are collected and their
attributes vectorized into six semicolon-delimited string columns appended to
the SNP row:

| Column | Source XML attribute |
|---|---|
| `codon_merge_linked_snp` | `<codon_merge linked_snp="…">` |
| `codon_merge_linkage` | `<codon_merge linkage="…">` |
| `codon_merge_spanning_depth` | `<codon_merge spanning_depth="…">` |
| `codon_merge_variant_bases` | `<combo type="variant" bases="…">` |
| `codon_merge_variant_count` | `<combo type="variant" count="…">` |
| `codon_merge_variant_percent` | `<combo type="variant" percent="…">` |

If a SNP links to N codon merges, each column holds N semicolon-separated
values. If no `<codon_merge>` elements are present, all six columns are
`NA_character_`. `<linked_snps>` XML data is **not** currently parsed into
the R data frame — that information is available only in the XML/JSON output.

**Small helper functions** (`_ASAP.get.{depth,nreads,proportions,quality_discards}.R`):
each unpacks a comma-delimited position array from the ASAP data frame into a
tall `(run, name, assay_name, position, value)` data frame, using
`foreach`/`doParallel` for parallelization.

**`process_xml.R`**: Per-sample Nextflow entry point. Invoked as
`Rscript process_xml.R <xml_file> <min_snp> <sample_id>`. Calls both
reader functions, applies numeric coercions, filters SNPs to those with
`snp_proportion ≥ min_snp`, unpacks position arrays via the four helper
functions, and saves `{sample_id}_XML_Data.Rdata` + `{sample_id}_Summary.csv`.

### 6.3 Combine and Downstream Scripts

**`process_combine_rdata.R`**: Invoked as
`Rscript process_combine_rdata.R <poi_csv|NULL> <file_name> <rdata1> …`.
Uses `foreach`/`doParallel` (SLURM-aware core detection via
`parallelly::availableCores()`) to load each Rdata in parallel, then
merges with `data.table::rbindlist()`. If a POI CSV is provided, filters
`array_info` to positions within gene intervals.

**`process_asaptools_cov_table.R`**: Coverage as percentage of positions
above a depth threshold, pivoted wide by sample, styled as a 10-bin
color gradient (red→green) in Excel via openxlsx.

**`process_asaptools_generate_figures.R`**: Five ggplot2 plots — coverage
depth, N-read proportion, breadth-of-coverage heatmap, alignment summary,
and read funnel — each exported as static JPG and interactive HTML via
plotly/ggplotly. Rolling-mean downsampling targets ≤10,000 points per facet
for manageable file size.

The **Alignment Summary** plot (counts + percentage panels) stacks
`mapped_reads` ("Aligned"), `unassigned_reads` ("Unassigned"),
`unmapped_reads − unassigned_reads` ("Unmapped"), and `total_reads −
trimmed_reads` ("Removed by FastP"). The subtraction matters: `unmapped_reads`
(pysam `.unmapped`) is the *total* unmapped-read count in the BAM, and
`unassigned_reads` (pysam `.nocoordinate`) is the subset of those with no
alignment coordinate at all (neither mate aligned anywhere) — so
`unassigned_reads` is already included inside `unmapped_reads`, not a
separate pool. Stacking both raw would double-count nearly the entire
non-aligned read pool. The remainder, `unmapped_reads − unassigned_reads`,
is just the reads flagged unmapped but placed at a coordinate because their
mate *did* align (mate-rescued placement) — usually a small sliver. See
§5.1 for where `mapped_reads`/`unmapped_reads`/`unassigned_reads` are
computed.

**`process_asaptools_fasta_export.R`**: Exports `consensus_seq` from
`final_asap` as one FASTA per assay, filtering to samples meeting the
`breadth_threshold`. Vectorized string construction avoids slow row-by-row
loops.

**`process_asaptools_snps_amino_acids.R`**: Invoked as
`Rscript process_asaptools_snps_amino_acids.R <rdata> <ref1> [<ref2> …]`.
Sources `_expand_codon_merges.R`, `_genome.snp.to.gene.snp.R`, and
`_snps.to.amino.R`. Calls `expand_codon_merges(SNPS)` first so that
combined codon-merge rows (e.g. `"T5118A|T5119A"`) are present in the SNP
table before amino acid translation; saves
`SNP_Amino_Acid_Table.Rdata` (contains `Amino_Acids` and `Gene_SNPS` data
frames, keyed by `(assay_name, SNP)`).

**`process_asaptools_snp_table.R`**: The most complex script (327 lines).
Also independently calls `expand_codon_merges(SNPS)` (see §6.4 below), then
left-joins `Amino_Acids` and `Gene_SNPS` by `(assay_name, SNP)`. Produces
wide-format and linelist CSV outputs with per-position coverage, SNP
proportion, gene annotation, amino acid change, and optional primer-region
annotation. Sample QC filtering by `--max_snp_count` and `--remove_names`.

### 6.4 Codon-Merge Expansion and AA Translation Walkthrough

This is the least-obvious part of the asaptools pipeline. Three components
must be understood together:

**`_expand_codon_merges.R`** (75 lines): Takes the SNPS data frame with
`codon_merge_*` semicolon-delimited columns and:

1. For SNPs where every `codon_merge_linkage` entry for that SNP is
   `"complete"`, marks the individual SNP row for removal (`drop_keys`).
2. Constructs a new combined row per pair — e.g., T5118A + T5119A → a new
   row with `SNP = "T5118A|T5119A"` in canonical (low-position|high-position)
   order, using the `bases` field from `<combo type="variant">` to get the
   mutation alleles, and `count`/`percent`/`spanning_depth` from the variant
   combo.
3. Each pair appears twice in the source (once from each SNP's perspective);
   `distinct(pair_key)` keeps one.
4. Returns SNPS with individual complete-pair rows replaced by combined rows;
   for `"partial"` pairs, keeps both the individual rows AND adds the combined
   row.

**`_genome.snp.to.gene.snp.R`** and **`_snps.to.amino.R`** both support
pipe-joined (multi-token) SNP strings. A SNP string is split on `"|"`, each
token parsed as `{ref}{position}{mut}`, and all mutations are applied to the
same `Observed_Seq` before a single `Biostrings::translate()` call. This
means a combined row `"T5118A|T5119A"` produces one amino acid annotation for
the **joint codon change**, not two separate single-nucleotide annotations —
which is the correct behavior.

**Critical join invariant**: Both `process_asaptools_snps_amino_acids.R` and
`process_asaptools_snp_table.R` call `expand_codon_merges(SNPS)` independently
on the same input. The SNP strings they produce must be identical for the
`left_join(by = c("assay_name", "SNP"))` in `process_asaptools_snp_table.R` to
correctly attach amino acid annotations. If the two scripts diverge in how they
expand codon merges (e.g. by using different sort orders for the pair), the
join will silently drop AA annotations for combined rows. Always verify both
scripts source the same `_expand_codon_merges.R` function with the same input.

**Fallback behavior**: When a multi-token SNP includes an insertion (base
length > 1) or deletion (`_`), `snps.to.amino()` returns
`"Insertions Not Supported"` or `"Deletions Not Supported"` rather than
attempting translation.

---

## 7. Testing

### 7.1 Python Unit Tests (`pytest`, `asap_test` conda env)

Run from the repo root:

```bash
conda run -n asap_test pytest tests/ -v
```

| File | Coverage |
|---|---|
| `test_codon_correction.py` | `_apply_codon_correction` and `_add_codon_merges_node`; also `_apply_discover_roi`'s exclusion of codon-merge partners from `linked_snps` |
| `test_discover_roi.py` | `_apply_discover_roi` and `_add_linked_snps_node`; focuses on the `discover_roi_min_snp_perc` pre-filter and the four-bucket denominator semantics (`spanning_depth` vs. `comparable_count` vs. `count_a`) |
| `test_genbank_cds.py` | `_parse_genbank_cds`: forward/reverse-strand CDS, mid-codon frame skipping, multi-file input, LRU cache, and an integration test against the real H37Rv GenBank fixture if available |
| `test_identity_filter.py` | `identityFilter._identity_filter` with both `filter_pairs=True` and `False`; also verifies scope of filtering to specified references only |
| `test_linkage.py` | `_tally_allele_linkage` and `_build_fragment_allele_table` using synthetic BAMs: same-read linkage, cross-mate linkage, unlinked SNPs, non-spanning reads, supplementary/secondary exclusion, and `min_reads` gating |
| `test_pileup_naming.py` | De-novo SNP naming in `_process_pileup`: nucleotide-notation naming and position-0 wildcard override |
| `test_smor_consensus.py` | `generateSMORbam._write_bam` and `generateSMORbam_correction._write_bam`: read accounting for matched pairs, singletons, and non-overlapping pairs |

### 7.2 Nextflow End-to-End Tests (nf-test, SLURM)

End-to-end tests live in `nextflow/tests/ASAP_EtE.nf.test`. See README.rst
(Running the Test Suite) for the full `nf-test` workflow, tag list, and how
to run individual test cases or the full parallel suite via
`run_tests_parallel.sh`. The test suite uses snapshot-based assertion for
stable output values. Run `--update-snapshot` on the first run (or after
intentional output changes) to establish the baseline.

---

## 8. Extending ASAP

The `<codon_merge>` feature is the template for adding a new XML annotation
element with a matching Nextflow parameter and R integration. The pattern is:

### Python side (`newBamProcessor.py`)

1. **Compute**: add a function `_apply_<feature>(snp_list, …)` that mutates
   SNP dicts in place (e.g. `snp['my_data'] = …`).
2. **Annotate**: add `_add_<feature>_node(snp_node, data)` that serializes
   the dict to XML children of `<snp>`.
3. **Wire up in `main()`**: add argparse argument, call `_apply_<feature>`
   after `_build_fragment_allele_table` (if it needs read-level data), and
   call `_add_<feature>_node` in the SNP serialization loop.
4. **Cast types**: add an `if '@my_attr' in e` block in
   `cast_json_output_types` for each JSON attribute that should be a number
   or list. Add a list-wrapping rule if the element can appear multiple times
   under a parent.

### Schema and config side

5. Add a parameter entry to `nextflow/nextflow_schema.json`
   (`bam_processor_options` group) with a precise description.
6. Add a row to the Step 6 parameter table in `README.rst` (regenerate with
   a script matching the existing column widths: `W1=34`, `W2=10`,
   `W3=60`).
7. Add the `--my-param` CLI argument to the `PROCESS_BAM` process command
   in `nextflow/modules/asap/main.nf`.

### R side (`asap_tools/`)

8. **XML parsing**: update `_read.ASAP.snps.individual.R` to extract the new
   XML attributes into columns of the SNPS data frame.
9. **Expansion helper** (if the new element introduces combined-row logic):
   add a `_expand_<feature>.R` function following the pattern of
   `_expand_codon_merges.R`, and source it in **both**
   `process_asaptools_snps_amino_acids.R` and `process_asaptools_snp_table.R`
   to keep the SNP join key consistent.
10. **SNP table**: update `process_asaptools_snp_table.R` to include the new
    columns in the wide and linelist outputs.

### Tests

11. Add a test file in `tests/` following the pattern of
    `test_codon_correction.py`: construct `pos_table`/SNP dicts by hand (no
    BAM needed for pure annotation logic), cover the main classification
    branches, and verify XML serialization via `_add_<feature>_node`.
12. Add a tagged nf-test case in `nextflow/tests/ASAP_EtE.nf.test` that
    exercises the new parameter end-to-end with a snapshot assertion.

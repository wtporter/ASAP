# Changelog

## [Unreleased] — `development` (2026-07-09)

SMOR consensus-correction parameters are now namespaced, and the way agreeing
base qualities are combined is configurable.

### Changed

- **Renamed `qual_diff_threshold` → `smor_correction_qual_diff_threshold`.** This
  parameter is consumed **only** by `SMOR_CORRECTION`
  (`generateSMORbam_correction.py -q`), so it now carries the `smor_correction_`
  prefix and lives under the **SMOR Processing Options** group. **Breaking:**
  update any command line / config that passes `--qual_diff_threshold`. The
  default (`10`) and the script's `-q` flag are unchanged.
- **Renamed `asap/newBamProcessor.py` → `asap/ASAPBamProcessor.py`** (and the
  `nextflow/bin/` symlink) for a consistent, self-descriptive name. Internal
  references updated: the `PROCESS_BAM` process invocation, the `asap.*` imports in
  the Python tests, and the docs. No behavioral change. The unintegrated
  `fasterBamProcessor.py` sibling is untouched.
- **Renamed the `whole_genome` parameter → `suppress_per_base`** (CLI
  `--suppress-per-base`). The flag never described the reference — it suppresses the
  per-position output arrays (`consensus`, `depths`, `proportions`, `n_reads`,
  `ref_positions`) so large/whole-genome references don't blow up XML size. The new
  name says what it does. **Breaking:** update any config/CLI passing `whole_genome`.

### Added

- **`smor_correction_agreement_method`** (`sum` | `max`, default `sum`) — controls
  how the two base qualities are combined when both reads of a pair **agree** at a
  position during SMOR correction (`generateSMORbam_correction.py -a`):
  `sum` = `min(q1+q2, 60)` (the previous fixed behavior), `max` = `max(q1, q2)`.
  The default preserves existing output.
- **`prune_per_base`** (`--prune-per-base`, default `false`) — a hybrid alternative
  to `suppress_per_base` for large/whole-genome references. Instead of dropping
  **all** per-base data, it retains the per-position numeric arrays (`depths`,
  `proportions`, `n_reads`, `quality_discards`) **only at positions with depth ≥
  `depth`**, and emits a matching sparse `ref_positions`; the `consensus_sequence`
  stays full-length so FASTA export is unaffected. This keeps the per-base data that
  asaptools coverage tables and QC/SNP figures need while shrinking output on
  references where most positions have ~zero coverage. Coexists with
  `suppress_per_base` (which still drops everything); the threshold reuses `depth`.
  **asaptools coordinate-awareness:** the R array extractors now derive genomic
  `position` from `ref_positions` (with a contiguous fallback) instead of assuming
  arrays start at position 1 — a correctness fix for offset assays as well — and the
  coverage table now uses the **exact reference length** (`nchar(consensus_sequence)`
  whole-reference; gene-range length for POI) as the denominator so pruned runs
  report true coverage instead of a constant 100%.

### Fixed

- Corrected an earlier changelog mislabel: `qual_diff_threshold` (now
  `smor_correction_qual_diff_threshold`) was listed under *Identity filtering*.
  It has always driven **SMOR consensus correction** (`generateSMORbam_correction.py`),
  not `identityFilter.py`.

## [Unreleased] — `ROI_Development` → `public` (2026-07-09)

**Scope:** 27 commits · 79 files · **+10,944 / −1,564** · 2026‑06‑09 → 2026‑07‑09

A large, multi-subsystem advance of the ASAP amplicon pipeline. It adds
read-level variant discovery and codon-aware allele linkage to the core engine, a
full **intron-aware / multi-contig amino-acid calling** path (works on segmented
viruses *and* eukaryotic multi-contig references), a substantially expanded
figure/report suite, automated primer-BED generation, and a real test suite
(Python + R) where there previously was little. A new 990-line `DEVELOPER.md`
documents the architecture and per-sample output schema.

### Highlights

- **Codon-aware allele linkage & read-level variant discovery** — new
  `asap/allele_linkage.py`, `asap/genbank_cds.py`; per-codon distributions,
  linked-SNP reporting, in-frame indel AA translation, and optional
  region-of-interest (ROI) discovery.
- **Intron-aware, multi-contig amino-acid calling** — the SNP→gene→AA path now
  handles spliced (`join()`) CDS and multi-record GenBank files, enabling AA
  calls on eukaryotic references (validated on a real fungal PMA1 gene) while
  staying byte-for-byte identical on single-exon viral references.
- **Pair-aware identity filtering** and **configurable SMOR consensus** in the
  core BAM processing chain.
- **Expanded, self-sizing figures** — new SNP genome-track / prevalence / density
  / strand / quality panels, a fastp QC panel, automatic figure sizing, and a
  corrected read-fate funnel.
- **Automated primer BED from a primer CSV** (degenerate-aware alignment) plus
  primer-masking QC, including an explicit `removed_reads` accounting column.
- **New test suites** — 7 Python test modules and ~15 R (testthat) test files
  with committed fixtures.

### Core engine — `asap/` (Python)

| Area | What changed |
|---|---|
| **Allele linkage** | New `allele_linkage.py`: codon-aware linkage of co-occurring variants; `codon_partner_names` tracking; linked-SNP sets. |
| **Codon correction** | Full 3-base codon distributions (incl. partial/full deletions `_`); per-CDS codon entries even for overlaps; richer `codon_merges` (codon depth, reference/observed codon, call %, full distribution); `excl_has_n`/`excl_no_span` exclusion reporting. Params: `codon_correction`, `codon_correction_error`, `codon_correction_min_reads`. |
| **GenBank CDS parsing** | New `genbank_cds.py` — CDS/allele-linkage logic extracted to a dedicated, testable module. |
| **Identity filtering** | `identityFilter.py`: pair-aware percent-identity filtering. Params: `filter_pairs`. |
| **SMOR consensus** | `generateSMORbam*.py`: configurable consensus-correction quality threshold. Param: `qual_diff_threshold` (renamed to `smor_correction_qual_diff_threshold` in the current dev cycle). |
| **BAM processing** | `ASAPBamProcessor.py` (formerly `newBamProcessor.py`) reworked to emit the richer per-amplicon funnel + SNP metrics. `fasterBamProcessor.py` is an **experimental, unintegrated** performance rewrite — committed for reference only; it is not wired into the pipeline (nothing imports or invokes it) because it was not reliably faster. `ASAPBamProcessor.py` remains the active processor. |
| **SNP metrics** | Per-base quality stats (mean/median/min/max) and read-strand distribution (R1/R2/SE) now in the XML and parsed downstream. |
| **ROI discovery** | Optional read-level variant / region-of-interest discovery. Params: `discover_roi`, `discover_roi_min_perc`, `discover_roi_min_reads`, `discover_roi_min_snp_perc`. |
| **Primer masking** | `maskPrimers.py`: now also emits an explicit **`removed_reads`** column (reads actually dropped — non-zero only with `--primer-only`), so read-loss is tracked, not inferred. |

### Multi-contig, intron-aware amino-acid calling — `asap_tools/` (R)

The SNP→gene→AA path was reimplemented so it is correct for **spliced genes and
multi-contig references**, not just single-exon viral genomes.

- **`_split_genbank_records.R`** (new) — splits a multi-record GenBank file into
  per-contig single-LOCUS files (genbankr cannot read multi-record input).
- **`_extract_gene_table.R`** (rewritten) — groups exon rows by transcript into
  one spliced-CDS record (introns excised); adds an `exons` list-column,
  `codon_start`, and a `genomic_to_spliced_pos()` coordinate mapper.
- **`_snps.to.amino.R` / `_genome.snp.to.gene.snp.R`** — exon-aware genomic→CDS
  mapping (intronic positions → non-coding); minus-strand handling via spliced
  length; a `ref_df` parameter to parse each contig once; and a fix for a crash
  on all-intergenic/all-intronic batches.
- **`process_asaptools_snps_amino_acids.R`** — globs `.gbf`/`.gbff`, iterates
  split contigs, matches SNPs by exact `assay_name`, and only parses
  SNP-bearing contigs.

Validated on a real fungal PMA1 gene (6 exons → correct intron-corrected codon)
and confirmed **byte-for-byte identical** output to the prior code on real RSV
(viral) data.

### Downstream R tools & figures — `asap_tools/`

- **New figure scripts:** `process_asaptools_snp_figures.R` (SNP genome track,
  position prevalence, proportion density, strand bias, base quality),
  `process_fastp_panel.R` (trim QC).
- **Genome track:** multi-contig aware (per-contig tracks) with SNP tiles
  **colored by predicted effect** (AA change / frameshift / synonymous /
  non-coding) sourced from the amino-acid table.
- **Automatic figure sizing** across QC and SNP figures, scaled to sample/facet
  counts and capped under ggplot2's 50-inch `ggsave` limit.
- **Read-fate funnel fix:** previously double-counted primer-masked-but-retained
  reads as "lost." Now the funnel reflects true conservation
  (`aligned = removed + identity + smor + final + other`) and uses the new
  `primer_removed_reads` column for the "No Primer" category (0 unless
  `--primer_only`).
- **`asaptools_interactive_plots` toggle** (default `false`) — gates the
  expensive self-contained interactive HTML widgets; JPGs always produced.
- **Primer BED from CSV:** `process_primers_to_bed.R` aligns primer sequences
  (degenerate-aware, `primer_max_mismatch`) to the reference to produce the BED.
- Supporting refactors: `_expand_codon_merges.R`, `_parse_snp_distribution.R`,
  `_shorten_sample_names.R`, `combine_masked_reads_wide.R`.

### Testing

- **Python (`tests/`, new):** `test_codon_correction.py`, `test_linkage.py`,
  `test_genbank_cds.py`, `test_identity_filter.py`, `test_smor_consensus.py`,
  `test_discover_roi.py`, `test_pileup_naming.py`.
- **R (`asap_tools/tests/`, new):** testthat suite + committed fixtures
  (`tiny_with_locus_tag.gb`, `tiny_cds_only.gb`, `tiny_spliced_cds.gb`,
  `tiny_multicontig.gb`, `sample_asap.xml`) covering gene-table extraction,
  SNP→AA translation (incl. spliced-CDS + multi-contig target specs), the
  read/parse helpers, and a real-RSV integration test.
- **Runner:** `asap_tools/tests/run_r_tests.sh` auto-discovers a prebuilt conda
  env (local run by default, `--sbatch` to submit).

### Nextflow orchestration, config & infra

- `main.nf` / modules rewired: amino-acid conversion hoisted to a shared channel
  feeding both SNP plots and the SNP table; primer/QC stats propagated.
- **Environment management:** `env_dir` param drives `conda.cacheDir` /
  `singularity.cacheDir`; `ASAP_CONDA_ENV` / `R_CONDA_ENV` overrides;
  `prebuild_envs.sh` (new); conda env YAMLs.
- **Resource strategy:** per-label CPU/memory/time with attempt-scaled memory and
  `retry` (maxRetries 3).
- `run_tests.sh` / `run_tests_parallel.sh` updates; `nextflow_schema.json` params.

### Documentation

- **`DEVELOPER.md`** (new, ~990 lines) — architecture, module responsibilities,
  per-sample output schema.
- **`README.rst`** substantially updated (+583/−372).

### New parameters

| Parameter | Default | Purpose |
|---|---|---|
| `skip_fastqc` / `skip_multiqc` | `true` | Skip FastQC / MultiQC steps |
| `filter_pairs` | `true` | Pair-aware identity filtering |
| `qual_diff_threshold` | `10` | SMOR correction quality delta (renamed to `smor_correction_qual_diff_threshold` in the current dev cycle) |
| `codon_correction` | `false` | Enable codon correction |
| `codon_correction_error` | `0.05` | Codon error rate |
| `codon_correction_min_reads` | `10` | Min reads for codon correction |
| `discover_roi` | `false` | Read-level ROI discovery |
| `discover_roi_min_perc` / `_min_reads` / `_min_snp_perc` | `0.1` / `10` / `0.05` | ROI thresholds |
| `primer_max_mismatch` | `2` | Mismatches when building BED from primer CSV |
| `asaptools_snp_plots` | `true` | Generate SNP figures |
| `asaptools_interactive_plots` | `false` | Export interactive HTML widgets |
| `asaptools_max_sample_snp_count` | `10000` | SNP-table cap |
| `env_dir` | `null` | Base dir for conda/singularity caches |

*(`primer_only`, `mask_bam`, `wiggle` already existed; `primer_only` still defaults `false`.)*

### Backward compatibility & migration notes

- **Viral / single-contig behavior is unchanged** — AA calling verified identical
  on real RSV data.
- **New XML/RData fields** (`primer_removed_reads`, SNP quality/strand metrics)
  are additive; the funnel and readers **default gracefully** when a field is
  absent from older outputs, so pre-existing RData still renders.
- New features are **opt-in** (`codon_correction`, `discover_roi`,
  `asaptools_interactive_plots` all default off/false; SMOR/identity via flags).
- `primer_removed_reads` and the richer SNP metrics only populate on **fresh
  pipeline runs** (produced at the mask/BAM-processing stage).

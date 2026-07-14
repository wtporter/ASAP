..
   Running the ASAP test suite. For an overview and quick start, see ../README.rst

This page covers how to run the ASAP end-to-end (``nf-test``) and R unit-test suites,
the available test tags, and a description of each end-to-end test. For an overview and
quick start, see the `main README <../README.rst>`_. For pipeline internals and how to
add tests, see `DEVELOPER.md <../DEVELOPER.md>`_ (§7 Testing).

.. contents:: On this page
   :local:
   :depth: 1

Running the Test Suite
======================

End-to-end tests are defined in ``tests/ASAP_EtE.nf.test`` and executed via
``nf-test``. A SLURM submission script is provided for convenience:

.. code-block:: bash

   cd /path/to/ASAP/nextflow
   mkdir -p logs/nf-test        # must exist before sbatch

   # Run all tests
   sbatch run_tests.sh

   # First run — generate snapshots for stable outputs
   sbatch run_tests.sh --update-snapshot

   # Run a single data type
   sbatch --job-name=ASAP_tb run_tests.sh --tag tb --keep-going

   # Run the whole suite in parallel — one SLURM job per test (~11 jobs at once)
   ./run_tests_parallel.sh

Available test tags: ``help``, ``rsv``, ``tb``, ``sc2``, ``bwa``, ``bowtie2``,
``minimap2``, ``multi_gb``, ``json_input``, ``excel_input``, ``fasta_input``,
``paired_end``, ``single_end``, ``ont``, ``primer_masking``, ``primer_only``,
``identity_filter``, ``smor``, ``ivar``, ``asaptools``, ``combine_output``,
``suppress_per_base``, ``prune_per_base``, ``sc2_se_bwa``, ``tb_json_snp_aa``.

``run_tests_parallel.sh`` runs the suite roughly in the time of its single
longest test rather than the sum of all of them, by submitting one
``sbatch run_tests.sh --tag <tag>`` job per test. It uses a curated list of
tags that each match **exactly one** test (``sc2_se_bwa`` and
``tb_json_snp_aa`` were added for tests that previously had no uniquely
matching tag) — see the script's header comment for the full tag → test
mapping and the uniqueness rule to follow when adding new tests. Note that
the other tag-scoped examples above overlap (e.g. ``tb``/``sc2`` each match
several tests), so they should be run one at a time, not concurrently.

Test Descriptions
-----------------

+------+------------------------------------------------------+-----------------------------------------------------+
| No.  | Test name                                            | What it validates                                   |
+======+======================================================+=====================================================+
| 1    | Help Documentation                                   | ``--help`` flag renders schema-driven help text     |
+------+------------------------------------------------------+-----------------------------------------------------+
| 2    | RSV – Multi-GenBank, BWA, Paired-End, Combine Output | Multiple ``.gb`` files → JSON; BWA aligner;         |
|      |                                                      | paired-end Illumina reads; OUTPUT_COMBINER;         |
|      |                                                      | FORMAT_OUTPUT (HTML report); MultiQC aggregation    |
+------+------------------------------------------------------+-----------------------------------------------------+
| 3    | SC2 – JSON Direct Input, Bowtie2, Combine Output     | JSON passed without conversion; Bowtie2 aligner;    |
|      |                                                      | paired-end reads; combined XML and HTML report      |
+------+------------------------------------------------------+-----------------------------------------------------+
| 4    | SC2 – Single-End Illumina, BWA                       | Single-end read detection; SE fastp trimming;       |
|      |                                                      | BWA alignment; ASAP BAM processing; combine output  |
+------+------------------------------------------------------+-----------------------------------------------------+
| 5    | SC2 – ONT Long Reads, Minimap2                       | ``--technology ont``; fastplong trimmer;            |
|      |                                                      | Minimap2 alignment; ASAP BAM processing             |
+------+------------------------------------------------------+-----------------------------------------------------+
| 6    | SC2 – Primer Masking + iVAR Variant Calling          | ``maskPrimers.py`` BED-based masking; iVAR trim,    |
|      |                                                      | variant calling, and consensus generation           |
+------+------------------------------------------------------+-----------------------------------------------------+
| 7    | TB – Full Feature Amplicon                           | Single GenBank reference; BWA; primer masking;      |
|      |                                                      | identity filter; SMOR correction; full asaptools    |
|      |                                                      | suite (XML→R, combine, cov table, QC plots, FASTA   |
|      |                                                      | export, SNP→AA conversion, SNP table);              |
|      |                                                      | positions of interest; ``combine_output=false``     |
+------+------------------------------------------------------+-----------------------------------------------------+
| 8    | TB – Excel Reference Input                           | Excel → JSON conversion via                         |
|      |                                                      | ``prepareJSONInput_nextflow.py``; BWA alignment;    |
|      |                                                      | per-sample XML generation                           |
+------+------------------------------------------------------+-----------------------------------------------------+
| 9    | TB – FASTA Reference Input, Whole Genome Mode        | FASTA → JSON conversion; ``--suppress_per_base``    |
|      |                                                      | flag (omit all per-position arrays); per-sample     |
|      |                                                      | XML generation                                      |
+------+------------------------------------------------------+-----------------------------------------------------+
| 10   | TB – FASTA Reference Input, Prune Per-Base Mode      | FASTA → JSON conversion; ``--prune_per_base`` flag  |
|      |                                                      | (retain per-position arrays only where depth ≥      |
|      |                                                      | ``--depth``); per-sample XML generation             |
+------+------------------------------------------------------+-----------------------------------------------------+
| 11   | SC2 – Primer-Only BAM Filtering                      | ``--primer-only`` flag in ``maskPrimers.py``;       |
|      |                                                      | only reads overlapping a primer region are kept;    |
|      |                                                      | downstream ASAP processing continues from the       |
|      |                                                      | filtered BAM; masking stats file validated          |
+------+------------------------------------------------------+-----------------------------------------------------+
| 12   | TB – JSON Input, asaptools + GenBank SNP→AA          | JSON reference with ``asaptools_genbank_location``  |
|      |                                                      | supplying a separate ``.gb`` file for               |
|      |                                                      | ``PROCESS_SNPS_TO_AMINOACIDS``; full asaptools      |
|      |                                                      | suite; validates ``snp_reports/``, ``plots/``, and  |
|      |                                                      | ``SNP_Amino_Acid_Table.Rdata`` output               |
+------+------------------------------------------------------+-----------------------------------------------------+

R Unit Tests
------------

Independent of the ``nf-test`` end-to-end suite above, the R post-processing functions
in ``asap_tools/`` have a dedicated ``testthat`` unit-test suite:

.. code-block:: bash

   cd /path/to/ASAP
   ./asap_tools/tests/run_r_tests.sh

``run_r_tests.sh`` syncs a dedicated R conda environment from
``nextflow/modules/asap_tools/r_env.yml`` (creating or updating it under
``nextflow/work/r_test_env``) and submits a SLURM job running
``testthat::test_dir()`` over ``asap_tools/tests/testthat/``. Logs are written to
``asap_tools/tests/logs/asap_r_tests_<JOBID>.log``.

Current coverage includes the depth/read-count/proportion accessor functions
(``_ASAP.get.*``), codon-merge expansion, SNP distribution parsing, XML parsing
(``_read.ASAP.*``), and SNP-to-amino-acid translation (``_snps.to.amino``) — a unit-test
layer for the R functions themselves, complementary to the full-pipeline ``nf-test`` suite
above.


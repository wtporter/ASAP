Installation & Setup
====================

Full setup guide for ASAP: prerequisites, environment creation, execution profiles,
shared environment caches for multi-user clusters, and troubleshooting. For the short
version see `Quick Start <../README.rst>`_ in the README; for supported platforms and
pinned tool versions see the README's Requirements section.

.. contents:: On this page
   :local:
   :depth: 1

Prerequisites
-------------

- **Conda / Mamba** — builds the Nextflow control environment and each per-step tool
  environment. Mamba is strongly recommended for speed.
- **Singularity / Apptainer** — runs the containerized alignment and QC tools under the
  ``slurm`` profile.
- **SLURM** — required only for the ``slurm`` execution profile (HPC production runs).
- A **reference file** in FASTA, GenBank, Excel (``.xlsx``), or JSON format.

Everything else is provided by the environment YAML below: Nextflow 25.10.4 and nf-test
0.9.5 come from ``ASAP_nextflow_env.yml``, and the **nf-schema** plugin
(``nf-schema@2.5.1``) is declared in ``nextflow.config`` and downloaded automatically on
first run.

Step-by-step
------------

.. code-block:: bash

   # 1. Clone the repository
   git clone https://github.com/wtporter/ASAP.git
   cd ASAP

   # 2. Create and activate the Nextflow control environment
   conda env create -f ASAP_nextflow_env.yml      # Nextflow 25.10.4 + nf-test 0.9.5
   conda activate ASAP_nextflow_env

   # 3. Move into the pipeline directory and verify
   cd nextflow
   nextflow -version
   nextflow run main.nf --help                    # full, schema-driven parameter list

The ``ASAP_nextflow_env.yml`` file at the repository root creates a conda environment
named ``ASAP_nextflow_env``. Singularity containers and all per-step Conda environments
are resolved automatically the first time each step runs — no manual per-tool setup is
required.

Execution profiles
------------------

Select a profile with ``-profile`` (a single dash — it is a Nextflow core option, not an
ASAP ``--parameter``):

+---------------+---------------------+---------------------+-----------------------------------+
| Profile       | Executor            | Containers          | Best For                          |
+===============+=====================+=====================+===================================+
| ``slurm``     | SLURM (child jobs)  | Singularity + Conda | Production HPC runs               |
+---------------+---------------------+---------------------+-----------------------------------+
| ``conda``     | Local (current node)| Conda               | Interactive ``srun`` or laptop    |
+---------------+---------------------+---------------------+-----------------------------------+

Under ``slurm`` each process is submitted as its own SLURM job, with memory and wall-time
doubling on each automatic retry. Under ``conda`` everything runs on the current node —
use it inside an interactive ``srun`` session or on a workstation.

Shared environment cache (``--env_dir``)
----------------------------------------

By default Conda/Singularity environments are cached under Nextflow's working directory,
so each new run rebuilds them. Set ``--env_dir <path>`` to redirect the caches to
``<path>/conda`` and ``<path>/singularity`` — point multiple runs (or all users on a
cluster) at one shared location to build each environment once and reuse it:

.. code-block:: bash

   nextflow run main.nf ... --env_dir /shared/asap_envs -profile slurm

For large parallel runs, pre-build the environments once up front with
``nextflow/prebuild_envs.sh`` so that many jobs do not race to build the same environment
(concurrent builds can contend on the shared mamba lock and corrupt a half-built env).

Troubleshooting
---------------

- **Environment build is slow or fails.** Make sure Mamba is available, and use a shared
  ``--env_dir`` so environments build once and are reused. For big parallel runs,
  pre-build with ``prebuild_envs.sh``.
- **A run died partway.** Re-run the same command with ``-resume`` to reuse cached results
  and continue from the failed step.
- **A SLURM job hit a memory or time limit.** Processes automatically double memory and
  wall-time on each retry; persistent failures usually mean the target partition needs
  more headroom. See the performance notes in `OUTPUTS.rst <OUTPUTS.rst>`_.
- **Tests fail to launch.** Create the log directory first (``mkdir -p logs/nf-test``)
  before submitting. See `TESTING.rst <TESTING.rst>`_.

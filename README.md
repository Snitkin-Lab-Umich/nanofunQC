# nanofunQC
A quality control workflow implemented in Snakemake to assemble long read data for Candida auris. It is a modified version of nanoQC (https://github.com/Snitkin-Lab-Umich/nanoQC).

## Summary

This pipeline is very similar to [funQCD](https://github.com/Snitkin-Lab-Umich/funQCD) in terms of installation, setup, and output. The steps unique to this pipeline are:

- Runs [Filtlong](https://github.com/rrwick/Filtlong) to remove low quality reads and discards reads less than 1000 Bp.
- Generates pre and post-Filtlong QC plots using [Nanoplot](https://github.com/wdecoster/NanoPlot).
- Assemble clean filtlong nanopore reads with [Flye](https://github.com/fenderglass/Flye) assembler.
- Flye assembly is then polished with long reads using [Medaka](https://github.com/nanoporetech/medaka)
- The medaka assembly is then passed through [Prokka](https://github.com/tseemann/prokka) for annotation, [skani](https://github.com/bluenote-1577/skani) to identify closest reference genome, and [MLST](https://github.com/tseemann/mlst) for determining sequence type based on sequences of housekeeping genes
- Flye and medaka assemblies are then run on [BUSCO](https://busco.ezlab.org/) for assembly completeness statistics and [QUAST](https://quast.sourceforge.net/) for assembly statistics.

## Installation

> Clone the github directory onto your system.

```

git clone https://github.com/Snitkin-Lab-Umich/nanofunQC

```

> Load Bioinformatics, snakemake, and singularity modules.

```

module load Bioinformatics snakemake singularity

```

## Setup

Similar to funQCD, you'll first need to run setup.py with the path to your raw reads:

```

python setup.py [path_to_reads]

```

Next, you'll need to update config/config.yaml with a new prefix and long_reads path. The long_reads field should be the same path to your raw data used above, and the prefix field should be a unique name for the current run. Further options can be adjusted in the profile/config.yaml file.

## Running nanofunQC

> First, perform a dry run of nanofunQC by running the following command. This will show the steps and commands that the pipeline will execute in the real run, without actually executing anything. Remove the --quiet flag for a more detailed view.

```

snakemake -s workflow/nfqc.smk -p --configfile config/config.yaml --profile ./profile/ -n --quiet

```

> The snakemake options present in profile/config.yaml should be visible in the detailed dry run (such as memory and runtime for each rule). By default, --slurm is enabled in these options, and snakemake will submit jobs to the cluster using the account in your profile. If everything looks correct, start the run using a job script with minimal CPUs, moderate memory, and a long runtime. An example job script is provided in `run_nfqc.job`.

```
#!/bin/bash

#SBATCH --mail-user=[your_email]
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=[your_account]
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=20:00:00

module load Bioinformatics snakemake singularity
snakemake -s workflow/nfqc.smk -p --configfile config/config.yaml --profile ./profile/

```




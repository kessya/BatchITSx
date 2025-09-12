# BatchITSx

![Version](https://img.shields.io/badge/version-0.0.1-blue)
![Lifecycle](https://img.shields.io/badge/status-experimental-orange)
![License](https://img.shields.io/github/license/vmikk/BatchBlaster)

BatchITSx is a bioinformatics pipeline designed to run the [ITSx](https://microbiology.se/software/itsx/) tool in batch mode.

It is built using the [Nextflow](https://www.nextflow.io/) workflow management system, ensuring portability and reproducibility across multiple computing environments. The pipeline is primarily intended for High Performance Computing (HPC) clusters and supports task submission via the SLURM job scheduler.

## Features

- Run [ITSx](https://microbiology.se/software/itsx/) on multiple sequences in batch mode
- Sequence quality filtering based on ITSx results
- Scalable and reproducible analysis powered by Nextflow
- Compatible with Linux, MacOS, and Windows environments

## Quick Start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

    ```bash
    curl -s https://get.nextflow.io | bash
    ```
2. Run BatchITSx

    ```bash
    nextflow run kessya/BatchITSx -r main \
        --input 'path/to/your/input' \
        --itsx_target "F,All" \
        --region "itsfull" \
        --itsx_bin "ITSx" \
        --vsearch_bin "vsearch" \
        --itsx_chunksize "100"
    ```

## Parameters

- `--input` : Path to the input file containing the sequences (Required)
- `--outdir` : Path to the output directory (Default: `./results`)
- `--itsx_target` : Organism groups targeted by ITSx (Default: `F,All`)
- `--region` : For which region (itsfull or its2) the pipeline is customised for (Default: `itsfull`)
- `--itsx_bin` : Path to the ITSx executable (Default: `ITSx`)
- `--vsearch_bin` : Path to the VSEARCH executable (Default: `vsearch`)
- `--itsx_chunksize` : Sequnece count to be used when splitting FASTA for running ITSx in batch mode (Default: `100`)
- ...

## Output

Results will be saved in the specified output directory (`./results`, by default). Output includes:

- A file with ITSx results (ITS1, 5.8S, ITS2) in a FASTA format
- Separate files per target group with ITSx results in a FASTA format
- Separate summary reports per target group

## Dependencies

- [Nextflow](https://www.nextflow.io/) (>=25.04.6)
- [Singularity](https://sylabs.io/singularity/) or [Docker](https://www.docker.com/)

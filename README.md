# Transcriptome Analysis Pipeline

This repository contains a bioinformatics pipeline for processing paired-end RNA sequencing data to assemble and evaluate a transcriptome. It is designed to run on a Slurm-based High-Performance Computing (HPC) system, with all dependencies managed using `conda`.

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Repository Setup](#repository-setup)
3. [Environment Configuration](#environment-configuration)
   - [Installing Conda](#installing-conda)
   - [Creating the Conda Environment](#creating-the-conda-environment)
4. [Pipeline Overview](#pipeline-overview)
5. [Directory Structure](#directory-structure)
6. [Running the Pipeline](#running-the-pipeline)
   - [Preparing Input Data](#preparing-input-data)
   - [Executing the Script](#executing-the-script)
   - [Debug Mode](#debug-mode)
7. [Output Files](#output-files)
8. [Troubleshooting](#troubleshooting)
9. [Contributing](#contributing)
10. [License](#license)

## Prerequisites

- Access to a Slurm-based HPC cluster.
- Git installed (`git --version` to check).
- A working internet connection for downloading dependencies.
- Basic familiarity with command-line interfaces and HPC job scheduling.

## Repository Setup

1. **Clone the Repository**
   Clone from GitHub to your HPC system:
   ```bash
   git clone https://github.com/<username>/<repository-name>.git
   cd <repository-name>
   ```
   Replace `<username>` and `<repository-name>` with your GitHub details.

2. **Verify Directory Structure**
   After cloning, ensure the structure matches:
   ```
   .
   └── transcriptome/
       ├── bin/
       │   └── main.py
       ├── config/
       │   ├── transcriptome_1.yml
       │   └── transcriptome_2.yml
       ├── data/
       │   ├── raw_reads/
       │   ├── trimmed_reads/
       │   └── transcriptome/
       ├── results/
       │   ├── trimming/
       │   ├── merging/
       │   ├── normalization/
       │   ├── assembly/
       │   └── busco/
       ├── logs/
       │   ├── main/
       │   ├── trimming/
       │   ├── merging/
       │   ├── normalization/
       │   ├── assembly/
       │   └── busco/
       ├── .gitignore
       └── README.md
   ```

## Environment Configuration

### Installing Conda
Conda manages the Python environment and all bioinformatics tools.

1. **Install Miniconda (if not already installed)**
   Download and install Miniconda:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
   Follow the prompts and initialize Conda:
   ```bash
   source ~/.bashrc
   ```

### Creating the Conda Environment

1. **Create the Environment**
   Use the provided transcriptome_1.yml file:
   ```bash
   conda env create -f config/transcriptome_1.yml
   ```
   This creates an environment named `transcriptome` with all required tools.

2. **Activate the Environment**
   Activate the environment before running the pipeline:
   ```bash
   conda activate transcriptome
   ```

3. **Verify Tools**
   Check that all tools are installed:
   ```bash
   fastp --version
   pigz --version
   repair.sh -h
   Trinity --version
   rnaspades.py --version
   busco --version
   ```

4. **Set Trinity Path**
   The normalization step requires TRINITY_HOME. Set it in your shell:
   ```bash
   export TRINITY_HOME=$(dirname $(which Trinity))
   echo 'export TRINITY_HOME=$(dirname $(which Trinity))' >> ~/.bashrc
   ```

## Pipeline Overview

The pipeline processes paired-end RNA-seq data through these steps:

1. **Trimming**: Uses fastp to trim adapters and low-quality bases.
2. **Merging**: Combines trimmed reads into single files with pigz.
3. **Checking**: Verifies paired-end integrity with repair.sh from BBtools.
4. **Normalization**: Normalizes read coverage using Trinity's insilico_read_normalization.pl.
5. **Assembly**: Assembles transcripts with rnaSPAdes and Trinity.
6. **Evaluation**: Assesses assembly quality with BUSCO.

Each step submits jobs to Slurm with appropriate dependencies.

## Directory Structure

- `bin/main.py`: Core script orchestrating the pipeline.
- `config/`: Conda environment YAML files.
- `data/raw_reads/`: Input directory for raw FASTQ files.
- `data/trimmed_reads/`: Output directory for trimmed reads.
- `data/transcriptome/`: Final transcriptome assemblies.
- `results/`: Intermediate and final results.
- `logs/`: Log files for each step.

## Running the Pipeline

### Preparing Input Data

1. **Place Raw Reads**
   Copy paired-end FASTQ files into `data/raw_reads/`. Files must follow the naming convention:
   - `sample_R1.fastq` (forward reads)
   - `sample_R2.fastq` (reverse reads)
   
   Example:
   ```bash
   cp /path/to/reads/*_R1.fastq data/raw_reads/
   cp /path/to/reads/*_R2.fastq data/raw_reads/
   ```
   
   **Note**: Make sure all paired end read files have R1.fastq or R2.fastq for code to work

2. **Compress Files (Optional)**
   If files are uncompressed, compress them with pigz:
   ```bash
   pigz data/raw_reads/*.fastq
   ```

### Executing the Script

1. **Navigate to bin/**
   ```bash
   cd bin
   ```

2. **Run the Pipeline**
   Ensure the transcriptome environment is active, then execute:
   ```bash
   python main.py
   ```
   The script submits Slurm jobs and prints job IDs for tracking.

3. **Monitor Jobs**
   Check job status with:
   ```bash
   squeue -u $USER
   ```

### Debug Mode

To skip completed steps (useful for testing or resuming):

```bash
python main.py -d
```

The `-d` flag checks for existing output files and skips corresponding steps if they are present and non-empty.

## Output Files

- **Trimming**: `data/trimmed_reads/sample_trimmed_R[1,2].fastq.gz`, `results/trimming/sample_fastp.[html,json]`
- **Merging**: `results/merging/merged_R[1,2].fastq.gz`, `results/merging/fixed_R[1,2].fastq.gz`
- **Normalization**: `results/normalization/left.norm.fq`, `results/normalization/right.norm.fq`
- **Assembly**:
  - rnaSPAdes: `results/assembly/rnaspades/transcripts.fasta`
  - Trinity: `results/assembly/trinity/Trinity.fasta`
- **BUSCO**: `results/busco/[rnaspades,trinity]/short_summary.*.txt`

Logs are stored in `logs/` subdirectories for each step.

## Troubleshooting

- **Job Submission Fails**: Check Slurm error logs (e.g., `logs/trimming/sample_fastp.err`) and ensure resource requests match HPC limits.
- **Tool Not Found**: Verify the environment is active (`conda activate transcriptome`) and tools are installed (`conda list`).
- **Missing Pairs**: Ensure raw read filenames include `_R1` and `_R2` suffixes.
- **Permission Issues**: Confirm write access to all directories.

For further assistance, consult Slurm logs or your HPC administrator.

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/your-feature`).
3. Commit changes (`git commit -m "Add your feature"`).
4. Push to the branch (`git push origin feature/your-feature`).
5. Open a Pull Request.

## License

[Add license information here]
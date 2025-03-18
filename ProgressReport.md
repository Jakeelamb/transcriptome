# Progress Report: De Novo Transcriptome Assembly Using Single-Cell and Bulk RNA-Seq Datasets

## Project Goal

The primary objective of this project is to construct a high-quality de novo transcriptome assembly by integrating single-cell and bulk RNA sequencing (RNA-Seq) datasets. By leveraging the complementary strengths of these data types—single-cell RNA-Seq for capturing cellular heterogeneity and bulk RNA-Seq for deeper sequencing coverage—the aim is to produce a comprehensive and accurate representation of the transcriptome. This assembly will serve as a foundational resource for downstream biological analyses, such as gene expression profiling and functional annotation, with an emphasis on achieving the highest possible quality as assessed by completeness and accuracy metrics.

## Workflow Outline and Progress

The pipeline for this project is designed to process large-scale RNA-Seq datasets efficiently on a Slurm-based High-Performance Computing (HPC) system. Below is an outline of the workflow, current progress, and key decisions made at each step.

### 1. Read Trimming with fastp

**Process:** Raw paired-end FASTQ files are trimmed using fastp with the following parameters: a quality threshold of 20 (`-q 20`), deduplication enabled (`--dedup`), and automatic detection of paired-end adapters (`--detect_adapter_for_pe`). These settings ensure the removal of low-quality bases, duplicate reads, and adapter sequences that could compromise downstream assembly.

**Progress:** 

**Rationale:** A quality threshold of 20 (Phred score) balances the removal of sequencing errors with the preservation of usable data, as it filters out bases with a >1% error rate. Deduplication reduces redundancy, which is critical for single-cell data where PCR amplification artifacts are common. Adapter detection for paired-end reads ensures compatibility with the diverse naming conventions in the dataset (e.g., `_R1.fastq.gz` and `_R2.fastq.gz`), as outlined in the README.md.

![Figure 1: Distribution of raw vs. trimmed read sizes]()

### 2. Read Merging and Normalization

**Process:**

* **Merging:** Trimmed FASTQ files are converted to FASTA format using `seqtk seq -A`, reducing file size by approximately half due to the exclusion of quality scores. These files are then concatenated into two large files: `left.fa` (1.88 TB) and `right.fa` (1.88 TB), representing all forward and reverse reads, respectively.
* **Normalization:** The merged FASTA files undergo k-mer normalization using `bbnorm.sh` from the BBTools suite, targeting a coverage depth of 100x with a k-mer size of 25 (`k=25`). Normalization reduces the dataset's computational complexity by capping high-abundance k-mers while retaining low-abundance ones.

**Progress:** Merging has been completed, with `left.fa` and `right.fa` generated as reported.

**Rationale:** Merging the trimmed reads before normalization allows for true normalization. Normalizing subsets of data could mask true signal. 

* **Merging:** Combining all trimmed reads into two files is necessary for normalization, as `bbnorm.sh` operates on paired-end datasets holistically. Converting to FASTA halves the storage footprint, reducing the computational footprint.
* **Normalization:** Normalization reduces sequencing depth variability, which is particularly beneficial for mixed single-cell and bulk datasets with uneven coverage. By targeting 100x coverage, the pipeline retains sufficient depth for rare transcripts while minimizing computational load for assembly. A k-mer size of 25 was chosen as it strikes a balance between specificity and sensitivity: smaller k-mers (e.g., 21) increase overlap detection but risk spurious matches in repetitive regions, while larger k-mers (e.g., 31) enhance specificity but may miss short or lowly expressed transcripts. Given the diverse dataset, 25 optimizes assembly quality across both single-cell and bulk RNA-Seq data.

![Figure 2: K-mer abundance histograms before and after normalization (to be updated with final k-mer data)]()

### 3. Transcriptome Assembly

**Plan:** The normalized FASTA files (`left.norm.fa` and `right.norm.fa`) will be assembled using rnaSPAdes, a de novo assembler optimized for RNA-Seq data. Depending on the size of the normalized dataset and computational resources, Trinity may also be employed as a complementary assembler to compare assembly outcomes.

**Progress:** Assembly has not yet started, as normalization is ongoing. 
**Rationale:** rnaSPAdes was selected for its efficiency with large datasets and ability to handle variable transcript expression levels, which is ideal for the mixed single-cell and bulk RNA-Seq data. Trinity is considered as a secondary option due to its robustness with complex transcriptomes, though its higher memory demands may necessitate evaluation post-normalization. Using both assemblers could provide a comparative analysis to maximize quality, contingent on resource availability.

### 4. Quality Assessment with BUSCO

**Plan:** The assembled transcriptomes will be evaluated using BUSCO (Benchmarking Universal Single-Copy Orthologs) with the eukaryota_odb10 lineage to assess completeness based on conserved eukaryotic genes.

**Progress:** An initial run (details not provided in the current report) yielded a BUSCO score of 37%, indicating significant room for improvement. The current assembly step has not yet produced outputs for re-evaluation.

**Rationale:** BUSCO provides a standardized metric for transcriptome completeness, critical for validating the "highest quality possible" goal. The low initial score (37%) suggests issues with data preprocessing or assembly parameters in the prior attempt, motivating the current optimized workflow with normalization and advanced assemblers.

![Figure 3: BUSCO completeness scores for initial and final assemblies]()

## Decision-Making Rationale Summary

The pipeline design reflects a strategic balance of quality, efficiency, and computational feasibility:

* **Trimming Parameters:** Quality 20, deduplication, and adapter detection ensure a clean dataset, addressing sequencing artifacts prevalent in single-cell RNA-Seq.
* **Normalization Strategy:** Merging and normalizing with `bbnorm.sh` at 100x coverage and k=25 reduces data complexity while preserving transcript diversity, tailored to the mixed dataset's needs.
* **Assembler Choice:** rnaSPAdes prioritizes efficiency, with Trinity as a potential backup to leverage its strengths if resources permit.
* **Quality Focus:** BUSCO evaluation drives iterative improvements, targeting a high completeness score to meet the project goal.

## Next Steps

2. Update the statistics report with k-mer data post-normalization by re-running `generate_stats.py --parse-only`.
3. Proceed with rnaSPAdes assembly, followed by optional Trinity assembly if feasible.
4. Reassess BUSCO scores and refine the pipeline as needed to enhance transcriptome quality.

This workflow positions the project to achieve a high-quality de novo transcriptome, leveraging the strengths of both single-cell and bulk RNA-Seq datasets.
# Project 1. Candida auris Fungal NGS Pipeline

![Python](https://img.shields.io/badge/Python-3.11-blue?logo=python&logoColor=white)
![Nextflow](https://img.shields.io/badge/Nextflow-Workflow-23AA62?logo=nextflow&logoColor=white)
![HISAT2](https://img.shields.io/badge/HISAT2-Alignment-orange)
![SAMtools](https://img.shields.io/badge/SAMtools-BAM%20Processing-4C8CBF)
![BEDTools](https://img.shields.io/badge/BEDTools-Coverage-purple)
![Linux](https://img.shields.io/badge/Linux-WSL2-FCC624?logo=linux&logoColor=black)
![Tests](https://img.shields.io/badge/Tests-6%20passed-brightgreen)
![Status](https://img.shields.io/badge/Status-Documentation%20in%20progress-yellow)

A workflow-based NGS preprocessing and alignment project for *Candida auris*, built with Nextflow, Linux/Bash, and Python. It covers read quality control, preprocessing, reference-genome alignment, BAM processing, alignment QC, and gene-level coverage summaries.

> **Scope note:** This repository contains multiple related workflows, not one unified workflow that processes all six cohort samples end to end. `workflow/main.nf` uses a separate technical sample (`SRR26185932`); `workflow/resistant_cohort.nf` processes the three resistant samples; existing results for the three sensitive samples were generated separately.

## Objectives

- Inspect raw sequencing-read quality and preprocessing outcomes.
- Trim adapters and preprocess reads.
- Align reads to the *C. auris* B8441 V2 reference genome using HISAT2.
- Sort and index BAM alignments and inspect alignment statistics with SAMtools.
- Calculate gene-level coverage with BEDTools.
- Summarize alignment QC, coverage breadth, and estimated mean sequencing depth using Python.
- Organize workflow code, scripts, tests, metadata, documentation, and results.

## Dataset and sample groups

The cohort contains six samples labelled `sensitive` or `resistant` in `metadata/samples.tsv`.

| GEO sample | SRA experiment | SRA run | Phenotype |
|---|---|---|---|
| GSM4057610 | SRX6793583 | SRR10059687 | sensitive |
| GSM4057611 | SRX6793584 | SRR10059688 | sensitive |
| GSM4057612 | SRX6793585 | SRR10059689 | sensitive |
| GSM4057613 | SRX6793586 | SRR10059690 | resistant |
| GSM4057614 | SRX6793587 | SRR10059691 | resistant |
| GSM4057615 | SRX6793588 | SRR10059692 | resistant |

`SRR26185932` is an additional technical/workflow sample used by `workflow/main.nf`; it is not a seventh member of the six-sample phenotype cohort.

## Reference genome

- Organism: *Candida auris*
- Reference strain: B8441, V2
- Assembly accession: `GCA_002759435.2`
- Reference FASTA: `data/reference/ncbi_dataset/data/GCA_002759435.2/GCA_002759435.2_Cand_auris_B8441_V2_genomic.fna`
- Annotation: `data/reference/ncbi_dataset/data/GCA_002759435.2/genomic.gff`
- Genome-length file: `data/reference/candida_auris.genome`
- Reported genome size: 12,365,959 bp across 15 contigs

The HISAT2 index is expected under `data/reference/hisat2_index_current/`, with prefix `candida_auris_B8441_v2`.

## Workflow overview

### `workflow/main.nf`

Processes the separate sample `SRR26185932` through raw-read QC, Cutadapt trimming, fastp preprocessing, post-processing QC, HISAT2 alignment, SAMtools BAM processing, and per-contig BEDTools coverage. Its default output directory is `results/nextflow/`.

### `workflow/resistant_cohort.nf`

Defines a workflow for the three resistant cohort samples (`SRR10059690`, `SRR10059691`, and `SRR10059692`). It includes SRA-to-FASTQ conversion, QC, read preprocessing, alignment/BAM processing, and gene-level coverage. Its default output directory is `results/resistant_cohort/`.

### `workflow/bedtools_coverage.nf`

Runs a separate BEDTools gene-coverage workflow on the existing BAM for `SRR26185932`. Its default output directory is `results/nextflow/bedtools_verified/`.

**Important:** The workflows have different scopes and are not a single end-to-end execution for the full six-sample cohort. The sensitive-sample results are present in project output directories, but are not produced by `workflow/resistant_cohort.nf`.

## Software

The project uses Nextflow, SRA Toolkit (`fasterq-dump`), FastQC, Cutadapt, fastp, HISAT2, SAMtools, BEDTools, and Python with NumPy and pandas.

Recorded environment versions are listed in [`docs/tool_versions.txt`](docs/tool_versions.txt). Python dependencies are listed in [`requirements.txt`](requirements.txt). The analysis was developed in a Linux/WSL environment using Conda.

## Python utilities and tests

- `scripts/alignment_qc.py` — extracts and summarizes alignment QC metrics.
- `scripts/coverage_summary.py` — summarizes gene-level coverage, including the fraction of each gene covered.
- `scripts/depth_summary.py` — estimates mean depth from mapped bases and reference genome length.

Unit tests are located in `test/test_pipeline.py`.

Run the tests from the project root:

```bash
python -m unittest discover -s test -p 'test_*.py' -v
```

Compile-check the Python scripts:

```bash
python -m py_compile scripts/alignment_qc.py scripts/coverage_summary.py scripts/depth_summary.py test/test_pipeline.py
```

## Selected results

The following metrics are taken from the existing summary tables in `results/summary/`.

### Alignment

| Sample | Phenotype | Primary mapping rate | Properly paired |
|---|---|---:|---:|
| SRR10059687 | sensitive | 91.78% | 86.85% |
| SRR10059688 | sensitive | 89.47% | 83.17% |
| SRR10059689 | sensitive | 85.78% | 78.68% |
| SRR10059690 | resistant | 91.41% | 86.92% |
| SRR10059691 | resistant | 89.75% | 83.80% |
| SRR10059692 | resistant | 91.30% | 87.71% |

Mean primary mapping rate across the six samples: approximately **89.92%**.

### Gene-level coverage breadth

The coverage summary contains 5,583 annotated gene features per sample. Mean gene fraction covered ranges from approximately **99.87% to 99.91%**; 3–5 genes per sample have zero covered bases according to the summary.

Coverage breadth is the fraction of a feature's bases with coverage; it is not the same as sequencing depth. Neither measure establishes gene expression or biological activity.

### Estimated mean sequencing depth

The current summary estimates mean depth from mapped bases divided by reference genome length. Estimates range from approximately **384× to 746×** across the six samples. These are derived estimates, not direct per-position depth measurements.

## Running workflow previews

Run these commands from the project root in an environment with the required tools installed:

```bash
nextflow run workflow/main.nf -preview
nextflow run workflow/resistant_cohort.nf -preview
nextflow run workflow/bedtools_coverage.nf -preview
```

`-preview` checks workflow planning and process resolution; it **does not execute the analysis**. A full run can require substantial time, disk space, and memory. Whole-BAM gene-coverage calculations can be resource-intensive. Use existing outputs for inspection unless you have reviewed resource requirements and intend to rerun the workflow.

## Project structure

```text
.
├── workflow/
│   ├── main.nf
│   ├── resistant_cohort.nf
│   └── bedtools_coverage.nf
├── scripts/
│   ├── alignment_qc.py
│   ├── coverage_summary.py
│   └── depth_summary.py
├── test/
│   └── test_pipeline.py
├── metadata/
│   └── samples.tsv
├── data/
│   └── reference/
├── results/
│   └── summary/
├── docs/
│   └── tool_versions.txt
├── nextflow.config
└── requirements.txt
```

Large raw sequencing files, BAM files, indexes, temporary work directories, and detailed intermediate outputs may be excluded from version control; see `.gitignore`. Choose the exact set of files included in a public repository deliberately to avoid committing large or redundant files.

## Limitations and interpretation

- The three workflows are separate and do not constitute a single unified run across all six cohort samples.
- The sensitive and resistant sample outputs were produced through separate processing paths; differences in preprocessing should be considered when making comparisons.
- Current summary metrics describe read quality/alignment and gene-level coverage. They do not by themselves identify differentially expressed genes or establish antifungal-resistance mechanisms.
- Estimated mean depth is calculated from mapped bases and genome length; it is not a direct per-position depth profile.
- Overrepresented-sequence results are documented in `results/summary/overrepresented_sequences_interpretation.md`. Their origin was not definitively established, so they should not be labelled contamination without further evidence.
- Workflow previews and passing unit tests do not substitute for a full reproducible rerun of every workflow.

## Data provenance and citation

Sample identifiers and phenotype labels are recorded in `metadata/samples.tsv`. The reference assembly used by the project is NCBI assembly `GCA_002759435.2` (*Candida auris* B8441 V2). Before publishing, add the primary sequencing-study citation and any required dataset/reference citations, and specify a license appropriate for the code and documentation.

## Author

Add your preferred author name or profile link here.

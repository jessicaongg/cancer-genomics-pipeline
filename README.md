# Cancer Genomics Pipeline

## Overview
This repository contains reproducible cancer genomics workflows developed for paired-end next-generation sequencing (NGS) analysis in Linux environments.

## Workflow Overview

![Pipeline Overview](workflow/pipeline_overview.png)

The pipeline includes:
- Read preprocessing
- Alignment
- BAM processing
- Variant calling
- Structural variant detection
- Functional annotation
- Quality control and downstream analysis

## Tools Used
- bwa
- samtools
- bcftools
- GATK
- Manta
- snpEff
- SnpSift
- bedtools

## Skills Demonstrated
- Bash scripting
- Variant calling workflows
- Structural variant analysis
- Linux command-line bioinformatics
- Reproducible computational pipelines
- NGS data processing

## Repository Structure

```text
scripts/        Bash workflows and pipeline scripts
results/        Example outputs and summary tables
workflow/       Pipeline diagrams and workflow figures
docs/           Additional methodology notes
environment/    Software versions and dependencies

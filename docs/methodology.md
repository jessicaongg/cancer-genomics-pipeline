# Methodology

This repository demonstrates a paired tumour-normal cancer genomics workflow for next-generation sequencing data.

## Pipeline Summary

1. **Reference alignment**
   - Paired-end FASTQ files are aligned to a human reference genome using `bwa mem`.

2. **BAM processing**
   - SAM/BAM files are sorted and indexed using `samtools`.

3. **Duplicate marking**
   - Duplicate reads are identified using `samtools markdup`.

4. **Read group assignment**
   - Read group metadata is added using `gatk AddOrReplaceReadGroups`.

5. **Variant calling**
   - Germline SNVs and INDELs are called using `gatk HaplotypeCaller`.

6. **Variant annotation**
   - Variants are annotated using `snpEff` and `SnpSift`.

7. **Structural variant calling**
   - Structural variants are identified using `Manta`.

8. **Summary generation**
   - Mapping statistics and annotated variant tables are generated for downstream interpretation.

## Notes

Raw sequencing data and large intermediate files are not included in this repository.

# Cancer Genomics Pipeline Workflow

```text
Raw FASTQ files
        ↓
Quality Control (FastQC)
        ↓
Adapter Trimming / Preprocessing
        ↓
Reference Alignment (bwa)
        ↓
SAM → BAM Conversion (samtools)
        ↓
BAM Sorting & Indexing
        ↓
Variant Calling (bcftools / GATK)
        ↓
Structural Variant Detection (Manta)
        ↓
Variant Annotation (snpEff / SnpSift)
        ↓
Summary Statistics & Downstream Analysis

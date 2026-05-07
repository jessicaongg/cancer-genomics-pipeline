#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# NGS Read Mapping and Variant Calling Pipeline
# ============================================================
# This script performs:
# 1. Reference genome indexing
# 2. Paired-end read alignment using BWA
# 3. SAM to BAM conversion
# 4. BAM sorting
# 5. Variant calling using bcftools
#
# Note:
# Replace file paths with your own input FASTQ and reference files.
# ============================================================

# -----------------------------
# Input files
# -----------------------------

REF="reference_genome.fasta"
DIR="directory"

SAMPLE_1="sample_1"
SAMPLE_2="sample_2"

SAMPLE_1_R1="sample_1_R1.fastq"
SAMPLE_1_R2="sample_1_R2.fastq"

SAMPLE_2_R1="sample_2_R1.fastq"
SAMPLE_2_R2="sample_2_R2.fastq"

READS=(
 "sample_1_R1.fastq"
 "sample_1_R2.fastq"
 "sample_2_R1.fastq"
 "sample_2_R2.fastq"
)                                 # Paired-end read files

log(){ printf "[%s] %s\n" "$(date +%T)" "$*"; } # log messages with timestamp
# ---------------------- Step 1: Create working directory ---------------------
cd
mkdir -p sai sam bam vcf logs

# ------------- Step 2: Make shortcuts to practice files & reference ----------
# using symbolic links to avoid copying large files --> save dsk space
ln -sf "${DIR}/.../"* . # Link all files into current directory
ln -sf "${DIR}/${REF}" . # Link the reference genome file

###############################################################################
# ------------------------------ Preprocessing step ---------------------------
###############################################################################
# ---------------------------- Step 3: Index reference ----------------------- 

echo "Indexing reference genome..."
bwa index "$REF"

# ----------------------- Step 4: Align reads with BWA ----------------------- 
echo "Aligning reads for ${SAMPLE_1}..."
bwa aln "$REF" "$SAMPLE_1_R1" > "sai/${SAMPLE_1}_R1.sai"       #forward
bwa aln "$REF" "$SAMPLE_1_R2" > "sai/${SAMPLE_1}_R2.sai"       #reverse

echo "Aligning reads for ${SAMPLE_2}..."
bwa aln "$REF" "$SAMPLE_2_R1" > "sai/${SAMPLE_2}_R1.sai"       #forward
bwa aln "$REF" "$SAMPLE_2_R2" > "sai/${SAMPLE_2}_R2.sai"      #reverse

###############################################################################
# --------------------------------- Read mapping ------------------------------
###############################################################################

# ------------------- Step 5: Generate SAM files ----------------------------- 
echo "Generating paired-end (SAM) files..."

bwa sampe "$REF" \
  "sai/${SAMPLE_1}_R1.sai" \
  "sai/${SAMPLE_1}_R2.sai" \
  "$SAMPLE_1_R1" "$SAMPLE_1_R2" \
  > "sam/${SAMPLE_1}.sam"

bwa sampe "$REF" \
  "sai/${SAMPLE_2}_R1.sai" \
  "sai/${SAMPLE_2}_R2.sai" \
  "$SAMPLE_2_R1" "$SAMPLE_2_R2" \
  > "sam/${SAMPLE_2}.sam"

###############################################################################
# ------------------------------- Convert and sort ----------------------------
###############################################################################

# -------------------- Step 6: Convert and sort SAM → BAM --------------------- 
echo "Converting SAM to BAM..."

samtools view -bS "sam/${SAMPLE_1}.sam" > "bam/${SAMPLE_1}.bam"
samtools view -bS "sam/${SAMPLE_2}.sam" > "bam/${SAMPLE_2}.bam"

echo "Sorting BAM files..."

samtools sort "bam/${SAMPLE_1}.bam" -o "bam/${SAMPLE_1}.sorted.bam"
samtools sort "bam/${SAMPLE_2}.bam" -o "bam/${SAMPLE_2}.sorted.bam"

# ---------------------------- Step 7: Index BAM files ------------------------- 

echo "Indexing BAM files..."

samtools index "bam/${SAMPLE_1}.sorted.bam"
samtools index "bam/${SAMPLE_2}.sorted.bam"

###############################################################################
# -------------------------- Cleanup and Organization -------------------------
###############################################################################

# ----------------------- Step 8: Remove unecessary files ---------------------

rm -f "sam/${SAMPLE_1}.sam" "sam/${SAMPLE_2}.sam"      # Remove large SAM files

# ----------------------------- Step 9: Rename files --------------------------

mv -f "bam/${SAMPLE_1}.sorted.bam" "bam/${SAMPLE_1}.bam" 
mv -f "bam/${SAMPLE_2}.sorted.bam" "bam/${SAMPLE_2}.bam"

###############################################################################
# ------------------------------- Variant calling -----------------------------
###############################################################################

# -------------------- Step 10: Variant calling with bcftools ----------------- 

echo "Calling variants..."

bcftools mpileup -f "$REF" "bam/${SAMPLE_1}.sorted.bam" \
  | bcftools call --ploidy 1 -vc \
  > "vcf/${SAMPLE_1}_variants.vcf"

bcftools mpileup -f "$REF" "bam/${SAMPLE_2}.sorted.bam" \
  | bcftools call --ploidy 1 -vc \
  > "vcf/${SAMPLE_2}_variants.vcf"

# -------------- Step 11: QC and compare variant calling results --------------
# Extract columns:
# f2 - position, f4 - reference base, f5 - alternate base, f6 - quality 

echo "Generating variant summaries..."

grep -v "^#" "vcf/${SAMPLE_1}_variants.vcf" | cut -f2,4-6 \
  > "vcf/${SAMPLE_1}_variant_summary.tsv"

grep -v "^#" "vcf/${SAMPLE_2}_variants.vcf" | cut -f2,4-6 \
  > "vcf/${SAMPLE_2}_variant_summary.tsv"

echo "Pipeline completed successfully."

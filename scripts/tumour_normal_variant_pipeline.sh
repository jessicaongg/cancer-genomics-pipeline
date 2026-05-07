#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Tumour-Normal Cancer Genomics Variant Calling Pipeline
# ============================================================
# This script demonstrates a paired tumour-normal NGS workflow:
# 1. Reference alignment
# 2. BAM sorting and indexing
# 3. Duplicate marking
# 4. Germline SNV/INDEL calling
# 5. Variant annotation
# 6. Structural variant calling
#
# Replace placeholder paths with your own reference, FASTQ and resource files.
# ============================================================

# -----------------------------
# Input files and directories
# -----------------------------
# project directories
PROJECT_ROOT="$(pwd)"

DATA_DIR="${PROJECT_ROOT}/data"
REF_DIR="${PROJECT_ROOT}/reference"

# Reference genome
REF="Homo_sapiens_assembly38.fasta"
REF_FA="${REF_DIR}/${REF}"
INTERVALS="${REF_DIR}/target_regions.bed"
GNOMAD_AF="${REF_DIR}/af-only-gnomad.hg38.vcf.gz"

# Define FASTQ inputs
NORMAL_R1="data/normal_R1.fastq.gz"
NORMAL_R2="data/normal_R2.fastq.gz"
TUMOUR_R1="data/tumour_R1.fastq.gz"
TUMOUR_R2="data/tumour_R2.fastq.gz"

THREADS=4
log() {
    printf "[%s] %s\n" "$(date +%T)" "$*"
}

mkdir -p bam vcf metrics annotation structural_variants summary logs tmp

#################################################################################
# -------------------- Preprocessing step and read mapping-----------------------
#################################################################################

# ------ Step 1 & 2: bwa mem alignment of normal and tumour to reference --------

log "Aligning normal sample..."

bwa mem -t "$THREADS" "$REF_FA" "$NORMAL_R1" "$NORMAL_R2" \
    2> "$LOG_DIR/normal_bwa.log" \
    | samtools sort -@ "$THREADS" -o "$BAM_DIR/normal.sorted.bam"

samtools index "$BAM_DIR/normal.sorted.bam"

samtools flagstat "$BAM_DIR/normal.sorted.bam" \
    > "$METRICS_DIR/normal_flagstat.txt"

log "Aligning tumour sample..."

bwa mem -t "$THREADS" "$REF_FA" "$TUMOUR_R1" "$TUMOUR_R2" \
    2> "$LOG_DIR/tumour_bwa.log" \
    | samtools sort -@ "$THREADS" -o "$BAM_DIR/tumour.sorted.bam"

samtools index "$BAM_DIR/tumour.sorted.bam"

samtools flagstat "$BAM_DIR/tumour.sorted.bam" \
    > "$METRICS_DIR/tumour_flagstat.txt"

# ------------------------ Step 3: marking duplicates ---------------------------

log "Marking duplicates for normal sample..."

samtools sort -n -@ "$THREADS" \
    -o "$TMP_DIR/normal.namesort.bam" \
    "$BAM_DIR/normal.sorted.bam"                 # sort read name

samtools fixmate -m \
    "$TMP_DIR/normal.namesort.bam" \
    "$TMP_DIR/normal.fixmate.bam"               # fixmate info

samtools sort -@ "$THREADS" \
    -o "$TMP_DIR/normal.possort.bam" \
    "$TMP_DIR/normal.fixmate.bam"                # sort by position

samtools markdup -s \
    "$TMP_DIR/normal.possort.bam" \
    "$BAM_DIR/normal.dedup.bam"                 # mark duplicates

samtools index "$BAM_DIR/normal.dedup.bam"       # create BAM index

# compute mapping stats and save metrics
samtools flagstat "$BAM_DIR/normal.dedup.bam" \
    > "$METRICS_DIR/normal_dedup_flagstat.txt"

log "Marking duplicates for tumour sample..."

samtools sort -n -@ "$THREADS" \
    -o "$TMP_DIR/tumour.namesort.bam" \
    "$BAM_DIR/tumour.sorted.bam"

samtools fixmate -m \
    "$TMP_DIR/tumour.namesort.bam" \
    "$TMP_DIR/tumour.fixmate.bam"

samtools sort -@ "$THREADS" \
    -o "$TMP_DIR/tumour.possort.bam" \
    "$TMP_DIR/tumour.fixmate.bam"

samtools markdup -s \
    "$TMP_DIR/tumour.possort.bam" \
    "$BAM_DIR/tumour.dedup.bam"

samtools index "$BAM_DIR/tumour.dedup.bam"

samtools flagstat "$BAM_DIR/tumour.dedup.bam" \
    > "$METRICS_DIR/tumour_dedup_flagstat.txt"
    
# --------- Step 4: Create sequence dictionary and Add read groups ------------------

DICT="${REF_FA%.fasta}.dict"

log "Checking sequence dictionary..."

if [[ ! -f "$DICT" ]]; then
    log "Creating sequence dictionary..."
    gatk CreateSequenceDictionary \
        -R "$REF_FA" \
        -O "$DICT"
fi

log "Adding read groups to normal BAM..."

gatk AddOrReplaceReadGroups \
    -I "$BAM_DIR/normal.dedup.bam" \
    -O "$BAM_DIR/normal.rg.bam" \
    -RGID normal \
    -RGLB lib1 \
    -RGPL ILLUMINA \
    -RGPU unit1 \
    -RGSM normal

samtools index "$BAM_DIR/normal.rg.bam"

log "Adding read groups to tumour BAM..."

gatk AddOrReplaceReadGroups \
    -I "$BAM_DIR/tumour.dedup.bam" \
    -O "$BAM_DIR/tumour.rg.bam" \
    -RGID tumour \
    -RGLB lib1 \
    -RGPL ILLUMINA \
    -RGPU unit1 \
    -RGSM tumour

samtools index "$BAM_DIR/tumour.rg.bam"


#################################################################################
# ------------------------------- Variant calling -------------------------------
#################################################################################

# ------------------------ Step 5:Germline SNV/INDEL calling ---------------------------

echo "Calling germline SNVs/INDELs using GATK HaplotypeCaller..."

gatk HaplotypeCaller \
    -R "$REF_FA" \
    -I "$BAM_DIR/normal.rg.bam" \
    -L "$INTERVALS" \
    -O "$VCF_DIR/normal.germline.vcf.gz"

# ------------------------ Step 6: Functional annotation ---------------------------
echo "Annotating variants with snpEff..."

snpEff GRCh38.105 "$VCF_DIR/normal.germline.vcf.gz" \
    > "$ANNOTATION_DIR/normal.germline.snpeff.vcf"

log "Annotating variants with gnomAD allele frequencies..."

SnpSift annotate "$GNOMAD_AF" \
    "$ANNOTATION_DIR/normal.germline.snpeff.vcf" \
    > "$ANNOTATION_DIR/normal.germline.annotated.vcf"

# ----------- Step 7: Extract annotated moderate/high impact variants -----------------
log "Extracting MODERATE/HIGH impact variants..."

SnpSift extractFields -s "," \
    "$ANNOTATION_DIR/normal.germline.annotated.vcf" \
    CHROM POS ID REF ALT "ANN[*].IMPACT" "ANN[*].GENE" \
    "ANN[*].EFFECT" "ANN[*].HGVS_P" AF \
    > "$SUMMARY_DIR/annotated_variants_raw.tsv"

awk -F'\t' 'BEGIN {
    OFS="\t";
    print "Chrom","Pos","ID","Ref","Alt","Impact","Gene","Effect","AA_change","gnomAD_AF"
}
NR > 1 && ($6 ~ /MODERATE|HIGH/) {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10
}' "$SUMMARY_DIR/annotated_variants_raw.tsv" \
    > "$SUMMARY_DIR/moderate_high_impact_variants.tsv"
    
# ------------------- Step 8: Structural variant calling ------------------------------
log "Configuring Manta structural variant workflow..."

configManta.py \
    --tumorBam "$BAM_DIR/tumour.rg.bam" \
    --normalBam "$BAM_DIR/normal.rg.bam" \
    --referenceFasta "$REF_FA" \
    --callRegions "$INTERVALS" \
    --runDir "$SV_DIR/manta_run"

log "Running Manta..."

python "$SV_DIR/manta_run/runWorkflow.py" \
    -m local \
    -j "$THREADS"

cp "$SV_DIR/manta_run/results/variants/diploidSV.vcf.gz" \
    "$SV_DIR/tumour_normal_structural_variants.vcf.gz"
    
#################################################################################
# ------------------------------ Generate summary -------------------------------
#################################################################################

log "Generating alignment summary table..."

{
    echo -e "Sample\tTotal_Reads\tMapped_Reads\tUnmapped_Reads\tDuplicates"

    for sample in normal tumour; do
        total=$(grep "in total" "$METRICS_DIR/${sample}_flagstat.txt" | awk '{print $1}')
        mapped=$(grep "mapped (" "$METRICS_DIR/${sample}_flagstat.txt" | head -1 | awk '{print $1}')
        unmapped=$((total - mapped))
        duplicates=$(grep "duplicates" "$METRICS_DIR/${sample}_dedup_flagstat.txt" | awk '{print $1}')

        echo -e "${sample}\t${total}\t${mapped}\t${unmapped}\t${duplicates}"
    done
} > "$SUMMARY_DIR/alignment_summary.tsv"

# ------------------- Step 9: Cleanup temporary files ------------------------------

log "Cleaning temporary files..."

rm -f "$TMP_DIR"/*.namesort.bam
rm -f "$TMP_DIR"/*.fixmate.bam
rm -f "$TMP_DIR"/*.possort.bam

log "Tumour-normal cancer genomics pipeline completed successfully."

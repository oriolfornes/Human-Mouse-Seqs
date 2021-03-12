#!/usr/bin/env bash

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From https://screen.encodeproject.org/
# SCREEN: Search Candidate cis-Regulatory Elements by ENCODE
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

URL="https://api.wenglab.org/screen_v13/fdownloads"
RANDOM_SEED=123

# For each genome...
for GENOME in "GRCh38" "mm10"; do

    # Get genome FASTA file
    GENOME_FA="../Genomes/${GENOME}/${GENOME}.fa"
    if [ ! -f ${GENOME_FA} ]; then
        GENOME_FA="../Genomes/hg38/hg38.fa"
    fi

    # Download BED file
    if ! [ -f "${GENOME}-ccREs.bed" ]; then
        curl -O "${URL}/${GENOME}-ccREs.bed"
    fi

    # For each file...
    for FILE in "PLS" "pELS" "dELS" "CTCF-only" "DNase-H3K4me3"; do

        # Initialize
        PREFIX="${GENOME}-ccREs.${FILE}"

        # Download BED file
        if ! [ -f "${PREFIX}.bed" ]; then
            curl -O "${URL}/${PREFIX}.bed"
        fi

        # Resize to 250bp
        if [ ! -f "${PREFIX}.250bp.bed" ]; then
            awk '{C=$2+(($3-$2)/2);printf $1"\t%.0f\t%.0f\n",C-125,C+125;}' \
                "${PREFIX}.bed" > "${PREFIX}.250bp.bed"
        fi

        # Get FASTA
        if [ ! -f "${PREFIX}.250bp.fa" ]; then
            bedtools getfasta -fi ${GENOME_FA} -fo "${PREFIX}.250bp.fa" \
                -bed "${PREFIX}.250bp.bed"
        fi

    done

    # Get non-overlapping, 250bp-long, random intervals
    if [ ! -f "${GENOME}-random.nov.250bp.bed" ]; then
        bedtools random -l 250 -seed ${RANDOM_SEED} -g "${GENOME_FA}.sizes" \
            | bedtools sort > "${GENOME}-random.250bp.bed"
        bedtools intersect -a "${GENOME}-random.250bp.bed" \
            -b "${GENOME}-ccREs.bed" -v > "${GENOME}-random.nov.250bp.bed.tmp"
        bedtools closest -t first -a "${GENOME}-random.nov.250bp.bed.tmp" \
            -b "${GENOME}-random.nov.250bp.bed.tmp" | cut -f 1-3 \
            | bedtools sort | uniq > "${GENOME}-random.nov.250bp.bed"
        rm "${GENOME}-random.nov.250bp.bed.tmp"
    fi

    # Get FASTA
    if [ ! -f "${GENOME}-random.nov.250bp.fa" ]; then
        bedtools getfasta -fo "${GENOME}-random.nov.250bp.fa" \
            -fi ${GENOME_FA} -bed "${GENOME}-random.nov.250bp.bed"
    fi

done
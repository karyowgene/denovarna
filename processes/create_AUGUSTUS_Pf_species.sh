#!/bin/bash

# This script trains Augustus on VSA genes from Plasmodium falciparum
# using user-specified FASTA, GFF, and species name.

set -euo pipefail

# Check for the correct number of arguments
if [[ "$#" -ne 3 ]]; then
    echo "Usage: $0 <FASTA_FILE> <GFF_FILE> <SPECIES_NAME>"
    exit 1
fi

# Assign input arguments to variables
FASTA="$1"
GFF="$2"
SPECIES_NAME="$3"

# Check input files
if [[ ! -f "$FASTA" || ! -f "$GFF" ]]; then
    echo "ERROR: Required FASTA or GFF file not found."
    exit 1
fi

# Create and move into a working directory
mkdir -p augustus_training
cd augustus_training

echo "Extracting VSA gene IDs..."
awk -F '\t' '($3 == "protein_coding_gene" || $3 == "pseudogene") && tolower($9) ~ /pfemp1|erythrocyte membrane protein 1|stevor|rifin/ {
    split($9, a, ";");
    for (i in a) {
        if (a[i] ~ /^ID=/) {
            split(a[i], b, "=");
            print b[2];
        }
    }
}' "$GFF" | sort -u > vsa_gene_ids.txt

echo "Filtering GFF for VSA gene models..."
grep "^#" "$GFF" > ${SPECIES_NAME}_vsa.gff
while read geneid; do
    grep -P "$geneid([.;]|$)" "$GFF"
done < vsa_gene_ids.txt >> ${SPECIES_NAME}_vsa.gff

echo "Converting GFF to GenBank format..."
gff2gbSmallDNA.pl ${SPECIES_NAME}_vsa.gff "$FASTA" 1000 ${SPECIES_NAME}.gb

echo "Splitting training and test sets..."
randomSplit.pl ${SPECIES_NAME}.gb 100

echo "Creating new Augustus species: $SPECIES_NAME"
new_species.pl --species=$SPECIES_NAME

echo "Running etraining..."
etraining --species=$SPECIES_NAME ${SPECIES_NAME}.gb.train

echo "Testing Augustus on test set..."
augustus --species=$SPECIES_NAME ${SPECIES_NAME}.gb.test > ${SPECIES_NAME}.test.prediction.gff

echo "Augustus training and testing complete."

#!/bin/bash

CLASSIFICATION="$1"

if [ -z "$CLASSIFICATION" ]; then
    echo "Usage:"
    echo "./scripts/genomeguard_analyze.sh <classification.tsv>"
    exit 1
fi

if [ ! -f "$CLASSIFICATION" ]; then
    echo "ERROR: File not found:"
    echo "$CLASSIFICATION"
    exit 1
fi

echo "=================================================="
echo "              GENOMEGUARD v0.1"
echo "          TAXONOMIC INTEGRITY ANALYSIS"
echo "=================================================="

echo
echo "Input:"
echo "$CLASSIFICATION"

echo
echo "### TOTAL CLASSIFICATION RECORDS"

TOTAL=$(awk 'NR>1 {count++} END {print count+0}' "$CLASSIFICATION")
echo "$TOTAL"

echo
echo "### UNCLASSIFIED READS"

UNCLASSIFIED=$(awk 'NR>1 && $3==0 {count++} END {print count+0}' "$CLASSIFICATION")
echo "$UNCLASSIFIED"

echo
echo "### CLASSIFIED READS"

CLASSIFIED=$(awk 'NR>1 && $3!=0 {count++} END {print count+0}' "$CLASSIFICATION")
echo "$CLASSIFIED"

echo
echo "### CLASSIFICATION RATE"

awk -v c="$CLASSIFIED" -v t="$TOTAL" \
'BEGIN {
    if(t>0)
        printf "%.4f%%\n", (c/t)*100
    else
        print "0%"
}'

echo
echo "### TOP 20 TAXA"

awk '
NR>1 && $3!=0 {
    count[$3]++
}
END {
    for(tax in count)
        print count[tax] "\t" tax
}' "$CLASSIFICATION" |
sort -k1,1nr |
head -20

echo
echo "=================================================="
echo "              GENOMEGUARD COMPLETE"
echo "=================================================="

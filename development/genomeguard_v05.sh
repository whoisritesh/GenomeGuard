#!/bin/bash

CLASSIFICATION="$1"
REPORT="$2"

if [ -z "$CLASSIFICATION" ] || [ -z "$REPORT" ]; then
    echo "Usage:"
    echo "./scripts/genomeguard_v05.sh classification.tsv report.tsv"
    exit 1
fi

if [ ! -f "$CLASSIFICATION" ]; then
    echo "ERROR: Classification file not found:"
    echo "$CLASSIFICATION"
    exit 1
fi

if [ ! -f "$REPORT" ]; then
    echo "ERROR: Report file not found:"
    echo "$REPORT"
    exit 1
fi

echo "=================================================="
echo "          GENOMEGUARD v0.5"
echo "       SAMPLE INTEGRITY ENGINE"
echo "=================================================="

# --------------------------------------------------
# 1. CLASSIFICATION SUMMARY
# --------------------------------------------------

TOTAL=$(awk 'NR>1 {n++} END {print n+0}' "$CLASSIFICATION")

UNCLASSIFIED=$(awk '
NR>1 && $3==0 {n++}
END {print n+0}
' "$CLASSIFICATION")

CLASSIFIED=$((TOTAL - UNCLASSIFIED))

if [ "$TOTAL" -gt 0 ]; then
    CLASS_RATE=$(awk -v c="$CLASSIFIED" -v t="$TOTAL" \
    'BEGIN {printf "%.4f", (c/t)*100}')

    UNCLASS_RATE=$(awk -v u="$UNCLASSIFIED" -v t="$TOTAL" \
    'BEGIN {printf "%.4f", (u/t)*100}')
else
    CLASS_RATE="0.0000"
    UNCLASS_RATE="0.0000"
fi

echo
echo "### READ-LEVEL SUMMARY"
echo "Total records       : $TOTAL"
echo "Classified records  : $CLASSIFIED"
echo "Unclassified        : $UNCLASSIFIED"
echo "Classification rate  : ${CLASS_RATE}%"
echo "Unclassified rate    : ${UNCLASS_RATE}%"

# --------------------------------------------------
# 2. TAXONOMIC DISTRIBUTION
# --------------------------------------------------

echo
echo "### TAXONOMIC DISTRIBUTION"

awk '
NR>1 {
    tax=$3
    count[tax]++
}
END {
    for(tax in count)
        print count[tax] "\t" tax
}
' "$CLASSIFICATION" |
sort -k1,1nr > /tmp/genomeguard_v05_taxa.txt

DOM_READS=$(awk 'NR==1 {print $1}' /tmp/genomeguard_v05_taxa.txt)
DOM_TAXID=$(awk 'NR==1 {print $2}' /tmp/genomeguard_v05_taxa.txt)

SECOND_READS=$(awk 'NR==2 {print $1}' /tmp/genomeguard_v05_taxa.txt)
SECOND_TAXID=$(awk 'NR==2 {print $2}' /tmp/genomeguard_v05_taxa.txt)

echo "Dominant taxID      : $DOM_TAXID"
echo "Dominant reads      : $DOM_READS"
echo "Second taxID        : $SECOND_TAXID"
echo "Second reads        : $SECOND_READS"

# --------------------------------------------------
# 3. DOMINANCE
# --------------------------------------------------

if [ "$CLASSIFIED" -gt 0 ]; then

    DOMINANCE=$(awk \
    -v d="$DOM_READS" \
    -v c="$CLASSIFIED" \
    'BEGIN {printf "%.4f", (d/c)*100}')

    SECOND_RATIO=$(awk \
    -v s="$SECOND_READS" \
    -v d="$DOM_READS" \
    'BEGIN {
        if(d>0) printf "%.6f", s/d;
        else print "0"
    }')

else
    DOMINANCE="0.0000"
    SECOND_RATIO="0.000000"
fi

echo
echo "### DOMINANCE METRICS"
echo "Dominance           : ${DOMINANCE}%"
echo "Secondary/Dominant  : $SECOND_RATIO"

# --------------------------------------------------
# 4. TOP 15 TAXA
# --------------------------------------------------

echo
echo "### TOP 15 TAXA"

head -15 /tmp/genomeguard_v05_taxa.txt

# --------------------------------------------------
# 5. SECONDARY BURDEN
# --------------------------------------------------

if [ "$CLASSIFIED" -gt 0 ]; then

    SECONDARY_READS=$((CLASSIFIED - DOM_READS))

    SECONDARY_BURDEN=$(awk \
    -v s="$SECONDARY_READS" \
    -v c="$CLASSIFIED" \
    'BEGIN {printf "%.4f", (s/c)*100}')

else
    SECONDARY_READS=0
    SECONDARY_BURDEN="0.0000"
fi

echo
echo "### SECONDARY SIGNAL"
echo "Secondary reads     : $SECONDARY_READS"
echo "Secondary burden    : ${SECONDARY_BURDEN}%"

# --------------------------------------------------
# 6. SIMPLE INTEGRITY CLASS
# --------------------------------------------------

INTEGRITY="REVIEW"

awk -v d="$DOMINANCE" -v s="$SECONDARY_BURDEN" '
BEGIN {
    if(d >= 90 && s <= 10)
        exit 10
    else if(d >= 70 && s <= 30)
        exit 20
    else
        exit 30
}
'

STATUS=$?

if [ "$STATUS" -eq 10 ]; then
    INTEGRITY="HIGH DOMINANCE"
elif [ "$STATUS" -eq 20 ]; then
    INTEGRITY="MODERATE DOMINANCE"
else
    INTEGRITY="COMPLEX / REVIEW"
fi

echo
echo "### SAMPLE INTEGRITY"
echo "Integrity category  : $INTEGRITY"

# --------------------------------------------------
# 7. FINAL INTERPRETATION
# --------------------------------------------------

echo
echo "=================================================="
echo "              GENOMEGUARD v0.5"
echo "             FINAL INTERPRETATION"
echo "=================================================="

echo
echo "Dominant taxID: $DOM_TAXID"
echo "Dominant reads: $DOM_READS"
echo "Dominance: ${DOMINANCE}%"
echo "Secondary burden: ${SECONDARY_BURDEN}%"

echo
echo "Interpretation:"
echo

if [ "$STATUS" -eq 10 ]; then

    echo "The sample shows strong dominance by a single taxonomic signal."
    echo "Secondary signals are comparatively limited."
    echo "No contamination conclusion should be made from dominance alone."

elif [ "$STATUS" -eq 20 ]; then

    echo "The sample shows a dominant taxonomic signal."
    echo "However, substantial secondary taxonomic evidence is present."
    echo "The secondary signals should be investigated before assigning"
    echo "a simple single-organism interpretation."

else

    echo "The sample contains substantial taxonomic complexity."
    echo "Further investigation is recommended before interpreting"
    echo "the sample as a simple single-organism dataset."

fi

echo
echo "=================================================="
echo "GenomeGuard v0.5 complete"
echo "=================================================="

rm -f /tmp/genomeguard_v05_taxa.txt


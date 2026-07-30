#!/bin/bash

CLASSIFICATION="$1"
REPORT="$2"

if [ -z "$CLASSIFICATION" ] || [ -z "$REPORT" ]; then
    echo "Usage:"
    echo "./scripts/genomeguard_v05_1.sh classification.tsv report.tsv"
    exit 1
fi

if [ ! -f "$CLASSIFICATION" ]; then
    echo "ERROR: Classification file not found"
    exit 1
fi

if [ ! -f "$REPORT" ]; then
    echo "ERROR: Report file not found"
    exit 1
fi

echo "=================================================="
echo "         GENOMEGUARD v0.5.1"
echo "       SAMPLE INTEGRITY ENGINE"
echo "          QC-CORRECTED VERSION"
echo "=================================================="

# ==================================================
# 1. READ SUMMARY
# ==================================================

TOTAL=$(awk 'NR>1 {n++} END {print n+0}' "$CLASSIFICATION")

UNCLASSIFIED=$(awk '
NR>1 && $3==0 {n++}
END {print n+0}
' "$CLASSIFICATION")

CLASSIFIED=$((TOTAL - UNCLASSIFIED))

CLASS_RATE=$(awk -v c="$CLASSIFIED" -v t="$TOTAL" \
'BEGIN {
    if(t>0) printf "%.4f",100*c/t;
    else print "0.0000"
}')

UNCLASS_RATE=$(awk -v u="$UNCLASSIFIED" -v t="$TOTAL" \
'BEGIN {
    if(t>0) printf "%.4f",100*u/t;
    else print "0.0000"
}')

echo
echo "### READ-LEVEL SUMMARY"
echo "Total records       : $TOTAL"
echo "Classified records  : $CLASSIFIED"
echo "Unclassified        : $UNCLASSIFIED"
echo "Classification rate : ${CLASS_RATE}%"
echo "Unclassified rate   : ${UNCLASS_RATE}%"

# ==================================================
# 2. CLASSIFIED TAXONOMIC DISTRIBUTION
#    TAXID 0 EXCLUDED
# ==================================================

awk '
NR>1 && $3!=0 {
    count[$3]++
}
END {
    for(tax in count)
        print count[tax] "\t" tax
}
' "$CLASSIFICATION" |
sort -k1,1nr > /tmp/genomeguard_v051_taxa.txt

DOM_READS=$(awk 'NR==1 {print $1}' /tmp/genomeguard_v051_taxa.txt)
DOM_TAXID=$(awk 'NR==1 {print $2}' /tmp/genomeguard_v051_taxa.txt)

SECOND_READS=$(awk 'NR==2 {print $1}' /tmp/genomeguard_v051_taxa.txt)
SECOND_TAXID=$(awk 'NR==2 {print $2}' /tmp/genomeguard_v051_taxa.txt)

echo
echo "### CLASSIFIED TAXONOMIC DISTRIBUTION"
echo "Dominant taxID     : $DOM_TAXID"
echo "Dominant reads     : $DOM_READS"
echo "Second taxID       : $SECOND_TAXID"
echo "Second reads       : $SECOND_READS"

# ==================================================
# 3. DOMINANCE
# ==================================================

DOMINANCE=$(awk \
-v d="$DOM_READS" \
-v c="$CLASSIFIED" \
'BEGIN {
    if(c>0) printf "%.4f",100*d/c;
    else print "0.0000"
}')

SECOND_RATIO=$(awk \
-v s="$SECOND_READS" \
-v d="$DOM_READS" \
'BEGIN {
    if(d>0) printf "%.6f",s/d;
    else print "0.000000"
}')

SECONDARY_READS=$((CLASSIFIED - DOM_READS))

SECONDARY_BURDEN=$(awk \
-v s="$SECONDARY_READS" \
-v c="$CLASSIFIED" \
'BEGIN {
    if(c>0) printf "%.4f",100*s/c;
    else print "0.0000"
}')

echo
echo "### DOMINANCE METRICS"
echo "Dominance           : ${DOMINANCE}%"
echo "Secondary/Dominant  : $SECOND_RATIO"
echo "Secondary reads     : $SECONDARY_READS"
echo "Secondary burden    : ${SECONDARY_BURDEN}%"

# ==================================================
# 4. TOP 15 CLASSIFIED TAXA
# ==================================================

echo
echo "### TOP 15 CLASSIFIED TAXA"
head -15 /tmp/genomeguard_v051_taxa.txt

# ==================================================
# 5. TAXONOMIC COMPLEXITY
# ==================================================

UNIQUE_TAXA=$(awk '
NR>1 && $3!=0 {
    taxa[$3]=1
}
END {
    n=0
    for(t in taxa) n++
    print n
}
' "$CLASSIFICATION")

echo
echo "### TAXONOMIC COMPLEXITY"
echo "Unique classified taxa : $UNIQUE_TAXA"

# ==================================================
# 6. SECONDARY COMMUNITY
# ==================================================

echo
echo "### SECONDARY COMMUNITY"

awk -v dominant="$DOM_TAXID" '
NR>1 && $3!=0 && $3!=dominant {
    count[$3]++
}
END {
    for(tax in count)
        print count[tax] "\t" tax
}
' "$CLASSIFICATION" |
sort -k1,1nr |
head -15

# ==================================================
# 7. RELATIVE CONTRIBUTION OF TOP TAXA
# ==================================================

echo
echo "### TOP TAXA RELATIVE CONTRIBUTION"

awk -v total="$CLASSIFIED" '
NR<=15 {
    pct=($1/total)*100
    printf "%-10s %-12s %.4f%%\n",$2,$1,pct
}
' /tmp/genomeguard_v051_taxa.txt

# ==================================================
# 8. SAMPLE INTEGRITY
# ==================================================

if awk -v d="$DOMINANCE" -v s="$SECONDARY_BURDEN" '
BEGIN {
    exit !(d >= 90 && s <= 10)
}'; then

    INTEGRITY="HIGH DOMINANCE"

elif awk -v d="$DOMINANCE" -v s="$SECONDARY_BURDEN" '
BEGIN {
    exit !(d >= 70 && s <= 30)
}'; then

    INTEGRITY="MODERATE DOMINANCE"

else

    INTEGRITY="COMPLEX / REVIEW"

fi

echo
echo "### SAMPLE INTEGRITY"
echo "Integrity category : $INTEGRITY"

# ==================================================
# 9. QC FLAGS
# ==================================================

echo
echo "### QUALITY-CONTROL FLAGS"

if awk -v x="$UNCLASS_RATE" 'BEGIN {exit !(x > 10)}'; then
    echo "WARNING: High unclassified fraction"
else
    echo "PASS: Unclassified fraction not excessive"
fi

if awk -v x="$DOMINANCE" 'BEGIN {exit !(x >= 70)}'; then
    echo "PASS: Strong dominant taxonomic signal"
else
    echo "REVIEW: Dominant taxonomic signal is weak"
fi

if awk -v x="$SECONDARY_BURDEN" 'BEGIN {exit !(x >= 20)}'; then
    echo "FLAG: Substantial secondary taxonomic burden"
else
    echo "PASS: Secondary burden relatively limited"
fi

if awk -v x="$UNIQUE_TAXA" 'BEGIN {exit !(x > 100)}'; then
    echo "FLAG: High taxonomic complexity"
else
    echo "INFO: Taxonomic complexity within current threshold"
fi

# ==================================================
# 10. FINAL INTERPRETATION
# ==================================================

echo
echo "=================================================="
echo "       GENOMEGUARD v0.5.1 INTERPRETATION"
echo "=================================================="

echo
echo "Primary taxonomic signal : TaxID $DOM_TAXID"
echo "Primary read fraction    : ${DOMINANCE}%"
echo "Secondary burden         : ${SECONDARY_BURDEN}%"
echo "Unclassified fraction    : ${UNCLASS_RATE}%"
echo "Unique classified taxa  : $UNIQUE_TAXA"
echo "Integrity category       : $INTEGRITY"

echo
echo "Interpretation:"

if [ "$INTEGRITY" = "HIGH DOMINANCE" ]; then

    echo "The classified read population is strongly dominated"
    echo "by a single taxonomic signal."

elif [ "$INTEGRITY" = "MODERATE DOMINANCE" ]; then

    echo "The sample has a clear dominant taxonomic signal,"
    echo "but a substantial secondary classified population"
    echo "is also present."

else

    echo "The classified read population shows substantial"
    echo "taxonomic complexity and should be investigated"
    echo "before assigning a simple single-organism interpretation."

fi

echo
echo "IMPORTANT:"
echo "Secondary taxonomic evidence does not by itself prove"
echo "contamination, co-infection, or a mixed biological sample."
echo "Additional evidence is required for such conclusions."

echo
echo "=================================================="
echo "GenomeGuard v0.5.1 complete"
echo "=================================================="

rm -f /tmp/genomeguard_v051_taxa.txt


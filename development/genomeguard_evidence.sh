#!/bin/bash

CLASSIFICATION="$1"
REPORT="$2"

if [ -z "$CLASSIFICATION" ] || [ -z "$REPORT" ]; then
    echo "Usage:"
    echo "./scripts/genomeguard_evidence.sh <classification.tsv> <report.tsv>"
    exit 1
fi

if [ ! -f "$CLASSIFICATION" ]; then
    echo "ERROR: Classification file not found."
    exit 1
fi

if [ ! -f "$REPORT" ]; then
    echo "ERROR: Centrifuge report not found."
    exit 1
fi

echo "=================================================="
echo "          GENOMEGUARD v0.2"
echo "           EVIDENCE ENGINE"
echo "=================================================="

TOTAL=$(awk 'NR>1 {n++} END {print n+0}' "$CLASSIFICATION")

CLASSIFIED=$(awk '
NR>1 && $3 != 0 {n++}
END {print n+0}
' "$CLASSIFICATION")

UNCLASSIFIED=$(awk '
NR>1 && $3 == 0 {n++}
END {print n+0}
' "$CLASSIFICATION")

echo
echo "### READ-LEVEL SUMMARY"
echo "Total classification records : $TOTAL"
echo "Classified records           : $CLASSIFIED"
echo "Unclassified records         : $UNCLASSIFIED"

CLASS_RATE=$(awk -v c="$CLASSIFIED" -v t="$TOTAL" \
'BEGIN {
    if(t>0) printf "%.4f", (c/t)*100;
    else print "0"
}')

UNCLASS_RATE=$(awk -v u="$UNCLASSIFIED" -v t="$TOTAL" \
'BEGIN {
    if(t>0) printf "%.4f", (u/t)*100;
    else print "0"
}')

echo "Classification rate          : ${CLASS_RATE}%"
echo "Unclassified fraction        : ${UNCLASS_RATE}%"

echo
echo "### TAXONOMIC READ DISTRIBUTION"

awk '
NR>1 && $3 != 0 {
    count[$3]++
}
END {
    for(tax in count)
        print count[tax] "\t" tax
}' "$CLASSIFICATION" |
sort -k1,1nr > /tmp/genomeguard_taxa.txt

DOMINANT_READS=$(head -1 /tmp/genomeguard_taxa.txt | awk '{print $1}')
DOMINANT_TAXID=$(head -1 /tmp/genomeguard_taxa.txt | awk '{print $2}')

SECOND_READS=$(sed -n '2p' /tmp/genomeguard_taxa.txt | awk '{print $1}')
SECOND_TAXID=$(sed -n '2p' /tmp/genomeguard_taxa.txt | awk '{print $2}')

echo "Dominant taxID                : $DOMINANT_TAXID"
echo "Dominant reads                : $DOMINANT_READS"
echo "Second taxID                  : $SECOND_TAXID"
echo "Second taxon reads            : $SECOND_READS"

DOMINANCE=$(awk -v d="$DOMINANT_READS" -v c="$CLASSIFIED" \
'BEGIN {
    if(c>0) printf "%.4f", (d/c)*100;
    else print "0"
}')

SECONDARY_RATIO=$(awk -v s="$SECOND_READS" -v d="$DOMINANT_READS" \
'BEGIN {
    if(d>0) printf "%.6f", s/d;
    else print "0"
}')

echo
echo "### EVIDENCE METRICS"

echo "Dominance                     : ${DOMINANCE}%"
echo "Secondary/Dominant ratio      : $SECONDARY_RATIO"

echo
echo "### TOP 10 TAXA"

head -10 /tmp/genomeguard_taxa.txt

echo
echo "### UNIQUE-READ SUPPORT"

awk '
NR>1 {
    for(i=1;i<=NF;i++) {
        if($i=="species") {

            reads=$(i+2)
            unique=$(i+3)
            abundance=$(i+4)

            name=""

            for(j=1;j<i;j++) {
                name=name (j>1 ? " " : "") $j
            }

            if(reads>0) {
                ratio=unique/reads
                printf "%.6f\t%d\t%d\t%s\n",
                       ratio, reads, unique, name
            }

            break
        }
    }
}' "$REPORT" |
sort -k1,1nr |
head -20

echo
echo "=================================================="
echo "              PRELIMINARY INTERPRETATION"
echo "=================================================="

if awk -v d="$DOMINANCE" 'BEGIN {exit !(d >= 90)}'; then
    echo "Dominant signal: STRONG"
elif awk -v d="$DOMINANCE" 'BEGIN {exit !(d >= 70)}'; then
    echo "Dominant signal: MODERATE"
else
    echo "Dominant signal: WEAK"
fi

if awk -v r="$SECONDARY_RATIO" 'BEGIN {exit !(r >= 0.05)}'; then
    echo "Secondary signal: SUBSTANTIAL"
else
    echo "Secondary signal: LOW"
fi

echo
echo "GenomeGuard interpretation:"
echo "Taxonomic evidence should be reviewed before"
echo "calling the sample contaminated or mixed."

echo "=================================================="

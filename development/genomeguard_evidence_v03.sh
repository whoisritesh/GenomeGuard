#!/bin/bash

CLASSIFICATION="$1"
REPORT="$2"

if [ -z "$CLASSIFICATION" ] || [ -z "$REPORT" ]; then
    echo "Usage:"
    echo "./scripts/genomeguard_evidence_v03.sh <classification.tsv> <report.tsv>"
    exit 1
fi

if [ ! -f "$CLASSIFICATION" ] || [ ! -f "$REPORT" ]; then
    echo "ERROR: Input file not found."
    exit 1
fi

echo "=================================================="
echo "             GENOMEGUARD v0.3"
echo "          EVIDENCE-WEIGHTED ANALYSIS"
echo "=================================================="

CLASSIFIED=$(awk 'NR>1 && $3!=0 {n++} END{print n+0}' "$CLASSIFICATION")

echo
echo "Classified records: $CLASSIFIED"

echo
echo "### TAXONOMIC EVIDENCE TABLE"
printf "%-12s %-12s %-12s %-12s %-12s %s\n" \
"TaxID" "Reads" "Unique" "Unique%" "Relative%" "Evidence"

awk -v total="$CLASSIFIED" '
NR>1 {
    for(i=1;i<=NF;i++) {

        if($i=="species") {

            name=""

            for(j=1;j<i;j++) {
                name=name (j>1 ? " " : "") $j
            }

            reads=$(i+2)
            unique=$(i+3)

            if(reads>0) {

                unique_ratio=(unique/reads)*100
                relative=(reads/total)*100

                if(reads>=10000 && unique>=1000)
                    evidence="HIGH"
                else if(reads>=1000 && unique>=100)
                    evidence="MODERATE"
                else if(reads>=100 && unique>=10)
                    evidence="LOW"
                else
                    evidence="TRACE"

                printf "%-12s %-12d %-12d %-12.3f %-12.4f %s\t%s\n",
                $(i+1), reads, unique, unique_ratio, relative, evidence, name
            }

            break
        }
    }
}' "$REPORT" |
sort -k5,5nr |
head -30

echo
echo "=================================================="
echo "Evidence categories:"
echo
echo "HIGH     = strong read + unique-read support"
echo "MODERATE = substantial support"
echo "LOW      = limited support"
echo "TRACE    = very small signal"
echo "=================================================="

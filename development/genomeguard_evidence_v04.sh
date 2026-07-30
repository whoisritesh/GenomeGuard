#!/bin/bash

CLASSIFICATION="$1"
REPORT="$2"

if [ -z "$CLASSIFICATION" ] || [ -z "$REPORT" ]; then
    echo "Usage:"
    echo "./scripts/genomeguard_evidence_v04.sh <classification.tsv> <report.tsv>"
    exit 1
fi

if [ ! -f "$CLASSIFICATION" ] || [ ! -f "$REPORT" ]; then
    echo "ERROR: Input file not found."
    exit 1
fi

echo "=================================================="
echo "             GENOMEGUARD v0.4"
echo "       NORMALIZED TAXONOMIC EVIDENCE"
echo "=================================================="

CLASSIFIED=$(awk '
NR>1 && $3 != 0 {n++}
END {print n+0}
' "$CLASSIFICATION")

echo
echo "Classified records: $CLASSIFIED"

echo
echo "### NORMALIZED SPECIES EVIDENCE"

printf "%-12s %-12s %-12s %-12s %-12s %-10s %s\n" \
"TAXID" "READS" "UNIQUE" "UNIQUE%" "RELATIVE%" "EVIDENCE" "SPECIES"

awk -v total="$CLASSIFIED" '
NR>1 {

    for(i=1;i<=NF;i++) {

        if($i=="species") {

            name=""

            for(j=1;j<i;j++) {
                name=name (j>1 ? " " : "") $j
            }

            taxid=$(i+1)
            reads=$(i+4)
            unique=$(i+5)

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

                printf "%-12s %-12d %-12d %-12.3f %-12.4f %-10s %s\n",
                taxid,
                reads,
                unique,
                unique_ratio,
                relative,
                evidence,
                name
            }

            break
        }
    }
}' "$REPORT" |
sort -k5,5nr |
head -30

echo
echo "=================================================="
echo "### STRONG SECONDARY SIGNALS"
echo "=================================================="

awk -v total="$CLASSIFIED" '
NR>1 {

    for(i=1;i<=NF;i++) {

        if($i=="species") {

            name=""

            for(j=1;j<i;j++) {
                name=name (j>1 ? " " : "") $j
            }

            taxid=$(i+1)
            reads=$(i+4)
            unique=$(i+5)

            relative=(reads/total)*100

            if(reads>=1000 && unique>=100)
                printf "%10d %10d %10.4f %s [%s]\n",
                reads,
                unique,
                relative,
                name,
                taxid

            break
        }
    }
}' "$REPORT" |
sort -k1,1nr

echo
echo "=================================================="
echo "GenomeGuard v0.4 complete"
echo "=================================================="

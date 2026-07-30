#!/bin/bash

CLASS="$1"
REPORT="$2"

if [ ! -f "$CLASS" ] || [ ! -f "$REPORT" ]; then
    echo "ERROR: Input file missing"
    exit 1
fi

echo "=================================================="
echo "          GENOMEGUARD v0.7"
echo "      READ-LEVEL EVIDENCE VALIDATION"
echo "=================================================="

TOTAL=$(awk 'NR>1 {n++} END{print n+0}' "$CLASS")
UNCLASSIFIED=$(awk 'NR>1 && $3==0 {n++} END{print n+0}' "$CLASS")
CLASSIFIED=$((TOTAL-UNCLASSIFIED))

echo
echo "### READ-LEVEL SUMMARY"
echo "Total records       : $TOTAL"
echo "Classified records  : $CLASSIFIED"
echo "Unclassified        : $UNCLASSIFIED"

RATE=$(awk -v c="$CLASSIFIED" -v t="$TOTAL" 'BEGIN{printf "%.4f",100*c/t}')
echo "Classification rate : ${RATE}%"

echo
echo "### TOP TAXA BY READ-LEVEL ASSIGNMENT"

awk '
NR>1 && $3!=0 {
    count[$3]++
}
END {
    for(t in count)
        print count[t] "\t" t
}' "$CLASS" | sort -k1,1nr | head -20

echo
echo "### READ-LEVEL SCORE SUMMARY"

awk '
NR>1 && $3!=0 {
    score=$4
    total++
    sum+=score

    if(score>0) positive++
    if(score>=10000) high++
    if(score>=20000) veryhigh++
}
END {
    if(total>0) {
        printf "Mean score          : %.2f\n",sum/total
        printf "Positive-score reads: %d\n",positive
        printf "Score >=10000       : %d\n",high
        printf "Score >=20000       : %d\n",veryhigh
    }
}' "$CLASS"

echo
echo "### TOP SECONDARY TAXA WITH READ SUPPORT"

awk '
NR>1 && $3!=0 {
    count[$3]++
}
END {
    for(t in count)
        print count[t] "\t" t
}' "$CLASS" |
sort -k1,1nr |
awk '$2 != 562' |
head -15

echo
echo "### READ-LEVEL VALIDATION OF MAJOR SECONDARY SIGNALS"

for TAXID in 1813821 624 621 623 622 564 28901 208962 573
do
    READS=$(awk -v t="$TAXID" 'NR>1 && $3==t {n++} END{print n+0}' "$CLASS")

    if [ "$READS" -gt 0 ]; then
        echo
        echo "TaxID: $TAXID"
        echo "Assigned reads: $READS"

        awk -v t="$TAXID" '
        NR>1 && $3==t {
            n++
            sum+=$4

            if($4>0) positive++
            if($4>=10000) high++
            if($4>=20000) veryhigh++

            if($5>0) {
                delta=$4-$5
                dsum+=delta
            }
        }
        END {
            if(n>0) {
                printf "Mean score       : %.2f\n",sum/n
                printf "Positive score   : %.2f%%\n",100*positive/n
                printf "Score >=10000    : %.2f%%\n",100*high/n
                printf "Score >=20000    : %.2f%%\n",100*veryhigh/n

                if(dsum>0)
                    printf "Mean score gap   : %.2f\n",dsum/n
            }
        }' "$CLASS"
    fi
done

echo
echo "=================================================="
echo "          GENOMEGUARD v0.7 INTERPRETATION"
echo "=================================================="

echo
echo "Primary signal: TaxID 562 (Escherichia coli)"
echo
echo "v0.7 evaluates whether secondary taxonomic signals"
echo "are supported at the individual-read level."
echo
echo "IMPORTANT:"
echo "Read-level taxonomic assignments are evidence signals."
echo "They do not independently establish contamination,"
echo "co-infection, or mixed biological origin."
echo
echo "=================================================="
echo "        GenomeGuard v0.7 complete"
echo "=================================================="

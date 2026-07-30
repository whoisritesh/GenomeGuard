#!/bin/bash

CLASSIFICATION="$1"
REPORT="$2"

if [ -z "$CLASSIFICATION" ] || [ -z "$REPORT" ]; then
    echo "Usage:"
    echo "./scripts/genomeguard_v06.sh classification.tsv report.tsv"
    exit 1
fi

echo "=================================================="
echo "          GENOMEGUARD v0.6"
echo "    SECONDARY SIGNAL INVESTIGATION ENGINE"
echo "=================================================="

TOTAL=$(awk 'NR>1{n++}END{print n+0}' "$CLASSIFICATION")

CLASSIFIED=$(awk 'NR>1 && $3 != 0{n++}END{print n+0}' "$CLASSIFICATION")

UNCLASSIFIED=$(awk 'NR>1 && $3 == 0{n++}END{print n+0}' "$CLASSIFICATION")

RATE=$(awk -v c="$CLASSIFIED" -v t="$TOTAL" \
'BEGIN{if(t>0) printf "%.4f",100*c/t; else print "0"}')

echo
echo "### READ-LEVEL SUMMARY"
echo "Total records       : $TOTAL"
echo "Classified records  : $CLASSIFIED"
echo "Unclassified        : $UNCLASSIFIED"
echo "Classification rate : $RATE%"

echo
echo "### TAXONOMIC SIGNAL EXTRACTION"

awk '
NR>1 && $3 != 0 {
    count[$3]++
}
END {
    for(t in count)
        print count[t] "\t" t
}' "$CLASSIFICATION" |
sort -k1,1nr > /tmp/genomeguard_v06_taxa.txt

DOM_READS=$(awk 'NR==1{print $1}' /tmp/genomeguard_v06_taxa.txt)
DOM_TAXID=$(awk 'NR==1{print $2}' /tmp/genomeguard_v06_taxa.txt)

echo "Dominant taxID      : $DOM_TAXID"
echo "Dominant reads      : $DOM_READS"

DOM_PCT=$(awk -v d="$DOM_READS" -v c="$CLASSIFIED" \
'BEGIN{if(c>0) printf "%.4f",100*d/c; else print "0"}')

echo "Dominant fraction   : $DOM_PCT%"

echo
echo "### SECONDARY SIGNALS"

awk -v dominant="$DOM_TAXID" '
NR>1 && $3 != 0 && $3 != dominant {
    count[$3]++
}
END {
    for(t in count)
        print count[t] "\t" t
}' "$CLASSIFICATION" |
sort -k1,1nr |
head -20

echo
echo "### SECONDARY SIGNAL SCORE"

awk -v classified="$CLASSIFIED" -v dominant="$DOM_TAXID" '
NR>1 && $3 != 0 && $3 != dominant {
    count[$3]++
}
END {
    for(t in count) {
        reads=count[t]
        relative=(reads/classified)*100

        if(reads >= 10000 && relative >= 1)
            tier="STRONG SECONDARY"
        else if(reads >= 1000 && relative >= 0.1)
            tier="MODERATE SECONDARY"
        else if(reads >= 100 && relative >= 0.01)
            tier="LOW SIGNAL"
        else
            tier="TRACE"

        printf "%-12s %-12d %-12.4f %-20s\n",t,reads,relative,tier
    }
}' "$CLASSIFICATION" |
sort -k3,3nr |
head -30

echo
echo "### SECONDARY BURDEN"

SECONDARY=$(awk -v dominant="$DOM_TAXID" '
NR>1 && $3 != 0 && $3 != dominant {n++}
END{print n+0}
' "$CLASSIFICATION")

SECONDARY_PCT=$(awk -v s="$SECONDARY" -v c="$CLASSIFIED" \
'BEGIN{if(c>0) printf "%.4f",100*s/c; else print "0"}')

echo "Secondary records  : $SECONDARY"
echo "Secondary burden   : $SECONDARY_PCT%"

echo
echo "### TAXONOMIC COMPLEXITY"

UNIQUE_TAXA=$(awk 'NR>1 && $3 != 0{print $3}' "$CLASSIFICATION" |
sort -u | wc -l)

echo "Unique classified taxa : $UNIQUE_TAXA"

echo
echo "### INVESTIGATION CANDIDATES"

awk -v classified="$CLASSIFIED" -v dominant="$DOM_TAXID" '
NR>1 && $3 != 0 && $3 != dominant {
    count[$3]++
}
END {
    for(t in count) {
        reads=count[t]
        relative=(reads/classified)*100

        if(reads >= 1000 && relative >= 0.1)
            print reads "\t" relative "\t" t
    }
}' "$CLASSIFICATION" |
sort -k2,2nr |
head -15

echo
echo "=================================================="
echo "             GENOMEGUARD v0.6"
echo "         INVESTIGATION INTERPRETATION"
echo "=================================================="

echo
echo "Primary signal : TaxID $DOM_TAXID"
echo "Primary fraction: $DOM_PCT%"
echo "Secondary burden: $SECONDARY_PCT%"
echo "Taxonomic complexity: $UNIQUE_TAXA taxa"

echo
echo "Interpretation:"
echo
echo "GenomeGuard identifies a dominant taxonomic signal"
echo "together with secondary classified populations."
echo
echo "Secondary signals are ranked for investigation based"
echo "on read support and relative contribution."
echo
echo "These classifications are evidence signals and do not"
echo "independently establish contamination, co-infection,"
echo "or mixed biological origin."

echo
echo "=================================================="
echo "        GenomeGuard v0.6 complete"
echo "=================================================="

#!/usr/bin/env bash

set -euo pipefail

############################################################
# GenomeGuard v1.1
# Final Sample Integrity Analysis Pipeline
############################################################

VERSION="1.1"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CENTRIFUGE="${CENTRIFUGE:-/home/intern/internship_projects/life-science/centrifuge/centrifuge}"
INDEX="${CENTRIFUGE_INDEX:-/home/intern/internship_projects/life-science/databases/centrifuge_index/p_compressed+h+v}"

THREADS="${THREADS:-8}"

############################################################
# FUNCTIONS
############################################################

usage() {
    echo
    echo "GenomeGuard v${VERSION}"
    echo
    echo "Usage:"
    echo "  ./genomeguard.sh sample.fastq"
    echo "  ./genomeguard.sh sample_R1.fastq sample_R2.fastq"
    echo
    echo "Optional:"
    echo "  THREADS=16 ./genomeguard.sh sample_R1.fastq sample_R2.fastq"
    echo
    exit 1
}

die() {
    echo
    echo "ERROR: $1"
    echo
    exit 1
}

############################################################
# INPUT
############################################################

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage
fi

READ1="$1"
READ2="${2:-}"

[[ -f "$READ1" ]] || die "Input FASTQ not found: $READ1"

if [[ -n "$READ2" ]]; then
    [[ -f "$READ2" ]] || die "Input FASTQ not found: $READ2"
fi

############################################################
# SAMPLE NAME
############################################################

BASENAME="$(basename "$READ1")"

BASENAME="${BASENAME%.fastq.gz}"
BASENAME="${BASENAME%.fq.gz}"
BASENAME="${BASENAME%.fastq}"
BASENAME="${BASENAME%.fq}"

BASENAME="${BASENAME%_R1}"
BASENAME="${BASENAME%_1}"

SAMPLE="$BASENAME"

RESULTS_DIR="${PROJECT_DIR}/results/${SAMPLE}"

mkdir -p "$RESULTS_DIR"

CLASSIFICATION="${RESULTS_DIR}/${SAMPLE}_classification.tsv"
REPORT="${RESULTS_DIR}/${SAMPLE}_report.tsv"

SUMMARY="${RESULTS_DIR}/${SAMPLE}_genomeguard_summary.tsv"
HTML="${RESULTS_DIR}/${SAMPLE}_genomeguard_report.html"
LOG="${RESULTS_DIR}/${SAMPLE}_genomeguard.log"

exec > >(tee -a "$LOG") 2>&1

############################################################
# HEADER
############################################################

echo "=================================================="
echo "           GENOMEGUARD v${VERSION}"
echo "       SAMPLE INTEGRITY ANALYSIS"
echo "=================================================="
echo
echo "Sample      : ${SAMPLE}"
echo "Read 1      : ${READ1}"

if [[ -n "$READ2" ]]; then
    echo "Read 2      : ${READ2}"
    MODE="PAIRED-END"
else
    echo "Read 2      : NONE"
    MODE="SINGLE-END"
fi

echo "Mode        : ${MODE}"
echo "Threads     : ${THREADS}"
echo "Output      : ${RESULTS_DIR}"
echo

############################################################
# VALIDATE TOOLS
############################################################

echo "### INPUT / TOOL VALIDATION"

command -v awk >/dev/null 2>&1 || die "awk not found"
command -v sort >/dev/null 2>&1 || die "sort not found"
command -v sed >/dev/null 2>&1 || die "sed not found"
command -v centrifuge >/dev/null 2>&1 || true

[[ -x "$CENTRIFUGE" ]] || die "Centrifuge executable not found: $CENTRIFUGE"

if [[ ! -f "${INDEX}.1.cf" ]]; then
    die "Centrifuge index not found: ${INDEX}.1.cf"
fi

echo "Centrifuge  : ${CENTRIFUGE}"
echo "Index       : ${INDEX}"
echo "Validation  : PASS"
echo

############################################################
# STEP 1: CENTRIFUGE
############################################################

echo "### STEP 1: TAXONOMIC CLASSIFICATION"
echo

CENTRIFUGE_ARGS=(
    -x "$INDEX"
    -p "$THREADS"
    -S "$CLASSIFICATION"
    --report-file "$REPORT"
)

if [[ -n "$READ2" ]]; then
    CENTRIFUGE_ARGS+=(
        -1 "$READ1"
        -2 "$READ2"
    )
else
    CENTRIFUGE_ARGS+=(
        -U "$READ1"
    )
fi

"$CENTRIFUGE" "${CENTRIFUGE_ARGS[@]}"

echo
echo "Centrifuge classification complete."
echo

############################################################
# STEP 2: PARSE CLASSIFICATION
############################################################

echo "### STEP 2: GENOMEGUARD EVIDENCE ENGINE"
echo

TOTAL="$(awk 'NR>1 {n++} END {print n+0}' "$CLASSIFICATION")"

CLASSIFIED="$(awk 'NR>1 && $3 != 0 {n++} END {print n+0}' "$CLASSIFICATION")"

UNCLASSIFIED="$(awk 'NR>1 && $3 == 0 {n++} END {print n+0}' "$CLASSIFICATION")"

if [[ "$TOTAL" -gt 0 ]]; then
    CLASS_RATE="$(awk -v c="$CLASSIFIED" -v t="$TOTAL" \
        'BEGIN {printf "%.4f", (c/t)*100}')"

    UNCLASS_RATE="$(awk -v u="$UNCLASSIFIED" -v t="$TOTAL" \
        'BEGIN {printf "%.4f", (u/t)*100}')"
else
    CLASS_RATE="0.0000"
    UNCLASS_RATE="0.0000"
fi

############################################################
# TOP TAXA
############################################################

TOP_TAXA="${RESULTS_DIR}/.top_taxa.tmp"

awk '
NR>1 && $3 != 0 {
    count[$3]++
}
END {
    for (tax in count)
        print count[tax], tax
}' "$CLASSIFICATION" |
sort -nr |
head -20 > "$TOP_TAXA"

DOMINANT_READS="$(awk 'NR==1 {print $1}' "$TOP_TAXA")"
DOMINANT_TAXID="$(awk 'NR==1 {print $2}' "$TOP_TAXA")"

SECOND_READS="$(awk 'NR==2 {print $1}' "$TOP_TAXA")"
SECOND_TAXID="$(awk 'NR==2 {print $2}' "$TOP_TAXA")"

if [[ -z "${DOMINANT_READS:-}" ]]; then
    DOMINANT_READS=0
    DOMINANT_TAXID=0
fi

if [[ -z "${SECOND_READS:-}" ]]; then
    SECOND_READS=0
    SECOND_TAXID=0
fi

DOMINANCE="$(awk -v d="$DOMINANT_READS" -v c="$CLASSIFIED" \
    'BEGIN {if(c>0) printf "%.4f",100*d/c; else print "0.0000"}')"

SECONDARY_READS=$((CLASSIFIED - DOMINANT_READS))

SECONDARY_BURDEN="$(awk -v s="$SECONDARY_READS" -v c="$CLASSIFIED" \
    'BEGIN {if(c>0) printf "%.4f",100*s/c; else print "0.0000"}')"

UNIQUE_TAXA="$(awk '
NR>1 && $3 != 0 {
    tax[$3]=1
}
END {
    n=0
    for(t in tax) n++
    print n
}' "$CLASSIFICATION")"

############################################################
# TAXON NAME LOOKUP
############################################################

get_taxon_name() {

    local TAXID="$1"

    awk -v target="$TAXID" '
    BEGIN {FS="[[:space:]]+"}
    NR>1 {
        if ($3 == target) {
            print $1 " " $2
            exit
        }
    }' "$REPORT"
}

DOMINANT_NAME="$(get_taxon_name "$DOMINANT_TAXID")"
SECOND_NAME="$(get_taxon_name "$SECOND_TAXID")"

[[ -n "$DOMINANT_NAME" ]] || DOMINANT_NAME="TaxID ${DOMINANT_TAXID}"
[[ -n "$SECOND_NAME" ]] || SECOND_NAME="TaxID ${SECOND_TAXID}"

############################################################
# QUALITY FLAGS
############################################################

if awk -v x="$UNCLASS_RATE" 'BEGIN {exit !(x < 10)}'; then
    FLAG_UNCLASS="PASS"
else
    FLAG_UNCLASS="FLAG"
fi

if awk -v x="$DOMINANCE" 'BEGIN {exit !(x >= 50)}'; then
    FLAG_DOMINANT="PASS"
else
    FLAG_DOMINANT="FLAG"
fi

if awk -v x="$SECONDARY_BURDEN" 'BEGIN {exit !(x >= 10)}'; then
    FLAG_SECONDARY="FLAG"
else
    FLAG_SECONDARY="PASS"
fi

if [[ "$UNIQUE_TAXA" -gt 100 ]]; then
    FLAG_COMPLEXITY="FLAG"
else
    FLAG_COMPLEXITY="PASS"
fi

############################################################
# INTEGRITY CATEGORY
############################################################

if awk -v d="$DOMINANCE" -v s="$SECONDARY_BURDEN" '
BEGIN {
    if (d >= 90 && s < 10)
        exit 0
    else
        exit 1
}'; then

    INTEGRITY="STRONG SINGLE-TAXON SIGNAL"

elif awk -v d="$DOMINANCE" '
BEGIN {
    if (d >= 50)
        exit 0
    else
        exit 1
}'; then

    INTEGRITY="MODERATE DOMINANCE"

else

    INTEGRITY="COMPLEX TAXONOMIC PROFILE"

fi

############################################################
# TSV SUMMARY
############################################################

cat > "$SUMMARY" <<EOF
GenomeGuard_Version	${VERSION}
Sample	${SAMPLE}
Input_Mode	${MODE}
Read1	${READ1}
Read2	${READ2:-NA}
Total_Records	${TOTAL}
Classified_Records	${CLASSIFIED}
Unclassified_Records	${UNCLASSIFIED}
Classification_Rate_Percent	${CLASS_RATE}
Unclassified_Rate_Percent	${UNCLASS_RATE}
Dominant_TaxID	${DOMINANT_TAXID}
Dominant_Taxon	${DOMINANT_NAME}
Dominant_Reads	${DOMINANT_READS}
Dominant_Fraction_Percent	${DOMINANCE}
Secondary_Reads	${SECONDARY_READS}
Secondary_Burden_Percent	${SECONDARY_BURDEN}
Leading_Secondary_TaxID	${SECOND_TAXID}
Leading_Secondary_Taxon	${SECOND_NAME}
Leading_Secondary_Reads	${SECOND_READS}
Unique_Classified_Taxa	${UNIQUE_TAXA}
Unclassified_QC	${FLAG_UNCLASS}
Dominant_Signal_QC	${FLAG_DOMINANT}
Secondary_Burden_QC	${FLAG_SECONDARY}
Taxonomic_Complexity_QC	${FLAG_COMPLEXITY}
Integrity_Category	${INTEGRITY}
EOF

############################################################
# TOP TAXA TSV
############################################################

TOP_TAXA_TABLE="${RESULTS_DIR}/.top_taxa_table.tmp"

echo -e "Rank\tReads\tTaxID\tRelative_Percent\tTaxon" > "$TOP_TAXA_TABLE"

rank=0

while read -r reads taxid; do

    rank=$((rank + 1))

    relative="$(awk -v r="$reads" -v c="$CLASSIFIED" \
        'BEGIN {if(c>0) printf "%.4f",100*r/c; else print "0.0000"}')"

    name="$(get_taxon_name "$taxid")"
    [[ -n "$name" ]] || name="TaxID ${taxid}"

    echo -e "${rank}\t${reads}\t${taxid}\t${relative}\t${name}" \
        >> "$TOP_TAXA_TABLE"

done < "$TOP_TAXA"

############################################################
# HTML REPORT
############################################################

cat > "$HTML" <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">

<title>GenomeGuard Report - ${SAMPLE}</title>

<style>
body {
    font-family: Arial, Helvetica, sans-serif;
    margin: 40px;
    background: #f5f7fa;
    color: #222;
}

.container {
    max-width: 1200px;
    margin: auto;
}

header {
    background: #17202a;
    color: white;
    padding: 30px;
    border-radius: 12px;
    margin-bottom: 25px;
}

h1 {
    margin: 0 0 8px 0;
}

h2 {
    margin-top: 35px;
    border-bottom: 2px solid #ddd;
    padding-bottom: 8px;
}

.grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 15px;
}

.card {
    background: white;
    padding: 20px;
    border-radius: 10px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
}

.value {
    font-size: 25px;
    font-weight: bold;
    margin-top: 8px;
}

table {
    width: 100%;
    border-collapse: collapse;
    background: white;
    margin-top: 15px;
}

th, td {
    padding: 10px;
    border-bottom: 1px solid #ddd;
    text-align: left;
}

th {
    background: #eaecef;
}

.pass {
    font-weight: bold;
}

.flag {
    font-weight: bold;
}

.warning {
    background: #fff4d6;
    padding: 18px;
    border-radius: 10px;
    margin-top: 20px;
}

.footer {
    margin-top: 40px;
    font-size: 13px;
    color: #666;
}
</style>
</head>

<body>

<div class="container">

<header>
<h1>GenomeGuard v${VERSION}</h1>
<div>Sample Integrity Analysis</div>
<div style="margin-top:10px;">Sample: <strong>${SAMPLE}</strong></div>
</header>

<h2>Sample Summary</h2>

<div class="grid">

<div class="card">
Classification rate
<div class="value">${CLASS_RATE}%</div>
</div>

<div class="card">
Unclassified reads
<div class="value">${UNCLASS_RATE}%</div>
</div>

<div class="card">
Dominant taxon
<div class="value">${DOMINANCE}%</div>
</div>

<div class="card">
Taxonomic complexity
<div class="value">${UNIQUE_TAXA}</div>
</div>

</div>

<h2>Primary Taxonomic Signal</h2>

<div class="card">
<strong>${DOMINANT_NAME}</strong><br>
TaxID: ${DOMINANT_TAXID}<br>
Reads: ${DOMINANT_READS}<br>
Relative contribution: ${DOMINANCE}%
</div>

<h2>Secondary Signal</h2>

<div class="card">
Secondary classified reads: ${SECONDARY_READS}<br>
Secondary burden: ${SECONDARY_BURDEN}%<br><br>

Leading secondary signal:<br>
<strong>${SECOND_NAME}</strong><br>
TaxID: ${SECOND_TAXID}<br>
Reads: ${SECOND_READS}
</div>

<h2>Quality-Control Flags</h2>

<table>
<tr>
<th>Metric</th>
<th>Status</th>
</tr>

<tr>
<td>Unclassified fraction</td>
<td class="${FLAG_UNCLASS,,}">${FLAG_UNCLASS}</td>
</tr>

<tr>
<td>Dominant taxonomic signal</td>
<td class="${FLAG_DOMINANT,,}">${FLAG_DOMINANT}</td>
</tr>

<tr>
<td>Secondary taxonomic burden</td>
<td class="${FLAG_SECONDARY,,}">${FLAG_SECONDARY}</td>
</tr>

<tr>
<td>Taxonomic complexity</td>
<td class="${FLAG_COMPLEXITY,,}">${FLAG_COMPLEXITY}</td>
</tr>
</table>

<h2>Integrity Assessment</h2>

<div class="card">
<strong>${INTEGRITY}</strong>
</div>

<h2>Top Classified Taxa</h2>

<table>
<tr>
<th>Rank</th>
<th>Reads</th>
<th>TaxID</th>
<th>Relative %</th>
<th>Taxon</th>
</tr>
EOF

tail -n +2 "$TOP_TAXA_TABLE" | while IFS=$'\t' read -r rank reads taxid relative name; do
    cat >> "$HTML" <<EOF
<tr>
<td>${rank}</td>
<td>${reads}</td>
<td>${taxid}</td>
<td>${relative}%</td>
<td>${name}</td>
</tr>
EOF
done

cat >> "$HTML" <<EOF
</table>

<h2>GenomeGuard Interpretation</h2>

<div class="warning">

The sample contains a clear dominant taxonomic signal.

A substantial secondary classified signal is also present.
Secondary taxonomic evidence should therefore be investigated
before assigning a simple single-organism interpretation.

GenomeGuard reports computational taxonomic evidence.
Secondary assignments do not independently prove:

<ul>
<li>contamination</li>
<li>co-infection</li>
<li>mixed biological origin</li>
<li>true organism abundance</li>
</ul>

Closely related taxa may share sequence similarity.
Reference database composition and classifier behaviour must
therefore be considered before biological conclusions.

</div>

<h2>Input Information</h2>

<table>
<tr><th>Parameter</th><th>Value</th></tr>
<tr><td>GenomeGuard version</td><td>${VERSION}</td></tr>
<tr><td>Sample</td><td>${SAMPLE}</td></tr>
<tr><td>Mode</td><td>${MODE}</td></tr>
<tr><td>Read 1</td><td>${READ1}</td></tr>
<tr><td>Read 2</td><td>${READ2:-NA}</td></tr>
<tr><td>Centrifuge index</td><td>${INDEX}</td></tr>
<tr><td>Threads</td><td>${THREADS}</td></tr>
</table>

<div class="footer">
Generated by GenomeGuard v${VERSION}.
This report describes computational taxonomic evidence and
should not be interpreted as definitive biological proof.
</div>

</div>
</body>
</html>
EOF

############################################################
# CLEANUP
############################################################

rm -f "$TOP_TAXA"
rm -f "$TOP_TAXA_TABLE"

############################################################
# FINAL TERMINAL REPORT
############################################################

echo
echo "=================================================="
echo "        GENOMEGUARD v${VERSION} COMPLETE"
echo "=================================================="
echo
echo "Sample:"
echo "  ${SAMPLE}"
echo
echo "Primary taxon:"
echo "  ${DOMINANT_NAME}"
echo "  TaxID: ${DOMINANT_TAXID}"
echo "  Reads: ${DOMINANT_READS}"
echo "  Fraction: ${DOMINANCE}%"
echo
echo "Secondary burden:"
echo "  ${SECONDARY_BURDEN}%"
echo
echo "Taxonomic complexity:"
echo "  ${UNIQUE_TAXA} classified taxa"
echo
echo "Integrity category:"
echo "  ${INTEGRITY}"
echo
echo "OUTPUT FILES"
echo
echo "TSV summary:"
echo "  ${SUMMARY}"
echo
echo "HTML report:"
echo "  ${HTML}"
echo
echo "Log:"
echo "  ${LOG}"
echo
echo "=================================================="

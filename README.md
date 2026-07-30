# GenomeGuard

**Evidence-Weighted Taxonomic Sample Integrity Analysis from FASTQ Reads**

GenomeGuard is a command-line bioinformatics pipeline built around **Centrifuge**. It accepts sequencing reads in FASTQ format, performs taxonomic classification against a pre-built Centrifuge reference index, and converts the raw classification output into an evidence-weighted sample-level report.

The central philosophy of GenomeGuard is:

> **Taxonomic evidence should guide investigation, not be mistaken for biological proof.**

GenomeGuard therefore does not simply report the organism with the largest number of reads. It examines dominant and secondary signals, read-level support, score separation, taxonomic complexity, and relationships between related taxa.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Why GenomeGuard Was Created](#why-genomeguard-was-created)
3. [What GenomeGuard Does](#what-genomeguard-does)
4. [Architecture](#architecture)
5. [Complete Workflow](#complete-workflow)
6. [Input Data](#input-data)
7. [The 5.48 GB Centrifuge Index](#the-548-gb-centrifuge-index)
8. [Why the Index Is Required](#why-the-index-is-required)
9. [What We Know About the Index](#what-we-know-about-the-index)
10. [GenomeGuard Development History](#genomeguard-development-history)
11. [Evidence Engine](#evidence-engine)
12. [Final v1.1 Pipeline](#final-v11-pipeline)
13. [Installation and Environment](#installation-and-environment)
14. [Running GenomeGuard](#running-genomeguard)
15. [SRA-to-FASTQ Workflow](#sra-to-fastq-workflow)
16. [Understanding the Output](#understanding-the-output)
17. [Demonstration Dataset](#demonstration-dataset)
18. [Interpretation of the Demonstration](#interpretation-of-the-demonstration)
19. [Scientific Limitations](#scientific-limitations)
20. [Reproducibility](#reproducibility)
21. [Recommended Repository Structure](#recommended-repository-structure)
22. [Future Roadmap](#future-roadmap)
23. [Project Significance](#project-significance)
24. [Final Statement](#final-statement)

---

# Project Overview

GenomeGuard was developed as a practical bioinformatics project for analyzing sequencing datasets at the taxonomic level.

A conventional taxonomic classifier can answer:

> "Which taxonomic IDs were assigned to the reads?"

GenomeGuard goes one step further and asks:

> "How strong and consistent are those taxonomic signals, and which signals deserve further investigation?"

The pipeline is designed to work with ordinary FASTQ sequencing data and a reusable Centrifuge index.

The final intended user experience is simple:

```bash
./genomeguard.sh sample_1.fastq sample_2.fastq
```

For paired-end data:

```bash
./genomeguard.sh SRR39879768_1.fastq SRR39879768_2.fastq
```

The pipeline automatically performs the taxonomic classification and produces structured results.

---

# Why GenomeGuard Was Created

Raw taxonomic classification output can contain hundreds of taxonomic assignments.

For example, a sequencing dataset may show:

- one dominant organism,
- several closely related taxa,
- many low-abundance taxa,
- unclassified reads,
- ambiguous read assignments,
- and signals that are difficult to interpret biologically.

Simply seeing a secondary taxon does **not** necessarily mean that the sample is contaminated or contains a second organism.

Possible explanations include:

- genuine biological mixtures,
- shared genomic sequences,
- sequence similarity between related organisms,
- reference database composition,
- classifier behaviour,
- ambiguous read assignments.

GenomeGuard was therefore developed incrementally to provide an evidence layer around Centrifuge.

---

# What GenomeGuard Does

The final pipeline performs the following major operations:

```text
FASTQ
  |
  v
Input validation
  |
  v
Centrifuge taxonomic classification
  |
  v
Read-level classification evidence
  |
  v
Taxonomic distribution analysis
  |
  +--> Dominant taxon
  |
  +--> Secondary burden
  |
  +--> Taxonomic complexity
  |
  +--> Read-level evidence
  |
  +--> Score discrimination
  |
  +--> Related taxonomic signals
  |
  v
GenomeGuard interpretation
  |
  +--> TSV summary
  +--> HTML report
  +--> Log
```

---

# Architecture

GenomeGuard has two major layers.

## Layer 1 — Taxonomic Classification

Centrifuge performs the actual read classification.

It uses a pre-built reference index to classify sequencing reads.

The classification stage produces information including:

- read ID,
- taxonomic ID,
- Centrifuge score,
- second-best score,
- hit length,
- query length,
- number of matches.

The taxonomic report provides aggregated information such as:

- taxon,
- taxonomic ID,
- taxonomic rank,
- read support,
- unique-read support,
- abundance-related information.

## Layer 2 — GenomeGuard Evidence Engine

GenomeGuard processes those Centrifuge results.

It calculates and evaluates:

- total records,
- classified records,
- unclassified records,
- classification rate,
- dominant taxon,
- dominant fraction,
- secondary burden,
- taxonomic complexity,
- unique-read support,
- relative contribution,
- score support,
- score separation,
- taxonomic relationships.

This separation is important:

> **Centrifuge performs classification. GenomeGuard performs evidence organization and interpretation.**

---

# Complete Workflow

The complete recommended workflow is:

```text
1. Obtain sequencing dataset
        |
        v
2. SRA accession
        |
        v
3. prefetch
        |
        v
4. SRA file
        |
        v
5. fasterq-dump --split-files
        |
        v
6. FASTQ files
        |
        +--------------------+
        |                    |
        v                    v
     FastQC             GenomeGuard
        |                    |
        |                    v
        |              Centrifuge
        |                    |
        |                    v
        |              Taxonomic evidence
        |                    |
        |                    v
        |              GenomeGuard engine
        |                    |
        +-----------> TSV + HTML + LOG
```

---

# Input Data

GenomeGuard accepts FASTQ sequencing reads.

For paired-end sequencing:

```text
sample_1.fastq
sample_2.fastq
```

For example:

```text
SRR39879768_1.fastq
SRR39879768_2.fastq
```

The final pipeline was designed so that users do not need to manually run every internal analysis stage.

---

# The 5.48 GB Centrifuge Index

One of the most important resources used by GenomeGuard is the Centrifuge reference index.

The local project environment contained an index with the prefix:

```text
p_compressed+h+v
```

The observed index files included:

```text
p_compressed+h+v.1.cf
p_compressed+h+v.2.cf
p_compressed+h+v.3.cf
```

The total local index size was approximately **5.48 GB**.

## What is this 5.48 GB file set?

It is not a sequencing sample.

It is not the FASTQ data.

It is not the GenomeGuard output.

It is a **pre-built taxonomic reference/index resource used by Centrifuge**.

The index represents the searchable reference information and indexing structures required for efficient taxonomic classification.

---

# Why the Index Is Required

Without a reference index, Centrifuge would not have the reference information required to classify sequencing reads taxonomically.

Conceptually:

```text
Sequencing read
      |
      v
Centrifuge
      |
      | compares against
      v
Reference index
      |
      v
Taxonomic assignment
```

The index is prepared once and can then be reused across many sequencing samples.

This is why GenomeGuard does not need to rebuild a 5.48 GB database every time a new FASTQ file is analyzed.

---

# What We Know About the Index

The project environment confirmed:

```text
Centrifuge version:
1.0.5
```

Observed executable:

```text
/home/intern/internship_projects/life-science/centrifuge/centrifuge
```

Observed index prefix:

```text
/home/intern/internship_projects/life-science/databases/centrifuge_index/p_compressed+h+v
```

The project demonstrated assignments including:

- Escherichia coli
- Shigella-associated taxa
- other Enterobacterales-associated taxa
- lower-abundance additional taxa.

However, the original project records did **not** establish a complete organism-by-organism inventory of the index or the exact reference release/build provenance.

Therefore this repository should not claim that the index contains a specific complete list of organisms unless that information is independently verified.

For publication-grade reproducibility, the database source, release, build date, taxonomy source, and checksums should be recorded.

---

# Why We Use a Reusable Reference Index

The index provides several practical advantages:

### 1. Speed

The reference has already been indexed.

### 2. Reusability

The same database can be used for multiple FASTQ datasets.

### 3. Standardization

Samples can be analyzed against the same reference resource.

### 4. Automation

GenomeGuard can call Centrifuge automatically.

### 5. Evidence generation

Centrifuge provides the read-level scores and taxonomic information that GenomeGuard needs.

---

# GenomeGuard Development History

GenomeGuard was developed incrementally.

## v0.1 — Taxonomic Integrity Analysis

Initial metrics:

- total records,
- classified records,
- unclassified reads,
- classification rate,
- top taxa.

Purpose:

Establish the basic taxonomic profile.

---

## v0.2 — Evidence Engine

Added:

- dominant taxon,
- secondary taxon,
- dominance,
- secondary/dominant ratio,
- unique-read support.

Purpose:

Move from simple counts toward evidence.

---

## v0.3 — Evidence-Weighted Analysis

Added:

- unique-read percentage,
- relative contribution,
- evidence categories.

Purpose:

Distinguish stronger and weaker taxonomic signals.

---

## v0.4 / v0.4.1 / v0.4.2 — Robust Parsing

The project discovered and corrected a report-field interpretation problem.

The Centrifuge report contained fields such as:

```text
Escherichia coli
562
species
12319210
493767
426861
0.798753
```

The corrected implementation properly mapped:

- TaxID,
- taxon,
- read count,
- unique-read count,
- report abundance-related field.

This was an important step because downstream evidence calculations depend on correctly understanding the source data.

---

## v0.5 / v0.5.1 — Sample Integrity Engine

Added:

- secondary burden,
- taxonomic complexity,
- QC flags,
- conservative interpretation.

Purpose:

Determine whether a sample has a dominant signal while explicitly reporting secondary complexity.

---

## v0.6 — Secondary Signal Investigation Engine

Added:

- ranked secondary signals,
- investigation candidates,
- secondary evidence ranking.

Purpose:

Identify which secondary signals deserve attention first.

---

## v0.7 — Read-Level Evidence Validation

Added:

- mean Centrifuge score,
- positive-score fraction,
- high-score support,
- very-high-score support,
- score gap.

Purpose:

Check whether secondary signals are supported at the individual-read level.

---

## v0.8 — Cross-Evidence Consistency

Combined several evidence dimensions:

1. read abundance,
2. relative contribution,
3. positive-score support,
4. high-score support,
5. very-high-score support.

Purpose:

Prioritize signals that are consistently supported across several dimensions.

---

## v0.9 — Taxonomic Conflict Analysis

Introduced related-group reasoning.

The demonstration sample contained:

```text
Escherichia-associated signal
Shigella-associated signal
Other Enterobacterales-associated signal
```

The purpose was to avoid automatically treating every taxonomic label as an independent organism.

Closely related taxa can share sequence similarity.

---

## v0.10 — Discriminative Evidence Engine

Added comparison of:

```text
Mean Centrifuge score
Mean second-best score
Mean score gap
```

A larger score gap can provide stronger computational separation.

A small or zero score gap can indicate ambiguity.

Evidence categories included:

```text
SUPPORTED
MODERATE
AMBIGUOUS
WEAK
```

---

# v1.0 — Integrated GenomeGuard Pipeline

The individual analysis stages were integrated into a single user-facing pipeline.

The intended workflow became:

```bash
./genomeguard.sh sample_1.fastq sample_2.fastq
```

The pipeline automatically:

1. validates inputs,
2. runs Centrifuge,
3. analyzes classification results,
4. creates a sample-specific output directory,
5. writes GenomeGuard results.

---

# v1.1 — Final Practical Pipeline

The final working version added:

- input validation,
- tool validation,
- Centrifuge index validation,
- paired-end handling,
- automatic sample identification,
- TSV summary,
- HTML report,
- execution log,
- integrated classification output.

The final output is designed to be useful both to a researcher reading the result and to a computational workflow consuming the TSV.

---

# Evidence Engine

GenomeGuard uses multiple evidence dimensions.

## 1. Read Support

How many reads were assigned to a taxon?

A large number of reads provides stronger quantitative support than a tiny number of reads.

---

## 2. Relative Contribution

What fraction of classified reads belongs to the taxon?

This provides context for raw read counts.

---

## 3. Unique-Read Support

How many reads provide unique support for the taxon according to the Centrifuge report?

This helps distinguish broad assignments from more specific support.

---

## 4. Positive Score Support

What fraction of assigned reads have a positive Centrifuge score?

---

## 5. High-Score Support

What fraction of reads have scores above the selected high-score threshold?

---

## 6. Very-High-Score Support

What fraction of reads have very high scores?

---

## 7. Score Separation

GenomeGuard compares the best assignment score with the second-best score.

Conceptually:

```text
Score gap = Best score - Second-best score
```

A larger gap can indicate stronger computational discrimination.

However, the score gap should not be treated as proof of biological presence.

---

# Final v1.1 Pipeline

The final command is:

```bash
./genomeguard.sh sample_1.fastq sample_2.fastq
```

For example:

```bash
./genomeguard.sh SRR39879768_1.fastq SRR39879768_2.fastq
```

The pipeline then performs:

```text
INPUT VALIDATION
       |
       v
CENTRIFUGE CLASSIFICATION
       |
       v
READ-LEVEL CLASSIFICATION
       |
       v
TAXONOMIC REPORT
       |
       v
GENOMEGUARD EVIDENCE ENGINE
       |
       +-------------------+
       |                   |
       v                   v
     TSV                 HTML
       |
       v
     LOG
```

---

# Installation and Environment

The project was developed in a Linux environment.

Observed Centrifuge version:

```text
centrifuge-class version 1.0.5
```

Observed executable:

```text
/home/intern/internship_projects/life-science/centrifuge/centrifuge
```

Example project location:

```text
~/internship_projects/life-science/GenomeGuard
```

The Centrifuge index was located under:

```text
~/internship_projects/life-science/databases/centrifuge_index/
```

The exact paths may differ on another machine.

GenomeGuard should therefore be configured with the appropriate Centrifuge executable and index location for the target environment.

---

# Running GenomeGuard

## Paired-End Data

```bash
./genomeguard.sh SRR12345678_1.fastq SRR12345678_2.fastq
```

## Output

The final pipeline creates a sample-specific result directory.

The exact structure depends on the current script configuration, but the intended output includes:

```text
sample/
├── sample_genomeguard_summary.tsv
├── sample_genomeguard_report.html
├── sample_genomeguard.log
├── sample_classification.tsv
└── sample_report.tsv
```

The TSV is intended for structured analysis.

The HTML is intended for human-readable review.

The log records execution information.

---

# SRA-to-FASTQ Workflow

If the dataset originates from NCBI SRA, the workflow used during development was:

## Step 1 — Download the SRA Run

```bash
prefetch SRR12345678
```

---

## Step 2 — Convert SRA to FASTQ

For paired-end sequencing:

```bash
fasterq-dump --split-files SRR12345678
```

This generates:

```text
SRR12345678_1.fastq
SRR12345678_2.fastq
```

---

## Step 3 — Quality Assessment

Run:

```bash
fastqc SRR12345678_1.fastq SRR12345678_2.fastq
```

FastQC is used for basic sequencing-read quality assessment.

GenomeGuard is not intended to replace FastQC.

---

## Step 4 — GenomeGuard

Run:

```bash
./genomeguard.sh SRR12345678_1.fastq SRR12345678_2.fastq
```

---

# Understanding the Output

## Total Records

The total number of classification records processed.

---

## Classified Records

The number of records receiving a taxonomic assignment.

---

## Unclassified Records

Reads for which Centrifuge did not provide a classified taxonomic assignment.

---

## Classification Rate

```text
classified records / total records
```

A high classification rate means most records received a taxonomic assignment.

It does not automatically mean the biological interpretation is correct.

---

## Dominant Taxon

The taxon with the largest number of classified reads.

---

## Dominant Fraction

The fraction of classified reads assigned to the dominant taxon.

---

## Secondary Burden

The fraction of classified reads assigned outside the dominant taxon.

For example:

```text
Dominant fraction = 76.0932%

Secondary burden = 23.9068%
```

These values describe the computational taxonomic distribution.

They do not automatically establish contamination.

---

## Taxonomic Complexity

The number of unique classified taxa.

A high number may indicate:

- genuine complexity,
- closely related sequence assignments,
- database effects,
- classifier behaviour,
- or other technical/biological factors.

---

# Demonstration Dataset

The project was demonstrated using:

```text
SRR39879768
```

The corresponding paired-end files were:

```text
SRR39879768_1.fastq
SRR39879768_2.fastq
```

---

# Demonstration Results

GenomeGuard produced:

```text
Total records       : 654997
Classified records  : 648898
Unclassified        : 6099
Classification rate : 99.0689%
```

Primary taxonomic signal:

```text
Taxon        : Escherichia coli
TaxID        : 562
Reads        : 493767
Fraction     : 76.0932%
```

Secondary burden:

```text
23.9068%
```

Taxonomic complexity:

```text
473 classified taxa
```

Final integrity category:

```text
MODERATE DOMINANCE
```

---

# Interpretation of the Demonstration

The sample contains a clear dominant computational taxonomic signal:

```text
Escherichia coli
```

However, a substantial secondary classified population is also present.

Therefore the appropriate interpretation is:

> The sample has a clear dominant taxonomic signal, but substantial secondary classified evidence is present and should be investigated.

It is **not** scientifically appropriate to conclude from this analysis alone that the sample is:

- contaminated,
- co-infected,
- a mixed biological sample,
- or composed of a particular organism at a specific biological abundance.

---

# Scientific Limitations

GenomeGuard is an evidence-reporting pipeline.

It does not independently prove:

- contamination,
- co-infection,
- mixed biological origin,
- true organism abundance.

Taxonomic assignments can be affected by:

### Shared Sequence Similarity

Related organisms can contain similar sequences.

### Reference Database Composition

What is represented in the reference database influences what can be classified.

### Classifier Behaviour

The classification algorithm determines how reads are assigned.

### Ambiguous Reads

Some reads can have similarly good matches to multiple references.

### Database Version

Changing the database can change classification results.

Therefore, GenomeGuard results should be interpreted in the context of:

- sequencing quality,
- experimental design,
- reference database,
- independent genomic evidence,
- and biological knowledge of the sample.

---

# Reproducibility

For a publication-quality analysis, record:

```text
GenomeGuard version
Centrifuge version
Centrifuge index prefix
Database source
Database release/version
Database build date
Taxonomy source/version
Input accession
Sequencing layout
Command used
Date of analysis
Output files
```

It is also recommended to record checksums for important database files.

This is particularly important because the same FASTQ data can produce different taxonomic profiles when analyzed against different reference databases.

---

# Recommended Repository Structure

A clean GitHub repository can be organized as:

```text
GenomeGuard/
│
├── README.md
│
├── genomeguard.sh
│
├── scripts/
│   └── development/
│       ├── genomeguard_analyze.sh
│       ├── genomeguard_evidence.sh
│       ├── genomeguard_evidence_v03.sh
│       ├── genomeguard_evidence_v04.sh
│       ├── genomeguard_evidence_v04_1.sh
│       ├── genomeguard_evidence_v04_2.sh
│       ├── genomeguard_v05.sh
│       ├── genomeguard_v05_1.sh
│       ├── genomeguard_v06.sh
│       ├── genomeguard_v07.sh
│       ├── genomeguard_v08.sh
│       ├── genomeguard_v09.sh
│       └── genomeguard_v10.sh
│
├── docs/
│   ├── GenomeGuard_Project.txt
│   ├── GenomeGuard_Workflow.txt
│   ├── GenomeGuard_Roadmap.txt
│   └── GenomeGuard_Index_Database.txt
│
├── examples/
│   └── README.md
│
├── results/
│   └── .gitkeep
│
└── LICENSE
```

Large FASTQ files and the 5.48 GB reference index should generally **not** be committed directly to the GitHub repository.

Instead, document how users can obtain or build the required resources.

---

# Future Roadmap

The current v1.1 pipeline should be treated as a foundation for further validation.

Potential future work:

## 1. Multi-Dataset Validation

Test GenomeGuard across many independent sequencing datasets.

---

## 2. Known Single-Organism Controls

Use appropriate benchmark samples where the expected dominant organism is known.

---

## 3. Known Mixture Controls

Use benchmark datasets with known mixtures to determine whether GenomeGuard can reliably prioritize secondary signals.

---

## 4. Independent Tool Comparison

Compare GenomeGuard/Centrifuge findings with other taxonomic classification approaches.

The purpose should be validation rather than simply increasing the number of tools.

---

## 5. Database Benchmarking

Evaluate how results change with different reference databases.

---

## 6. Visualization

Future releases could add:

- taxonomic composition plots,
- evidence-strength plots,
- secondary-signal ranking charts,
- score-gap distributions,
- HTML interactive visualizations.

---

## 7. Automated Testing

Create test datasets and regression tests so that every new GenomeGuard version can be checked against expected results.

---

## 8. Database Provenance

Add automatic reporting of:

- database version,
- database build date,
- reference source,
- taxonomy version,
- checksums.

---

# Project Significance

GenomeGuard demonstrates a broader bioinformatics principle:

> **A classifier output is not the same thing as a biological conclusion.**

A raw taxonomic table can be technically correct while still being difficult to interpret.

GenomeGuard adds an evidence-oriented layer that asks:

```text
How much support exists?
        |
        v
How consistent is that support?
        |
        v
How discriminative is the assignment?
        |
        v
Are related taxa being counted independently?
        |
        v
Does the signal deserve further investigation?
```

This makes the workflow more transparent and more scientifically cautious.

---

# Final Statement

GenomeGuard is an evidence-weighted taxonomic sample-integrity pipeline built around Centrifuge.

It transforms:

```text
FASTQ reads
```

into:

```text
Taxonomic classification
        +
Read-level evidence
        +
Taxonomic distribution
        +
Secondary-signal analysis
        +
Discriminative evidence
        |
        v
Structured sample-level report
```

The final project emphasizes:

- reproducibility,
- automation,
- transparent evidence,
- conservative interpretation,
- practical bioinformatics,
- machine-readable output,
- human-readable reporting.

The most important principle of the project is:

> **GenomeGuard identifies evidence that deserves investigation; it does not turn taxonomic assignments into biological proof.**

---

## Project Status

**Current project version:** GenomeGuard v1.1

**Core classifier:** Centrifuge 1.0.5

**Primary input:** FASTQ

**Primary outputs:** TSV + HTML + log

**Demonstration dataset:** SRR39879768

**Demonstration dominant signal:** Escherichia coli (TaxID 562)

**Demonstration integrity category:** MODERATE DOMINANCE

---

## License

Add the license selected for the GitHub repository before publication.

---

## Citation

If this project is eventually used in academic work, cite the GenomeGuard repository together with the relevant software/database resources and document the exact database version used for the analysis.

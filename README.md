# Denovo RNA Var Assembly Pipeline (denovarna)

A Nextflow pipeline for de novo assembly and functional annotation of VAR genes from _Plasmodium falciparum_ bulk RNAseq data.

## Table of Contents

- [Overview](#overview)
- [Pipeline Phases](#pipeline-phases)
- [Input Requirements](#input-requirements)
- [Running the Pipeline](#running-the-pipeline)

## Overview

The **denovar** pipeline performs de novo assembly of VAR genes from Illumina paired-end short reads originating from _Plasmodium falciparum_ bulk RNAseq. It produces assembled VAR gene contigs and their corresponding domain annotations.

### What the Pipeline Does

This pipeline processes RNAseq data from _Plasmodium falciparum_ samples to:

1. **Extract** polymorphic VAR gene reads from sequencing data by filtering out human and core genomic sequences
2. **Assemble** VAR genes de novo using RNASPAdes and remove duplicate sequences
3. **Correct** and validate assembled sequences using read mapping and error correction
4. **Predict** genes and translate to proteins using Augustus
5. **Annotate** functional protein domains (DBL, CIDR, NTS, ATS) using HMMER
6. **Classify** assembled genes by comparing against the P. falciparum 3D7 reference and transfer functional annotations

## Pipeline Phases

The denovar pipeline consists of 4 major phases. Each phase is broken down into detailed steps:

### Phase 1: Preprocessing and Polymorphic Read Extraction

The goal is to isolate VAR gene reads from the bulk RNAseq data, filtering out human contamination and core genomic sequences.

#### Step 1.1: Prepare Modified Reference Genome (RemovePfEMP1fromGenome)

**Purpose:** Avoid mapping polymorphic VAR reads to the core reference by masking known PfEMP1 and Stevor gene regions for later specific capture.

**Process:**

1. Parse the P. falciparum GFF file to identify genomic coordinates of pfemp1 and stevor genes
2. Create a BED file (pfemp1_regions.bed) with these coordinates, extended by 500bp on each side
3. Use bedtools maskfasta to soft-mask these regions in the original genome
4. Output: `PlasmoDB-68_Pfalciparum3D7_noPfEMP1.fasta` (masked reference genome)

#### Step 1.2: Index Reference Databases (IndexGenomes)

**Purpose:** Create index files required for rapid read alignment in subsequent steps.

**Process:**

1. Index the human genome (hg38)
2. Index the masked P. falciparum genome
3. Index the varDB exon 1 database (collection of diverse var exon 1 sequences)

#### Step 1.3: Extract Polymorphic Reads (ExtractPolymorphicPfemp1Reads)

**Purpose:** Isolate reads likely derived from non-core, polymorphic VAR genes.

**Process per sample:**

1. **Map to human genome:** Align reads to hg38 and extract read pairs that are UNMAPPED (not human)
2. **Map to masked core genome:** Align reads to the masked P. falciparum genome. Extract read pairs where BOTH mates are UNMAPPED (non-core reads)
3. **Map to varDB exon 1:** Align reads to the varDB exon 1 database. Extract read pairs where BOTH mates are MAPPED (confirmed VAR-like reads)
4. **Combine read sets:** The final "polymorphic reads" = (unmapped to human) ∪ (unmapped to core) ∩ (mapped to varDB)
5. **Extract FASTQs:** Use seqtk subseq to pull identified read pairs from original FASTQ files
6. **Output:** `*_polymorphic_reads_R1.fq` and `*_polymorphic_reads_R2.fq`

---

### Phase 2: VAR Gene Assembly & Refinement

The goal is to assemble the extracted polymorphic reads into contigs and correct any assembly errors.

#### Step 2.1: De novo Assembly and Dereplication (AssembleDereplicate)

**Purpose:** Assemble polymorphic reads into contigs and reduce redundancy.

**Process per sample:**

1. **Assemble with RNASPAdes:** Run rnaspades.py assembler on the polymorphic reads
2. **Dereplicate with CD-HIT:** Use cd-hit to cluster and remove highly similar contigs (≥99% identity)
3. **Size selection:** Filter unique contigs, retaining only those ≥300 bp
4. **Output:** `*_unique_contigs_300bp.fasta`

#### Step 2.2: Read Mapping for Validation (MapPolyReadsContigs)

**Purpose:** Map original polymorphic reads back to assembled contigs to create a BAM file for assembly correction.

**Process per sample:**

1. Build a Bowtie2 index from the unique contigs
2. Align polymorphic reads to the contigs
3. Sort and index the resulting BAM file
4. **Output:** `*_mapped_to_unique_contigs.bam`

#### Step 2.3: Assembly Correction (AssemblyCorrection)

**Purpose:** Improve accuracy of assembled contigs using read alignment information to fix base errors and fill gaps.

**Process per sample:**

1. Run Pilon using assembled contigs and BAM file from Step 2.2 to generate corrected genome
2. Clean FASTA headers in Pilon output
3. Generate basic assembly statistics
4. **Output:** `*_genome.fasta` (corrected VAR gene assembly)

---

### Phase 3: Gene Prediction and Functional Annotation

The goal is to predict genes in the contigs and identify their functional domains.

#### Step 3.1: Gene Prediction with AUGUSTUS (VargeneAnnotation)

**Purpose:** Predict open reading frames (ORFs) and protein-coding genes in assembled contigs.

**Process per sample:**

1. Run Augustus with Plasmodium-specific model (--species=Pf_VSA) to predict genes
2. Generate FASTA files:
   - `*.genes.aa.fasta` (predicted proteins)
   - `*.exons.na.fasta` (exon nucleotide sequences)
   - `*.genes.na.fasta` (gene nucleotide sequences)
3. Convert AUGUSTUS output to standard GTF and GFF formats
4. Use featureCounts to count reads mapping to each predicted gene
5. Search predicted proteins for key VAR gene motifs (LARSFADIG and DYVPQYLRW)
6. **Output:** Gene predictions in GTF/GFF, protein and nucleotide sequences, motif search results

#### Step 3.2: Protein Domain Annotation (DomainsAnnotation)

**Purpose:** Identify and classify major protein domains (DBL-alpha, CIDR-alpha, NTS, ATS) within predicted genes.

**Process per sample:**

1. Run hmmscan on predicted protein sequences against custom HMM database of VAR domains
2. Parse domain hits into GFF file using custom Perl scripts
3. Generate domain association files summarizing domain architecture (e.g., DBLa-CIDRa-)
4. **Output:** `*.Domain.gff` (domain annotations in GFF format), domain association summaries

#### Step 3.3: Protein Subdomain Annotation (SubDomainsAnnotation)

**Purpose:** Identify finer-scale subdomains within major domains identified in Step 3.2.

**Process per sample:**

1. Run hmmscan against separate HMM database of VAR subdomains
2. Parse results into subdomain GFF file
3. Generate subdomain association summaries
4. **Output:** `*.Subdomain.gff` (subdomain annotations in GFF format)

---

### Phase 4: Classification and Final Annotation

The goal is to transfer functional annotations from the well-annotated P. falciparum 3D7 reference.

#### Step 4.1: BLAST Annotation against 3D7 (BlastAnnotation)

**Purpose:** Transfer functional annotations from P. falciparum 3D7 reference strain to newly assembled genes.

**Process per sample:**

1. Create BLAST database from 3D7 reference protein sequences
2. Run blastp to find best 3D7 match for each predicted gene
3. Transfer annotation from best hit to sequences in FASTA files
4. Extract final set of confidently annotated VAR genes based on domain presence
5. **Output:**
   - `*.3D7ref.genes.aa.fasta` (3D7-annotated proteins)
   - `*.VARgenes.aa.fasta` (finalized VAR protein sequences)
   - `*.VARgenes.na.fasta` (finalized VAR nucleotide sequences)
   - BLAST results file

### Prerequisites

- Linux/Unix environment
- Nextflow (v25.10.4 or higher)
- Conda

### Step 1: Clone the Repository

If this is a new installation:

```bash
cd /path/to/workspace
git clone <repository-url>
cd denovarna
```

### Step 2: Install Dependencies

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate denovo-rna-var
```

### Step 3: Check Reference Data

Ensure all required reference data exists in the `data/` directory:

```bash
ls -la data/Reference/           # Reference genomes
ls -la data/HmmerDomains/        # HMM profiles for domains
ls -la data/Proteins_3D7/        # 3D7 reference proteins
ls -la data/varDB/               # VAR gene database
```

Required reference files include:

- `PlasmoDB-67_Pfalciparum3D7_Genome.fasta` - P. falciparum 3D7 genome
- `PlasmoDB-67_Pfalciparum3D7.gff3` - Genome annotation
- `var3kb_exon1.fasta` or similar - VAR exon 1 sequences
- `hg38.fasta` - Human reference genome
- HMMER domain profiles (`.hmm` files)
- 3D7 protein database

---

## Input Requirements

### Sample Data

Place your **trimmed, paired-end FASTQ files** in the `data/trimmed_reads/` directory:

**Important:** Reads should be:

- **Paired-end** (two files per sample: R1 and R2)
- **Trimmed** (adapter sequences removed, low-quality bases trimmed)
- **Named consistently:** `{sample_name}_R1.fastq.gz` and `{sample_name}_R2.fastq.gz`

If your files aren't trimmed yet, use Trimmomatic or similar before running this pipeline.

### Sample List

Create a sample list file (`data/sample_list.txt`) with one sample name per line (without the `_R1`/`_R2` suffix):

```
sample1
sample2
sample3
FTSSSKL
```

### Directory Structure

Your working directory should look like:

```
denovar/
├── main.nf
├── nextflow.config
├── environment.yml
├── processes/
│   ├── identify_non_core_reads.nf
│   ├── Var_assembly.nf
│   ├── annotation_FunctionalAnnotatio.nf
│   └── *.pl/py (custom scripts)
├── data/
│   ├── trimmed_reads/          ← Your input FASTQ files
│   │   ├── sample1_R1.fastq.gz
│   │   ├── sample1_R2.fastq.gz
│   │   └── ...
│   ├── sample_list.txt         ← List of sample names
│   ├── Reference/
│   │   ├── PlasmoDB-67_Pfalciparum3D7_Genome.fasta
│   │   ├── PlasmoDB-67_Pfalciparum3D7.gff3
│   │   ├── hg38.fasta
│   │   ├── var3kb_exon1.fasta
│   │   └── ...
│   ├── HmmerDomains/
│   │   ├── DBL_CIDR_NTS_ATS.hmmdb
│   │   ├── Subdomain.hmm
│   │   └── ...
│   └── Proteins_3D7/
│       ├── Pf3D7.proteins.fasta
│       └── ...
└── results/                     ← Output directory (created by pipeline)
```

---

## Running the Pipeline

### Configuration

Before running, check `nextflow.config` and customize parameters if needed:

### Basic Execution

#### Local Machine (Single CPU):

```bash
nextflow run main.nf -profile local
```

#### HPC/Cluster (SLURM):

```bash
nextflow run main.nf -profile cluster
```

### Important Output Files Explained

| File                                   | Phase | Description                                                  |
| -------------------------------------- | ----- | ------------------------------------------------------------ |
| `*_polymorphic_reads_R1/R2.fq`         | 1     | Extracted VAR gene reads (input for Phase 2)                 |
| `*_unique_contigs_300bp.fasta`         | 2     | Assembled, deduplicated contigs (≥300 bp)                    |
| `*_genome.fasta`                       | 2     | **Main assembly output** - corrected VAR gene sequences      |
| `*_mapped_to_unique_contigs.bam`       | 2     | Read mapping for validation                                  |
| `*.genes.aa.fasta`                     | 3     | Predicted protein sequences from Augustus                    |
| `*Motif_counts.txt`                    | 3     | LARSFADIG and DYVPQYLRW motif detection results              |
| `*Domain.association.structured.txt`   | 3     | Summary of domain                                            |
| `*Subomain.association.structured.txt` | 3     | Summary of subdomain                                         |
| `*.VARgenes.aa.fasta`                  | 4     | **Final annotated VAR proteins** - high-confidence output    |
| `*.VARgenes.na.fasta`                  | 4     | **Final annotated VAR nucleotides** - high-confidence output |

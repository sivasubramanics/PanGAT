# PanG2P

**PanG2P** is a pangenome-informed genome-to-phenome Java framework for trait mapping and prediction. It integrates high-performance tools such as **KMC** and **MMseqs2** to extract single-copy anchor regions, quantify sample k-mer coverage, compile cohort matrices, and generate filtered genotype matrices for downstream GWAS and phenotypic analyses.

---

## Prerequisites

Ensure the following dependencies are installed and available in your environment (`PATH`):

- **Java Runtime Environment (JRE)**: Java 11 or higher
- **KMC**: For k-mer counting and single-copy region extraction
- **MMseqs2**: For sequence clustering and coverage analysis
- **R** *(Optional)*: Required if using the `--gapit` flag during matrix export

---

## Installation
```bash
# Create a new conda environment and install dependencies
conda create -n pang2p openjdk kmc mmseqs2 r-base
# Activate the environment
conda activate pang2p
# Run the downloaded PanG2P JAR file to verify installation
java -jar pang2p-0.0.1.jar -h
```

## Quick Start Pipeline

A standard PanG2P workflow proceeds in four main steps:

```bash
# 1. Build Reference Anchors
java -jar pang2p.jar build-ref -g genomes.csv -b ./pang2p_output -t 16

# 2. Genotype Samples
java -jar pang2p.jar genotype -b ./pang2p_output -s samples.csv -t 16

# 3. Compile Cohort Matrix
java -jar pang2p.jar cohort -b ./pang2p_output

# 4. Export & Filter Genotype Matrix
java -jar pang2p.jar export-matrix \
  -p ./pang2p_output/cohort.pan_matrix.tsv.gz \
  -o ./pang2p_output/matrix_export \
  --gapit

```
---
## Input File Formats

### 1. Genomes Mapping File (`genomes.csv`)

Used by `build-ref` to map genome identifiers to FASTA file locations.

```csv
id,path
genome_01,/path/to/genomes/genome_01.fasta
genome_02,/path/to/genomes/genome_02.fasta
genome_03,/path/to/genomes/genome_03.fasta
```

> **Note:**
> * Header row must be exactly `id,path`.
> * Each genome FASTA must contain the same number of chromosomes or contigs, ordered consistently across all files.


### 2. Samples Mapping File (`samples.csv`)

Used by `genotype` to map sample identifiers to sequencing read files (FASTQ/FASTA, plain or `.gz`).

```csv
id,path
sample_01,/path/to/reads/sample_01.fastq.gz
sample_02,/path/to/reads/sample_02.fastq.gz
```

---

## Command Reference

### Global Options

```bash
java -jar pang2p.jar [-hV] [COMMAND]
```
```bash
  -h, --help        : Show help message and exit.
  -V, --version     : Print version information and exit.
```
---

### Pipeline Commands

#### 1. `build-ref`

Prepares and builds pan-anchor reference sequences across input genomes.

```bash
java -jar pang2p.jar build-ref -g genomes.csv -b ./pang2p_output [options]

```

| Option | Short | Type | Default | Description |
| --- | --- | --- | --- | --- |
| `--genomes` | `-g` | File | *Required* | Path to input genome mapping CSV (`id,path`) |
| `--basedir` | `-b` | Dir | *Required* | Base output directory |
| `--threads` | `-t` | Int | `2` | Number of CPU threads |
| `--kmer` | `-k` | Int | `31` | K-mer size for KMC analysis |
| `--min_count` | `-m` | Int | `1` | Minimum k-mer count threshold for single-copy region extraction |
| `--min_seq_identity` | — | Float | `0.8` | Minimum sequence identity threshold for MMseqs2 clustering |
| `--min_coverage` | — | Float | `0.8` | Minimum alignment coverage threshold for MMseqs2 clustering |
| `--cov_mode` | — | Int | `0` | MMseqs2 coverage mode selection |
| `--min_length` | — | Int | `100` | Minimum sequence length for extracted single-copy regions |
| `--force` | — | Flag | `false` | Overwrite existing output directory |

---

#### 2. `genotype`

Genotypes sample reads against reference pan-anchors.

```bash
java -jar pang2p.jar genotype -b ./pang2p_output -s samples.csv [options]

```

| Option | Short | Type | Default | Description |
| --- | --- | --- | --- | --- |
| `--basedir` | `-b` | Dir | *Required* | Base directory created during `build-ref` |
| `--samples` | `-s` | File | *Required* | Path to sample mapping CSV (`id,path`) |
| `--threads` | `-t` | Int | `2` | Number of threads to use |
| `--kmer` | `-k` | Int | `31` | K-mer size for KMC |
| `--min_count` | `-m` | Int | `1` | Minimum k-mer count for coverage calculation |
| `--memory` | — | Flag | `false` | Load KMC database into memory |
| `--keep-kmc` | — | Flag | `false` | Keep KMC database files after processing |

----

#### 3. `export-matrix`

Exports normalized coverage matrices and discrete genotype matrices from the mapped k-mer matrix.

```bash
java -jar pang2p.jar export-matrix -p <panMatrixFile> -o <outPrefix> [options]

```

| Option | Short | Type | Default | Description |
| --- | --- | --- | --- | --- |
| `--panmatrix` | `-p` | File | *Required* | Path to input pan-genome matrix file (`cohort.pan_matrix.tsv.gz`) |
| `--outprefix` | `-o` | String | *Required* | Prefix path and filename for output files |
| `--include` | `-i` | File | Optional | File containing list of genotype IDs to include |
| `--exclude` | `-e` | File | Optional | File containing list of genotype IDs to exclude |
| `--min-maf` | — | Float | `0.05` | Minimum Minor Allele Frequency filter |
| `--max-missing-marker` | — | Float | `0.1` | Maximum missing data proportion per marker |
| `--max-heterozygosity` | — | Float | `0.5` | Maximum heterozygosity threshold |
| `--cv` | — | Float | `0.0` | Coefficient of Variation cutoff for filtering markers |
| `--gapit` | — | Flag | `false` | Convert output to GAPIT input format (requires R) |

---

### Utility Commands

* **`cohort`**: Add and compile per-sample k-mer coverage into a single cohort matrix..
* **`anchors`**: Extract single-copy regions (anchors) for a single genome.
* **`kmer-coverage`**: Calculate k-mer coverage for a single sample read set.

---

## Output Directory Structure

```text
pang2p_output/
├── 00_prep_ref/                    # Reference preparation files
│   ├── 00_genomes/                 # Input genome FASTA copies/formatting
│   └── 01_pan_anchors/             # Extracted pangenome anchor sequences
├── 01_genotyping/                  # Per-sample k-mer coverage output
│   └── cohort.pan_matrix.tsv.gz    # Compiled matrix generated by `cohort`
└── genomesMap.csv                  # Renamed chromosome and contig mapping file
```
---
## License

See the [LICENSE](https://www.google.com/search?q=LICENSE) file for details on licensing.

---

## GenAI Involvement

This project was developed with the assistance of GenAI tools (Copilot), which helped only with drafting documentation and refactoring code. The final content has been reviewed and verified.

---

## Contact
For questions, issues, or feature requests, please open an issue on the [GitHub repository](https://github.com/sivasubramanics/pang2p)
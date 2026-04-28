# Session 02 – Searching and Downloading Sequencing Data from NCBI SRA

> **Course:** Bioinformatics B2S SS2026
> 
> **Topics:** NCBI SRA database · SRA Explorer · SRA Toolkit · wget · FASTQ format
> 
> **Prerequisites:** Session 01 (Linux basics, connecting to the server)

---

## Table of Contents

1. [Background: What is the SRA?](#1-background-what-is-the-sra)
2. [The FASTQ Format](#2-the-fastq-format)
3. [Part A – Searching the SRA via the Web (SRA Explorer)](#3-part-a--searching-the-sra-via-the-web-sra-explorer)
4. [Part B – Downloading with wget](#4-part-b--downloading-with-wget)
5. [Part C – Downloading with SRA Toolkit](#5-part-c--downloading-with-sra-toolkit)
6. [Part D – Inspecting your downloaded data](#6-part-d--inspecting-your-downloaded-data)
7. [Practical Exercises](#7-practical-exercises)
8. [Quick Reference Cheat Sheet](#8-quick-reference-cheat-sheet)

---
## Before start wirte the following command

```bash
tmux
```

## 1. Background: What is the SRA?

The **Sequence Read Archive (SRA)** is a public database maintained by the NCBI (National Center for Biotechnology Information). It is the largest repository of raw sequencing data in the world.

Whenever researchers publish a study that involves sequencing (e.g. 16S amplicon, metagenomics, RNA-seq, whole genome), they are required to deposit their raw sequencing files here so anyone can reproduce or build on their work.

### Key concepts

| Term | Meaning |
|------|---------|
| **BioProject** | The umbrella entry for a study (e.g. `PRJNA…`) |
| **BioSample** | Metadata for a single biological sample |
| **SRR / ERR / DRR** | A run accession – one actual sequencing file or pair of files |
| **FASTQ** | The file format storing raw sequencing reads + quality scores |

### Today's demo datasets

| Dataset | Accession | Description |
|---------|-----------|-------------|
| PhiX control | `SRR13457500` | Illumina PhiX174 spike-in control – tiny, fast to download |
| 16S microbiome | `SRR2726675` | Human gut 16S rRNA amplicon (V4 region) – Illumina MiSeq |
| 16S microbiome | `SRR27375990` | Human gut comparison Nanopore Illumina |

> 💡 **Why PhiX?** PhiX is a bacteriophage genome commonly used as a sequencing control. Its reads are short, the dataset is tiny (~50 MB), and it downloads in seconds — perfect for testing your tools before working with larger datasets.

---

## 2. The FASTQ Format

Before downloading anything, it helps to know what you are downloading.

A FASTQ file stores sequencing reads. Each read takes up exactly **4 lines**:

```
@SRR13457500.1 1 length=151
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC...
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF...
```

| Line | Content |
|------|---------|
| Line 1 | `@` + read name / identifier |
| Line 2 | DNA sequence |
| Line 3 | `+` (separator) |
| Line 4 | Quality scores (one character per base, ASCII encoded) |

Quality scores use the **Phred scale**: a score of 30 means a 1-in-1000 chance of a sequencing error (99.9% accuracy). Higher = better.

For **paired-end** sequencing there are two files:
- `*_1.fastq` — forward reads (Read 1)
- `*_2.fastq` — reverse reads (Read 2)

---

## 3. Part A – Searching the SRA via the Web (SRA Explorer)

**SRA Explorer** is a user-friendly website that makes finding and downloading SRA data much easier than the NCBI interface.

🔗 [https://sra-explorer.info](https://sra-explorer.info)

### Step-by-step: Find a 16S microbiome dataset

1. Open [https://sra-explorer.info](https://sra-explorer.info) in your browser.
2. In the search box, type a BioProject accession or a keyword, for example:

   ```
   PRJNA310034
   PRJNA1058165 
   ```

   or search by topic:

   ```
   human gut 16S microbiome MiSeq
   ```

3. Press **Enter** or click the search icon.
4. You will see a list of SRR runs. Click the **shopping cart** icon next to the ones you want.
5. Click **"Saved datasets"** (top right) to see your selection.
6. Click **"Raw FastQ Download URLs"** to get direct download links — these are plain `https://` links you can use with `wget`.

> 💡 SRA Explorer also generates ready-made `wget` and `curl` commands you can copy directly into your terminal.

### Step-by-step: Look up a known accession on NCBI

If you already have an accession number (e.g. from a paper's supplementary data):

1. Go to [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)
2. Paste the accession (e.g. `PRJNA1058165 or SRR27375990`) into the search bar and press Enter.
3. Click on the result to see full metadata: organism, instrument, library strategy, number of reads, file size.
4. Note the accession — you will use it in the terminal in the next steps.

---

## 4. Part B – Downloading with `wget`

`wget` is a simple Linux command for downloading files from the internet using a URL. It is the fastest approach when you already have a direct download link (e.g. from SRA Explorer).

### Setup: Create your working directory

Log in to the server, then:

```bash
mkdir -p ~/session01/database
cd ~/session01/database
```

### Download the PhiX dataset

SRA stores files on Amazon S3 and ENA (European Nucleotide Archive) servers. SRA Explorer gives you these direct links. For our PhiX example:

```bash
# Download Read 1 and Read 2 of the PhiX run
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR37888883/SRR37888883

# Or use the ENA FTP mirror (often faster in Europe):
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR246/071/SRR24695071/SRR24695071_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR246/071/SRR24695071/SRR24695071_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/DR000/SRR37888883/SRR37888883_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/DR000/SRR37888883/SRR37888883_2.fastq.gz
```

> 🌍 **Tip for European servers:** The ENA FTP mirror is usually faster than the NCBI servers when you are working from a European institution.

### Download the 16S microbiome dataset

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR272/005/SRR2726675/SRR2726675_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR272/005/SRR2726675/SRR2726675_2.fastq.gz
```

### Check what was downloaded

```bash
ls -lh
```

Example output:

```
-rw-r--r-- 1 user user  47M Apr  4 10:12 SRR13457500_1.fastq.gz
-rw-r--r-- 1 user user  49M Apr  4 10:12 SRR13457500_2.fastq.gz
-rw-r--r-- 1 user user 120M Apr  4 10:15 SRR2726675_1.fastq.gz
-rw-r--r-- 1 user user 118M Apr  4 10:15 SRR2726675_2.fastq.gz
```

The `.gz` extension means the files are **gzip-compressed** — this is normal and expected. Most bioinformatics tools accept `.fastq.gz` directly without needing to decompress them.

---

## 5. Part C – Downloading with SRA Toolkit

The **SRA Toolkit** is the official NCBI software for downloading SRA data. It gives you more control and works entirely from the accession number — no URL needed.

### Install SRA Toolkit on Ubuntu 24.04

```bash
# Download the latest SRA Toolkit for Ubuntu
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

# Extract it
tar -xzf sratoolkit.current-ubuntu64.tar.gz

# Move it to a permanent location
mv sratoolkit.*-ubuntu64 ~/sratoolkit

# Add the binaries to your PATH for this session
export PATH=$PATH:~/sratoolkit/bin

# Verify the installation
fastq-dump --version
```

Expected output:

```
fastq-dump : 3.x.x
```

> 💡 To make the PATH change permanent, add the `export` line to your `~/.bashrc` file using `nano ~/.bashrc`.

### Configure SRA Toolkit (first time only)

```bash
vdb-config --interactive
```

This opens a configuration menu. Press `X` to exit with default settings — this creates the necessary config files.

### Step 1: `prefetch` — Download the SRA file

`prefetch` downloads the compressed `.sra` file to your local cache:

```bash
cd ~/session01/database

# Download PhiX
 prefetch SRR37888883

# Download 16S microbiome run
prefetch SRR2726675

# downalod multi files
prefetch --option-file SraAccList.txt --output-directory ./sra --max-size 100G

```

You will see a progress bar. The `.sra` files are saved in `~/ncbi/public/sra/` by default.

### Step 2: `fasterq-dump` — Convert to FASTQ

`fasterq-dump` extracts the reads from the `.sra` file and writes them as FASTQ:

```bash
# Convert PhiX
 /courses/master_b2s/software/sratoolkit/bin/fasterq-dump SRR37888883 --outdir ~/session01/data/ --split-files

# Convert 16S microbiome
fasterq-dump SRR2726675 --outdir ~/session02/database/ --split-files

```

```bash

# Convert multi files
## flatten the directory (because prefetch nests them):
find ./sra -name "*.sra" -exec mv {} ./sra \;
find ./sra -type d -empty -delete
## Now all .sra files are in ./sra.

# Use fasterq-dump (recommended over fastq-dump):
fasterq-dump ./sra/*.sra -O ./fastq -e 8

# You can use seqtk (fast and simple) to convert sra to fastq
for f in ./fastq/*.fastq; do
    seqtk seq -a "$f" > "${f%.fastq}.fasta"
done



# OTher direct solution
## Use fasterq-dump directly from accession list:
cat SraAccList.txt | xargs -n 1 -P 8 fasterq-dump -O ./fastq

## -P 8 = parallel downloads (huge speedup)

# OTher direct solution
while read acc; do
    fastq-dump --split-files "$acc"
done < accession.txt

```

| Flag | Meaning |
|------|---------|
| `--outdir` | Where to write the output FASTQ files |
| `--split-files` | Write Read 1 and Read 2 into separate files (paired-end) |

### (Optional) Compress the output

`fasterq-dump` writes uncompressed files. To save disk space, compress them:

```bash
gzip ~/session02/data/SRR13457500_1.fastq
gzip ~/session02/data/SRR13457500_2.fastq
gzip ~/session02/data/SRR2726675_1.fastq
gzip ~/session02/data/SRR2726675_2.fastq
```

---

## 6. Part D – Inspecting your downloaded data

Once your files are downloaded, always do a quick sanity check before running any analysis.

### Check file sizes

```bash
ls -lh ~/session02/data/
```

### Count the number of reads

Each read in a FASTQ file takes exactly 4 lines. So the number of reads = total lines ÷ 4.

```bash
# For a compressed file:
zcat SRR13457500_1.fastq.gz | wc -l
zcat SRR37888883_1.fastq.gz | wc -l | awk '{print $1 / 4}'
# Divide the result by 4 to get the number of reads
```

Or use a one-liner:

```bash
echo $(( $(zcat SRR37888883_1.fastq.gz | wc -l) / 4 )) reads
```

### Peek at the first few reads

```bash
# View the first 12 lines (= first 3 reads)
zcat SRR37888883_1.fastq.gz | head -12
```

Example output:

```
@SRR13457500.1 1 length=151
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF
@SRR13457500.2 2 length=151
GCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SRR13457500.3 3 length=151
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```

### Check read length distribution

```bash
zcat SRR37888883_1.fastq.gz | awk 'NR%4==2 {print length($0)}' | sort | uniq -c
```

This prints how many reads have each length — useful for spotting truncated reads or adapter contamination.

### Confirm paired-end files have the same number of reads

```bash
echo "R1: $(( $(zcat SRR37888883_1.fastq.gz | wc -l) / 4 )) reads"
echo "R2: $(( $(zcat SRR37888883_2.fastq.gz | wc -l) / 4 )) reads"
```

Both numbers must be identical. If they differ, the file may be corrupted.

---

## 7. Practical Exercises

Work through these exercises on the server. Use the commands above as reference.

---

### Exercise 1 – Set up your workspace

```bash
mkdir -p ~/session02/data
cd ~/session02/data
pwd
```

---

### Exercise 2 – Explore SRA Explorer

1. Open [https://sra-explorer.info](https://sra-explorer.info) in your browser.
2. Search for the BioProject **`PRJNA310034`**.
3. Find the run **`SRR2726675`** and note:
   - The organism
   - The sequencing instrument
   - The library strategy
   - The number of reads and file size

> ❓ **Question:** Is this single-end or paired-end sequencing? How can you tell?

---

### Exercise 3 – Download PhiX data with wget

```bash
cd ~/session02/data

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/000/SRR13457500/SRR13457500_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/000/SRR13457500/SRR13457500_2.fastq.gz

ls -lh
```

> ❓ **Question:** How large are the files? Are they compressed?

---

### Exercise 4 – Inspect the PhiX data

```bash
# Look at the first 8 lines
zcat SRR13457500_1.fastq.gz | head -8

# Count the reads
echo $(( $(zcat SRR13457500_1.fastq.gz | wc -l) / 4 )) reads

# Confirm R1 and R2 have equal read counts
echo "R1: $(( $(zcat SRR13457500_1.fastq.gz | wc -l) / 4 )) reads"
echo "R2: $(( $(zcat SRR13457500_2.fastq.gz | wc -l) / 4 )) reads"
```

> ❓ **Question:** What is the read length? What do the quality score characters mean?

---

### Exercise 5 – Download the 16S microbiome dataset with SRA Toolkit

```bash
# Make sure your PATH is set
export PATH=$PATH:~/sratoolkit/bin

cd ~/session02/data

prefetch SRR2726675
fasterq-dump SRR2726675 --outdir . --split-files

ls -lh
```

---

### Exercise 6 – Compare the two datasets

Fill in the table below for your report:

| Property | PhiX (`SRR13457500`) | 16S microbiome (`SRR2726675`) |
|----------|----------------------|-------------------------------|
| Number of reads (R1) | | |
| Read length | | |
| File size (compressed) | | |
| Organism / source | | |
| Single-end or paired-end? | | |

---

## 8. Quick Reference Cheat Sheet

| Command | What it does |
|---------|-------------|
| `wget <URL>` | Download a file from a URL |
| `prefetch <SRR>` | Download an SRA file by accession |
| `fasterq-dump <SRR> --split-files` | Convert SRA to FASTQ (paired-end) |
| `ls -lh` | List files with human-readable sizes |
| `zcat file.fastq.gz \| head -8` | Preview first 2 reads of a compressed FASTQ |
| `zcat file.fastq.gz \| wc -l` | Count lines in a compressed file |
| `gzip file.fastq` | Compress a FASTQ file |
| `zcat f.fastq.gz \| awk 'NR%4==2{print length}' \| sort \| uniq -c` | Read length distribution |

### Useful links

| Resource | URL |
|----------|-----|
| SRA Explorer | https://sra-explorer.info |
| NCBI SRA | https://www.ncbi.nlm.nih.gov/sra |
| ENA Browser | https://www.ebi.ac.uk/ena/browser/home |
| SRA Toolkit docs | https://github.com/ncbi/sra-tools/wiki |

---

*Next session: Quality control with FastQC and Trimmomatic 🔬*

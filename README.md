# Recovery and Phylogenetic Analysis of Mitochondrial Cytochrome b (cytb) Sequences from Rodent Samples

A reproducible command-line pipeline for the targeted recovery, assembly, and phylogenetic placement of mitochondrial *cytochrome b* (*cytb*) sequences from shotgun or capture-based Illumina paired-end sequencing of rodent samples.

---

## Overview

This pipeline takes raw paired-end Illumina FASTQ files and a reference *cytb* sequence, and proceeds through four logical stages:

1. **Read quality assessment** — FastQC inspection of raw reads.
2. **Reference-guided read capture** — BWA-MEM alignment to a *cytb* reference, followed by stringent extraction of on-target reads with SAMtools.
3. **De novo assembly and contig validation** — MEGAHIT assembly of mapped reads, QUAST-based assembly assessment, and BLASTn confirmation of *cytb* identity.
4. **Phylogenetic reconstruction** — Extraction of the assembled *cytb* contig, concatenation with reference taxa, MAFFT alignment with strand-aware adjustment, trimAl trimming, and maximum-likelihood inference with IQ-TREE.

The design philosophy is deliberately conservative: reads are first restricted to the mitochondrial target by mapping, assembled *de novo* from that enriched subset, and the resulting contig is only accepted into the phylogeny after sequence-similarity confirmation against the reference. This minimises the inclusion of nuclear mitochondrial pseudogenes (*numts*) or off-target contigs — a well-known pitfall in short-read mitochondrial recovery from rodent material.

---

## Workflow

```
Raw PE FASTQ
     │
     ▼
FastQC (quality check)
     │
     ▼
bwa mem ──► samtools view/sort ──► flagstat (mapping QC)
     │
     ▼
samtools view -F 4 (keep mapped reads only)
     │
     ▼
samtools fastq (extract mapped PE reads)
     │
     ▼
MEGAHIT (de novo assembly)
     │
     ▼
QUAST (assembly metrics)          BLASTn vs cytb reference (identity check)
     │
     ▼
seqkit grep (extract target contig)
     │
     ▼
Concatenate with reference taxa
     │
     ▼
MAFFT --adjustdirection ──► trimAl ──► IQ-TREE (-m MFP -B 1000)
     │
     ▼
Maximum-likelihood phylogeny
```

---

## Prerequisites

| Tool | Purpose | Suggested version |
|---|---|---|
| FastQC | Raw read quality assessment | ≥ 0.11.9 |
| BWA | Read alignment (BWA-MEM) | ≥ 0.7.17 |
| SAMtools | BAM filtering, sorting, indexing, read extraction | ≥ 1.15 |
| MEGAHIT | Metagenomic de novo assembler | ≥ 1.2.9 |
| QUAST | Assembly quality metrics | 5.3.0 |
| NCBI BLAST+ | Sequence identity confirmation | ≥ 2.12 |
| Prodigal | Open reading frame prediction (optional; see caveats) | ≥ 2.6.3 |
| seqkit | FASTA/Q manipulation | ≥ 2.3 |
| MAFFT | Multiple sequence alignment | ≥ 7.475 |
| trimAl | Alignment trimming | ≥ 1.4 |
| IQ-TREE | Maximum-likelihood phylogenetic inference | ≥ 2.0 |

A Conda environment (`cytb_phylo`) is used to isolate dependencies:

```bash
conda activate cytb_phylo
```

---

## Input files

| File | Description |
|---|---|
| `ZA7L_S11_L001_R1_001.fastq.gz` / `R2_001.fastq.gz` | Paired-end Illumina reads for one sample (per-sample run; loop or parallelise across samples) |
| `cytb_db.fa` | Reference *cytb* sequence used for read mapping |
| `cytb_12.fa` | Reference *cytb* sequences from related taxa (12 sequences) used in the final phylogeny |

---

## Pipeline steps

### Stage 1 — Quality control

```bash
fastqc   # inspect per-base quality, adapter content, duplication
```

### Stage 2 — Reference-guided capture of mitochondrial reads

```bash
bwa index cytb_db.fa

bwa mem -t 12 cytb_db.fa ZA7L_S11_L001_R1_001.fastq.gz ZA7L_S11_L001_R2_001.fastq.gz \
  | samtools view -b - \
  | samtools sort -@ 12 -m 512M -o mapped.sorted.bam

samtools flagstat mapped.sorted.bam                      # overall mapping statistics

samtools view -b -F 4 mapped.sorted.bam > mapped_only.bam  # retain mapped reads
samtools flagstat mapped_only.bam                          # post-filter statistics
samtools sort mapped_only.bam -o mapped_only.sorted.bam
samtools index mapped_only.sorted.bam
samtools idxstats mapped_only.sorted.bam > idxstats.txt    # per-reference read counts

samtools fastq mapped_only.sorted.bam \
  -1 mapped_R1.fastq.gz -2 mapped_R2.fastq.gz \
  -0 /dev/null -s /dev/null -n                             # export mapped PE reads
```

The `-F 4` flag discards unmapped reads; `-n` in `samtools fastq` prevents `/1` `/2` suffixes from being appended, keeping read names compatible with the assembler.

### Stage 3 — De novo assembly and validation

```bash
mkdir assembly && cd assembly
mv ../mapped_R1.fastq.gz ../mapped_R2.fastq.gz .

megahit -1 mapped_R1.fastq.gz -2 mapped_R2.fastq.gz -o megahit_out
cd megahit_out

cat final.contigs.fa                                       # inspect contigs

python3 ~/tools/quast-5.3.0/quast.py final.contigs.fa -o quast_output

mkdir blast && cd blast
mv ../final.contigs.fa .
cp /home/vmsowoya/mbewa_study/Version_2/ZA7L/cytb_12.fa .

makeblastdb -in cytb_12.fa -dbtype nucl -out cytb_db
blastn -query final.contigs.fa -db cytb_db -out cytb_vs_contigs.tsv \
  -outfmt 6 -evalue 1e-10 -num_threads 8
```

The BLASTn tabular output (`-outfmt 6`) identifies which assembled contig represents *cytb* and provides the percent identity and coverage needed to accept or reject it.

### Stage 4 — Phylogenetic reconstruction

```bash
prodigal -i final.contigs.fa -a proteins.faa -d genes.fna -o prodigal.gff -p meta   # optional ORF annotation

mkdir phylogeny && cd phylogeny
cp ../final.contigs.fa .

seqkit grep -r -p "^k79_1\b" final.contigs.fa > ZA7L_cytb_extract.fa
# Rename the extracted contig header to the rodent sample ID (e.g. ZA7L)

cp /home/vmsowoya/mbewa_study/Version_2/ZA7L/cytb_12.fa .
cat ZA7L_cytb_extract.fa cytb_12.fa > cytb_13.fa

mafft --adjustdirection --auto cytb_13.fa > aligned.fa
trimal -in aligned.fa -out trimmed.fa
sed 's/^>_R_/>/' cleaned_13.fa < /dev/null 2>/dev/null || sed 's/^>_R_/>/' trimmed.fa > cleaned_13.fa

iqtree -s cleaned_13.fa -m MFP -B 1000 -nt AUTO
```

Key choices:

- **`mafft --adjustdirection`** automatically reverses sequences oriented in the opposite direction relative to the majority — essential when contigs may assemble on either strand.
- **`-m MFP`** lets IQ-TREE select the best-fit nucleotide substitution model (e.g. GTR+F+I+G4 for *cytb*).
- **`-B 1000`** applies 1,000 ultrafast bootstrap replicates for node support.

### Outputs

| File | Description |
|---|---|
| `mapped.sorted.bam` / `flagstat` | Mapping statistics for the full library |
| `idxstats.txt` | Reads per reference sequence |
| `megahit_out/final.contigs.fa` | Assembled contigs from mapped reads |
| `quast_output/` | Assembly metrics (N50, total length, contig counts) |
| `cytb_vs_contigs.tsv` | BLASTn hits of contigs against the *cytb* reference |
| `cleaned_13.fa` | Final trimmed alignment (sample + 12 references) |
| `IQ-TREE` outputs (`*.treefile`, `*.iqtree`, `*.contree`) | Best-scoring ML tree, model summary, and consensus tree |

---

## Running multiple samples

The commands above are written for a single sample (`ZA7L`). For a cohort, wrap Stages 2–4 in a loop over sample IDs, renaming the extracted contig header to the sample ID before concatenation so that all samples can be combined into a single alignment prior to `mafft`.

---

## Important caveats

1. **Prodigal is not appropriate for mitochondrial genomes.** Prodigal assumes bacterial genetic code and codon usage; vertebrate mitochondria use an alternative genetic code (e.g. TGA = tryptophan, AGA/AGG = stop). For *cytb* phylogenetics, gene calling is unnecessary — the BLASTn-confirmed nucleotide contig is sufficient. If protein-level annotation is required, use a mitochondria-aware annotator (MFannot, MITOS, or `prodigal` alternatives supporting `-g 2`-style mito codes).
2. **Hard-coded contig names are fragile.** `seqkit grep -p "^k79_1\b"` depends on MEGAHIT's contig naming, which changes between runs and samples. Prefer selecting the longest contig with a significant BLASTn hit to the *cytb* reference.
3. **Read trimming is not included.** Consider adding adapter/quality trimming (fastp or Trimmomatic) after FastQC if per-base quality warrants it, and report this in methods.
4. **trimAl default settings.** Consider `-automated1` for a reproducible, publication-defensible trimming heuristic, and report the final alignment length.
5. **Outgroup rooting.** IQ-TREE will root the tree at the first sequence by default; specify an explicit outgroup with `-o` if a rooted depiction is needed.
6. **numt caution.** Short-read data can recover nuclear mitochondrial pseudogenes. BLASTn confirmation plus inspection of coverage evenness and stop-codon/frameshift checks in the aligned *cytb* reading frame are recommended before phylogenetic placement.

---

## Citation

If you use this pipeline, please cite the underlying tools: BWA (Li & Durbin, 2009), SAMtools (Danecek et al., 2021), MEGAHIT (Li et al., 2015), QUAST (Mikheenko et al., 2018), BLAST+ (Camacho et al., 2009), seqkit (Shen et al., 2016), MAFFT (Katoh & Standley, 2013), trimAl (Capella-Gutiérrez et al., 2009), and IQ-TREE 2 (Minh et al., 2020).

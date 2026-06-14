# PacBio Assembly and Annotation of *Helicops angulatus*

This repository documents the complete bioinformatic pipeline used to assemble, scaffold, and annotate the reference genome of the South American water snake, *Helicops angulatus* (Specimen IAvH-10017).

---

## 1. K-mer Profiling & Genome Characterization

### Objectives
Prior to running the main assembly graph algorithms, we performed an alignment-free **$k$-mer analysis** ($k=21$) using raw PacBio HiFi reads. The primary goals were to:
* Empirically estimate sequencing depth (coverage).
* Quantify the core **heterozygosity rate** of the individual.
* Preview the genomic repeat landscape to calibrate downstream assembly parameters.

---

### Workflow & Scripts

Canonical 21-mers were counted using **Jellyfish v2**, and the frequency distribution was exported as a histogram. To resolve the high-coverage profile of modern PacBio HiFi data, the maximum $k$-mer parameter was truncated to `1000` to ensure proper mathematical convergence.

#### SLURM Job Script: Jellyfish Count & Histo

```bash
#!/bin/bash
#SBATCH --job-name=jellyfish_angulatus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/jellyfish/jellyfish.%j.out
#SBATCH --error=/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/jellyfish/jellyfish.%j.err

source ~/.bash_profile
conda activate assembly

############################
# USER PARAMETERS
############################
READS="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/Data-X202SC26050352-Z01-F001/RP_10017/Revio/FPAC260238538-1A/RP_10017.hifi_reads.fastq.gz"
SAMPLE="Helicops_angulatus_10017"
KMER=21
HASH_SIZE="10G"
THREADS=${SLURM_CPUS_PER_TASK}
OUTDIR="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/jellyfish"

# K-mer counting
jellyfish count -C -m ${KMER} -s ${HASH_SIZE} -t ${THREADS} <(zcat ${READS}) -o ${OUTDIR}/${SAMPLE}.k${KMER}.jf

# Histogram generation
jellyfish histo -t ${THREADS} ${OUTDIR}/${SAMPLE}.k${KMER}.jf > ${OUTDIR}/${SAMPLE}.k${KMER}.histo
```


<img width="1651" height="844" alt="Screenshot 2026-06-12 at 4 28 59 AM" src="https://github.com/user-attachments/assets/cc714e00-18cf-401e-928a-0b6c7a4bcc00" />

<img width="1633" height="818" alt="Screenshot 2026-06-12 at 10 37 05 AM" src="https://github.com/user-attachments/assets/c2c99c08-8554-450e-84f7-2f6d26499655" />

### Summary of Results

* Homozygous Profile (aa):	85.2%	Proportion of 21-mer blocks sharing identical homologous alleles.

| Genomic Metric | Value from Plot | Notes / Interpretation |
| :--- | :--- | :--- |
| **Homozygous Profile (`aa`)** | **85.2%** | Proportion of $k$-mer blocks sharing identical homologous alleles. |
| **Heterozygous Profile (`ab`)** | **14.8%** | Proportion of 21-mer blocks containing at least one variable site. |
| **K-mer Coverage Peak (`kcov`)**| **57.1$\times$** | High-quality sequencing depth, providing robust assembly power. |
| **Model Sequencing Error Rate** | **0.462%** | Estimated base-calling error rate, typical for highly accurate HiFi data. |
| **Estimated Heterozygosity Rate** | **~0.76%** | Mathematically derived from the $ab$ fraction ($1 - (1 - 0.148)^{1/21}$). |

* The Heterozygous Peak (First Blue Mountain): Distinctly located at half-coverage (~28$\times$), this prominent peak represents polymorphic alleles. Its substantial height relative to the main peak confirms a biologically significant degree of genetic variation (~0.76%).
* The Homozygous Main Peak (Highest Blue Mountain): Centered at $57.1\times$, this dominant peak captures the genomic regions where maternal and paternal chromosomes match identically, establishing our true sequencing coverage depth.
* The estimated model length (~727 Mb) represents a severe underestimation of the true haploid size (~2.16 Gb). This artifact is caused by GenomeScope filtering out high-frequency repetitive strings to run its calculations.


## 2. Contig Assembly with Hifiasm

### Objectives
Following $k$-mer characterization, we initiated the primary *de novo* assembly using **hifiasm**. The goal of this step was to construct highly contiguous, partially phased primary and alternate pseudohaplotypes directly from the PacBio HiFi reads, taking advantage of the empirical sequence variation to resolve homologous chromosomes.

---

### Workflow & Scripts

The assembly was executed on the SLURM cluster. To handle the highly repetitive nature of the genome, we aggressively purged haplotig duplications (using the `-l3` parameter) to prevent homologous regions from assembling as separate, false-duplicate contigs.

The -l flag in hifiasm controls the Level of Purge Duplication. There are four operational levels available:

* l0 (Disable): Turns off purging entirely. Used only for highly inbred/homozygous lines or when providing external trio-binning data (parental short reads).
* l1 (Light): Only purges contained haplotigs. If a smaller contig is completely identical to a sub-region of a larger contig, it gets removed. If a shorter piece of DNA starts and ends entirely inside a larger, almost identical piece of DNA, it assumes the short piece is just the redundant maternal/paternal copy.
* l2 (Moderate): Purges all types of haplotigs, including both contained and overlapping sequences. Expands its search. It now looks for both perfectly contained bubbles AND these dangling, overlapping alternate alleles.
* l3 (Aggressive): Purges all types of haplotigs in the most aggressive mathematical way, lowering the similarity threshold to force the collapse of highly divergent alleles. It forces the algorithm to assume that even if two paths look quite divergent or heavily mutated, as long as their surrounding graph topology looks like they are alleles, they should be collapsed.

If you don't explicitly specify the -l flag in your hifiasm command, the default behavior is -l3 (Aggressive purging).

#### SLURM Job Script: Hifiasm Assembly
```bash
#!/bin/bash
#SBATCH --job-name=hifiasm_angulatus
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=3-00:00:00
#SBATCH --output=hifiasm_H.angulatus_l3_%j.out
#SBATCH --error=hifiasm_H.angulatus_l3%j.err

source ~/.bash_profile
conda activate assembly

cd /home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/hifiasm_assembly/hifiasm_assembly_l3

############################
# USER PARAMETERS
############################
READS="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/Data-X202SC26050352-Z01-F001/RP_10017/Revio/FPAC260238538-1A/RP_10017.hifi_reads.fastq.gz"

############################
# DE NOVO ASSEMBLY
############################
hifiasm \
  -o Helicops_angulatus_IAvH_10017.l3 \
  -t 32 \
  -l3 \
  ${READS} \
  2> Helicops_angulatus_IAvH_10017.hifiasm_l3.log
```

### Summary of Assembly Results (Primary Contigs)

The primary pseudohaplotype (`.p_ctg.fasta`) generated by `hifiasm` using the `-l3` parameter was evaluated for spatial contiguity and gene-space completeness. 

#### 1. Contiguity Metrics
Standard assembly statistics were computed for the primary contigs:

| Metric | Value | 
| :--- | :--- | 
| **Total Assembly Length** | **2,167,390,149 bp (~2.16 Gb)** | 
| **Total Contigs** | **140** | 
| **Largest Contig** | **215,930,598 bp (~215.9 Mb)** | 
| **Contig N50** | **131,633,754 bp (~131.6 Mb)** | 
| **L50** | **7** | 
| **Gaps / N_count** | **0** | 

#### 2. Gene-Space Completeness (BUSCO)
To assess the biological completeness of the assembly, we ran BUSCO against the **Sauropsida (`sauropsida_odb10`)** lineage dataset (n = 7,480 orthologs):

| BUSCO Category | Percentage | Count |
| :--- | :--- | :--- |
| **Complete (C)** | **93.2%** | 6,970 |
| Complete and Single-copy (S) | 91.6% | 6,855 |
| Complete and Duplicated (D) | 1.5% | 115 |
| **Fragmented (F)** | **1.2%** | 87 |
| **Missing (M)** | **5.7%** | 423 |

---

### Result Interpretation

1. **Contiguity:** Achieving a **Contig N50 of ~131.6 Mb** straight out of the assembler (prior to any Hi-C scaffolding) is a great result. An L50 of 7 means that half of the entire 2.16 Gb genome is contained within just 7 massive DNA segments. Given that the largest contig is ~215.9 Mb, `hifiasm` has successfully assembled entire chromosomal arms—or even complete macrochromosomes—telomere-to-telomere without gaps.
2. **Haploid Size Resolution:** The final assembly length is **~2.16 Gb**. This recaptures the massive highly-repetitive fraction of the genome (see below) that GenomeScope 1.0 struggled to model mathematically, confirming that the true physical size of the *Helicops angulatus* haploid genome is >2 Gb.
3. **Validation of the `-l3` Purge Parameter:** The BUSCO duplication rate is low (**1.5%**). 
4. **High Gene Completeness:** Recovering **93.2%** of the core Sauropsida orthologs completely intact indicates a highly accurate and biologically viable reference genome.


## 3. Reference-Free Quality Evaluation (Merqury)

### Objectives
To rigorously evaluate the base-level accuracy and structural completeness of the primary assembly without introducing reference biases. We utilized **Merqury**, a $k$-mer-based evaluation tool, to cross-reference the raw PacBio HiFi sequencing reads against the assembled contigs.

---

### Workflow & Scripts

The evaluation was performed in two steps:
1. **Meryl** was used to build a comprehensive database of all valid 21-mers present in the raw sequencing data.
2. **Merqury** compared this empirical database against the `.p_ctg.fasta` assembly to calculate the Consensus Quality Value (QV) and $k$-mer completeness.

#### SLURM Job Script: Meryl & Merqury
```bash
#!/bin/bash
#SBATCH --job-name=merqury_angulatus
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --output=merqury_%j.out
#SBATCH --error=merqury_%j.err

source ~/.bash_profile
mamba activate merqury_env

WORKDIR="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/hifiasm_assembly/hifiasm_assembly_l3/merqury"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

READS="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/Data-X202SC26050352-Z01-F001/RP_10017/Revio/FPAC260238538-1A/RP_10017.hifi_reads.fastq.gz"
ASM="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/hifiasm_assembly/hifiasm_assembly_l3/Helicops_angulatus_IAvH_10017.l3.bp.p_ctg.fasta"

PREFIX="H_angulatus_primary"
KMER=21

# 1. Build k-mer database
meryl count k=${KMER} threads=${SLURM_CPUS_PER_TASK} memory=140 ${READS} output ${WORKDIR}/reads.meryl

# 2. Run Merqury evaluation
merqury.sh ${WORKDIR}/reads.meryl ${ASM} ${PREFIX}
```

### Results

| Metric | Value | 
| :--- | :--- | 
| **K-mer Completeness** | **95.23%** | 
| **Consensus Quality Value (QV)** | **71.58** | 
| **Calculated Error Rate** | **6.95e-08** | 

* A Consensus QV of 71.58 is a remarkably high accuracy score, highlighting the immense power of PacBio HiFi's circular consensus sequencing (CCS). It guarantees that the nucleotide sequence exceeds 99.99999% accuracy (less than 1 sequence error per 10 million base pairs)
* The score of 95.23% means that the vast majority of the biological information sequenced from the snake was successfully stitched into your primary contigs.


### Primary Assembly Comparison: `-l3` vs. `-l2` assemblies: 

| Metric | `-l3` (Aggressive Purge) | `-l2` (Moderate Purge) | The Difference |
| :--- | :--- | :--- | :--- |
| **Total Length** | 2,167,390,149 bp (~2.16 Gb) | 2,183,292,827 bp (~2.18 Gb) | **+ 15.9 Mb** in `-l2` |
| **Total Contigs** | 140 | 153 | **+ 13 contigs** in `-l2` |
| **Largest Contig** | 215,930,598 bp (~215.9 Mb) | 227,656,864 bp (~227.6 Mb) | **+ 11.7 Mb** in `-l2` |
| **Contig N50** | 131,633,754 bp (~131.6 Mb) | 147,694,758 bp (~147.6 Mb) | **+ 16.0 Mb** in `-l2` |
| **L50** | 7 | 6 | **1 fewer contig** needed in `-l2` |
| **BUSCO Complete** | 93.2% | 93.2% | Identical |
| **BUSCO Single (S)**| 91.6% | 91.5% | - 0.1% in `-l2` |
| **BUSCO Duplicated (D)**| 1.5% | 1.7% | **+ 0.2%** in `-l2` |
| **Merqury Completeness**| 95.23% | 95.28% | **+ 0.05%** in `-l2` |
| **Merqury QV** | 71.57 | 71.51 | Effectively Identical |

## 4. Chromosome-Level Scaffolding (Hi-C)

### Objectives
To elevate the highly contiguous primary assembly to a full chromosome-level reference. By analyzing the 3D physical interactions of the DNA within the nucleus using Hi-C sequencing, we aimed to anchor, order, and orient the massive `.p_ctg` contigs into complete macro- and micro-chromosomes.

### Workflow & Scripts
Hi-C reads were mapped against the primary assembly to generate spatial contact maps. A scaffolding algorithm was then applied to mathematically link the contigs. Finally, the automated scaffolds were visualized as a contact heatmap for manual curation to correct any structural misassemblies or inversions.

#### SLURM Job Script: Hi-C Alignment & Scaffolding
```bash
#!/bin/bash
#SBATCH --job-name=hic_scaffolding
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --output=/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/hic_scaffolding/hic_scaffolding_%j.out
#SBATCH --error=/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/hic_scaffolding/hic_scaffolding_%j.err

source ~/.bash_profile
conda activate scaffolding_env

WORKDIR="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/hic_scaffolding"
cd ${WORKDIR}

############################
# USER PARAMETERS
############################
ASM="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/hifiasm_assembly/hifiasm_assembly_l3/Helicops_angulatus_IAvH_10017.l3.bp.p_ctg.fasta"
HIC_R1="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/HiC/HiC_R1.fastq.gz"
HIC_R2="/home/dgarcia/mendel-nas1/PacBio/H.angulatus_IAvH_10017/HiC/HiC_R2.fastq.gz"
PREFIX="H_angulatus_scaffolded"

############################
# ALIGNMENT & SCAFFOLDING
############################
# 1. Align Hi-C reads to the primary assembly
chromap -i -r ${ASM} -o ${WORKDIR}/asm.index
chromap --preset hic -r ${ASM} -x ${WORKDIR}/asm.index -1 ${HIC_R1} -2 ${HIC_R2} --SAM -o ${WORKDIR}/aligned_hic.sam -t ${SLURM_CPUS_PER_TASK}

# 2. Run Scaffolding 
yahs ${ASM} ${WORKDIR}/aligned_hic.sam -o ${PREFIX}
```

### Summary of Scaffolding Results

The Hi-C scaffolding was performed using **YaHS**. The tables below demonstrate the structural refinement from the primary `hifiasm` (-l3) contigs to the final Hi-C scaffolds.

| Metric | Primary Contigs (Hifiasm -l3) | Final Scaffolds (YaHS) | Delta / Notes |
| :--- | :--- | :--- | :--- |
| **Total Length** | 2,167,390,149 bp | 2,167,399,349 bp | **+ 9,200 bp** (Spacer 'N's) |
| **Total Sequences** | 140 | 178 | **+ 38** (Structural breaks/edits) |
| **Largest Sequence** | 215,930,598 bp | 215,868,000 bp | Minor length reduction |
| **N50** | 131,633,754 bp | 131,826,854 bp | **+ 193.1 kb** |
| **L50** | 7 | 7 | Maintained macro-structure |
| **Gaps** | 0 | 92 | Representing 9,200 `N`s |








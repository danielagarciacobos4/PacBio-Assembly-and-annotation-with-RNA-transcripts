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








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





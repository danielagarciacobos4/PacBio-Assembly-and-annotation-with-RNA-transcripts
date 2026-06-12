<img width="1544" height="820" alt="Screenshot 2026-06-12 at 10 30 37 AM" src="https://github.com/user-attachments/assets/479def66-1378-454b-9ba5-210a115fd58c" />
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



<img width="1651" height="844" alt="Screenshot 2026-06-12 at 4 28 59 AM" src="https://github.com/user-attachments/assets/cc714e00-18cf-401e-928a-0b6c7a4bcc00" />

<img width="1633" height="818" alt="Screenshot 2026-06-12 at 10 37 05 AM" src="https://github.com/user-attachments/assets/c2c99c08-8554-450e-84f7-2f6d26499655" />





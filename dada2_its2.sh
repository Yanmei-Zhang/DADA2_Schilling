#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhan4170@umn.edu
#SBATCH -p ram1t
#SBATCH -o dada2_its2.out
#SBATCH -e dada2_its2.err

cd ~/dada2_tutorial

# Load the primer trimming tool that will remove the primers
module load cutadapt
# Check the cutadapt path
which cutadpat

# Load R
module load R/4.1.0

R --slave --vanilla<dada2_pipeline_its2.R

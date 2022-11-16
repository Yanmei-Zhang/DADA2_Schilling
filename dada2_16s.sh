#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhan4170@umn.edu
#SBATCH -p ram1t
#SBATCH --job-name="dada2_16s_job"
#SBATCH -o dada2_16s.out
#SBATCH -e dada2_16s.err

cd ~/dada2_tutorial

# Load R
module load R/4.1.0

R --slave --vanilla<dada2_pipeline_16s.R

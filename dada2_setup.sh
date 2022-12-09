#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhan4170@umn.edu
#SBATCH -p small
#SBATCH -o dada2_setup.out
#SBATCH -e dada2_setup.err

mkdir dada2_tutorial
cd dada2_tutorial

module load R/4.1.0

R --slave --vanilla<dada2_pipeline_setup.R

mkdir db_files
cd db_files

# Download the latest version of  General Fasta release files from the UNITE ITS database and used as a reference
wget https://files.plutof.ut.ee/public/orig/AB/5F/AB5F19B2369481FFC6025F7F95A58998FAE7648C9E20A8A33FE3226E724509B0.tgz # Change to the LATEST version
tar -xzvf AB5F19B2369481FFC6025F7F95A58998FAE7648C9E20A8A33FE3226E724509B0.tgz
rm AB5F19B2369481FFC6025F7F95A58998FAE7648C9E20A8A33FE3226E724509B0.tgz

# Download FUNGuild
wget https://github.com/UMNFuN/FUNGuild/archive/refs/heads/master.zip
unzip master.zip
rm master.zip

# Download the Silva reference database for 16S

wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

cd ..

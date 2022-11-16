# dada2_SchillingLab

*This is a tutorial for dada2 workdlow used in Schilling Lab, created by Yanmei Zhang*

The original pipeline on which this tutorial is based can be found [here](https://benjjneb.github.io/dada2/ITS_workflow.html). Please refer to this original online workflow to understand more about each step.

This pipeline runs the dada2 workflow for amplicon dataset from R on the HPC server at U of M.

## Starting Point - A Checklist

1. Check the primer set you used for sequencing and make sure the length of your target "Amplicon". This can vary between primer sets, as well as within primer sets. For example, ITS is higly variable in length even amplified using the same primer sers, which has a significant influence on the slection of sequencing platform and parameters for dada2 filtering and trimming steps.
   
   In this tutorial as one of the example, we will use the primers 5.8S and ITS4 to amplify the ITS2 region, which commnonly range from 100-500 bp.

   *5.8SR*: "TCGATGAAGAACGCAGCG" # Change it to YOUR primer

   *ITS4*: "TCCTCCGCTTATTGATATGC" # Change it to YOUR primer

2. Check the sequencing platform you used and make sure your read length. When you choose paired-end sequencing and try to merge the paired end reads, make sure the read length is long enough to allow for substantial overlap between the forward and reverse read. For example, with a MiSeq 2 x 300 bp run, if we assume a minimum of 20 bp overlap, you will get a maximum fragment length of 542 bp (300+300-17-20-20=542). This will work for ITS2.

   In this tutorial as an example, the sequence data were from Illumina MiSeq 2 x 300 bp run, which will generate bidirectional reads of 300 bp.

3. Check the sequencing dataset and make sure they meet certain criteria:

   - Samples have been demultiplexed, i.e., split into individual per-sample fastq files.
   - If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order.
  
## Set Up the Working Environment

**First, you need log into the wasabi ot mangi server through Schilling Lab:**

`ssh yourMSIusername@resourcename.msi.umn.edu`

Replace yourMSIusername with your MSI username (x500) and resourcename with the MSI resource that you wish to connect to (such as mesabi, mangi, etc). If you don't have access to MSI, you need your PI or gdesignated Group Administrator to add you the group.

 > **_NOTE:_** You need to be connected to the eduroam network or the UMN VPN for this command to be successful.\
 > Please refer to [MSI](https://www.msi.umn.edu/) for guidance to use their resources.
 
**Second, set up your working folder:**

Let's assume you have a working folder under this path:

~/dada2_tutorial/

```bash
mkdir dada2_tutorial
cd dada2_tutorial
```

You can move your demultiplexed data to this folder

**Third, download this tutorial from github:**

```bash
wget https://github.com/Yanmei-Zhang/dada2_SchillingLab/archive/refs/heads/main.zip  
unzip main.zip
rm main.zip

mv dada2_SchillingLab-main/* ./
```

**Fourth, download the dada2-formatted reference database of your choice:**

 > **_NOTE:_** DADA2 requires databases be in a custom format! Using the link here to [download](https://benjjneb.github.io/dada2/training.html):

```bash
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
```

**Last, set up the R enviroment:**

If this is your first time to run this pipeline, I recommend starting an [interactive job](https://www.msi.umn.edu/content/interactive-queue-use-srun) on any of MSI’s clusters that support job submission. After you have no problem with this pipeline, you can submit your batch job using [Slurm](https://www.msi.umn.edu/content/job-submission-and-scheduling-slurm) to schdule your job. This is truely effeciency if your have big data.

*Option 1:* In an interactive job run, first we load R:

```bash
module load R/4.1.0 # version 4.1.0 is required to use Bioconductor version '3.13' 
R
```

Install dada2 & other necessary packages. If this is your first time using R in HPC, when you install a package you might get a prompt asking if you want to create your own library. Answer 'yes' twice in the console to continue.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(lib = "~/R/x86_64-pc-linux-gnu-library/4.1")
BiocManager::install("dada2", lib = "~/R/x86_64-pc-linux-gnu-library/4.1")
BiocManager::install("ShortRead", lib = "~/R/x86_64-pc-linux-gnu-library/4.1")
```

> **_NOTE:_** This installation may take a long time, so only run this code if these packages are not already installed!

```r
# load packages and check the version
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
```

*Option 2:* In a batch job run, you submite the script using Slurm.

```bash
# In this srcipt, all the set up scripts is included
sbatch dada2_setup.sh
```

## Run dada2 workflow

### Run dada2 workflow for ITS2 amplicon data

In this tutorial I have included an [example data](https://github.com/Yanmei-Zhang/DADA2_SchillingLab/tree/main/example/ITS2) of ITS2 amplicon from Illumina Miseq 2x300. This data set contains six decayed wood samples from Cloquet. 

*Option 1:* In an interactive job run, [here](https://github.com/Yanmei-Zhang/DADA2_SchillingLab/wiki/The-dada2-ITS-Pipeline-for-Pair_ended-Reads) is the step-by-step illustration of the dada2 ITS2 pipeline in Rstudio. Pay attention to dada2 filter and trim parameters for ITS amplicon processing. 

*Option 2:* In a batch job run, submit the script using Slurm. If you have learned well of each step, optimized the parameters of each step, this is really efficient especially if you have large dataset and can save your a lot of time! 

```bash
sbatch dada2_its2.sh
```

### Run dada2 workflow for 16s amplicon data

In this tutorial I have included an [example data](https://github.com/Yanmei-Zhang/DADA2_SchillingLab/tree/main/example/16S) of 16s amplicon from Illumina Miseq 2x300. This data set contains six decayed wood samples from Cloquet. 

*Option 1:* In an interactive job run, [here](https://github.com/Yanmei-Zhang/DADA2_SchillingLab/wiki/The-dada2-Pipeline-for-16s-Amplicon-Pair_ended-Reads) is the step-by-step illustration of the dada2 pipeline in Rstudio. 

*Option 2:* In a batch job run, submit the script using Slurm. If you have learned well of each step, optimized the parameters of each step, this is really efficient especially if you have large dataset and can save your a lot of time! 

```bash
sbatch dada2_16s.sh
```

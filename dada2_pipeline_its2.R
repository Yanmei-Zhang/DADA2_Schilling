# DADA2 pipeline for ITS used by R/4.1.0

## 0. Get ready

# Load dada2 and other required packages
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(ggplot2)

# Define data path
data.path <- "~/dada2_tutorial/example/ITS2" # Change it to YOUR data folder
list.files(data.path)

# Generate matched lists of the forward and reverse read files, as well as parsing out the sample name. 
fnFs <- sort(list.files(data.path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(data.path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Set up work path in YOUR directory where you want data; you can also use the same path as the data is in, but it was nice to keep the data folder as original and store the processed data in another folder. 
work.path <- "~/dada2_tutorial/dada2_ITS2" # Change it to YOUR work path

## 1. Pre-processing data for dada2 -  remove sequences with Ns, remove primers with cutadapt

# Set up names of sub directories to stay organized
preprocess.fp <- file.path(work.path, "01_preprocess")
filtN.fp <- file.path(preprocess.fp, "filtN")
cut.fp <- file.path(preprocess.fp, "cutadapt")
filter.fp <- file.path(work.path, "02_filter") 
table.fp <- file.path(work.path, "03_tabletax") 

# Identify primers
FWD <- "TCGATGAAGAACGCAGCG" #5.8SR # Change it to YOUR sequencing primer
REV <- "TCCTCCGCTTATTGATATGC" #ITS4 # Change it to YOUR sequencing primer

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# “Pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))
trimN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
head(trimN)

# Count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. 
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Check the primers for the first sample before running the cutadapt. 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Remove Primers
cutadapt <- "/panfs/roc/msisoft/cutadapt/1.18/bin/cutadapt"# CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

if(!dir.exists(cut.fp)) dir.create(cut.fp)
fnFs.cut <- file.path(cut.fp, basename(fnFs))
fnRs.cut <- file.path(cut.fp, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Count the presence of primers in the first cutadapt-ed sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(cut.fp, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(cut.fp, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

sample.namesR <- unname(sapply(cutRs, get.sample.name))
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

## 2. Run dada2 pipeline 

### 2.1 Filter and Trim

# Inspect read quality profiles
# If the number of samples is 10 or less, plot them all, otherwise, just plot 10 randomly selected samples
if( length(cutFs) <= 10) {
  fwd_qual_plots <- plotQualityProfile(cutFs)
  rev_qual_plots <- plotQualityProfile(cutRs)
} else {
  rand_samples <- sample(size = 10, 1:length(cutFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(paste0(cutFs[rand_samples]))
  rev_qual_plots <- plotQualityProfile(paste0(cutRs[rand_samples]))
}

if(!dir.exists(filter.fp)) dir.create(filter.fp)
ggsave(paste0(filter.fp, "/its2_fwd_qual_10.png"), fwd_qual_plots)
ggsave(paste0(filter.fp, "/its2_rev_qual_10.png"), rev_qual_plots)

# Put filtered reads into separate sub-directories for big data workflow
subF.fp <- file.path(filter.fp, "filt_F") 
subR.fp <- file.path(filter.fp, "filt_R") 
dir.create(subF.fp)
dir.create(subR.fp)

# Filter and Trim
filtFs <- file.path(subF.fp, basename(cutFs))
filtRs <- file.path(subR.fp, basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 5), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

### 2.2 Infer Sequence Variants and Merge

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

errF_plot <- plotErrors(errF, nominalQ = TRUE) 
errR_plot <- plotErrors(errR, nominalQ=TRUE)

ggsave(paste0(filter.fp, "/its2_errF_plot.png"), errF_plot)
ggsave(paste0(filter.fp, "/its2_errR_plot.png"), errR_plot)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = TRUE)
# Inspecting the returned dada-class object:
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#Inspect the merger data.frame from the first sample
head(mergers[[2]])

# Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

### 2.3 Remove Chimeras and Summary of Reads

# Remove Chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
# Print percentage of our seqences that were not chimeric.
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Save file
if(!dir.exists(table.fp)) dir.create(table.fp)
write.table(t(seqtab.nochim), paste0(table.fp, "/dada2_its2_counts.txt"), sep="\t", quote=F, row.names=T)
write.table(track , paste0(table.fp, "/dada2_its2_track.txt"), sep="\t", quote=F, row.names = T)

### 2.4 Assign Taxonomy

# Assign taxonomy
unite.ref <- "~/dada2_tutorial/db_files/sh_general_release_dynamic_16.10.2022.fasta"  # Change it to location on YOUR device
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, paste0(table.fp, "/dada2_its2_taxa.txt"), sep="\t", quote=F, row.names = T)

seqtable.taxa<- cbind('#seq'=rownames(taxa), t(seqtab.nochim), taxa)
write.table(seqtable.taxa, paste0(table.fp, "/dada2_its2_counts.taxon.species.txt"), sep="\t", quote=F, row.names=F)

# Save image
save.image(file= paste0(work.path, "/dada2_its2.RData"))

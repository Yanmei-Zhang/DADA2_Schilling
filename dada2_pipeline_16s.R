# DADA2 pipeline for 16s used by R/4.1.0

## 0. Get ready

# Load dada2 and other required packages
library(dada2)
packageVersion("dada2")
library(ggplot2)

# Define data path
data.path <- "~/dada2_tutorial/example/16S" # Change it to YOUR data folder
list.files(data.path)

# Generate matched lists of the forward and reverse read files, as well as parsing out the sample name. 
fnFs <- sort(list.files(data.path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(data.path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

sample.namesR <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

# Set up work path in YOUR directory where you want data; you can also use the same path as the data is in, but it was nice to keep the data folder as original and store the processed data in another folder. 
work.path <- "~/dada2_tutorial/dada2_16S" # Change it to YOUR work path

# Set up names of sub directories to stay organized
filter.fp <- file.path(work.path, "01_filter") 
table.fp <- file.path(work.path, "02_tabletax") 

## 1. Run dada2 pipeline 

### 1.1 Inspect Read Quality Profiles

# Inspect read quality profiles
# If the number of samples is 10 or less, plot them all, otherwise, just plot 10 randomly selected samples
if( length(fnFs) <= 10) {
  fwd_qual_plots <- plotQualityProfile(fnFs)
  rev_qual_plots <- plotQualityProfile(fnRs)
} else {
  rand_samples <- sample(size = 10, 1:length(fnFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(paste0(fnFs[rand_samples]))
  rev_qual_plots <- plotQualityProfile(paste0(fnRs[rand_samples]))
}

if(!dir.exists(work.path)) dir.create(work.path)
if(!dir.exists(filter.fp)) dir.create(filter.fp)
ggsave(paste0(filter.fp, "/16s_fwd_qual_10.png"), fwd_qual_plots)
ggsave(paste0(filter.fp, "/16s_rev_qual_10.png"), rev_qual_plots)

# Put filtered reads into separate sub-directories for big data workflow
subF.fp <- file.path(filter.fp, "filt_F") 
subR.fp <- file.path(filter.fp, "filt_R") 
dir.create(subF.fp)
dir.create(subR.fp)

### 1.2 Filter and Trim

# Filter and Trim
filtFs <- file.path(subF.fp, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(subR.fp, paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), truncLen=c(240,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


### 1.3 Infer Sequence Variants and Merge

#### Learn the Error Rates

# Learn the Error Rates
# Set seed to ensure that randomized steps can be replicated
set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

errF_plot <- plotErrors(errF, nominalQ = TRUE) 
errR_plot <- plotErrors(errR, nominalQ=TRUE)

ggsave(paste0(filter.fp, "/16s_errF_plot.png"), errF_plot)
ggsave(paste0(filter.fp, "/16s_errR_plot.png"), errR_plot)

#### Dereplication, Sequence Inference, and Merging of paired-end Reads

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
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
dim(seqtab2)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

### 1.4 Remove Chimeras and Summary of Reads

# Remove Chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
# Print percentage of our seqences that were not chimeric.
sum(seqtab.nochim)/sum(seqtab2)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Save file
if(!dir.exists(table.fp)) dir.create(table.fp)
write.table(t(seqtab.nochim), paste0(table.fp, "/dada2_16s_counts.txt"), sep="\t", quote=F, row.names=T)
write.table(track , paste0(table.fp, "/dada2_16s_track.txt"), sep="\t", quote=F, row.names = T)

### 1.5 Assign Taxonomy

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/dada2_tutorial/db_files/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)

# Let’s inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Make species level assignments based on exact matching
taxa.plus <- addSpecies(taxa, "~/dada2_tutorial/db_files/silva_species_assignment_v138.1.fa.gz", verbose=TRUE)

# Let’s inspect the taxonomic assignments:
taxa.plus.print <- taxa.plus # Removing sequence rownames for display only
rownames(taxa.plus.print) <- NULL
head(taxa.plus.print)
write.table(taxa.plus, paste0(table.fp, "/dada2_16s_taxa.txt"), sep="\t", quote=F, row.names = T)

seqtable.taxa<- cbind('#seq'=rownames(taxa), t(seqtab.nochim), taxa)
write.table(seqtable.taxa, paste0(table.fp, "/dada2_16s_counts.taxon.txt"), sep="\t", quote=F, row.names=F)

seqtable.taxa.plus<- cbind('#seq'=rownames(taxa.plus), t(seqtab.nochim), taxa.plus)
write.table(seqtable.taxa.plus, paste0(table.fp, "/dada2_16s_counts.taxon.species.txt"), sep="\t", quote=F, row.names=F)

save.image(file= paste0(work.path, "/dada2_16s.RData"))

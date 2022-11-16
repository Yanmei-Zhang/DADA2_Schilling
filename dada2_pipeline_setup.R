# Install dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(lib = "~/R/x86_64-pc-linux-gnu-library/4.1")
BiocManager::install("dada2", lib = "~/R/x86_64-pc-linux-gnu-library/4.1")
BiocManager::install("ShortRead", lib = "~/R/x86_64-pc-linux-gnu-library/4.1") 

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

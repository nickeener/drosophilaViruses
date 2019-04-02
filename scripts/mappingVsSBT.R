# Creates graph of kmer match proportion vs fraction of reads that map to Twyford and colored by source

library(dplyr)
library(ggplot2)

# Get kmer match proportion data, and sort by run number
#setwd("~/projects/drosophilaViruses/mapping/SRP133370/")
setwd("~/projects/drosophilaViruses/mapping/ERP012119/")
file = Sys.glob("bloomTrees/corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
match <- read.csv(file = file, sep = '\t')$Decimal.Fraction
data <- data.frame(runs = runs, match = match)
data <- data[order(runs),]
match <- data$match

# Get mapped read data, read count data, and source data
file = Sys.glob("bowtie/individual/mapped_reads_per_run.txt")
mapped_counts <- as.double(read.csv(file = file, sep = '\t')$Mapped.Reads)
file = Sys.glob("readcounts.txt")
read_counts <- as.double(read.csv(file = file, sep = '\t')$Read.Count)
file = Sys.glob("sources.txt")
sources <- as.character((read.csv(file = file, sep = '\t'))$Source)

# Create new vector with the proportion of mapped reads out of total reads
percent_mapped <- (mapped_counts/read_counts)*100

# Create new dataframe from match proportion, source, and percent of reads mapped data
#data <- data.frame(match = match, percent_mapped = percent_mapped, Sources = sources)
data <- data.frame(match = match, percent_mapped = percent_mapped)

# Plot match proportion vs percent of reads mapped and color by source
#data %>% ggplot(aes(x=match, y=percent_mapped, color=sources)) +
data %>% ggplot(aes(x=match, y=percent_mapped)) +
  geom_point() +
  geom_vline(xintercept = .075, size = 0.5) +
  theme_bw() +
  xlab("Kmer Match Proportion") +
  ylab("Percent of Reads Mapped") +
  ggtitle("Relation of Percentage of Twyford Mapped\n Reads to SBT Match Proportion") +
  labs(color="Sources") +
  theme(plot.title = element_text(hjust = 0.5))

# Creates graph of kmer match proportion vs fraction of reads that map to Twyford and colored by source

library(dplyr)
library(ggplot2)

# Get kmer match proportion data, and sort by run number
setwd("~/projects/drosophilaViruses/mapping/SRP133370/")
file = Sys.glob("bloomTrees/corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
run_matches <- read.csv(file = file, sep = '\t')$Decimal.Fraction
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- c()
match <- append(match, run_matches)

# Get mapped read data, read count data, and source data
file = Sys.glob("bowtie/individual/mapped_reads_per_run.txt")
mapped_counts <- as.double(read.csv(file = file, sep = '\t')$Mapped.Reads)
file = Sys.glob("readcounts.txt")
read_counts <- as.double(read.csv(file = file, sep = '\t')$Read.Count)
file = Sys.glob("sources.txt")
sources <- as.character((read.csv(file = file, sep = '\t'))$Source)

# Create new vector with the proportion of mapped reads out of total reads
percent_mapped <- (mapped_counts/read_counts)*100

# Do the same for every other study
# ERP012119
setwd("~/projects/drosophilaViruses/mapping/ERP012119/")
file = Sys.glob("bloomTrees/corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
run_matches <- read.csv(file = file, sep = '\t')$Decimal.Fraction
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

file = Sys.glob("bowtie/individual/mapped_reads_per_run.txt")
mapped_counts <- as.double(read.csv(file = file, sep = '\t')$Mapped.Reads)
file = Sys.glob("readcounts.txt")
read_counts <- as.double(read.csv(file = file, sep = '\t')$Read.Count)
for (run in runs) {
  sources <- append(sources, "ERP012119")
}

percent_mapped <- append(percent_mapped, (mapped_counts/read_counts)*100)

#SRP053390
setwd("/media/nickeener/External_Drive/SRP053390")
file = Sys.glob("bloomTrees/corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
run_matches <- read.csv(file = file, sep = '\t')$Decimal.Fraction
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

file = Sys.glob("bowtie/individual/mapped_reads_per_run.txt")
mapped_counts <- as.double(read.csv(file = file, sep = '\t')$Mapped.Reads)
file = Sys.glob("readcounts.txt")
read_counts <- as.double(read.csv(file = file, sep = '\t')$Read.Count)
for (run in runs) {
  sources <- append(sources, "SRP053390")
}

percent_mapped <- append(percent_mapped, (mapped_counts/read_counts)*100)

# SRP070549
setwd("/media/nickeener/External_Drive/SRP070549")
file = Sys.glob("bloomTrees/corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
run_matches <- read.csv(file = file, sep = '\t')$Decimal.Fraction
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

file = Sys.glob("bowtie/individual/mapped_reads_per_run.txt")
mapped_counts <- as.double(read.csv(file = file, sep = '\t')$Mapped.Reads)
file = Sys.glob("readcounts.txt")
read_counts <- as.double(read.csv(file = file, sep = '\t')$Read.Count)
for (run in runs) {
  sources <- append(sources, "SRP070549")
}

percent_mapped <- append(percent_mapped, (mapped_counts/read_counts)*100)

# SRP106416
setwd("/media/nickeener/External_Drive/SRP106416")
file = Sys.glob("bloomTrees/corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
run_matches <- read.csv(file = file, sep = '\t')$Decimal.Fraction
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

file = Sys.glob("bowtie/individual/mapped_reads_per_run.txt")
mapped_counts <- as.double(read.csv(file = file, sep = '\t')$Mapped.Reads)
file = Sys.glob("readcounts.txt")
read_counts <- as.double(read.csv(file = file, sep = '\t')$Read.Count)
for (run in runs) {
  sources <- append(sources, "SRP106416")
}

percent_mapped <- append(percent_mapped, (mapped_counts/read_counts)*100)

# SRP117133
setwd("/media/nickeener/External_Drive/SRP117133")
file = Sys.glob("bloomTrees/corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
run_matches <- read.csv(file = file, sep = '\t')$Decimal.Fraction
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

file = Sys.glob("bowtie/individual/mapped_reads_per_run.txt")
mapped_counts <- as.double(read.csv(file = file, sep = '\t')$Mapped.Reads)
file = Sys.glob("readcounts.txt")
read_counts <- as.double(read.csv(file = file, sep = '\t')$Read.Count)
for (run in runs) {
  sources <- append(sources, "SRP117133")
}

percent_mapped <- append(percent_mapped, (mapped_counts/read_counts)*100)

# Create new dataframe from match proportion, source, and percent of reads mapped data
data <- data.frame(match = match, percent_mapped = percent_mapped, sources = sources)

# Plot match proportion vs percent of reads mapped and color by source
data %>% ggplot(aes(x=match, y=percent_mapped, color=sources)) +
  geom_point() +
  #geom_vline(xintercept = .075, size = 0.5) +
  theme_bw() +
  xlab("Kmer Match Proportion") +
  ylab("Percent of Reads Mapped") +
  ggtitle("Relation of Percentage of Twyford Mapped\n Reads to SBT Match Proportion") +
  labs(color="Sources") +
  theme(plot.title = element_text(hjust = 0.5))

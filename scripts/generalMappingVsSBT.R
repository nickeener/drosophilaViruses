# Creates graph of kmer match proportion vs fraction of reads that map to DCV or Kallithea Virus and colored by source

library(dplyr)
library(ggplot2)

# Get kmer match proportion data, and sort by run number for SRP119720
setwd("~/projects/drosophilaViruses/mapping/SRP119720/")
file = Sys.glob("bloomTrees/queries_dcv.dat")
runs <- as.character(read.csv(file = file, sep = ' ')$X.NC_001834.1)
run_matches <- read.csv(file = file, sep = ' ')$C
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- c()
match <- append(match, run_matches)

# Get read count data for SRP119720
file = Sys.glob("bowtie/individual/readcounts.txt")
mapped_counts <- c()
total_counts <- c()
mapped_counts <- read.csv(file = file, sep = '\t')$Mapped.Read.Count
total_counts <- read.csv(file = file, sep = '\t')$Total.Read.Count
percent_mapped <- (mapped_counts/total_counts)*100

# Count number of runs in study and add to source list
source <- c()
for (i in runs) {
  source <- append(source, 'SRP119720')
}

# Get kmer match proportion data, and sort by run number for ERP111111 (dcv control)
setwd("~/projects/drosophilaViruses/mapping/ERP111111/")
file = Sys.glob("bloomTrees/queries_dcv.dat")
runs <- as.character(read.csv(file = file, sep = ' ')$X.NC_001834.1)
run_matches <- read.csv(file = file, sep = ' ')$C
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

# Get read count data for ERP111111 dcv
file = Sys.glob("bowtie/individual_dcv/readcounts.txt")
mapped_counts <- read.csv(file = file, sep = '\t')$Mapped.Read.Count
total_counts <- read.csv(file = file, sep = '\t')$Total.Read.Count
percent_mapped <- append(percent_mapped, (mapped_counts/total_counts)*100)

# Count number of runs in study and add to source list
for (i in runs) {
  source <- append(source, 'ERP023609 Control DCV')
}

# Get kmer match proportion data, and sort by run number for ERP111111 (kallithea control)
setwd("~/projects/drosophilaViruses/mapping/ERP111111/")
file = Sys.glob("bloomTrees/queries_kallithea.dat")
runs <- as.character(read.csv(file = file, sep = ' ')$X.Kallithea)
run_matches <- read.csv(file = file, sep = ' ')$X22
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

# Get read count data for ERP111111 
file = Sys.glob("bowtie/individual_kallithea/readcounts.txt")
mapped_counts <- read.csv(file = file, sep = '\t')$Mapped.Read.Count
total_counts <- read.csv(file = file, sep = '\t')$Total.Read.Count
percent_mapped <- append(percent_mapped, (mapped_counts/total_counts)*100)


# Count number of runs in study and add to source list
for (i in runs) {
  source <- append(source, 'ERP023609 Control Kallithea')
}

# Get kmer match proportion data, and sort by run number for ERP222222 (DCV)
setwd("~/projects/drosophilaViruses/mapping/ERP222222/")
file = Sys.glob("bloomTrees/queries_dcv.dat")
runs <- as.character(read.csv(file = file, sep = ' ')$X.NC_001834.1)
run_matches <- read.csv(file = file, sep = ' ')$C
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

# Get read count data for ER222222
file = Sys.glob("bowtie/individual/readcounts.txt")
mapped_counts <- read.csv(file = file, sep = '\t')$Mapped.Read.Count
total_counts <- read.csv(file = file, sep = '\t')$Total.Read.Count
percent_mapped <- append(percent_mapped, (mapped_counts/total_counts)*100)


# Count number of runs in study and add to source list
for (i in runs) {
  source <- append(source, 'ERP023609 DCV')
}

# Get kmer match proportion data, and sort by run number for ERP333333 (kallithea)
setwd("~/projects/drosophilaViruses/mapping/ERP333333/")
file = Sys.glob("bloomTrees/queries_kallithea.dat")
runs <- as.character(read.csv(file = file, sep = ' ')$X.Kallithea)
run_matches <- read.csv(file = file, sep = ' ')$X24
data <- data.frame(runs = runs, run_matches = run_matches)
data <- data[order(runs),]
run_matches <- data$run_matches
match <- append(match, run_matches)

# Get read count data for ERP333333
file = Sys.glob("bowtie/individual/readcounts.txt")
mapped_counts <- read.csv(file = file, sep = '\t')$Mapped.Read.Count
total_counts <- read.csv(file = file, sep = '\t')$Total.Read.Count
percent_mapped <- append(percent_mapped, (mapped_counts/total_counts)*100)


# Count number of runs in study and add to source list
for (i in runs) {
  source <- append(source, 'ERP023609 Kallithea')
}

# Create new dataframe from match proportion, source, and percent of reads mapped data
data <- data.frame(match = match, percent_mapped = percent_mapped, sources = source)

# Plot match proportion vs percent of reads mapped and color by source
data %>% ggplot(aes(x=match, y=percent_mapped, color=source)) +
  geom_point() +
  #geom_vline(xintercept = .075, size = 0.5) +
  theme_bw() +
  xlab("Kmer Match Proportion") +
  ylab("Percent of Reads Mapped") +
  ggtitle("Relation of Percentage of Twyford Mapped\n Reads to SBT Match Proportion") +
  labs(color="Sources") +
  theme(plot.title = element_text(hjust = 0.5))
# Libraries
library(dplyr)
library(ggplot2)

# Get run accession numbers and read count data from file
setwd("~/projects/drosophilaViruses/mapping")
data <- read.csv(file = "SRP133370/readcounts.txt", sep = '\t')
run_acc <- data$Run.Accession.Number
counts <- data$Read.Count
read_counts <- c()
subsets <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 1)
for (set in subsets) {
  read_counts <- c(read_counts, counts*set)
}

# Get source data
file <- Sys.glob("SRP133370/bloomTrees/test_corrected_query_twyford_full.dat")
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
sources <- as.character((read.csv(file = file, sep = '\t'))$Source)
data <- data.frame(runs, sources)
newdata <- data[order(runs),]
sources <- newdata$sources
all_sources <- c()
for (set in subsets) {
  all_sources <- c(all_sources, sources)
}

# Get kmer match proportion data from uncorrected data
files <- Sys.glob("SRP0000*/bloomTrees/queries_twyford_full.dat")
files <- files[2:11] # Remove first entry (not the correct file)
all_proportions <- c()
for (i in seq(length(run_acc))) {
  all_proportions <- append(all_proportions, 0)
}
for (file in files) {
  if (file == "SRP000005/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000010/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000015/bloomTrees/queries_twyford_full.dat") {
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000020/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000025/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000030/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000040/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000050/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
   if (file == "SRP000070/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
   }
  if (file == "SRP000090/bloomTrees/queries_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
    kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    all_proportions <- c(all_proportions, newdata$kmer_proportion)
  }
}
file = "SRP133370/bloomTrees/queries_twyford_full.dat"
runs <- as.character(read.csv(file = file, sep = ' ', skip = 1, header = F)$V1)
kmer_proportion <- read.csv(file = file, sep= ' ', skip = 1, header = F)$V3
data <- data.frame(runs, kmer_proportion)
newdata <- data[order(runs),]
all_proportions <- c(all_proportions, newdata$kmer_proportion)
all_proportions[698] <- 0.044232


# Get kmer match proportion data from corrected data
files <- Sys.glob("SRP0000*/bloomTrees/corrected_query_twyford_full.dat")
corr_all_proportions <- c()
for (i in seq(length(run_acc))) {
  corr_all_proportions <- append(corr_all_proportions, 0)
}
for (file in files) {
  if (file == "SRP000005/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
  }
  if (file == "SRP000010/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000015/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000020/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000025/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000030/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000040/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000050/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000070/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
  if (file == "SRP000090/bloomTrees/corrected_query_twyford_full.dat") {
    runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
    kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
    data <- data.frame(runs, kmer_proportion)
    newdata <- data[order(runs),]
    corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
    }
}
file = "SRP133370/bloomTrees/corrected_query_twyford_full.dat"
runs <- as.character(read.csv(file = file, sep = '\t')$Run.Accession.Number)
kmer_proportion <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
data <- data.frame(runs, kmer_proportion)
newdata <- data[order(runs),]
corr_all_proportions <- c(corr_all_proportions, newdata$kmer_proportion)
corr_all_proportions[698] <- 0.0158784166950

# Create dataframes and convert run_acc into a factor vector
data <- data.frame(run_acc, read_counts, all_proportions)
corr_data <- data.frame(run_acc, read_counts, corr_all_proportions)

# Create plot for uncorrected data
data %>% ggplot(aes(x=read_counts, y=all_proportions, color=run_acc)) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("Number of Reads") +
  ylab("Kmer Match Proportion") +
  ggtitle("Background Kmers Included") +
  geom_line() 

# Create plot for corrected data
corr_data %>% ggplot(aes(x=read_counts, y=corr_all_proportions, color=run_acc)) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("Number of Reads") +
  ylab("Kmer Match Proportion") +
  ggtitle("Background Kmers Excluded") +
  geom_line()
  
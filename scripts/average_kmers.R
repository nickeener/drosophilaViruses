# Calculates average kmer match proportion for runs containing Twyford and for runs that don't
setwd("~/projects/drosophilaViruses/mapping")

# Get kmer match proportion data from Exposed Whole Fly SRP133370 runs
file <- Sys.glob("SRP133370/bloomTrees/test_corrected_query_twyford_full.dat")
kmer_proportions <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
sources <- as.character((read.csv(file = file, sep = '\t'))$Source)
data <- data.frame(sources, kmer_proportions)
data <- subset(data, sources == "Exposed Whole Fly")
kmer_proportions <- data$kmer_proportions

# Get kmer match proportion data from ERP012119
file <- Sys.glob("ERP012119/bloomTrees/corrected_query_twyford_full.dat")
kmer_proportions <- c(kmer_proportions, as.double(as.character((read.csv(file = file, sep = '\t'))$Decimal.Fraction)))

# Get kmer match proportion dataf from 
file = "/media/nickeener/External_Drive/SRP106416/bloomTrees/corrected_query_twyford_full.dat"
kmer_proportions <- append(kmer_proportions, as.double(as.character((read.csv(file = file, sep = '\t'))$Decimal.Fraction[1])))
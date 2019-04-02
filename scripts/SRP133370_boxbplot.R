library(dplyr)
library(ggplot2)

# Get data from SRP133370
file <- Sys.glob("SRP133370/bloomTrees/test_corrected_query_twyford_full.dat")
kmer_proportions <- as.double(as.character(read.csv(file = file, sep = '\t')$Decimal.Fraction))
sources <- as.character((read.csv(file = file, sep = '\t'))$Source)
data <- data.frame(kmer_proportions, sources)
data1 <- subset(data, sources == "Exposed Whole Fly" & kmer_proportions > 0.04)
sources1 <- as.character(data1$sources)
kmer_proportions1 <- data1$kmer_proportions
data2 <- subset(data, sources == "Unexposed Whole Fly")
sources2 <- as.character(data2$sources)
kmer_proportions2 <- data2$kmer_proportions
data3 <- subset(data, sources == "Infected Cadaver")
sources3 <- as.character(data3$sources)
kmer_proportions3 <- data3$kmer_proportions
data4 <- subset(data, sources == "Exposed Brain")
sources4 <- as.character(data4$sources)
kmer_proportions4 <- data4$kmer_proportions
data5 <- subset(data, sources == "Unexposed Brain")
sources5 <- as.character(data5$sources)
kmer_proportions5 <- data5$kmer_proportions
counts <- c(counts1,counts2,counts3,counts4,counts5)
kmer_proportions <- c(kmer_proportions1,kmer_proportions2,kmer_proportions3,kmer_proportions4,kmer_proportions5)
sources <- c(sources1,sources2,sources3,sources4,sources5)

# Get data from ERP012119
file <- Sys.glob("ERP012119/bloomTrees/corrected_query_twyford_full.dat")
kmer_proportions <- c(kmer_proportions, as.double(as.character((read.csv(file = file, sep = '\t'))$Decimal.Fraction)))
for

data <- data.frame(kmer_proportions, sources)

data %>% ggplot(aes(y=kmer_proportions, x=sources, color=sources)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Library Source") +
  ylab("Kmer Match Proportion") +
  ggtitle("SRP133370 SBT Matches to Twyford Virus")





# Get Bloom input sequences from studies that have match

#Libraries
library(SRAdb)
library(seqinr)

# Check if SRA metadata file is in working directory, if it exists:
# set it to sqlfile, if not: download it and then set it to sqlfile
setwd("~/scripts/drosphilaviruses")
if(!file.exists('SRAmetadb.sqlite')) {
  sqlfile <- getSRAdbFile()
} else {
  sqlfile <- 'SRAmetadb.sqlite'
}

# Set up connection to SRA database
sra_con <- dbConnect(SQLite(),sqlfile)

# Read each file in the specified directory and store accession numbers in a vector
files <- Sys.glob("bloomOutput/twyford/*.tabular")
run_acc = c()
for (file in files) {
  run_acc <- append(run_acc, as.character(read.table(file)$V1))
}

# Create vector with run all run accession numbers associated with study of interest
study <- "SRP067567"
allruns <- sraConvert(study, sra_con = sra_con)$run

# Create vector with only run accession numbers from the study of interest that had matches
studyruns <- intersect(allruns, run_acc)

# Create vector of Bloom input ids that had matches with a particular study
# First load all ids from each BloomOutput file and check if they are in 
# studyruns then remove duplicates
fileids <- c()
for (file in files) {
  runs <- as.character(read.table(file)$V1)
  for (run in runs) {
    if (run %in% studyruns) {
      fileids <- append(fileids, file)
    }
  }
}
fileids <- unique(fileids)
# Split filename into character vector and extract the actual id and put
# each in the ids vector
ids <- c()
for (id in fileids) {
  splitid <- strsplit(id, "")[[1]]
  for (i in seq(length(splitid))) {
    if (splitid[i] == "/") {
      break
    }
  }
  for (j in seq(i, length(splitid))) {
    if (splitid[j] == ".") {
      break
    }
  }
  i <- i + 1
  j <- j - 1
  ids <- append(ids, paste(splitid[i:j], collapse=""))
}

# Extract sequences with matching ids from BloomInput files
file <- Sys.glob("bloomInput/twyford.txt")
fileids <- as.character(read.table(file)$V1)
fileseqs <- as.character(read.table(file)$V2)
matchseqs <- fileseqs[match(ids,fileids)]

# Write multi-fasta file of matching IDs and sequences
write.fasta(as.list(matchseqs), ids, paste(study, ".fasta", sep=""), open="w", nbchar=70, as.string=FALSE)

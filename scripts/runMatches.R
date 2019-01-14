# Given a study accession number, creates a file with all run accession numbers
# from that study that had BloomTree matches

#Libraries
library(SRAdb)

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
all_runs <- c()
for (file in files) {
  run_acc <- append(run_acc, as.character(read.table(file)$V1))
}

# Create vector with run all run accession numbers associated with study of interest
study <- "SRP083991"
allruns <- sraConvert(study, sra_con = sra_con)$run

# Create vector with only run accession numbers from the study of interest that had matches
studyruns <- intersect(allruns, run_acc)

# Write to file
write(studyruns, file = "twyford.SRP083991.txt", sep ="\t")

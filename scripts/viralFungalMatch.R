#Find studies where there are BloomTree matches from both the viral genome and the fungal sequences

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

# Read each file in the specified directories and add them to a different vectors
files1 <- Sys.glob("bloomOutput/Emuscae/*/*.tabular")
run_acc1 = c()
for (file in files1) {
  run_acc1 <- append(run_acc1, as.character(read.table(file)$V1))
}
files2 <- Sys.glob("bloomOutput/twyford/*.tabular")
run_acc2 = c()
for (file in files2) {
  run_acc2 <- append(run_acc2, as.character(read.table(file)$V1))
}

# Remove duplicate accession numbers and calculate the intersection
run_acc1 <- unique(run_acc1)
run_acc2 <- unique(run_acc2)
shared_acc <- intersect(run_acc1, run_acc2)

# Find the study accession numbers of the shared run accession numbers
shared_study <- sraConvert(shared_acc, sra_con = sra_con)$study
shared_study <- unique(shared_study)

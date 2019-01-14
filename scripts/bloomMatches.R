# Given output from BloomTreeFilter, creates a spreadsheet with studies that contained 
# matches that includes the study title, abstract, library construction protocol,
# number of runs that had a BloomTree match, total runs, and fraction of runs with a match

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
for (file in files) {
  run_acc <- append(run_acc, as.character(read.table(file)$V1))
}

# Convert run accession numbers to study accession numbers
study_acc <- sraConvert(run_acc, sra_con = sra_con)$study

# Count number of duplicates of each study_acc and store in vector
study_acc <- sort(study_acc)
counts <- c()
duplicates <- duplicated(study_acc)
i <- 1
while (i < length(duplicates)) {
  if (duplicates[i] == TRUE) {
    count <- 0
    while (i+count != length(duplicates)+1 && duplicates[i+count] == TRUE){
      count <- count+1
    }
    counts <- append(counts, count+1)
    i <- i+count
  }
  else if (duplicates[i] == FALSE && duplicates[i+1] == TRUE) {
    i <- i+1
  }
  else if (duplicates[i] == FALSE && i == length(duplicates)) {
    counts <- append(counts, 1)
  }
  else {
    counts <- append(counts, 1)
    i <- i+1
  }
}

# Remove duplicates
study_acc <- unique(study_acc)

# Create vector with all run accession numbers in each study that had a 
# BloomTree match
all_match_runs <- c()
for (study in study_acc) {
  match_runs <- c()
  runs <- sraConvert(study, sra_con = sra_con)$run
  match_runs <- intersect(runs, run_acc)
  all_match_runs <- append(all_match_runs, paste(match_runs, collapse = ', '))
}

# Count number of runs for each study_acc and store in vector
runcounts <- c()
for (study in study_acc) {
  runcounts <- append(runcounts, length(sraConvert(study, sra_con = sra_con)$run))
}

# Calculate percentage of runs from each study_acc that had a match and store in vector
percent <- formatC(counts/runcounts, digits = 3, format = "f")

# Query SRA metadata database for each study accession number and retrive and store titles, abstracts, and library construction protocols
titles <- c()
abstracts <- c()
protocols <- c()
for (study in study_acc) {
  titles <- append(titles, dbGetQuery(sra_con, paste0("SELECT study_title FROM study WHERE study_accession = '", study,"'")))
  abstracts <- append(abstracts, dbGetQuery(sra_con, paste0("SELECT study_abstract FROM study WHERE study_accession = '", study,"'")))
  protocols <- append(protocols, dbGetQuery(sra_con, paste0("SELECT library_construction_protocol FROM experiment WHERE study_accession = '", study, "' GROUP BY library_construction_protocol")))
}

# Convert to character vectors
titles <- as.character(titles)
abstracts <- as.character(abstracts)
protocols <- as.character(protocols)

# Create dataframe from data and write to file
data <- data.frame("Study Accession Numbers"=study_acc, "Title"=titles, "Abstract"=abstracts, "Library Construction Protocol"=protocols, "Matched Runs"=counts, "Total Runs"=runcounts, "Fraction of Runs"=percent, "Matched Run Accession Numbers"=all_match_runs)
write.csv(data, "resuts/output.csv")


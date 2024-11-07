setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (dir.exists("output_file")) {
  print("The output file exists. ")
  quit(save = "no")}
source("functions.R")
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse short flags into a named list
parse_args <- function(args) {
  args_list <- list(
    path = NULL,
    rRNA = NULL,
    plot = FALSE  # Default value for -p is FALSE
  )
  
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    
    if (key == "-d") {
      args_list$path <- args[i + 1]  # Directory path as a string
      i <- i + 2  # Move to the next flag
    } else if (key == "-r") {
      args_list$rRNA <- args[i + 1]  # rRNA type (16S or ITS)
      i <- i + 2  # Move to the next flag
    } else if (key == "-p") {
      args_list$plot <- TRUE  # Set plot to TRUE if -p is present
      i <- i + 1  # Move to the next flag
    } else {
      i <- i + 1  # Move to the next flag if key is unrecognized
    }
  }
  
  return(args_list)
}

# Parse the input arguments
params <- parse_args(args)

# Access the parameters
path <- params$path
rRNA <- params$rRNA
plot <- params$plot


if (rRNA == "16S"){
  load("SILVA_SSU_r138_2019.RData")
} else if (rRNA == "ITS"){
  load("UNITE_v2023_July2023.RData")
} else{
  stop("Invalid rRNA argument. Please use '16S' or 'ITS'.")
}

# Load the data
in_fastq <- input_fastq(path)
fnFs <- in_fastq$fnFs
fnRs <- in_fastq$fnRs
sample.names <- in_fastq$sample.names

# Filter the data
out <- filter(fnFs, fnRs, path, sample.names)
filtFs <- out$filtFs
filtRs <- out$filtRs

# Learn the errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Run dada2
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Assign taxonomy
taxa <- taxa(seqtab.nochim)

# Create phyloseq object
ps_ <- ps(seqtab.nochim, taxa)

# Convert phyloseq object to biom format
otu2biom(ps_)

# Run PICRUSt2
# convert the phyloseq sequences to fasta
seq2fna(ps_)
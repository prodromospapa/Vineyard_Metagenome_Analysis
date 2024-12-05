setwd(dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])))
source("functions.R")
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse short flags into a named list
parse_args <- function(args) {
  args_list <- list(
    path = NULL,
    rRNA = NULL,
    method = "dada2"
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
    } else if(key == "-m"){
      args_list$method <- args[i + 1]  # method (dada2 or decipher)
      i <- i + 2  # Move to the next flag
    }
    else {
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
method <- params$method


# Load the data
in_fastq <- input_fastq(path)
fnFs <- in_fastq$fnFs
fnRs <- in_fastq$fnRs
sample.names <- in_fastq$sample.names

# Filter the data
out <- filter_names(rRNA, sample.names)
filtFs <- out$filtFs
filtRs <- out$filtRs
filter_data(fnFs, filtFs, fnRs, filtRs, path ,sample.names)

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
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

# Assign taxonomy
if (rRNA == "16S"){
  if (method == "dada2"){
    ref = "silva_nr99_v138.1_wSpecies_train_set.fa" # https://zenodo.org/records/4587955
  } else if (method == "decipher"){
    load("SILVA_SSU_r138_2019.RData") # https://figshare.com/ndownloader/files/46245217
  } else {
    stop("Invalid method for 16S. Please use 'dada2' or 'decipher'.")
  }
} else if (rRNA == "ITS"){
  if (method == "dada2"){
    ref = "Full_UNITE+INSD.fa" # https://doi.plutof.ut.ee/doi/10.15156/BIO/2959330
  } else if (method == "decipher"){
    load("UNITE_v2023_July2023.RData") # https://figshare.com/articles/dataset/UNITE_v2023_July2023/27018373?file=49181545
  } else {
    stop("Invalid method for ITS. Please use 'dada2' or 'decipher'.")
  }
} else {
  stop("Invalid rRNA argument. Please use '16S' or 'ITS'.")
}

if (method == "dada2"){
  taxa <- assignTaxonomy(seqtab.nochim, ref, multithread=TRUE, tryRC = TRUE)
} else if (method == "decipher"){
  taxa <- taxa(seqtab.nochim)
}

# Create phyloseq object
ps_ <- ps(seqtab.nochim, taxa)

# remove 
if (rRNA=="16S"){
  ps_ = subset_taxa(ps_, 
                   Family  != "Mitochondria" &
                     Class   != "Chloroplast")
} else if (rRNA=="ITS"){
  ps_ = subset_taxa(ps_,
                    Class != "Ascomycota_cls_Incertae_sedis")
}

# export percentages
write.csv(transform_sample_counts(otu_table(ps_), function(x) round(x / sum(x) * 100,2)),file=paste0("output/otu_",rRNA,".csv"))
write.csv(tax_table(ps_),file=paste0("output/taxa_",rRNA,".csv"))

# create abundance plots
system(paste0("mkdir -p output/plots_",rRNA))
ranks <- rank_names(ps_)
for (group in ranks){
  try({ggsave(filename = paste0("output/plots_",rRNA,"/plot_",group,".pdf"),
    plot = plot_tax(ps_,group,TRUE),width=13)
  },silent=TRUE)
  }

# create abundance dendrogram
pdf(paste0("output/plots_",rRNA,"/microbiome_cluster.pdf"))
print(create_dendrogram(ps_,rRNA))
dev.off()

# Convert phyloseq object to biom format
otu2biom(ps_,rRNA)

# PICRUSt2
# convert the phyloseq sequences to fasta
seq2fna(ps_,rRNA)

# run picrust 2
system(paste0("bash functional.sh ", rRNA))

# plot the results
df <- tsv2df(rRNA)

# heatmap
pdf(paste0("output/plots_",rRNA,"/functional_heatmap.pdf"))
print(plot_df(df))
dev.off()

# correlation
pdf(paste0("output/plots_",rRNA,"/correlation.pdf"))
plot(corr_coef(df))
dev.off()

# clustering
dist_matrix <- dist(t(df), method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "ward.D2")
pdf(paste0("output/plots_",rRNA,"/functional_cluster.pdf"))
plot(hclust_result, main = "Hierarchical Clustering of Wine Varieties", xlab = "", sub = "")
dev.off()
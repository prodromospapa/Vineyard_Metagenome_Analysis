# List of required packages
required_packages <- c(
  "dada2",
  "phyloseq",
  "Biostrings",
  "DECIPHER",
  "ape",
  "tools",
  "ggplot2",
  "tidyr",
  "dplyr",
  "vegan",
  "ggpubr",
  "reshape2",
  "biomformat",
  "phangorn"
)

# Install missing packages
installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)

if (length(missing_packages) > 0) {
  BiocManager::install(missing_packages, ask = FALSE)
} else {
  message("All required packages are already installed!")
}

library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(ape)
library(tools)
library(ggplot2)  
library(tidyr)
library(dplyr)
library(vegan)
library(ggpubr)
library(reshape2)
library(biomformat)
library(phangorn)



# input_fastq: get the forward and reverse fastq files and sample names
input_fastq <- function(path){
    fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE, recursive = TRUE))
    fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE, recursive = TRUE))
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    results <-list(
        fnFs=fnFs,
        fnRs=fnRs,
        sample.names=sample.names)
    return(results)
}

# filter_names: get the filtered forward and reverse fastq files
filter_names <- function(path, sample.names){
  filtFs <- file.path(path, "output/filtered" ,paste0(sample.names, ".filt.1.fastq"))
  filtRs <- file.path(path, "output/filtered", paste0(sample.names, ".filt.2.fastq"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  results <- list(
    filtFs=filtFs,
    filtRs=filtRs)
  return(results)
}

# filter: filter the forward and reverse fastq files
filter_ <- function(fnFs,fnRs, path ,sample.names){
  filter_names <- filter_names(path, sample.names)
  filtFs <- filter_names$filtFs
  filtRs <- filter_names$filtRs
  output<-filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
        maxN=0, maxEE=c(2,2),truncQ = 2,
        compress=FALSE, multithread=TRUE)
  return(output)
}

# assign_taxonomy: assign taxonomy to the sequences
taxa <- function(seqtab.nochim){
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  dna <- DNAStringSet(getSequences(seqtab.nochim))
  ids <- IdTaxa(dna, trainingSet, processors=NULL, verbose=FALSE, strand="both")
  taxa_ <- t(sapply(ids, function(x) {
          m <- match(ranks, x$rank)
          taxa1 <- x$taxon[m]
          taxa1[startsWith(taxa1, "unclassified_")] <- NA
          taxa1
  }))
  colnames(taxa_) <- ranks; rownames(taxa_) <- getSequences(seqtab.nochim)
  colnames(taxa_) <- tools::toTitleCase(colnames(taxa_))
  return(taxa_)
}

# phyloseq: create a phyloseq object
ps <- function(seqtab.nochim,taxa_){
  ps_ <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                tax_table(taxa_))
  dna <- Biostrings::DNAStringSet(taxa_names(ps_))
  names(dna) <- taxa_names(ps_)
  ps_ <- merge_phyloseq(ps_, dna)

  # add thee phy_tree
  alignment <- AlignSeqs(dna) 
  phylo_tree <- upgma(dist.dna(as.DNAbin(alignment)))
  phylo_tree <- ladderize(phylo_tree)  # Arrange tree branches
  ps_ <- merge_phyloseq(ps_, phy_tree(phylo_tree))

  taxa_names(ps_) <- paste0("ASV", seq(ntaxa(ps_)))
  return(ps_)
}

# make_biom: create a BIOM object from the OTU data
otu2biom <- function(ps_,rRNA){
# Define the output file path
output_file_path <- paste0("output/output_file",rRNA,".biom")
# Check if the output file already exists
if (!file.exists(output_file_path)) {
  system("mkdir -p output")
  # Convert the otu_table to a matrix (required for make_biom)
  otu <- as(t(otu_table(ps_)), "matrix")
  
  # Create the BIOM object from the OTU data
  biom_data <- make_biom(data = otu)
  
  # Add taxonomy data manually if desired
  tax <- as.data.frame(as(tax_table(ps_), "matrix"))
  biom_data$rows <- lapply(seq_len(nrow(tax)), function(i) {
    biom_data$rows[[i]]$metadata <- list(taxonomy = tax[i, ])
    biom_data$rows[[i]]
  })
  
  # Write the BIOM object to a file
  write_biom(biom_data, output_file_path)
}
}

# plot_tax: plot the taxonomic composition at a specified taxonomic level
plot_tax <- function(ps_, group,NArm=FALSE) {
  # Aggregate at the specified taxonomic level
  ps_genus <- tax_glom(ps_, taxrank = group,NArm=NArm)
  
  # Transform counts to relative abundance
  ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x) * 100)
  
  # Convert to data frame for ggplot
  taxa_df <- psmelt(ps_genus_rel)
  
  # Create the plot 
  plot <- ggplot(taxa_df, aes(x = Sample, y = Abundance, fill = !!sym(group))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste0("Taxonomic Composition at ", group, " Level"), 
         x = "Sample", 
         y = "Relative Abundance (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",          # Position legend at the bottom
    legend.box = "horizontal")
  
  # Save the plot
  return(plot)
}

seq2fna <- function(ps_,rRNA){
  if (!file.exists(paste0("output/study_seqs_",rRNA,".fna"))) {
  ps_ %>%
        refseq() %>%
        Biostrings::writeXStringSet(paste0("output/study_seqs_",rRNA,".fna"), append=FALSE,
                                    compress=FALSE, compression_level=NA, format="fasta")
}
}
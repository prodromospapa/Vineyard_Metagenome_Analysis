Requirements:

1) Download 
    - SILVA dataset: 
        1) https://figshare.com/articles/dataset/SILVA_SSU_r138_2019_RData/25799980?file=46245217
        2) https://zenodo.org/records/4587955
    - UNITE dataset: 
        1) https://figshare.com/articles/dataset/UNITE_v2023_July2023/27018373?file=49181545
        2) https://doi.plutof.ut.ee/doi/10.15156/BIO/2959330


2) Make the picrust2 environment and install picrust2 and all the other necessary libraries in it(https://huttenhower.sph.harvard.edu/picrust/):

```bash
conda create -n picrust2 -c bioconda -c conda-forge picrust2 r-igraph r-curl
```

also download this [picrust2 reference files folder](https://github.com/picrust/picrust2/tree/master/picrust2/default_files), as it is needed for ITS, and put it in the main directory (do not rename the file).

RUN:

1) Enter the picrust2 environment:

```bash
conda activate picrust2
```

2) Run the `main.R` to export the microbial populations abundance and function, with the respectively plots:

```bash
Rscript main.R -d < directory of fastq files > -r < rRNA (16S or ITS) > -m < set method to assign taxonomy (dada2 or decipher), default: dada2 >
```
example:

```bash
Rscript main.R -d wine_varieties -r 16S -m dada2
```

Results:

You can find the results in the `output` folder

```
Dataframe:
    functional_(rRNA)_EC: This file contains the percentages of different enzyme commission (EC) numbers related to rRNA functions across all samples
    functional_(rRNA)_KO: This file contains the percentages of KEGG Orthologs (KOs) associated with rRNA across each sample
    functional_(rRNA)_pathways: This file contains the percentages of various metabolic or signaling pathways related to rRNA functions in each sample
    otu_(rRNA): This file contains the OTU (Operational Taxonomic Unit) table for rRNA data
    percentage_removed: This file reports the percentage of data that was removed during preprocessing (e.g. Mitochondria, Chloroplast etc.)
    taxa_(rRNA): This file contains taxonomic information related to rRNA sequences
Plots:
    alpha: This file contains the results of alpha diversity metrics for the microbiome data.
    correlation_EC: This file contains correlation data between samples according to EC (Enzyme Commission) value.
    correlation_KO: This file contains correlation data between samples according to KEGG Orthologs (KOs) value.
    correlation_pathways: This file contains correlation data between samples according to functional pathways value.
    functional_cluster_EC: This file contains the clustering results for EC numbers based on their relative abundances across samples.
    functional_cluster_KO: This file contains the clustering results for KEGG Orthologs (KOs).
    functional_cluster_pathways: This file contains the clustering results for metabolic or signaling pathways.
    functional_heatmap_EC: This file contains a heatmap visualization of the EC numbers across samples.
    functional_heatmap_KO: This file contains a heatmap visualization of KEGG Orthologs (KOs) across samples.
    functional_heatmap_pathways: This file contains a heatmap visualization of functional pathways across samples.
    microbiome_cluster: This file contains clustering results for microbiome communities based on Genus level value.
    plot_Class: This file contains a plot showing the distribution of microbiome samples at the class level.
    plot_Family: This file contains a plot showing the distribution of microbiome samples at the family level.
    plot_Genus: This file contains a plot showing the distribution of microbiome samples at the genus level.
    plot_Kingdom: This file contains a plot showing the distribution of microbiome samples at the kingdom level.
    plot_Order: This file contains a plot showing the distribution of microbiome samples at the order level.
    plot_Phylum: This file contains a plot showing the distribution of microbiome samples at the phylum level.
    plot_Species: This file contains a plot showing the distribution of microbiome samples at the species level.

```
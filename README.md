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
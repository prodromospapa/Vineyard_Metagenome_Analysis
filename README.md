Requirements:

1) Download 
    - SILVA dataset: https://figshare.com/articles/dataset/SILVA_SSU_r138_2019_RData/25799980?file=46245217
    - UNITE dataset: https://figshare.com/articles/dataset/UNITE_v2023_July2023/27018373?file=49181545

2) Install picrust2 using (https://huttenhower.sph.harvard.edu/picrust/):

    conda create -n picrust2 -c bioconda -c conda-forge picrust2

also download this [picrust2](https://github.com/picrust/picrust2/tree/master/picrust2) folder, as it is needed for ITS, and put it in the main directory.

RUN:

1) 
    Rscript main.R -d < directory of fastq files > -r < rRNA (16S or ITS) > -p < optional to print plots >

2) 
    source abundance.sh

Read Results:

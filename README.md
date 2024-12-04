Requirements:

1) Download 
    - SILVA dataset: 
        1) https://figshare.com/articles/dataset/SILVA_SSU_r138_2019_RData/25799980?file=46245217
        2) https://zenodo.org/records/4587955
    - UNITE dataset: 
        1) https://figshare.com/articles/dataset/UNITE_v2023_July2023/27018373?file=49181545
        2) https://doi.plutof.ut.ee/doi/10.15156/BIO/2959330


2) Install picrust2 using (https://huttenhower.sph.harvard.edu/picrust/):

'''ruby
    conda install -n picrust2 -c bioconda -c conda-forge picrust2
'''

also download this [picrust2](https://github.com/picrust/picrust2/tree/master/picrust2) folder, as it is needed for ITS, and put it in the main directory.

RUN:

1) This exports the results for abundance:

'''ruby
    Rscript main.R -d < directory of fastq files > -r < rRNA (16S or ITS) > -p < optional to print plots, default: FALSE > -m < set method to assign taxonomy (dada2 or decipher), default: dada2 >
'''
example:

'''ruby
    Rscript main.R -d wine_varieties -r 16S -p -m dada2
'''

2) This calculates the total microbiome functional abundances (you must have installed the picrust2 environment above):

'''ruby
    source functional.sh
'''

Read Results:

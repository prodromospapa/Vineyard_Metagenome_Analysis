conda activate picrust2

if [ $# -ne 1 ]; then
    echo "Usage: $0 <your_argument>"
    exit 1
fi

if [[ "$1" != "16S" && "$1" != "ITS" ]]; then
    echo "Error: Argument must be '16S' or 'ITS'."
    exit 1
fi

# Execute the command based on the input argument
if [ "$1" == "16S" ]; then
    picrust2_pipeline.py -s output/study_seqs_16S.fna -i output/output_file_$1.biom -o output/picrust_$1 -p 12
elif [ "$1" == "ITS" ]; then
    picrust2_pipeline.py -s output/study_seqs_ITS.fna -i output/output_file_$1.biom -o output/picrust_$1 -p 12 \
        --ref_dir picrust2/default_files/fungi/fungi_ITS \
        --custom_trait_tables picrust2/default_files/fungi/ec_ITS_counts.txt.gz \
        --marker_gene_table picrust2/default_files/fungi/ITS_counts.txt.gz \
        --reaction_func picrust2/default_files/fungi/ec_ITS_counts.txt.gz
fi

conda deactivate
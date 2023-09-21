# try running IQ-Tree

conda create --name iqtree

screen -S iqtree

conda activate iqtree

mamba install -c bioconda iqtree

# create new folder for output
cd ~/dada2_pipeline_myxorpoB/data/raxml

# try and run iqtree

# first look at which model is best
# only use models available in raxml
iqtree -s alignment/alignment_91percent.fasta -m MF -g constraint_trees/constraint_tree_91percent.tre --prefix output/myxo_91percent/myxo91 -T 5 -mset raxml

iqtree -s alignment/alignment_91percent.fasta -m GTR+F+R13 -g constraint_trees/constraint_tree_91percent.tre --prefix output/myxo_91percent/iq_tree_myxo91 -T 5

# test based on the best model
iqtree -s alignment/alignment_91percent.fasta -m GTR+G -g constraint_trees/constraint_tree_91percent.tre --prefix output/myxo_91percent/iq_tree_myxo91_v2 -T 5 -pers 0.2 -nstop 500

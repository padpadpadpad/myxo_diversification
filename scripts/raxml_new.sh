# if not done already - create environment for running RAXML
# install raxml-ng into it

# check raxml env
conda env list

# create an environment for raxml-ng
mamba create -n raxml_env 

# open the environment
conda activate raxml_env

mamba install -c bioconda raxml-ng=1.1.0

# create screen
screen -S raxml

conda activate raxml_env

# create new folder for output
cd ~/dada2_pipeline_myxorpoB/data/output/run_myxo_gtdbr202
mkdir -p raxml/myxo_97.7_constraint

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignments/new/alignment_myxo_97.7percent_new.fasta --model GTR+G --prefix raxml/myxo_97.7_constraint --threads 20 --seed 2 --tree-constraint alignments/new/contree.tre --outgroup otu_outgroup1

raxml-ng --msa alignments/new/alignment_myxo_97.7percent_new.fasta --model GTR+G --prefix raxml/myxo_97.7_constraint --threads 20 --seed 2 --tree-constraint alignments/new/contree.tre --outgroup otu_outgroup1
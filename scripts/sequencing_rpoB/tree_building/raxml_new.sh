# if not done already - create environment for running RAXML
# install raxml-ng into it

# create screen
screen -S raxml

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_asv

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_asv.fasta --model GTR+G --prefix trees/myxo_asv/myxo_asv --threads 20 --seed 2 --tree-constraint constraint_trees/constraint_tree_asv.tre

raxml-ng --msa alignment/alignment_asv.fasta --model GTR+G --prefix trees/myxo_asv/myxo_asv --threads 20 --seed 2 --tree-constraint constraint_trees/constraint_tree_asv.tre
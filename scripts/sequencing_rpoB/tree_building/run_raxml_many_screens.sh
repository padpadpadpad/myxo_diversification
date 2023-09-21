# fit phylogenetic tree with RaxML for each otu resolution

# if not done already - create environment for running RAXML
# install raxml-ng into it

# create an environment for raxml-ng
mamba create -n raxml_env 

# open the environment
conda activate raxml_env

mamba install -c bioconda raxml-ng=1.1.0

raxml-ng -v
# current version is here: https://github.com/amkozlov/raxml-ng

conda deactivate

ls

# 1. check the number of cores available
cat /proc/cpuinfo | grep processor | wc -l

# 2. move to the folder where the trees are to be made
cd ~/raxml_new

# 2. start a screen and run the alignment

#--------------#
# ASV level ####
#--------------#

screen -S raxml_asv

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_asv

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_asv.fasta --model GTR+G --prefix trees/myxo_asv/myxo_asv --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_asv.tre

raxml-ng --msa alignment/alignment_asv.fasta --model GTR+G --prefix trees/myxo_asv/myxo_asv --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_asv.tre

#-------------------#
# 99% similarity ####
#-------------------#

screen -S raxml_99

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_99

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_99percent.fasta --model GTR+G --prefix trees/myxo_99/myxo_99 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_99percent.tre

raxml-ng --msa alignment/alignment_99percent.fasta --model GTR+G --prefix trees/myxo_99/myxo_99 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_99percent.tre
# exit screen

#-------------------#
# 98% similarity ####
#-------------------#

screen -S raxml_98

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_98

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_98percent.fasta --model GTR+G --prefix trees/myxo_98/myxo_98 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_98percent.tre

raxml-ng --msa alignment/alignment_98percent.fasta --model GTR+G --prefix trees/myxo_98/myxo_98 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_98percent.tre

# exit screen

#-------------------#
# 97% similarity ####
#-------------------#

screen -S raxml_97

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_97

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_97percent.fasta --model GTR+G --prefix trees/myxo_97/myxo_97 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_97percent.tre

raxml-ng --msa alignment/alignment_97percent.fasta --model GTR+G --prefix trees/myxo_97/myxo_97 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_97percent.tre

# exit screen

#-------------------#
# 96% similarity ####
#-------------------#

screen -S raxml_96

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_96

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_96percent.fasta --model GTR+G --prefix trees/myxo_96/myxo_96 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_96percent.tre

raxml-ng --msa alignment/alignment_96percent.fasta --model GTR+G --prefix trees/myxo_96/myxo_96 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_96percent.tre

#-------------------#
# 95% similarity ####
#-------------------#

screen -S raxml_95

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_95

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_95percent.fasta --model GTR+G --prefix trees/myxo_95/myxo_95 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_95percent.tre

raxml-ng --msa alignment/alignment_95percent.fasta --model GTR+G --prefix trees/myxo_95/myxo_95 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_95percent.tre

# exit screen

#-------------------#
# 94% similarity ####
#-------------------#

screen -S raxml_94

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_94

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_94percent.fasta --model GTR+G --prefix trees/myxo_94/myxo_94 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_94percent.tre

raxml-ng --msa alignment/alignment_94percent.fasta --model GTR+G --prefix trees/myxo_94/myxo_94 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_94percent.tre

#-------------------#
# 93% similarity ####
#-------------------#

screen -S raxml_93

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_93

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_93percent.fasta --model GTR+G --prefix trees/myxo_93/myxo_93 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_93percent.tre

raxml-ng --msa alignment/alignment_93percent.fasta --model GTR+G --prefix trees/myxo_93/myxo_93 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_93percent.tre

# exit screen

#-------------------#
# 92% similarity ####
#-------------------#

screen -S raxml_92

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_92

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_92percent.fasta --model GTR+G --prefix trees/myxo_92/myxo_92 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_92percent.tre

raxml-ng --msa alignment/alignment_92percent.fasta --model GTR+G --prefix trees/myxo_92/myxo_92 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_92percent.tre

# exit screen

#-------------------#
# 91% similarity ####
#-------------------#

screen -S raxml_91

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_91

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_91percent.fasta --model GTR+G --prefix trees/myxo_91/myxo_91 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_91percent.tre

raxml-ng --msa alignment/alignment_91percent.fasta --model GTR+G --prefix trees/myxo_91/myxo_91 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_91percent.tre

# exit screen

#-------------------#
# 90% similarity ####
#-------------------#

screen -S raxml_90

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_90

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_90percent.fasta --model GTR+G --prefix trees/myxo_90/myxo_90 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_90percent.tre

raxml-ng --msa alignment/alignment_90percent.fasta --model GTR+G --prefix trees/myxo_90/myxo_90 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_90percent.tre

# exit screen

# exit screen

#---------------------#
# 97.7% similarity ####
#---------------------#

screen -S raxml_97.7

conda activate raxml_env

# create new folder for output
cd ~/raxml_new
mkdir -p trees/myxo_97.7

# check raxml - run it with a check first and with constraints
raxml-ng --check --msa alignment/alignment_97.7percent.fasta --model GTR+G --prefix trees/myxo_97.7/myxo_97.7 --threads 10 --seed 2 --tree-constraint constraint_trees/constraint_tree_97.7percent.tre

raxml-ng --msa alignment/alignment_97.7percent.fasta --model GTR+G --prefix trees/myxo_97.7/myxo_97.7 --threads 8 --seed 2 --tree-constraint constraint_trees/constraint_tree_97.7percent.tre

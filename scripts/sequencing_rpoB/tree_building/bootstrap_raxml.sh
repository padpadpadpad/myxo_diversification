# list environments using mamba
mamba env list

# move to raxml_mbe_response folder
cd raxml_mbe_response

# activate raxml_env using mamba
mamba activate raxml_env

# check folders
ls

# move into ASV
cd asv

ls

# running standard RaxML
raxml-ng --msa alignment_asv.fasta --model GTR+G --prefix boots --threads 25 --bootstrap --bs-trees autoMRE{500} --seed 42 --tree-constraint constraint_tree_asv.tre

cd myxo_97.7

# running standard RaxML
raxml-ng --msa alignment_97.7percent.fasta --model GTR+G --prefix boots --threads 25 --bootstrap --bs-trees autoMRE{500} --seed 42 --tree-constraint constraint_tree_97.7percent.tre

# running standard RaxML
raxml-ng --msa alignment_95percent.fasta --model GTR+G --prefix boots --threads 25 --bootstrap --bs-trees autoMRE{500} --seed 42 --tree-constraint constraint_tree_95percent.tre
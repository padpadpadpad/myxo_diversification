# create environment to run secsse on a cluster

mamba create -n r_env

mamba activate r_env

# download most recent r
mamba install -c conda-forge r-base=4.2.2

# install packages here
mamba install -c conda-forge r-curl
mamba install -c conda-forge r-ragg
mamba install -c conda-forge r-pak

R

pak::pkg_install("tidyverse")

conda deactivate
# run bamm

screen -S bamm

cd data/sequencing_rpoB/bamm

# run bamm
bamm -c bamm_control.txt

# run bamm for sampling fraction of .4
screen -S bamm_percent40

cd data/sequencing_rpoB/bamm

# run bamm
bamm -c bamm_control_0.4frac.txt

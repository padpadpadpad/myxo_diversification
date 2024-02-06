# run bamm on multiple screens

# open screen
screen -S bamm_sf1

# move to correct folder on cluster
cd ~/raxml_mbe_response/asv/bamm

bamm -c bamm_control_sf1.txt

# run sf0.5
screen -S bamm_sf0.5

# move to correct folder on cluster
cd ~/raxml_mbe_response/asv/bamm

bamm -c bamm_control_sf0.5.txt

# run sf0.25
screen -S bamm_sf0.25

# move to correct folder on cluster
cd ~/raxml_mbe_response/asv/bamm

bamm -c bamm_control_sf0.25.txt

# run sf0.125
screen -S bamm_sf0.125

# move to correct folder on cluster
cd ~/raxml_mbe_response/asv/bamm

bamm -c bamm_control_sf0.125.txt

# run sf0.0625
screen -S bamm_sf0.0625

# move to correct folder on cluster
cd ~/raxml_mbe_response/asv/bamm

bamm -c bamm_control_sf0.0625.txt

# open screen
screen -S bamm_sf1

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_95/bamm

bamm -c bamm_control_sf1.txt

# run sf0.5
screen -S bamm_sf0.5

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_95/bamm

bamm -c bamm_control_sf0.5.txt

# run sf0.25
screen -S bamm_sf0.25

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_95/bamm

bamm -c bamm_control_sf0.25.txt

# run sf0.125
screen -S bamm_sf0.125

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_95/bamm

bamm -c bamm_control_sf0.125.txt

# run sf0.0625
screen -S bamm_sf0.0625

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_95/bamm

bamm -c bamm_control_sf0.0625.txt

# open screen
screen -S bamm_sf1

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_97.7/bamm

bamm -c bamm_control_sf1.txt

# run sf0.5
screen -S bamm_sf0.5

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_97.7/bamm

bamm -c bamm_control_sf0.5.txt

# run sf0.25
screen -S bamm_sf0.25

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_97.7/bamm

bamm -c bamm_control_sf0.25.txt

# run sf0.125
screen -S bamm_sf0.125

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_97.7/bamm

bamm -c bamm_control_sf0.125.txt

# run sf0.0625
screen -S bamm_sf0.0625

# move to correct folder on cluster
cd ~/raxml_mbe_response/myxo_97.7/bamm

bamm -c bamm_control_sf0.0625.txt
# DPIM2 - The next-generation Drosophila Protein Interaction Map
DPIM2 - established from affinity purification-mass spectrometry of 5,805 baits covering the largest fraction yet of the Drosophila proteome. It contains 32,668 interactions among 3,644 proteins, organized into 632 clusters representing putative functional modules

Access DPIM2 data at https://artavanis-tsakonas.hms.harvard.edu/dpim2/

Sigma.js output network files are located in the <b>output/</b> folder

## scripts
#1 setup conda environment and environmental variables
setup.sh 

#2 run batch1 and batch2 DPiM2 network generation
scripts/only_batch2_year_adjustedOrthologHandling_run.sh
scripts/without_batch2_only_batch2_year_adjustedOrthologHandling_run.sh


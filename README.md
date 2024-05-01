# DPIM2 - The next-generation Drosophila Protein Interaction Map
DPIM2 - established from affinity purification-mass spectrometry of 5,805 baits covering the largest fraction yet of the Drosophila proteome. It contains 32,668 interactions among 3,644 proteins, organized into 632 clusters representing putative functional modules. It serves as a resource for predicting protein co-complex memberships and functional associations, as well as for generating functional hypotheses regarding poorly characterized proteins.

* DPIM2 is the next generation Drosophila Protein Interaction Map.
* It provides annotation for poorly studied genes and postulates novel protein-protein interaction relationships.
* Genetic validation confirms network predicted interactors can also modulate Notch signaling.
* Comparison with yeast and human networks identify conserved functional interactions.

Access DPIM2 data at https://artavanis-tsakonas.hms.harvard.edu/dpim2/

Sigma.js output network files are located in the <b>output/</b> folder

## scripts
1. setup conda environment and environmental variables<br>
> setup.sh 

2. run batch1 and batch2 DPiM2 network generation<br>
> ./scripts/only_batch2_year_adjustedOrthologHandling_run.sh<br>
> ./scripts/without_batch2_only_batch2_year_adjustedOrthologHandling_run.sh


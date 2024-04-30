source ../setup.sh

##Raw data 
##apmsData/

# skip failed steps and continue run
$DPSCR/reprocessed/preFilterMaster.pl -in /home/glocke/DPiM/augRemap/apmsData/raw.list -out filtered.out -donotremake

# run hyperG code
mkdir -p hyper/json
mkdir -p hyper/qdir

$DPSCR/runHyper.pl -in *.nrBait

# wait for all hper jobs are done then move to the next step
while true;
do
        test=$(qstat|grep hyper)
        if [ ${#test} -lt 1 ];
        then
                break
        fi

done

# go to hyper/qdir folder
cd hyper/qdir
ls $PWD/*sim*.o* > sim.list
# run compute FDR cutoff code
$DPSCR/compute_fdr_GL.pl -real hyper.o* -sim sim.list > hyper_cutoff.txt
cutoff=$(grep hyper hyper_cutoff.txt |cut -f2) # keep track of the cutoff

#go back to run folder
cd ../..
# generate final network
#cutoff=95.209
$DPSCR/scoreCutoff.pl -in hyper/qdir/hyper.o* -out hyper/nrBait.net -cutoff $cutoff


### clustering and annotation ###

# go into hyper folder
cd hyper/
# mcl iterative clustering plus breakup
$DPSCR/iterateCluster.pl -in nrBait.net -dir working/ -i 1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,4,5,6,7,8,9,10
cd working
ls $PWD/*txt > txt.list

# tidy up the clustering results, generating .tidy files
$UTILSCR/applyScriptToList.pl -l txt.list -s $DPSCR/secondary/tidyIterativeClusterNames.pl -ext tidy -in in
ls $PWD/*tidy > tidy.list

# get modularity and term annotation, generating .mod files
$UTILSCR/applyScriptToList.pl -l tidy.list -s $DPSCR/secondary/modularity.pl -ext mod -in mod -extra "-net ../nrBait.net -mode iter"


mkdir -p ./hyper/working
cd ./hyper
cp /home/kli3/HGTestRun/run11_reordered_prefiltering_removedbadyear_08082018/hyper/nrBait.net ./
cp -r /home/kli3/HGTestRun/run11_reordered_prefiltering_removedbadyear_08082018/hyper/qdir ./
cd ./working
cp /home/kli3/HGTestRun/run11_reordered_prefiltering_removedbadyear_08082018/hyper/working/mcl.clusters.i*.txt ./
cp /home/kli3/HGTestRun/run11_reordered_prefiltering_removedbadyear_08082018/hyper/working/*.mod ./
cp /home/kli3/HGTestRun/run11_reordered_prefiltering_removedbadyear_08082018/hyper/working/*.tidy ./
cp /home/kli3/HGTestRun/run11_reordered_prefiltering_removedbadyear_08082018/hyper/working/*.hypGeo.termTest ./
ls $PWD/*tidy > tidy.list

# annotate clusters and generate .hypGeo.termTest files
$DPSCR/rscripts/annotateClusters.pl -in tidy.list -ext hypGeo.termTest

# generate cluster statistics for different i's
$DPSCR/secondary/tertiary/clustStats.pl -in tidy.list > ../mcl_i_differences.txt #the list of file names of tidy which is prefix before hypGeo.termTest, because the script is looking for .hypGeo.termTest files

# annotate nodes, generate .node.annot files
ls $PWD/*termTest > termTest.list
/home/glocke/DPiM/scripts/rscripts/annotateNodes.R mcl.clusters.i2.txt r.node.annot # any cluster file should generate the same output file, but column 2 is different because cluster ids are different
/home/glocke/utilScripts/applyScriptToList.pl -l tidy.list -s /home/glocke/DPiM/scripts/secondary/annotateNodes.pl -ext node.annot -in in -l2 termTest.list -l2arg clusterannot -extra "-terms r.node.annot"

# annotate network, generate net.annot files
ls $PWD/*node.annot |grep -v r.node > node.annot.list
/home/glocke/utilScripts/applyScriptToList.pl -l node.annot.list -s /home/glocke/DPiM/scripts/secondary/annotateNetworkMore.pl -ext net.annot -in annot -extra "-net ../nrBait.net"

# novelty metrics annotation (Dec 7 update pptx slide #6)
/home/glocke/DPiM/scripts/metrics/novelAnnotation.pl -out novelty_annotation.txt -newcluster mcl.clusters.i1.8.txt.tidy -newnet ../nrBait.net -newnodes mcl.clusters.i1.8.txt.tidy.node.annot
cut -f1,2,3,4,7,8,9 novelty_annotation.txt > ../novelty_annotation.txt

# assembly the deliverables
cd ../.. #back to run folder
mkdir deliverables
cd deliverables
ln -s ../hyper/nrBait.net
ln -s ../hyper/mcl_i_differences.txt
ln -s ../hyper/novelty_annotation.txt
ln -s ../hyper/qdir/hyper.o* all_Edges_nrBait_HGScore.txt
ln -s ../*rebrand.nrBait all_Edges_nrBait_SpecCount.txt
cd ..

## very end would be pick the best i value, then bring in the .tidy.node.annot and .tidy.node.annot.net.annot files here

## FDR 0.1 version ##
cd hyper/qdir
/home/glocke/DPiM/scripts/compute_fdr_GL.pl -real hyper.o* -sim sim.list -fdr 0.1 > hyper_cutoff_fdr_10_prct.txt
cutoff_fdr_10_prct=$(grep hyper hyper_cutoff_fdr_10_prct.txt |cut -f2) # keep track of the cutoff
#go back to run folder
cd ../..
# generate final network

$DPSCR/scoreCutoff.pl -in hyper/qdir/hyper.o* -out hyper/nrBait_fdr_10_prct.net -cutoff $cutoff_fdr_10_prct

# go into hyper folder
cd hyper/
# mcl iterative clustering plus breakup
$DPSCR/iterateCluster.pl -in nrBait_fdr_10_prct.net -dir working_fdr_10_prct/ -i 1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,4,5,6,7,8,9,10
cd working_fdr_10_prct
ls $PWD/*txt > txt.list

# tidy up the clustering results, generating .tidy files
$UTILSCR/applyScriptToList.pl -l txt.list -s $DPSCR/secondary/tidyIterativeClusterNames.pl -ext tidy -in in
ls $PWD/*tidy > tidy.list

# get modularity and term annotation, generating .mod files
$UTILSCR/applyScriptToList.pl -l tidy.list -s $DPSCR/secondary/modularity.pl -ext mod -in mod -extra "-net ../nrBait_fdr_10_prct.net -mode iter"

# annotate clusters and generate .hypGeo.termTest files
$DPSCR/rscripts/annotateClusters.pl -in tidy.list -ext hypGeo.termTest

# generate cluster statistics for different i's
$DPSCR/secondary/tertiary/clustStats.pl -in tidy.list > ./mcl_i_differences.txt #the list of file names of tidy which is prefix before hypGeo.termTest, because the script is looking for .hypGeo.termTest files

# annotate nodes, generate .node.annot files
ls $PWD/*termTest > termTest.list
/home/glocke/DPiM/scripts/rscripts/annotateNodes.R mcl.clusters.i2.txt r.node.annot # any cluster file should generate the same output file, but column 2 is different because cluster ids are different
/home/glocke/utilScripts/applyScriptToList.pl -l tidy.list -s /home/glocke/DPiM/scripts/secondary/annotateNodes.pl -ext node.annot -in in -l2 termTest.list -l2arg clusterannot -extra "-terms r.node.annot"

# annotate network, generate net.annot files
ls $PWD/*node.annot |grep -v r.node > node.annot.list
/home/glocke/utilScripts/applyScriptToList.pl -l node.annot.list -s /home/glocke/DPiM/scripts/secondary/annotateNetworkMore.pl -ext net.annot -in annot -extra "-net ../nrBait_fdr_10_prct.net"

# novelty metrics annotation (Dec 7 update pptx slide #6)
/home/glocke/DPiM/scripts/metrics/novelAnnotation.pl -out novelty_annotation.txt -newcluster mcl.clusters.i1.8.txt.tidy -newnet ../nrBait_fdr_10_prct.net -newnodes mcl.clusters.i1.8.txt.tidy.node.annot
cut -f1,2,3,4,7,8,9 novelty_annotation.txt > ./novelty_annotation.txt


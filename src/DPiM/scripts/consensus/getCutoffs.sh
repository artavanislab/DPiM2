#!/bin/bash

#perl compute_fdr_KJ_consensus.pl ../George_consensus_network_try1/qdir/HGConsensus000???.sh.o* > cutoffs.txt


real="/home/glocke/DPiM/dpim4/withInstr/consensus3_0.3-10-2016/qdir/real.hyperspec.list"
job="buildNet"
cutFile="cutoffs5.txt"
outFile="/home/glocke/DPiM/dpim4/withInstr/consensus3_0.3-10-2016/support.net"

#clean temp folder and previous job output
rm nets/*
#rm $job.*


# get all networks with cutoffs
/home/glocke/DPiM/scripts/qsubBuildNetwork.pl -real $real -job $job -outdir $PWD/nets

size=0
while [ $size -lt 2000 ]
do
	#echo $size
	size="$(ls -l $job.o*|wc -l)"
done

echo "done with all permutations"

# create file with cutoffs
cat *.$job.o* > cutoffs1.txt

# clean up the dir
rm *.$job.*

# merge all networks
cat nets/*fdr.out > nets/all.txt
echo "All network files merged"

# only keep the interactors
awk '{print $1, "\t", $2}' nets/all.txt > nets/all_interactions.txt

# count the frequency
sort nets/all_interactions.txt | uniq -c > nets/all_interactions_counts.txt
echo "consensus network file generated"

# mv the file to current folder
perl -E 'say join "\t", qw(prey1 prey2 support); while (<>) { chomp; @spl = split; say join "\t", $spl[1], $spl[2], $spl[0];} ' nets/all_interactions_counts.txt > $outFile

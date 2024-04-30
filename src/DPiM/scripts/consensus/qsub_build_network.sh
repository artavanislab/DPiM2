for i in $@
do
echo "perl compute_fdr_KJ_consensus_dir2.pl $i" | qsub -l h_rt=00:01:00 -cwd 
done

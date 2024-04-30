module add gsl
module add gcc # necessary for hyperspec

DP3DATA="DPiM_Data_Summary_2014_edited_KJ.txt"
DP2DATA="dpim2_nrtap.120123"
FBGN2SYMB="/home/glocke/DPiM/fbgn_id2_2col.txt"

THISDATA=$DP3DATA
NEXTDATA=$THISDATA.reformatted.out

echo "data_plumbing.sh"
$DPSCR/data_plumbing.sh $THISDATA $NEXTDATA

THISDATA=$NEXTDATA
NEXTDATA="dpim3.rmDupes"

echo rmDupes
$DPSCR/removeDuplicateRuns.pl -in $THISDATA -out $NEXTDATA

THISDATA=$NEXTDATA
NEXTDATA=$THISDATA.sumIso

echo sumIso
$DPSCR/sumIsoforms.pl -in $THISDATA -out $NEXTDATA

THISDATA=$NEXTDATA
COMMCONTAM=$THISDATA.commonContam

echo calcCommonContam
$DPSCR/calcCommonContam_KJ_GL.pl -in $THISDATA -out $COMMCONTAM

SEQCONTAM=$THISDATA.seqContam

echo sequentialContam
$DPSCR/find_sequential_contam_KJ_GL.pl $FBGN2SYMB $THISDATA $COMMCONTAM $SEQCONTAM

NEXTDATA=$THISDATA.applyLC
echo "apply_lc_results_KJ.pl"
$DPSCR/apply_lc_results_KJ.pl $SEQCONTAM $THISDATA  > $NEXTDATA

THISDATA=$NEXTDATA
NEXTDATA=dpim3.09-25-2015

cat $DP2DATA $THISDATA > $NEXTDATA

# produce the nrtap file by using all the data

THISDATA=$NEXTDATA
NEXTDATA=$THISDATA.nrBait

echo "produce_nrtap.pl"
perl $DPSCR/produce_nrtap_9-23-2015.pl $THISDATA > $NEXTDATA


## show some statistics about dpim2 and dpim3
#THISDATA=$NEXTDATA
#echo "check DPIM2 vs DPIM3"
#Rscript $DPSCR/check_on_DPIM2_vs_3.r $DP2DATA $THISDATA

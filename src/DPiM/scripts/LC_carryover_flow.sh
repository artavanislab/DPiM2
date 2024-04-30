module add gsl
module add gcc # necessary for hyperspec

FBGN2SYMB="/home/glocke/DPiM/fbgn_id2_2col.txt"

echo extractDateFromExcelCSV.pl
$DPSCR/extractDateFromExcelCSV.pl

THISDATA="/home/glocke/DPiM/dpim4/withInstr/all.combined.01-30-2016"

echo combine234.pl
$DPSCR/combine234.pl -out $THISDATA

NEXTDATA=$THISDATA.plumb

echo "data_plumbing.sh"
$DPSCR/data_plumbing.sh $THISDATA $NEXTDATA

THISDATA=$NEXTDATA
NEXTDATA=$THISDATA.rmDupes

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
$DPSCR/find_sequential_contam_KJ_GL.pl $FBGN2SYMB $THISDATA $COMMCONTAM > $SEQCONTAM

NEXTDATA=$THISDATA.applyLC
echo "apply_lc_results_KJ.pl"
$DPSCR/apply_lc_results_KJ.pl $SEQCONTAM $THISDATA  > $NEXTDATA

THISDATA=$NEXTDATA
NEXTDATA=$THISDATA.0_05.tscFilter
    my $tscCut = $sumIso.".0_05.tscFilter";
    $cmd = "$scrDir/tscCutoff.pl -in $sumIso -out $tscCut";
    say $cmd;
    system($cmd);

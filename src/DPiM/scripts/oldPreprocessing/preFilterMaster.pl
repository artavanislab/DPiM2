#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Path qw(make_path);
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4APMS);

# coordinate pre-filtering of DPiM data
#

my %opts = getCommandLineOptions();

{
    my $scrDir = $opts{scrdir}; # location of scripts
    #my $dataDir = $opts{data};
    my $baseOut = $opts{out};
    my $id2symb = $opts{id2symb};
    my $remakeFlag = ! $opts{donotremake}; # DO remake if this flag is true
    
    # update FBgn ids to their most recent versions
    my $updateFB = "$baseOut.newFBgn";
    if (!$remakeFlag && foundFile($updateFB)) {
	say "\tusing previous $updateFB";
    } else {
	my $cmd = "$scrDir/updateFBgn.pl -in $baseOut -out $updateFB -ref $opts{fbgnstandard}";
	say $cmd;
	system($cmd);
	$remakeFlag = 1;
    }
    
    my $rmDupes = "$updateFB.rmDupes";
    if (!$remakeFlag && foundFile($rmDupes)) {
	say "\tusing previous $rmDupes";
    } else {
	my $cmd = "$scrDir/removeDuplicateRuns.pl -in $updateFB -out $rmDupes";
	say $cmd;
	system($cmd);
	$remakeFlag = 1;
    }

	
    my $sumIso = $rmDupes.".sumIso";
    if (!$remakeFlag && foundFile($sumIso)) {
	say "\tusing previous $sumIso";
    } else {
	my $cmd = "$scrDir/sumIsoforms.pl -in $rmDupes -out $sumIso";
	say $cmd;
	system($cmd);
    }

    if (0) {
	my $oneHitWonder = $sumIso.".ohw";
	if (!$remakeFlag && foundFile($oneHitWonder)) {
	    say "\tusing previous $oneHitWonder";
	} else {
	    my $cmd = "$scrDir/oneHitWonders.pl -in $sumIso -out $oneHitWonder";
	    say $cmd;
	    system($cmd);
	}
    }

    my $untrans = "$sumIso.untrans.txt";
    if (!$remakeFlag && foundFile($untrans)) {
	say "\tusing previous $untrans";
    } else {
	my $cmd = "$scrDir/findUntranslatedFBgn.pl -in $sumIso -out $untrans -avlen $opts{lenfile}";
	say $cmd;
	system($cmd);
    }
   
    my $trans = "$sumIso.trans";
    if (!$remakeFlag && foundFile($trans)) {
	say "\tusing previous $trans";
    } else {
	my $cmd = "$scrDir/removeUntranslated.pl -in $sumIso -out $trans -remove $untrans";
	say $cmd;
	system($cmd);
    }
    
    my $baitStats = "$trans.statsBySID";
    if (!$remakeFlag && foundFile($baitStats)) {
	say "\tusing previous $baitStats";
    } else {
	my $cmd = "$scrDir/baitPullDown.pl -in $trans -out $baitStats";
	say $cmd;
	system($cmd);
    }
    
    my $withBait = "$trans.withBait";
    my $noBait = "$trans.noBait";
    if (!$remakeFlag && foundFile($withBait)) {
	say "\tusing previous $withBait";
    } else {
	my $cmd = "$scrDir/selectNoBait.pl -apms $trans -stats $baitStats -clean $withBait -nobait $noBait";
	say $cmd;
	system($cmd);
	$remakeFlag = 1;
    }

    my $rebrand = "$trans.rebrandedBait";
    if (!$remakeFlag && foundFile($rebrand)) {
	say "\tusing previous $rebrand";
    } else {
	my $cmd = "$scrDir/rebrandBait.pl -nobait $noBait -withbait $withBait -out $rebrand";
	say $cmd;
	system($cmd);
	$remakeFlag = 1;
    }

    my $tscCut = $rebrand.".0_05.tscFilter";
    if (!$remakeFlag && foundFile($tscCut)) {
	say "\tusing previous $tscCut";
    } else {
	my $cmd = "$scrDir/tscCutoff.pl -in $rebrand -out $tscCut -statsscr $scrDir/baitPullDown.pl";
	say $cmd;
	system($cmd);
    }
    
    my $commContam = $tscCut.".commContam";
    if (!$remakeFlag && foundFile($commContam)) {
	say "\tusing previous $commContam";
    } else {
	my $cmd = "$scrDir/calcCommonContam_KJ_GL.pl -in $tscCut -out $commContam";
	say $cmd;
	system($cmd);
	$remakeFlag = 1;
    }

    my $seqContam = $tscCut.".seqContam";
    if (!$remakeFlag && foundFile($seqContam)) {
	say "\tusing previous $seqContam";
    } else {
	my $cmd = "$scrDir/find_sequential_contam_KJ_GL.pl $id2symb $tscCut $commContam > $seqContam";
	say $cmd;
	say "\t...this should take about an hour";
	system($cmd);
	$remakeFlag = 1;
    }

    my $applyLC = $tscCut.".applyLC";
    if (!$remakeFlag && foundFile($applyLC)) {
	say "\tusing previous $applyLC";
    } else {
	my $cmd = "$scrDir/apply_lc_results_KJ.pl $seqContam $tscCut  > $applyLC";
	say $cmd;
	system($cmd);
	$remakeFlag = 1;
    }

    my $storable = $baseOut."toscore.storable";
    if (!$remakeFlag && foundFile($storable)) {
	say "I'm a little confused, because the end product of this script, $storable, already exists, so, apparently, this script has done nothing...";
    } else {
	my $cmd = "$scrDir/sortAPMSByBait.pl -in $applyLC -out $storable";
	say $cmd;
	system($cmd);
	$remakeFlag = 1;
    }
    say "\n\n\nFinished!";
    say "perl-readable APMS data are found at $storable";
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	dp2 => '/home/glocke/DPiM/dpim4/dpim2_all.120123.plusUniqMsInst.cp',
	dp3 => '/home/glocke/DPiM/dpim4/DPiM_Data_Summary_2014_edited.reformatted.out',
	dp4 => '/home/glocke/DPiM/dpim4/DPiM3_Raw_data_2015_text.colFix',
	dp4date => '/home/glocke/DPiM/dpim4/excelCsv/sid2DateInstrID.tsv',
	scrdir => $ENV{DPSCR} // "must_supply",
	removetsc => 0.05, # remove bottom fraction of TSC/experiment
	id2symb => '/home/glocke/DPiM/fbgn_id2_2col_03-17-2016.txt',
	fbgnstandard => '/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2016_01.tsv',
	lenfile => '/home/glocke/DPiM/nsaf/dmel-all-translation-r6.09.aveLen.tsv',
	#data => '/home/glocke/DPiM/dpim4',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in apms.in -out filtered.out < $defaultString ".
	"-donotremake > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "dp2=s", "dp3=s", "dp4=s", "dp4date=s", 
	       "scrdir=s", "dpscr=s", "removetsc=f", "id2symb=s", 
	       "fbgnstandard=s", "lenfile=s", "donotremake");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    die "must set a valid -scrdir where scripts are found" 
	if $opts{scrdir} eq "must_supply";
    
    checkExist('f', $opts{id2symb});
    checkExist('f', $opts{fbgnstandard});
    checkExist('f', $opts{lenfile});
    checkExist('d', $opts{scrdir});
    checkExist('f', $opts{dp2});
    checkExist('f', $opts{dp3});
    checkExist('f', $opts{dp4});
    checkExist('f', $opts{dp4date});

    my @scripts = qw(
	extractDateFromExcelCSV.pl
	combine234.pl
	removeDuplicateRuns.pl
	sumIsoforms.pl
	findUntranslatedFBgn.pl
	removeUntranslated.pl
	baitPullDown.pl
	selectNoBait.pl
	tscCutoff.pl
	calcCommonContam_KJ_GL.pl
	find_sequential_contam_KJ_GL.pl
	apply_lc_results_KJ.pl
	sortAPMSByBait.pl
	);
    #oneHitWonders.pl
    for my $s (@scripts) {
	my $scrDir = $opts{scrdir};
	die "can't find $scrDir/$s" unless -e "$scrDir/$s";
    }
    
    return %opts;
}

# return true if this file does not need to be made/remade
# that is, return true if 
#  * this file is readable
#  * this file has non-zero size
#  **and we are not intending to remake it anyway**
sub foundFile {
    my $f = shift;
    return defined $opts{donotremake} && -r $f && -s $f;
}


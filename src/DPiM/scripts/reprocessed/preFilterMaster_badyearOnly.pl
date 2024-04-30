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
# OLD steps: 
# 1. convert data to dp4 format and collect
# 2. assign bait
# 3. enforce peptide FDR using Julians identity score
# 4. collect peptides for a given FBgn within each experiment (sum isoforms)
# 5. remove duplicate runs
# 6. LC carryover
#    6.1.  find common contaminants 
#    6.2.  find sequential contaminants (lengthiest step in whole process)
#    6.3.  apply LC carryover
# 7. remove FBgn0000000 and reverse peptides
# 8. remove runs with particularly high or low peptide counts
# 9. Rebrand runs where the bait does not appear
#      note: there is no procedure to keep track of this for nrBait!!!
#my @order = qw(dp4 newBait rmDupes pepFDR sumIso applyLC trans tscCut rebrand nrBait);

# May 25 2018 redefine order
my @order = qw(dp4 newBait rebrand pepFDR applyLC sumIso rmDupes trans tscCut nrBait);
# NEW steps:
# # 1. convert data to dp4 format and collect
# # 1.1 removing badyear experiment (2009-04-01 till 2010-06-01) Added by Kejie, 08082018
# # 2. assign bait
# # 3. Rebrand runs where the bait does not appear
# #      note: there is no procedure to keep track of this for nrBait!!!
# # 4. enforce peptide FDR using Julians identity score
# # 5. LC carryover
# # #    5.1.  find common contaminants
# # #    5.2.  find sequential contaminants (lengthiest step in whole process)
# # #    5.3.  apply LC carryover
# # 6. collect peptides for a given FBgn within each experiment (sum isoforms)
# # 7. remove duplicate runs
# # 8. remove FBgn0000000 and reverse peptides
# # 9. remove runs with particularly high or low peptide counts
# # 10. feed the filtered data to HG-Scoring pipeline which generates nrBait net

my %scripts = (
    #dp4 => 'reprocessed/toDP4List.pl',
	#dp4 => 'reprocessed/toDP4List_rmBadYear.pl', #removing badyear experiment (2009-04-01 till 2010-06-01) Added by Kejie, 08082018
	dp4 => 'reprocessed/toDP4List_keepBadYearOnly.pl', #run badyear experiment (2009-04-01 till 2010-06-01) only Added by Kejie, 09122018
    newBait => 'reprocessed/identifyBait10-31-2016.pl',
    rmDupes => 'reprocessed/removeDuplicateRuns2.pl',
    pepFDR => 'reprocessed/enforceFDR.pl',
    sumIso => 'sumIsoforms.pl',
    commContam => 'calcCommonContam_KJ_GL.pl',
    seqContam => 'find_sequential_contam_KJ_GL.pl',
    applyLC => 'wrapLCCorrection.pl',
    qsub => 'qsubWrap.pl', # called by applyLC
    trans => 'removeUntranslated2.pl',
    tscCut => 'tscCutoff.pl',
    #rebrand => 'wrapRebrand.pl',
    rebrand => 'wrapRebrand_rare.pl',
    statsBySID => 'baitPullDown.pl',  # called by wrapRebrand.pl
    selectNoBait => 'selectNoBait.pl',
    nrBait => 'nrBait.pl',
    );

#/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.11-14-2016_reject.tsv
my %opts = getCommandLineOptions();

{
    #$scripts{newBait}.= " -bait /home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.11-14-2016_reject.tsv";
    #KJ modified this on 04/05/2017 to replace the old tsv by the new tsv file (decided by Harsha and Kejie)
    $scripts{newBait}.= " -bait /home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.04-05-2017_reject.tsv";
    
    my $in = $opts{in};
    my $scrDir = $opts{scrdir}; # location of scripts
    #my $dataDir = $opts{data};
    my $baseOut = $opts{out};
    my $id2symb = $opts{id2symb};
    my $remakeFlag = ! $opts{donotremake}; # DO remake if this flag is true

    my $thisIn = $in;
    my $thisOut = $baseOut;

    my $baitlessStatsFile; ## nrBait musn't prefer a rebranded experiment
    
    for my $scr (@order) {
	$thisOut.= ".$scr";
	if (!$remakeFlag && foundFile($thisOut)) {
	    say "\tusing previous $thisOut";
	} else {
		if ($scr eq 'rebrand') {
        # due to the  May 25 2018 redefine order change, the newBait out file needs to be modified a little
			my $cmd = 'cut -f1-6 '.$thisIn.' > temp; perl -pi -e \'s/search_id/\nsearch_id/\' temp; mv temp '.$thisIn;
			say $cmd;
			system($cmd);
		}
		if ($scr eq 'pepFDR') {
		# due to the  May 25 2018 redefine order change, the rebrand out file needs to be modified a little
			my $cmd = 'awk \'FNR==NR{a[$1$3$4$5$6]=$7;next} ($1$3$4$5$6 in a) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"a[$1$3$4$5$6]}\' filtered.out.dp4 '.$thisIn.'> temp; mv '.$thisIn.' '.$thisIn.'.bk; mv temp '.$thisIn;
			say $cmd;
			system($cmd);
		}
	    my $cmd = "$scrDir/$scripts{$scr} -in $thisIn -out $thisOut";
	    if ($scr eq 'nrBait') {
			$cmd .= " -baitless $baitlessStatsFile -log nrBait.log";
	    }
	    say $cmd;
	    system($cmd);
	}
	$thisIn = $thisOut;
	if ($scr eq 'rebrand') {
	    $baitlessStatsFile = "$thisOut.working.input.stats";
	}
    }

}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	scrdir => $ENV{DPSCR} // "must_supply",
	#data => '/home/glocke/DPiM/dpim4',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in raw.apms.list -out filtered.out < ".
	"$defaultString -donotremake > \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", 
	       "scrdir=s", "dpscr=s", "removetsc=f", "id2symb=s", 
	       "fbgnstandard=s", "lenfile=s", "donotremake");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    die "must set a valid -scrdir where scripts are found" 
	if $opts{scrdir} eq "must_supply";
    
    checkExist('d', $opts{scrdir});

    #oneHitWonders.pl
    for my $s (values %scripts) {
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


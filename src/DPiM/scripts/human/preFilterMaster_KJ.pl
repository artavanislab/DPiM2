#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Spec;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# 1. convert to biop format
#   1.2. separate by instrument # test effect of switching the order
# 2. LC carryover
#   2.2. common contaminants    # of these two
#   2.3. sequential contaminants
#   2.4. applyLC
#   2.5. cat the results back together
# 3. require min 2 peptides
# 4. Entropy filter
# 5. per-plate filter
# 6. require known length
# 7. ...
# 8. Profit!

my %opts = getCommandLineOptions();

my $SCR = '/home/glocke/DPiM/scripts';

#my @order = qw(biop sepInst commContam seqContam applyLC catInst minPeptides entropy plate);
#my @order = qw(biop applyLC minPepEntropy plate); # 01/22/2018 Harsha and KJ decided to take out applyLC step simply rely on entropy would be sufficient
my @order = qw(biop minPepEntropy plate);

my %scripts = (
    biop => "$SCR/human/toBiop_KJ.pl",
    sepInst => "$SCR/human/separateInstrument.pl",
    commContam => "$SCR/calcCommonContam_KJ_GL.pl -mode human",
    seqContam => "$SCR/find_sequential_contam_KJ_GL.pl",
	applyLC => "$SCR/wrapLCCorrection.pl",
    minPepEntropy => "$SCR/human/minPepEntropy.pl",
    plate => "$SCR/human/perPlate.pl",
    );
{
    #my $in = $opts{in};
    my $baseOut = File::Spec->rel2abs($opts{out});

    my $remakeFlag = ! $opts{donotremake}; # DO remake if this flag is true

    my $thisIn = $baseOut;
    my $thisOut = $baseOut;

    for my $scr (@order) {
    $thisOut.= ".$scr";
    if (!$remakeFlag && foundFile($thisOut)) {
        say "\tusing previous $thisOut";
    } else {
        my $cmd = "$scripts{$scr} -in $thisIn -out $thisOut";
        #if ($scr eq 'nrBait') {
        #$cmd .= " -baitless $baitlessStatsFile -log nrBait.log";
        #}
        if ($scr eq 'applyLC') {
			$cmd .= " -human"; # requires the -human tag
		}
        say $cmd;
        system($cmd);
    }
    $thisIn = $thisOut;
    #if ($scr eq 'rebrand') {
    #    $baitlessStatsFile = "$thisOut.working.input.stats";
    #}
    }
    
}


exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   


sub getCommandLineOptions {

    my %defaults = (
    output => '/home/kli3/proj/Interactome/data/BioPlex2.0Nature/gpsRes4_output.tsv',
    summary => '/home/kli3/proj/Interactome/data/BioPlex2.0Nature/gpsRes4_runSummary.tsv'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString -donotremake > \n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "output=s", "summary=s", "-donotremake");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{output});
    checkExist('f', $opts{summary});

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


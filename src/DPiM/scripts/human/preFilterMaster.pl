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

my %scripts = (
    biop => "$SCR/human/toBiop.pl",
    sepInst => "$SCR/human/separateInstrument.pl",
    commContam => "$SCR/calcCommonContam_KJ_GL.pl",
    seqContam => "$SCR/find_sequential_contam_KJ_GL.pl",
    applyLC => "$SCR/apply_lc_results_KJ.pl",
    catInst => "$SCR/human/concatInstrument.pl",
    minPeptides => "$SCR/human/minPep.pl",
    entropy => "$SCR/human/entropyFilter.pl",
    plate => "$SCR/human/perPlate.pl",
    );
{
    my $in = $opts{in};
    my $baseOut = File::Spec->rel2abs($opts{out});

    
}


exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   


sub getCommandLineOptions {

    my %defaults = (
	output => '/home/kli3/Interactome/data/BioPlex_8400_03282016/gpsRes6_output.tsv',
	summary => '/home/kli3/Interactome/data/BioPlex_8400_03282016/gpsRes6_runSummary.tsv'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "output=s", "summary=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{output});
    checkExist('f', $opts{summary});

    return %opts;
}

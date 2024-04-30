#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use List::MoreUtils qw(each_array);
use HomeBrew::IO qw(checkExist readList);
#use DpimLib qw(getLineDP4APMS);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $realList = $opts{real};
    my $simList = $opts{sim};
    my $job = $opts{job};
    my $outDir = $opts{outdir};
    
    my @reals = readList($realList);
    my @sims = readList($simList);

    checkExist('f', $_) for @reals, @sims;

    my $it = each_array(@reals, @sims);	    
    while (my ($r, $s) = $it->()) {
	my $out = basename($r);
	$out =~ s/\.o\d+$//;
	$out = "$outDir/$out.fdr.out";
	my $cmd = qq{echo "perl compute_fdr_KJ_consensus_dir2.pl $r $s $out" | qsub -l h_rt=00:01:00 -cwd -N $job};
	system($cmd);
    }
    
    
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	job => 'buildnet',
	outdir => "$ENV{PWD}/nets",
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -real real.list -sim sim.list";

    my %opts = ();
    GetOptions(\%opts, "real=s", "sim=s", "job=s", "outdir=s");
    die $usage unless exists $opts{real} && exists $opts{sim};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{real});
    checkExist('f', $opts{sim});

    return %opts;
}

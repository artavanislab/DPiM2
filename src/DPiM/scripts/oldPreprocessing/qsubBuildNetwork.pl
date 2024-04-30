#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use File::Path qw(make_path);
use List::MoreUtils qw(each_array);
use HomeBrew::IO qw(checkExist readList);
#use DpimLib qw(getLineDP4APMS);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $realList = $opts{real};
    my $jobExt = $opts{jobext};
    my $outDir = $opts{outdir};
    my $simList = $opts{sim};
    
    unless (-d $outDir) {
	make_path($outDir) or die "can't make directory $outDir";
    }

    my @reals = readList($realList);
    my @sims;
    @sims = readList($simList) if defined $simList;

    checkExist('f', $_) for @reals, @sims;

    my $qsss = "/home/glocke/utilScripts/qsubWrap.pl";
    my $compFdr = "/home/glocke/DPiM/scripts/compute_fdr_KJ_consensus_dir2.pl";
    my $rt = '00:01:00';
    
    for my $i (0..$#reals) {
	my $r = $reals[$i];
	my $job = basename($r);
	$job =~ s/\.o\d+$//;
	my $s;
	if (defined $simList) {
	    $s = $sims[$i];
	} else {
	    my $dir = dirname($r);
	    my @sims = `ls $dir/$job.sim*.o*`;
	    chomp @sims;
	    die "found ", (0+ @sims), " simulations for $job" if @sims != 1;
	    $s = $sims[0];
	} 

	my $out = "$outDir/$job.fdr.out";
	my $cmd = "$qsss -cmd '$compFdr $r $s $out' -job $job.$jobExt -rt $rt";
	say $cmd;
	system($cmd) unless exists $opts{nosub};
    }
    
    
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	jobext => 'fdr',
	outdir => "$ENV{PWD}/nets",
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -real real.list < $defaultString -sim sim.list ".
	"-nosub >";

    my %opts = ();
    GetOptions(\%opts, "real=s", "jobext=s", "outdir=s", "sim=s", "nosub");
    die $usage unless exists $opts{real};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{real});
    checkExist('f', $opts{sim}) if exists $opts{sim};

    return %opts;
}

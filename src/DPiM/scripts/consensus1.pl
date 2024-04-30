#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Spec;
use File::Path qw(make_path);
use File::Temp qw(tempfile);
use Storable;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw();

# this version of this script is aimed at evaluating convergence

# coordinate the parallel runs of hyperspec
# here are the various functions
# - create the input data permutation for each job
# - send jobs to the cluster
# - monitor the queue for when jobs are complete and collect job output
# - keep a tally of the scores for each edge and monitor convergence

my %opts = getCommandLineOptions();

{
    my $outDir = File::Spec->rel2abs($opts{out});
    my $nTest = $opts{ntest};
    my $apmsFile = File::Spec->rel2abs($opts{apms});
    my $baseName = $opts{name};
    my $nSimulations = $opts{nsim};
    
    unless (-d $outDir) {
	make_path($outDir) or die "can't make directory $outDir";
    }
    my $dataDir = "$outDir/data";
    unless (-d $dataDir) {
	make_path($dataDir) or die "can't make directory $dataDir";
    }
    my $jsonDir = "$outDir/json";
    unless (-d $jsonDir) {
	make_path($jsonDir) or die "can't make directory $jsonDir";
    }
    my $scriptDir = "$outDir/qdir";
    unless (-d $scriptDir) {
	make_path($scriptDir) or die "can't make directory $scriptDir";
    }
    chdir($scriptDir);

    say "retrieve";
    my $apms = retrieve($apmsFile);

    say "make jobs";
    #my @jobs;
    # job[] = { script => $path, name => $name };
    my $seed = time();
    for my $jobN (1..$nTest) {
	my $name = sprintf("%s_%06d", $baseName, $jobN);
	my $script = makeJob($name, $jsonDir, $dataDir, $scriptDir, $apms);
	system("qsub -N $name $script");
	next unless $nSimulations > 0;
	for my $n (1..$nSimulations) {	    
	    my $simName = sprintf("%s.sim%03d", $name, $n);
	    my $script = makeJob($name, $jsonDir, $dataDir, $scriptDir, $apms,
				 $simName, $seed++);
	    system("qsub -N $simName $script");
	}
    }
    exit;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	apms => '/home/glocke/DPiM/dpim4/all.combined.10-21-2015.rmDupes.sumIso.0_05.tscFilter.withBait.storable',
	ntest => 1000,
	name => 'HGConsensus',
	nsim => 0,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out out_dir < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "apms=s", "ntest=i", "name=s", "nsim=i");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{apms});

    return %opts;
}

# if $seed arg is defined, make a simulation
sub makeJob {
    my ($name, $outDir, $dataDir, $scriptDir, $apms, $simName, $seed) = @_;

    my $dataFile = "$dataDir/$name.apms.tsv";
    my @search_id;
    my $jsonFile = "$outDir/$name.json";
    my $scriptFile = "$scriptDir/$name.sh";
    if (defined $simName) {
	$jsonFile = "$outDir/$simName.json";
	$scriptFile = "$scriptDir/$simName.sh";
    } else {
	randomSubset($dataFile, $apms, \@search_id);
	my $listOut = "$outDir/$name.repl.list";
	open my $OUT, ">", $listOut or die "can't write to $listOut. $!";
	say $OUT $_ for @search_id;
	close $OUT;
    }

    my $hyper = '/home/glocke/DPiM/cpp/hyperspec';
    
    my $cmd = "$hyper --in $dataFile --out $jsonFile";
    if (defined $seed) {
	$cmd.= " --sim 1 --seed $seed";
    }
    $cmd.="\n";
    #$cmd.="gzip $dataFile";
    $cmd.="gzip $jsonFile";
    writeScript($scriptFile, $cmd);
    return $scriptFile;
}

# randomly select a run for each bait
sub randomSubset {
    my ($outFile, $data, $search_id) = @_;

    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    
    # output format for each line:
    # "$search_id\t$bait_ref\t$prey_ref\t$total_peptides"
    for my $bait (sort keys %$data) {
	my @runs = keys %{ $data->{$bait} };
	my $sid = $runs[rand(0+ @runs)];
	push @$search_id, $sid;
	for my $row (@{ $data->{$bait}{$sid}}) {
	    say $OUT join "\t" , $sid, $bait, $row->{prey_ref}, 
		$row->{total_peptides};
	}
    }
    close $OUT;

    return;
}

sub writeScript {
    my ($outFile, $cmd) = @_;

    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";

    print $OUT '#!/bin/bash

#$ -l rt=0:10:00
#$ -V
#$ -cwd 
#$ -b y

module add gcc
module add gsl
module add boost

';
    print $OUT $cmd;
    close $OUT;

    system("chmod 755 $outFile");
    
    # let "SEED='.$baseSeed.'+$SGE_TASK_ID"
    # echo seed = $SEED > hyperspec.out.$SGE_TASK_ID.s

    
    return;

}

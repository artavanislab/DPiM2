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
#use DpimLib qw(getLineDP4APMS);

# make a network consisting only of experiments with a single replicate
# then for each bait that has multiple replicates, make a network adding each
#   replicate.  
# *do the comparisons elsewhere*
#   But basically, the idea is to make heatmaps with z = network difference

my %opts = getCommandLineOptions();

{
    my $outDir = $opts{out};
    my $apmsFile = $opts{apms};
    my $baseName = $opts{name};
    my $nSim = $opts{nsim};

    $outDir = File::Spec->rel2abs($outDir);
    
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
    
    my %nRep = map { $_ => 0+ keys %{ $apms->{$_} } } keys %$apms;
    my @oneRepl = grep {$nRep{$_} == 1} keys %nRep;
    my @twoRepl = grep {$nRep{$_} == 2} keys %nRep;
    my @multiRepl = grep {$nRep{$_} > 2} keys %nRep;
        
    reportNRep("$outDir/replicateTable.tsv", \%nRep, 
	       "# ".(0+ @oneRepl)." single-replicate baits; ".
	       (0+ @twoRepl)." double-replicate baits; ".
	       (0+ @multiRepl)." multi-replicate baits; ");
    
    my $stableDataFile = "$dataDir/stableData.tsv";
    writeStableData($stableDataFile, \@oneRepl, $apms);
    my $seed = time();
    
    say "stable";
    makeAndQ($baseName."_stable", $stableDataFile, $jsonDir, $scriptDir, $nSim);
    exit if exists $opts{stableonly};
    for my $bait (@twoRepl, @multiRepl) {
	for my $sid (keys %{$apms->{$bait}}) {
	    my $dataFile = "$dataDir/$bait"."_$sid.tsv";
	    writeOneRepl($dataFile, $stableDataFile, $bait, $sid, 
			 $apms->{$bait}{$sid});
	    
	    my $job = join "_", $baseName, $bait, $sid;
	    makeAndQ($job, $dataFile, $jsonDir, $scriptDir, $nSim);
	}
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   


sub getCommandLineOptions {

    my %defaults = (
	apms => '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.01-30-2016.storable',
	name => 'HGConsensus',
	nsim => 11,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out out_dir < $defaultString -stableonly >\n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "apms=s", "name=s", "nsim=i", "stableonly");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{apms});

    return %opts;
}

# make a table reporting the number of replicates for each bait
sub reportNRep {
    my ($outFile, $repHash, $comments) = @_;

    my @keys = sort {$repHash->{$b} <=> $repHash->{$a}} keys %$repHash;

    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";

    say $OUT $comments if $comments;
    say $OUT join "\t", qw(bait replicates);
    for my $k (@keys) {
	say $OUT join "\t", $k, $repHash->{$k};
    }
    close $OUT;
}

sub writeStableData {
    my ($outFile, $baits, $data) = @_;

    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    
    # output format for each line:
    # "$search_id\t$bait_ref\t$prey_ref\t$total_peptides"
    for my $bait (@$baits) {
	my @runs = keys %{ $data->{$bait} };
	my $sid = $runs[0];
	for my $row (@{ $data->{$bait}{$sid}}) {
	    say $OUT join "\t" , $sid, $bait, $row->{prey_ref}, 
		$row->{total_peptides};
	}
    }
    close $OUT;
    
    return;
}

sub writeOneRepl {
    my ($outFile, $stableFile, $bait, $sid, $repl) = @_;

    system("cp $stableFile $outFile");
    open my $OUT, ">>", $outFile or die "can't write to $outFile. $!";
    for my $row (@$repl) {
	say $OUT join "\t" , $sid, $bait, $row->{prey_ref}, 
	$row->{total_peptides};
    }
    close $OUT;
}

sub makeAndQ{
    my ($jobName, $dataFile, $jsonDir, $scriptDir, $nSim, $seed) = @_;

    # base job
    { 
	my $script = makeJob($jobName, $dataFile, $jsonDir, $scriptDir);
	my $cmd = "qsub -N $jobName $script";
	say $cmd;
	system($cmd);
    }

    # simulations
    for my $i (1..$nSim) {
	my $simJob = sprintf("$jobName.sim%03d", $i);
	my $script = makeJob($simJob, $dataFile, $jsonDir, $scriptDir, $seed++);
	my $cmd = "qsub -N $simJob $script";
	say $cmd;
	system($cmd);
    }
}


# if $seed arg is defined, make a simulation
sub makeJob {
    my ($name, $dataFile, $jsonDir, $scriptDir, $seed) = @_;

    my @search_id;
    my $jsonFile = "$jsonDir/$name.json";
    my $scriptFile = "$scriptDir/$name.sh";

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



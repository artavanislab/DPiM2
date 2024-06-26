#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);
use feature qw(:5.10);

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist);


my @modes = qw(all analyze report);
my %modes = map {$_ => 1} @modes;
my %opts = getCommandLineOptions();

my $dpim_file = $opts{in};
my $num_sims = $opts{nsim};
my $mode = $opts{mode};

#my $DPIMDIR = $ENV{DPIMDIR};
#my $HYPDIR = "$DPIMDIR/HyperSpec";
# perl compute_fdr.pl dpim2_nrtap.120123 10
my $hyperspec_cmd = "/home/glocke/DPiM/cpp/hyperspec";
my $length_corr_file = "/home/glocke/DPiM/interfly_fbgn_avgTryp.out";
my $callingDir = $ENV{PWD};
my $baseSeed = time();

my $outfile = $dpim_file.".sim".$num_sims.".run";

if ($mode eq 'all') {
    runHyperSpec();
}
if ($mode eq 'all' || $mode eq 'analyze') {
    parseHyperSpec();
}
if ($mode eq 'report') { 
    die "can't find '$outfile'" unless -e $outfile;
    # if the run file exists, only need to read in it to calculate the 
    # HGScore cutoff

    my ($hg_cutoff, $num_sims_new) = compHGScoreCutoff($outfile); 
    # overwrite the input num of simulations
    print "HGScore cutoff is $hg_cutoff (median of $num_sims_new simulation)\n";
}
exit;

sub runHyperSpec {
    my $out = "hyperspec.out.r"; 
    if (-s $out) {
	say "found previous run of hyperspec.  Using that...";
    } else {
	my $cmd = "$hyperspec_cmd $dpim_file $length_corr_file 0 $baseSeed > $out";
	say $cmd;
	system($cmd);
	die "Fatal hyperspec failure (missing modules?)" if -z $out;
    }
    exit;

    
    if (glob("hyperspec.out*s")) {
	# remove all simulation results produced previously	
	my $tempDir = tempdir();
	warn "moving previous hyperspec output to $tempDir";
	system("mv hyperspec.out*s $tempDir/.");
    }

    my $i=0;
    my $hsScript = 
	writeHyperSpecScript($callingDir, $num_sims, $hyperspec_cmd, $dpim_file,
			     $length_corr_file, $baseSeed);
    system("chmod 755 $hsScript");
    system("qsub -N compfdr $hsScript");
    
    
    while (1) { # check utill all simulation files are ready
	# want correct number of files, 
	# and these files should not be changing in size
	my @files = glob 'hyperspec.out.*.s';
	if (@files == $num_sims) {
	    my @size1;
	    my @size2;
	    my $size_flag = 0; # this would control the simulation output size to be reasonable to consider as done
	    foreach my $file (@files) {
		my $size1 = (stat($file))[7];
		if ($size1 < 1000000) { # average file size is 9241416
		    $size_flag = 1;
		    last;
		}
		push @size1, $size1;
	    }
	    next if $size_flag; # jump to next while loop if the size does not meet the requirement
	    sleep(2); #sleep for 2 seconds to allow the file size change
	    foreach my $file (@files) {
                my $size1 = (stat($file))[7];
		if ($size1 < 1000000) {
		    $size_flag = 1;
		    last;
		}
                push @size2, $size1;
            }
	    next if $size_flag; # jump to next while loop if the size does not meet the requirement
	    my $sums = 0;
	    foreach my $i (0..$#size1) {
		$sums += $size2[$i] - $size1[$i];
		if ($sums > 0) {
		    last;
		}
	    }
	    if($sums==0) {
		last;
	    }
	} else {
	    sleep(10);
	}
    } 
}

sub parseHyperSpec {

    my @real_edges = ();
    open(IN, "grep FBgn hyperspec.out.r | ");
    while (my $line = <IN>) {
	chomp($line);
	my @field = split(/\t/, $line);
	push @real_edges, $field[2];
    }
    close(IN);
    


    my @simu_FDR012345_cutoffs;

    for my $i (1..$num_sims) {
	my @simulated_edges = ();
	open(IN, "grep FBgn hyperspec.out.$i.s | ");
	while (my $line = <IN>) {
	    chomp($line);
	    my @field = split(/\t/, $line);
	    push @simulated_edges, $field[2];
	}
	close(IN);
	
	my @FDR012345_cutoff = ();
	my $FDR = 0;
	my $s = 0;
	my $r;
	while ($FDR < 0.05 && $s < $#simulated_edges) {
	    $r = 0;
	    die "real_edges[$r] doesn't exist" unless defined $real_edges[$r];
	    die "simulated_edges[$r] doesn't exist" unless defined $simulated_edges[$r];
	    while ($real_edges[$r] >= $simulated_edges[$s]) {
		$r ++;
	    }
	    $FDR = ($s+1)/($r+1);
	    
	    if ($FDR >= 0.01 && scalar(@FDR012345_cutoff) == 0) {
		push @FDR012345_cutoff, ($r+1);
		#push @FDR012345_cutoff, $real_edges[$r];
	    } elsif ($FDR >= 0.02 && scalar(@FDR012345_cutoff) < 2) {
		push @FDR012345_cutoff, ($r+1);
		#push @FDR012345_cutoff, $real_edges[$r];
	    } elsif ($FDR >= 0.03 && scalar(@FDR012345_cutoff) < 3) {
		push @FDR012345_cutoff, ($r+1);
		#push @FDR012345_cutoff, $real_edges[$r];
	    } elsif ($FDR >= 0.04 && scalar(@FDR012345_cutoff) < 4) {
		push @FDR012345_cutoff, ($r+1);
		#push @FDR012345_cutoff, $real_edges[$r];
	    }
	    #print $simulated_edges[$s]."\t".($s+1)."\t".($r+1)."\t".$FDR."\n";
	    $s ++;
	}
	push @FDR012345_cutoff, ($r+1);
	#push @FDR012345_cutoff, $real_edges[$r];
	push @simu_FDR012345_cutoffs, \@FDR012345_cutoff;
	print join("\t", @FDR012345_cutoff)."\tedges passing at 0.01-0.05 FDR\n";
    }
    
    open(OUT, ">$outfile");
    foreach my $i (0..$num_sims-1) {
	print OUT join("\t", @{$simu_FDR012345_cutoffs[$i]})."\tedges passing at 0.01-0.05 FDR\n";
    }
    close OUT;

    system("cp hyperspec.out.r $dpim_file.hyperspec"); # rename the output file

    my ($hg_cutoff, $num_sims) = compHGScoreCutoff($outfile); # identify the HGScore cutoff
    print "HGScore cutoff is $hg_cutoff (median of $num_sims simulation)\n";
}




sub getCommandLineOptions {
    my %defaults = (
	nsim => 101,
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;

    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = join("/", @arr);
    my $usage = "usage: $0 -in dpim.nrBait < $defaultString -mode $modeString>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "nsim=i", "mode=s");
    die $usage unless exists $opts{in};

    checkExist('f', $opts{in});

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    return %opts;
}





sub median {
    my @values = @_;
    my $median;
    my $mid = int @values/2;
    my @sorted_values = sort { $a <=> $b} @values;
    if (@values % 2) {
	$median = $sorted_values[ $mid ];
    } else {
	$median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
    }
    return $median;
}

sub compHGScoreCutoff {
    my ($file, $fdr) = @_;
    if (! defined $fdr) {
	$fdr = 0.05; # by default look for cutoff for 0.05 FDR
    }
    $fdr = 100*$fdr - 1; #conver the FDR to column index from the file read in 

    # readin the output .run file to get a list of line numbers (cutoff) from each of the simulation runs
    my @line_cutoffs;
    open my $FILE, "<", $file or die "can't read from $file. $!";
    my $num = 0;
    foreach my $line (<$FILE>){
	chomp $line;
	my @cont = split("\t", $line);
	my $line_cutoff = $cont[$fdr];
	push @line_cutoffs, $line_cutoff;
	$num++;
    }
    close $FILE;
    
    # build a hash to store number of times a line is picked
    my %line_cutoffs;
    $line_cutoffs{$_}++ for @line_cutoffs;

    my $inFile = "$dpim_file.hyperspec";
    open my $IN, "<", $inFile or die "Can't read from $inFile. $!";
    my $count = 1;
    my @res;
    while (my $line = readHS($IN)) {
	if(exists $line_cutoffs{$count}) {
	    my $i = 0;
	    while ($i<$line_cutoffs{$count}) {
		my @cont = split("\t",$line);
		push @res, $cont[2];
		$i++;
	    }
	}
	$count++;
    }
    my $HGScoreCutoff = median(@res);
    return ($HGScoreCutoff, $#res+1);
}

sub readHS {
    my ($IN) = @_;

    return undef if eof($IN);

    my $goodLine = sub {
	return $_[0] =~ /^FBgn\d+\s+FBgn/;
    };
    
    my $line;
    do {
	$line = <$IN>;
    } while (!eof($IN) && ! $goodLine->($line));
    return undef if ! $goodLine->($line);

    $_ = $line;
    chomp;
    return $line;
}


# write a bash script that calls hyperspec X times
sub writeHyperSpecScript {
    my ($dir, $numSims, $hyperspec_cmd, $dpim_file, $length_corr_file, 
	$baseSeed) = @_;

    my ($OUT, $shFile) = tempfile(DIR=>$dir, SUFFIX=>".sh");
    print $OUT '#!/bin/bash

#$ -l rt=0:40:00
#$ -V
#$ -cwd 
#$ -b y
#$ -t 1:101

let "SEED='.$baseSeed.'+$SGE_TASK_ID"
echo seed = $SEED > hyperspec.out.$SGE_TASK_ID.s

module add gcc
module add gsl
';

    print $OUT "$hyperspec_cmd $dpim_file $length_corr_file 1 \$SEED >> hyperspec.out.\$SGE_TASK_ID.s\n";
    close $OUT;
    return $shFile;
}

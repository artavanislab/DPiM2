#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Path qw(make_path);
use Math::GSL::CDF qw(gsl_cdf_binomial_P);
use HomeBrew::IO qw(checkExist readCols);
use DpimLib qw(getLineBiopAPMS);

# remove proteins that are over-abundant on a per-plate basis

# algorithm
# 1. find the fraction of experiments in which every protein appears
# 2. find the number of experiments in which every protein appears on each plate
# 3. do binomial tests to see which proteins are significantly overabundant
# 4. do FDR correction
# 5. for every protein with q < 0.01, find its average tsc on the plate
# 5.1. omit proteins with fewer than $minAppear=10 appearances on the plate
# 6. remove it from experiments on this plate where it had less than this average


my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $minAppear = $opts{minappear};
    my $rScript = $opts{rscript};
    my $workingDir = $opts{workingdir};
    
    if (! -d $workingDir) {
	make_path($workingDir) or die "can't make $workingDir";
    }
    my $baseOut = "$workingDir/overAbundant";
    if (! exists $opts{skipr}) {
	my $cmd = "Rscript $rScript $in $baseOut";
	say $cmd;
	say "This should take about half an hour";
	system($cmd);
    }
    my @testFiles = `ls $baseOut*`;
    die "Found ".(0+ @testFiles)." files like $baseOut*. Quitting." 
	unless 4 < @testFiles; ## expect testFiles==163
    chomp @testFiles;

    my %offenders; # offenders{plate}{prey} = average
    findOffenders(\%offenders, \@testFiles, $minAppear);

    my @cols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date    ms_inst_run_id   rep     plate);
    open my $IN, "<", $in or die "can't read $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 removed over-abundant prey on a per-plate basis";
    say $OUT "# minAppear = $minAppear";
    say $OUT join "\t", @cols;
    my %row;
    while (getLineBiopAPMS(\%row, $IN, 'line')) {
	my ($prey, $plate, $tsc) = ($row{prey_ref}, $row{plate}, 
				    $row{total_peptides});
	next if exists $offenders{$plate}{$prey} && 
	    $tsc < $offenders{$plate}{$prey};
	print $OUT $row{line};
    }
    close $IN; 
    close $OUT; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	minappear=>10,
	rscript => '/home/glocke/DPiM/scripts/human/perPlate.R',
	workingdir => "$ENV{PWD}/perPlateFiles",
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -skipr >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "minappear=i", "rscript=s", 
	       "workingdir=s", "skipr");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{rscript});

    return %opts;
}

# ret{$plate}{$search_id}{$prey} = {tsc = 1234, line = $line}}
# where $prey is a prey protein ID
sub parseAPMS {
    my ($ret, $in) = @_;

    open my $IN, "<", $in or die "can't read $in. $!";
    my %row;
    while (getLineBiopAPMS(\%row, $IN, 'line')) {
	$ret->{$row{plate}}{$row{search_id}}{$row{prey_ref}} = {
	    tsc => $row{total_peptides},
	    line => $row{line}
	};
	#die Dumper($ret);
    }

    return;
}

#  $ret{$prey} = fraction of all experiments including $prey
sub globalStats {
    my ($ret, $apms) = @_;

    my $nExpts = 0;
    for my $plateHash (values %$apms) {
	for my $expt (values %$plateHash) {
	    $nExpts++;
	    for my $prey (keys %$expt) {
		$ret->{$prey}++;
	    }
	}
    }

    for my $prey (keys %$ret) {
	$ret->{$prey} /= $nExpts;
    }

    return;
}

if (0) {
    ## making sure we understand the binomial test
    
    ## double gsl_cdf_binomial_P (unsigned int k, double p, unsigned int n)
    ## p(k) = {n! \over k! (n-k)! } p^k (1-p)^{n-k}
    ## chance of exactly k successes in n bernoulli trials with probability p

    my $trials = 10;
    my $prob = 0.4;
    my $found = $trials-1;
    my $pVal = 1 - gsl_cdf_binomial_P($found-1, $prob, $trials);
    say "$pVal = 1 - gsl_cdf_binomial_P($found-1, $prob, $trials)";
    exit;
}

# %ret{plate}{prey} = average tsc on plate
sub findOffenders {
    my ($ret, $testFiles, $minAppear) = @_;

    my @cols = qw(gene appearHere avg);
    for my $f (@$testFiles) {
	my @spl = split /\./, $f;
	my $plate = pop (@spl);
	next if $plate eq 'exptsPerPlate';
	my @read = readCols($f, \@cols);
	for my $row (@read) {
	    next if $row->{appearHere} < $minAppear;
	    $ret->{$plate}{$row->{gene}} = $row->{avg};
	}
    }

    return;
}

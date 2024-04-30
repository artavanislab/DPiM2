#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4_1APMS);

# read in dp4.1 formatted apms data (contains logp)
# remove rows with the lowest scoring peptide ids until the ratio of
# reverse peptides to non-reverse peptides

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $targetFDR = $opts{fdr};
    
    my (%logp, %preyHash); # logp{$x}{$fbgn}=1 if $fbgn is found with logp=$x
    #                  # prey{$fbgn}{$x}=1 if $fbgn is found with logp=$x
    my ($bait) = readAPMS(\%logp, \%preyHash, $in, $opts{ignorebait});
    # bait = number of unique baits
    my $reverse = 0+ grep /reverse/, keys %preyHash;
    my $total = 0+ keys %preyHash;
    my ($initReverse, $initTotal) = ($reverse, $total);

    my $fdr = calcFDR($reverse, $total);
    my $logPBound; # require logP strictly greater than ~
    
    for my $p (sort {$a <=> $b} keys %logp) {
	$logPBound = $p;
	for my $prey (keys %{ $logp{$p} }) {
	    $preyHash{$prey}{$p}--;
	    if ($preyHash{$prey}{$p} == 0) {
		delete $preyHash{$prey}{$p};
	    } elsif ($preyHash{$prey}{$p} < 0) {
		die "found preyHash{$prey}{$p} = $preyHash{$prey}{$p}";
	    }
	    if (0 == keys %{ $preyHash{$prey} }) {
		$total--;
		if ($prey =~ /reverse/) {
		    $reverse--;
		}
	    }
	}
	$fdr = calcFDR($reverse, $total);
	last if $fdr <= $targetFDR;
    }
    printf "found fdr %.3e at logp %.3e\n", $fdr, $logPBound;
    say "\tinitial reverse, total = $initReverse, $initTotal";
    say "\tfinal reverse, total = $reverse, $total";

    sanitizeAPMS($in, $out, $logPBound);
    
    
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	fdr => 0.01,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -ignorebait ".
	">\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "fdr=f", "ignorebait");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

sub readAPMS {
    my ($logp, $prey, $inFile, $ignoreBait) = @_;

    my %bait = ();
    my $unknown = 'FBgn0000000';
    
    open my $IN, "<", $inFile or die "can't read $inFile. $!";
    my %row;
    while(getLineDP4_1APMS(\%row, $IN)) {
	next if (!$ignoreBait) && $row{bait_ref} eq $unknown;
	$bait{$row{bait_ref}}=1;
	my $p = $row{logp};
	my $pry = $row{prey_ref};
	$logp->{$p}{$pry}++;
	$prey->{$pry}{$p}++;
    }
    close $IN; 

    return 0+ keys %bait;
}

sub calcFDR {
    my ($reverse, $total) = @_;

    return $reverse/($total - $reverse);
}

sub sanitizeAPMS {
    my ($inFile, $outFile, $logPBound) = @_;

    my @cols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date     ms_inst_run_id);

    open my $OUT, ">", $outFile or die "can't write to $outFile. $!";
    say $OUT "# $0 filtered $inFile";
    say $OUT "# set open logPBound at $logPBound";
    say $OUT join "\t", @cols;
    
    open my $IN, "<", $inFile or die "can't read $inFile. $!";
    my %row;
    while(getLineDP4_1APMS(\%row, $IN)) {
	next if $row{logp} <= $logPBound;
	say $OUT join "\t", map { $row{$_} } @cols;
    }    
    close $OUT; 
    close $IN; 
}

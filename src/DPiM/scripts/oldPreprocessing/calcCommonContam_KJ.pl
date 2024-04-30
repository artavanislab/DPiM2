#!/usr/bin/perl -w

#use strict;   EH?? no STRICT????  Is that why I smell brimstone?

sub process_cmd_args();
process_cmd_args();
$input_file = $ARGV[0];
$output_file = $ARGV[1];

open INPUT, $input_file or die "Cannot open $input_file\n";
$buf = <INPUT>; # Strip the header

# Set up an enum for later reference
$EXPERIMENTS = 0;
$TOTALPEPSUM = 1;
$UNIQUESUM = 2;

my %peptideMap;
my %totalExpts;
# Chew through the file, accumulating information as we go
while($buf = <INPUT>){

    chomp $buf;
    @tokens = split("\t", $buf);
    my $peptide = $tokens[8];
    my $experiment = $tokens[7];
    my $total_peptides = eval($tokens[5]);
    my $unique_peptides = eval($tokens[6]);
    
    # Keep a track of total Experiments
    if( ! exists $totalExpts{$experiment} ){
	$totalExpts{$experiment} = '';
    }

    #if( exists $peptideMap{$peptide} ){
	if( defined $peptideMap{$peptide} ){
	my @results = @ { $peptideMap{$peptide} };
	$results[$TOTALPEPSUM] = ($total_peptides + $results[$TOTALPEPSUM]);
	$results[$UNIQUESUM] = ($unique_peptides + $results[$UNIQUESUM]);
	my %experiments = %{ $results[$EXPERIMENTS] };
	if (exists $experiments{$experiment}){
	    ;
	}
	else{
	    $experiments{$experiment} = "";
	}
	$results[$EXPERIMENTS] = \%experiments;
	$peptideMap{$peptide} = \@results;
    }
    else{

	my @results;
	$results[$TOTALPEPSUM] = $total_peptides;

	my %experiments;
	$experiments{$experiment} = "";

	$results[$EXPERIMENTS] = \%experiments;
	$results[$UNIQUESUM] = $unique_peptides;
	$peptideMap{$peptide} = \@results;
    }
}
close INPUT;
print("Finished reading file\n");

$totalExpts = scalar keys %totalExpts;
open OUTPUT, ">$output_file";
# Now summarize the information collected to produce the output file
print "Total experiments seen was: $totalExpts\n";
print OUTPUT "fbgn\tCount\tFraction\tavg_tot_pep\n";
my @out;
foreach $peptide(keys(%peptideMap)){
    @summary = @{ $peptideMap{$peptide} };
    %experiments = %{ $summary[$EXPERIMENTS] };
    @expts = keys %experiments;
    $numExpts = scalar @expts;

    #print "$peptide\t$totalExpts\t$numExpts\t$summary[$TOTALPEPSUM]\t$summary[$UNIQUESUM]\n";
    $fraction = eval($numExpts / $totalExpts);
    $avg_tot_pep = eval(($summary[$TOTALPEPSUM]) / $numExpts);
    #print OUTPUT "$peptide\t$numExpts\t$fraction\t$avg_tot_pep\n";
	push @out, [$peptide,$numExpts,$fraction,$avg_tot_pep];
}

# sort the output by Count
@sorted = sort { $b->[1] <=> $a->[1] } @out;

foreach my $o (@sorted) {
	print OUTPUT join("\t",@$o)."\n";
}
close(OUTPUT);

sub process_cmd_args(){
    if($#ARGV != 1){
	print "Usage: calcCommonContam.pl [dpim file] [output file]\n";
	die;
    }
}


#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);

my %opts = getCommandLineOptions();
my $input_file = $opts{in};
my $output_file = $opts{out};

open my $IN, "<", $input_file or die "Cannot open $input_file. $!";
<$IN>; # Strip the header

my %peptideMap; # peptideMap{$fbgn} = { totalPepSum = TSC,        emacs wants->}
#	{	#		        expt => { $search_id=> 1, ...} }
my %totalExpts; # totalExpts{$search_id} = 1 for each sID

my $reader = \&getLineDP4APMS;
$reader = \&getLineBiopAPMS if $opts{mode} eq 'human';

# Chew through the file, accumulating information as we go
my %row;
while($reader->(\%row, $IN)){
    my $peptide = $row{prey_ref};
    my $experiment = $row{bait_ref};
    my $total_peptides = $row{total_peptides};
    #my $unique_peptides = eval($tokens[6]);
    
    # Keep a track of total Experiments
    $totalExpts{$experiment} = 1;

    $peptideMap{$peptide}{totalPepSum} += $total_peptides;
    $peptideMap{$peptide}{totalPepSqr} += $total_peptides * $total_peptides;
    $peptideMap{$peptide}{expt}{$experiment} = 1;
    #$peptideMap{$peptide}{uniqSum} += $unique_peptides;
}
close $IN;

print("Finished reading file\n");

my $totalExpts = scalar keys %totalExpts;
print "Total experiments seen was: $totalExpts\n";

# Now summarize the information collected to produce the output file
my @out; # out[i] = [ $fbgn, $numExpts, $numExpts/$totalExpts, $avg_tot_pep ]
foreach my $peptide (keys(%peptideMap)){
    my %summary = %{ $peptideMap{$peptide} };
    my $numExpts = 0+ keys %{ $summary{expt} };

    my $fraction = $numExpts / $totalExpts;
    my $avg_tot_pep = $summary{totalPepSum} / $numExpts;
    push @out, [$peptide, $numExpts, $fraction, $avg_tot_pep];
    #die Dumper($peptide, $numExpts, \%summary) 
    #	if 0 > ($summary{totalPepSqr}/$numExpts - $avg_tot_pep**2);
    #my $sd_tot_pep = sqrt($summary{totalPepSqr}/$numExpts - $avg_tot_pep**2);
    #push @out, [$peptide, $numExpts, $fraction, $avg_tot_pep, $sd_tot_pep];
    #print $OUTPUT "$peptide\t$numExpts\t$fraction\t$avg_tot_pep\n";
}

# sort the output by Count
@out = sort { $a->[0] cmp $b->[0] } @out;
my @sorted = sort { $b->[1] <=> $a->[1] } @out;

open my $OUT, ">", $output_file or die "Can't write to $output_file. $!";

say $OUT "# totalExpts = $totalExpts";
#say $OUT join "\t", qw(fbgn Count Fraction avg_tot_pep sd_tot_pep);
# 01/15/2018 KJ adding human opt
if ($opts{mode} eq 'human') {
	say $OUT join "\t", qw(entrez Count Fraction avg_tot_pep);
} else {
	say $OUT join "\t", qw(fbgn Count Fraction avg_tot_pep);
}
foreach my $o (@sorted) {
    say $OUT join("\t",@$o);
}
close($OUT);

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    my @modes = qw(fly human);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: "
	, join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};
    
	    

    checkExist('f', $opts{in});

    return %opts;
}


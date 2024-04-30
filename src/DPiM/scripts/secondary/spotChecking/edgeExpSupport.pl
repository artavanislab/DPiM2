#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineHyperspecAPMS readHS);

# describe the purpose for this script

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $apmsFile = $opts{nrb};
    my $out = $opts{out};

    # apms{$prey}{$sid} = 1 iff this prey appeared in this experiment
    my %apms;
    say "read apms";
    my %reader = (
	"4cols" => \&getLineHyperspecAPMS,
	"6cols" => \&getLineDP4APMS,
    );
    parseAPMS(\%apms, $apmsFile, \&getLineHyperspecAPMS);


    say "find support";

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 added a support column to $in based on $apmsFile";
    say $OUT "# \"support\" tells you the number of experiments in which this pair of proteins appears";
    say $OUT join "\t", qw(protein1 protein2 score support);
    
    open my $IN, "<", $in or die "can't read $in. $!";
    while (my $line = readHS($IN)) {
	$_ = $line;
	my ($prot1, $prot2) = split;
	chomp $line;
	say $OUT "$line\t".support($prot1, $prot2, \%apms);
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
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in in.net -nrb nrBait -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "nrb=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{nrb} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{nrb});

    return %opts;
}


# read a whole APMS file into a hash
# ret{$prey}{$sid} = 1 iff this prey appeared in this experiment
sub parseAPMS {
    my ($ret, $inFile, $reader, %options) = @_;

    open my $IN, "<", $inFile or die "Can't read from $inFile. $!";
    my %row;
    while (getLineHyperspecAPMS(\%row, $IN)) {
	$ret->{$row{prey_ref}}{$row{search_id}} = 1
    }

    return;
}


# find the number of experiments in which both proteins appear
sub support {
    my ($prot1, $prot2, $apms) = @_;

    my $ret=0;
    for my $sid (keys %{ $apms->{$prot1} }) {
	$ret++ if exists $apms->{$prot2}{$sid};
    }

    return $ret;
}

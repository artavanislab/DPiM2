#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4APMS);

# remove prey that appear in only one experiment with a single peptide

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $nHit = $opts{nhit};

    say "find one-hit wonders";
    my %oneHit = findOHW($in, $nHit);

    say "\tnow filter them out";
    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# removed one-hit-wonders from $in";
    say $OUT join "\t", qw(search_id       bait_ref        prey_ref        total_peptides  sample_date    ms_inst_run_id);
    my %row;
    while(getLineDP4APMS(\%row, $IN, 'line')) {
	print $OUT $row{line} unless exists $oneHit{$row{prey_ref}};
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	nhit=>1,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "nhit=i");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# ret{$fbgn} = 1 iff $fbgn identifies a one-hit-wonder prey
sub findOHW {
    my ($in, $nHit) = @_;

    if (defined $nHit && $nHit!=1) {
	die "nhit != 1 is not implemented";
    }
    
    # hash{prey_ref} = number of appearances
    my (%singles, %multis);

    open my $IN, "<", $in or die "Can't read from $in. $!";

    my %row;
    while(getLineDP4APMS(\%row, $IN)) {
	my $prey = $row{prey_ref};
	next if $prey eq $row{bait_ref};
	my $tsc = $row{total_peptides};
	if ($tsc == 1) {
	    $singles{$prey}++;
	} else {
	    $multis{$prey}++;
	}
    }

    my %ret;
    for my $fb (keys %singles) {
	next if $singles{$fb} > 1;
	next if exists $multis{$fb};
	$ret{$fb} = 1;
    }
    return %ret;
}

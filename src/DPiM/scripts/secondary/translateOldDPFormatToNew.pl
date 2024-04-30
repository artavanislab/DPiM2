#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineRawAPMS getLineAPMS getLine2015DP3APMS );

# translate an apms data file from a previous format into the current format

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    
    my %readModes = (
	raw => \&getLineRawAPMS,
	dp2 => \&getLineAPMS,
	dp3 => \&getLine2015DP3APMS
    );
    my $reader = $readModes{$mode} // die "don't understand mode '$mode'";
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";

    my @cols = qw(search_id bait_ref prey_ref total_peptides sample_date ms_inst_run_id );
    say $OUT join "\t", @cols;
    my %row;
    while($reader->(\%row, $IN, 'line')) {
	#die Dumper(\%row);
	say $OUT join "\t", map { $row{$_} } @cols;
    }
    close $IN;
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    # mode:
    my @modes = qw(raw dp2 dp3);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;
    my $usage = "usage: $0 -in input -out output < $modeString>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    return %opts;
}


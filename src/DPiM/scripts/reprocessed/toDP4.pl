#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCols);
use DpimLib qw(getLineMarchAPMS getLineMayAPMS);

# convert reprocessed apms files from Julian to dp4 format

# do not attempt to identify the correct bait.  just report dummy baits.

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    
    my @cols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date     ms_inst_run_id);

    my $reader;
    if ($mode eq 'may') {
	$reader = \&getLineMayAPMS;
	push @cols, 'logp';
    } elsif ($mode eq 'march') {
	$reader = \&getLineMayAPMS;
    } else {
	die "unknown mode '$mode'\n";
    }
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";

    say $OUT "# converted $in to dpim4 format";
    say $OUT join "\t", @cols;

    my %row;
    my $dummyNum = 1;
    my %dummyBaits;
    while ($reader->(\%row, $IN)) {
	#die Dumper(\%row) if $row{search_id} eq '243487';
	if ($row{contam}) {
	    next unless $opts{yescontam};
	    die "yescontam not yet implemented; must parse prey_string, basically";
	}
	next if $row{prey_ref} eq 'TAG' && ! exists $opts{yestag};
	next if $row{prey_ref} eq 'FLAG-HA-tag' && ! exists $opts{yesflag};
	if ($row{reverse}) {
	    next if $opts{noreverse};
	    $row{prey_ref} = 'reverse_'.$row{prey_ref};
	}
	
	## convert date from MMDDyy to YYYY-MM-DD
	if (length($row{sample_date})==5) {
	    # in a few cases, the date string will omit the leading '0' for
	    # months prior to october.  Grr...
	    $row{sample_date} = '0'.$row{sample_date};
	}
	my @spl = $row{sample_date} =~ /../g;
	$spl[-1] = "20".$spl[-1];
	unshift @spl, pop @spl;
	$row{sample_date} = join "-", @spl;

	$row{bait_ref} = 'FBgn'.$row{search_id}.'_dummy';
	say $OUT join "\t", map { 
	    $row{$_} // die "undefined row{$_}".Dumper(\%row) 
	} @cols;
	
    }
    
    close $IN;
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(may march);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString -noreverse ".
	"-yescontam -yesflag -yestag>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "noreverse", "yescontam", 
	       "yesflag", "yestag");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    checkExist('f', $opts{in});

    return %opts;
}

# hash mapping from run id to (updated) fbgn
# ret{ms_inst_run_id} = fbgn
sub readBaitMap {
    my ($baitMapFile) = @_;

    my @cols = qw(ms_inst_run_id  newfbgn);
    my @read = readCols($baitMapFile, \@cols);

    return map { $_->{ms_inst_run_id} => $_->{newfbgn} } @read;
}

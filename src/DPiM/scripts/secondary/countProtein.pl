#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist writeCols);
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);

# count the appearances of each bait/prey

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $countTSC = $opts{scalebytsc};

    my %count;
    my $counted;
    my $col;
    if ($mode eq 'net') {
	%count = countNet($in);
	$counted = "node degree";
	$col = "all FBgns";
    } elsif ($mode eq 'rawfbgn') {
	%count = countRawFBgn($in);
	$counted = "all FBgns";
	$col = "every line";
    } elsif ($mode eq 'baitprey') {
	$col = "bait_ref and prey_ref";
	%count = countAppearances($in, 'prey_ref');
	my %bait = countAppearances($in, 'bait_ref');
	$count{$_}+= $bait{$_} for keys %bait;
	$counted = 'appearances';	
    } else {
	$col = ($mode eq 'bait')?'bait_ref':'prey_ref';
	if ($countTSC) {
	    %count = countTSC($in, $col);
	    $counted = 'tsc';
	} else {
	    %count = countAppearances($in, $col);
	    $counted = 'appearances';
	}
    }

    my @fbgn = sort { $count{$b} <=> $count{$a} } keys %count;
    my @count = map { $count{$_} } @fbgn;
    my $preComments = "# counted $counted for $col in $in";
    my $header = join "\t", qw(protein count);
    my $format = join "\t", qw(%s %d);
    writeCols($out, [\@fbgn, \@count], $header, $format, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    my %defaults = (
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;

    my @modes = qw(bait prey baitprey rawfbgn net);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);
    my $usage = "usage: $0 -in input -out output < $modeString -scalebytsc ".
	"-human >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "scalebytsc", "human");
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

sub countNet {
    my ($in) = @_;

    open my $IN, "<", $in or die "Can't read from $in. $!";

    my %ret;
    while (<$IN>) {
	next unless /^FBgn\d/;
	my @a = split;
	$ret{$a[0]}++;
	$ret{$a[1]}++;
    }

    return %ret;
}

## seek all "FBgn1234" and all "reverse_FBgn5678" (don't double count these)
sub countRawFBgn {
    my ($in) = @_;

    open my $IN, "<", $in or die "Can't read from $in. $!";

    my %ret;
    while (my $line = <$IN>) {
	while ($line =~ /(reverse_FBgn\d+)/g) {
	    $ret{$1}++;
	}
	$line =~ s/reverse_FBgn\d+//g;
	while ($line =~ /(FBgn\d++)/g) {
	    $ret{$1}++;
	}
    }

    return %ret;
}

sub countTSC {
    my ($in, $col) = @_;

    my %ret; # ret{$fbgn} = $count

    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    my $reader = \&getLineDP4APMS;
    $reader = \&getLineBiopAPMS if defined $opts{human};
    while ($reader->(\%row, $IN)) {
	$ret{$row{$col}} += $row{total_peptides};
    }

    return %ret;
}



sub countAppearances {
    my ($in, $col) = @_;

    my %mediate; # ret{$fbgn} = $count

    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    my $reader = \&getLineDP4APMS;
    $reader = \&getLineBiopAPMS if defined $opts{human};
    while ($reader->(\%row, $IN)) {
	$mediate{$row{$col}}{$row{search_id}} = 1;
    }

    my %ret;
    for my $prot (keys %mediate) {
	$ret{$prot} = 0+ keys %{$mediate{$prot}};
    }
    
    return %ret;
}



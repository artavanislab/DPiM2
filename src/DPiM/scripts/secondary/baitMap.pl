#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCols);
use DpimLib qw(getLineDP4APMS);

# look at existing (pre-julian remap) dpim apms data
# extract a table mapping from ms_inst_run_id to fbgn and gene symbol

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $updateFile = $opts{updatemap};
    my $f2sFile = $opts{fbgn2symb};
    my $oldf2sFile = $opts{oldfbgn2symb};
    my $ambigOut = $opts{ambig};
    
    open my $IN, "<", $in or die "Can't read from $in. $!";

    my %naiveMap; # map{$ms_inst_run_id} = $bait_ref
    my %row;
    say "parse";
    while (getLineDP4APMS(\%row, $IN)) {
	$naiveMap{$row{ms_inst_run_id}} = $row{bait_ref};
    }

    say "translate";
    
    my %f2s = getF2s($f2sFile);
    my %oldf2s;
    if ($oldf2sFile =~ /fbgn_annotation_ID_fb/) {
	flybaseF2S(\%oldf2s, $oldf2sFile);
    } else {
	%oldf2s = getF2s($oldf2sFile);
    }
    
    my %multiUpdate = multiRefMap($updateFile); 
    # multimap{ms_inst_run_id} = [updatedbait1, updatedbait2, ...]
    my %multiMap = map { $_ => $multiUpdate{$naiveMap{$_}} } keys %naiveMap;

    my %bestMap;
    {
	my @updated = grep { defined $multiMap{$_} } keys %multiMap;
	my @notUpdated = grep { ! defined $multiMap{$_} } keys %multiMap;
	$bestMap{$_} = findBestMap($_, $naiveMap{$_}, $multiMap{$_}, \%f2s, 
				   \%oldf2s) for @updated;
	$bestMap{$_} = $naiveMap{$_} for @notUpdated;
    }
    
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# updated $in according to $updateFile, $f2sFile, and $oldf2sFile";
    say $OUT join "\t", qw(ms_inst_run_id oldfbgn newfbgn symb oldsymb);
    for my $run (sort keys %bestMap) {
	my $naiveF = $naiveMap{$run};
	my $newF = $bestMap{$run};
	say $OUT join "\t", $run, $naiveF, $newF, $f2s{$newF} // 'notinfile'
	    , $oldf2s{$naiveF} // 'notinfile';
    }
    close $OUT;
    
    if (defined $ambigOut) {
	my @k = grep { defined $multiMap{$_} } keys %multiMap;
	@k = grep { @{ $multiMap{$_} } > 1 } @k;
	if (@k == 0) {
	    warn "no ambiguous results!!\n";
	    next;
	}
	
	open my $OUT, ">", $ambigOut or die "Can't write to $ambigOut. $!";
	say $OUT "# ambigous maps elided in $out";
	say $OUT join "\t", qw(ms_inst_run_id oldfbgn newfbgn symb oldsymb);
	for my $run (sort @k) {
	    my $naiveF = $naiveMap{$run};
	    for my $newF (@{ $multiMap{$run} }) {
		say $OUT join "\t", $run, $naiveF, $newF
		    , $f2s{$newF} // 'notinfile'
		    , $oldf2s{$naiveF} // 'notinfile';
	    }
	    say $OUT "" if exists $opts{addspace};
	}
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	updatemap => '/home/glocke/DPiM/prevDPIM/FlyBase_IDs_FB2016_01.txt',
	fbgn2symb => '/home/glocke/DPiM/fbgn_id2_2col_03-17-2016.txt',
	oldfbgn2symb => '/home/glocke/DPiM/fbgn_id2_2col.txt'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -ambig ".
	"write.list -addspace> \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "updatemap=s", "fbgn2symb=s", 
	       "oldfbgn2symb=s", "ambig=s", "addspace");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{fbgn2symb});
    checkExist('f', $opts{oldfbgn2symb});

    return %opts;
}


# ret{old_fbgn} = new_fbgn
sub singleRefMap {
    my ($inFile) = @_;

    my @read = readCols($inFile, [qw(submitted_id current_id converted_id)]);
    my %ret;
    for my $row (@read) {
	die "$row->{current_id} doesn't match $row->{converted_id}"
	    unless $row->{current_id} eq $row->{converted_id};
	next if $row->{submitted_id} eq $row->{current_id};
	if (exists $ret{$row->{submitted_id}}) {
	    $ret{$row->{submitted_id}} = $row->{current_id}
	        if $ret{$row->{submitted_id}} =~ m/^CG\d+$/;
	} else {
	    $ret{$row->{submitted_id}} = $row->{current_id};
	}
    }

    return %ret;
}

# ret{old_fbgn} = [all matches];
sub multiRefMap {
    my ($inFile) = @_;

    my @read = readCols($inFile, [qw(submitted_id current_id converted_id)]);
    my %ret;
    for my $row (@read) {
	die "$row->{current_id} doesn't match $row->{converted_id}"
	    unless $row->{current_id} eq $row->{converted_id};
	next if $row->{submitted_id} eq $row->{current_id};
	$ret{$row->{submitted_id}} //= [];
	push @{$ret{$row->{submitted_id}}}, $row->{current_id};
    }

    return %ret;
}

sub getF2s {
    my ($in) = @_;

    my @read = readCols($in, [qw(fbgn symbol)]);
    return map {$_->{fbgn} => $_->{symbol}} @read;
}

sub flybaseF2S {
    my ($ret, $mapFile) = @_;

    open my $IN, "<", $mapFile or die "Can't read from $mapFile. $!";
    while (<$IN>) {
	next if /^#/;
	next if length($_) < 2;
	chomp;
	my @spl = split;
	$ret->{$spl[1]} = $spl[0];
	if (@spl > 2) {
	    my @spl2 = split ',', $spl[2];
	    $ret->{$_} = $spl[0] for @spl2;
	}
	#die Dumper($ret);
    }
    return;
}

sub findBestMap {
    my ($run, $prev, $possibles, $f2s, $oldf2s) = @_;

    if (1 == @$possibles) {
	return $possibles->[0];
    }

    my $oldSymbol = $oldf2s->{$prev};
    my @newSymbols = map { $f2s->{$_} } @$possibles;
    my @match = grep { $newSymbols[$_] eq $oldSymbol } 0..$#newSymbols;
    die "can't find good match for $run; $prev -> ", (join ", ", @$possibles)
	if @match != 1;
    return $possibles->[$match[0]];
}

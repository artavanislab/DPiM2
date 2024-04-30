#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Storable qw(retrieve);
use Data::Dumper;
use Statistics::Basic qw(correlation);
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineDP4APMS);

# find the experiment in one storable that corresponds most closely to 
#   experiments in the dpim1

# output:
# bait dp1ID dp4ID extraTSC missingTSC sumTSC extraProt missingProt totalProt cor
#  * dp1ID is assigned by me
#  * dp4ID is search_id in received data for the closest experiment (least sumTSC)
#  * extraTSC is the count of prey peptides in dp1 not present in dp4
#  * missingTSC is the count of prey peptides in dp4 not present in dp1
#  * sumTSC sum of the two
#  * extraProt is the count of prey in dp1 not present in dp4
#  * missingProt is the count of prey in dp4 not present in dp1
#  * totalProt is the number of unique preys identified in either dp4 or dp1
#  * cor is the correlation of TSC vectors (absent proteins set to zero)

my %opts = getCommandLineOptions();

{
    my $dp1File = $opts{dp1};
    my $dp4File = $opts{dp4};
    my $out = $opts{out};

    
    my %dp1;
    # dp1{bait}{search_id}{prey_ref} = total_peptides
    readDP1 (\%dp1, $dp1File);

    #die Dumper($dp1{FBgn0000039});
    
    my %dp4;
    {
	my $read = retrieve($dp4File);
	reconfigureAPMS(\%dp4, $read);	
    }

    my @diff;
    # diff[$i] = {bait, dp1ID, dp4ID, extraTSC, missingTSC, sumTSC, extraProt, 
    #             missingProt};
    for my $bait (sort keys %dp1) {
	for my $dp1ID (sort {$a <=> $b} keys %{ $dp1{$bait} }) {
	    #say "$bait\t$dp1ID";
	    my $ld = leastDiff($bait, $dp1ID, $dp1{$bait}{$dp1ID}, \%dp4);
	    if (! ref($ld)) {
		die "$bait $dp1ID is giving a wonky leastDiff\n", Dumper($ld);
	    }
	    push @diff, $ld;
	}
    }

    my @cols = qw(bait dp1ID dp4ID cor totalProt sumTSC extraTSC missingTSC extraProt missingProt );
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# searching for experiments in $dp4File that most closely match those in $dp1File";
    say $OUT join "\t", @cols;
    for my $row (@diff) {
	if (! ref($row) ) {
	    die Dumper('not ref', $row);
	}
	say $OUT join "\t", map { $row->{$_} } @cols;
    }
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

    my $usage = "usage: $0 -dp1 dp1.tsv -dp4 dp4.storable -out output\n";

    my %opts = ();
    GetOptions(\%opts, "dp1=s", "dp4=s", "out=s");
    die $usage unless exists $opts{dp1} && exists $opts{dp4} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{dp1});
    checkExist('f', $opts{dp4});

    return %opts;
}

sub readDP1 {
    my ($sort, $in) = @_;
        
    my @cols = qw(prey_ref total_peptides);
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my %row;
    while (getLineDP1ExcelAPMS(\%row, $IN)) {
	my $bait = $row{bait_ref};
	my $id = $row{search_id};
	my $prey = $row{prey_ref};
	my $tsc = $row{total_peptides};
	$sort->{$bait}{$id}{$prey} = $tsc;
    }

    return;
}

sub getLineDP1ExcelAPMS {
    my ($ret, $IN, $wholeLine) = @_;

    return undef if eof($IN);

    my @cols = qw( search_id bait_ref prey_ref total_peptides  );
    
    my @spl;
    my $line;
    do {{ # perl's next statement doesn't work inside a do block.
	$line = <$IN>;
	next if $line =~ /^#/;
	next if $line =~ /$cols[0]/;
	next if length($line) < 4;
	$_ = $line;
	chomp;
	@spl = split;
	}
    } while (!eof($IN) && @spl != @cols);
    return undef if @spl != @cols; # checks for eof

    $ret->{$cols[$_]} = $spl[$_] for 0..$#cols ;
    $ret->{$wholeLine} = $line if defined $wholeLine;
    return 1;
}

# in:  in {bait}{search_id} = [ {prey_ref=>'fbgn', 'total_peptides'=>n},...]
# ret: ret{bait}{search_id}{prey_ref} = total_peptides
sub reconfigureAPMS {
    my ($ret, $in) = @_;

    for my $bait (keys %$in) {
	for my $sid (keys %{ $in->{$bait} }) {
	    for my $row (@{ $in->{$bait}{$sid} }) {
		my $prey = $row->{prey_ref};
		my $tsc = $row->{total_peptides};
		$ret->{$bait}{$sid}{$prey} = $tsc;
	    }
	}
    }

    return;
}
    
# %ret = (bait, dp1ID, dp4ID, extraTSC, missingTSC, sumTSC, extraProt, 
#         missingProt, totalProt, cor)
sub leastDiff {
    my ($bait, $dp1ID, $dp1, $dp4) = @_;

    my @cols = qw(dp4ID extraTSC missingTSC sumTSC extraProt missingProt totalProt cor);
    my %ret = (bait => $bait, dp1ID => $dp1ID);

    if (!exists $dp4->{$bait}) {
	$ret{$_} = -1 for @cols;
	return \%ret;
    }

    my %diffs; # diffs{dp4SID} = (@cols => x);
    my $leastDiff = 1234567890;
    my $leastID;
    my $maxCor = 0;
    my $maxID;
    for my $dp4SID (keys %{$dp4->{$bait}}) {
	$diffs{$dp4SID} = thisDiff($dp1, $dp4->{$bait}{$dp4SID});
	$diffs{$dp4SID}{dp4SID} = $dp4SID;
	if ($diffs{$dp4SID}{sumTSC} < $leastDiff) {
	    $leastDiff = $diffs{$dp4SID}{sumTSC};
	    $leastID = $dp4SID;
	}
	if ($diffs{$dp4SID}{cor} > $maxCor) {
	    $maxCor = $diffs{$dp4SID}{cor};
	    $maxID = $dp4SID;
	}
    }

    warn Dumper(\%diffs, $bait, $dp1ID, $leastID) if $maxID ne $leastID;
    
    $ret{$_} = $diffs{$leastID}{$_} for @cols;
    $ret{dp4ID} = $leastID;
    #die Dumper(\%ret);
    return \%ret;
}

# ret = {extraTSC, missingTSC, sumTSC, extraProt,  missingProt, totalProt, cor}
sub thisDiff {
    my ($dp1, $dp4) = @_;

    my %ret = map { $_ => 0 } qw(extraTSC missingTSC sumTSC extraProt missingProt totalProt cor);
    for my $prey (keys %$dp1) {
	$ret{totalProt}++;
	if (! exists $dp4->{$prey}) {
	    $ret{extraTSC} += $dp1->{$prey};
	    $ret{extraProt}++;
	} else {
	    if ($dp1->{$prey} > $dp4->{$prey}) {
		$ret{extraTSC} += $dp1->{$prey} - $dp4->{$prey};
	    } else { 
		$ret{missingTSC} += $dp4->{$prey} - $dp1->{$prey};
	    }
	}
    }
    for my $prey (keys %$dp4) {
	if (! exists $dp1->{$prey}) {
	    $ret{totalProt}++;
	    $ret{missingTSC} += $dp4->{$prey};
	    $ret{missingProt}++;
	} # else we handled this prey in the previous loop
    }

    $ret{sumTSC} = $ret{extraTSC} + $ret{missingTSC};

    {
	my %tmpHash = map {$_ => 1} ((keys %$dp1), (keys %$dp4));
	my @allProt = keys %tmpHash;
	my @dp1Vec = map { $dp1->{$_} // 0 } @allProt;
	my @dp4Vec = map { $dp4->{$_} // 0 } @allProt;
	$ret{cor} = sprintf "%.3f", correlation(\@dp1Vec, \@dp4Vec);
    }
    
    return \%ret;
}

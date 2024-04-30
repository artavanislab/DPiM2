#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsHashRef readColsRef);
use DpimLib qw(getLineDP4APMS getLineDP4_1APMS);

# identify the bait for this run according to all the available information
#   sources of information include
#   * previous records from DPIM2-4
#   * extensive human inspection, primarily by Bob Obar

# 1. Identify the FBgn previously associated with this experiment
# 2. Identify potential updates for that FBgn (“targets”)
# 3. If it hasn’t been updated, great.  If it has only one target, great.  Stop.  Otherwise, I’ll have to pick one of multiple targets.  Proceed.
# 4. Check to see if there are peptides for any of the targets recovered in this experiment.
# 5. If exactly one target has peptides recovered, great.  Stop.  If more than one, pick the greatest.  Stop. Otherwise, proceed.
# 6. At this point I have almost no basis to prefer one target over another.  I’ll first prefer an FBgn whose associated gene symbol is not something like CG12345
#### 6.1 For example, the targets of FBgn0020439 are FBgn0266446, FBgn0266448, and FBgn0266451, which identify CG45076, CG45078, and fau, respectively.  I’ll prefer FBgn0266451/fau.
# 7. If I’m still picking from more than one target, I’ll just choose the alphanumerically first option.
# 8. Every time I update an FBgn, make a record of the change along with the criterion used to select it (equivalently, record which step I stopped on).

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $apmsOut = $opts{out};
    my $baitKeyOut = $opts{keyout};
    my $baitIDFile = $opts{baitid};
    my $flybaseFile = $opts{flybase};
    my $flFile = $opts{fltab};
    
    say "reading previous baits";
    # sid2Bait{$search_id} = FBgnXXX
    my (%sid2Bait, %rejectMap);
    {
	my @read;
	my @cols = qw(search_id bait_ref tap_id retain notes);
	readColsRef(\@read, $baitIDFile, \@cols);
	for my $row (@read) {
	    my $rid = $row->{search_id};
	    $sid2Bait{$rid} = $row->{bait_ref};
	    if ($row->{retain} eq 'no') {
		$rejectMap{$rid} = $row->{notes};
		$rejectMap{$rid} = 'humanCantAnnotate_-_'.$row->{notes}
		    if $row->{notes} eq 'DISAGREE';
		$rejectMap{$rid} = 'FL_-_'.$row->{notes} if 
		    $row->{tap_id} =~ /^FL/;
		$rejectMap{$rid} = 'FH0000_-_'.$row->{notes} if 
		    $row->{tap_id} eq 'FH0000';
	    }
	}
    }
    my %rejectCodes = map {$_ => 1} values %rejectMap;
    
    say "reading fbgn udate map";
    # updateFBgn{prev}{new} = new gene symbol
    my %updateFBgn;
    makeUpdateMap(\%updateFBgn, $flybaseFile);

    say "reading AP-MS data";
    my (%apms, %msID);
    readAPMS(\%apms, $in, \%msID);
    die "found no experiments - maybe you want to use -nologp?"
	if 0 == keys %apms;
    #die Dumper($apms{249877});

    my %updateCodes = (
	unknown => 'unknown', # don't recognize this run id
	nochange => 'nochange', # didn't update
	single => 'singleUpdate', # single possible update
	singlePep => 'singlePeptideMatch', # only one update matches peptides
	mostPep => 'mostPeptideMatches', # multiple peptide matches
	nonCG => 'nonCG', # the only symbol that's not CGXXXX 
	alphaCG => 'alphaGC', # alphanumerically first CG
	alphaNonCG => 'alphaNonCG', # alphanumerically first not CG
	);

    # newBait{$search_id} = FBgnXXX
    # *this is the work product of this script*
    my %newBait;
    my $dummyFBgn = 'FBgn0000000';

    say "updateFBgn";
    my %update; # update{search_id} = $updateCode
    my %targets; # targets{search_id} = number of targets
    for my $id (keys %apms) {
	my $prev = $sid2Bait{$id} // $dummyFBgn;
	$prev = $dummyFBgn if exists $rejectMap{$id};
	my @targets = keys %{ $updateFBgn{$prev} };
	$targets{$id} = \@targets;
	if ($prev eq $dummyFBgn) {
	    my $new = $prev;
	    $newBait{$id} = $new;
	    if (exists $rejectMap{$id}) {
		$update{$id} = $rejectMap{$id};
	    } else {
		$update{$id} = $updateCodes{unknown};
	    }
	} elsif (@targets == 0) {
	    warn "don't recognize $prev";
	    $update{$id} = $updateCodes{unknown};
	    $newBait{$id} = $prev;
	} elsif (1 == @targets) {
	    # one possible FBgn id
	    my $new = $targets[0];
	    $newBait{$id} = $new;
	    if ($new eq $prev) {
		$update{$id} = $updateCodes{nochange};
	    } else {
		$update{$id} = $updateCodes{single};
	    }
	} else {
	    # there are many possible updates for this FBgn id
	    # seek peptides for all possibilities
	    my $nonZero = 0; # how many targets have peptides
	    my $mostPep = 'fb'; # which target has the most peptides
	    my $maxPep = 0; # what's the max tsc count for any target
	    for my $t (@targets) {
		my $tsc = $apms{$id}{$t} // 0;
		#say "tsc{$t} = $tsc" if ($id eq 'f05783');
		$nonZero++ if $tsc > 0;
		if ($tsc > $maxPep) {
		    $maxPep = $tsc;
		    $mostPep = $t;
		}
	    }
	    #die Dumper($nonZero, $mostPep, $maxPep, \@targets, $apms{f05783}{FBgn0265189}) if ($id eq '249877');
	    if ($nonZero == 0) {
		
		# find no peptides for any possible updates
		my @nonCG = grep { $updateFBgn{$prev}{$_} !~ /^CG\d+/ } @targets;
		#my @tSymb = map  { $updateFBgn{$prev}{$_} } @targets;
		#die Dumper($prev, \@nonCG, \@tSymb);
		if (@nonCG == 1) {
		    $newBait{$id} = $nonCG[0];
		    $update{$id} = $updateCodes{nonCG};
		} elsif (@nonCG > 1) {
		    $newBait{$id} = (sort @nonCG)[0];
		    $update{$id} = $updateCodes{alphaNonCG};
		} else {
		    $newBait{$id} = (sort @targets)[0];
		    $update{$id} = $updateCodes{alphaCG};		    
		}
	    } else {
		$newBait{$id} = $mostPep;
		if ($nonZero == 1) {
		    $update{$id} = $updateCodes{singlePep};
		} else {
		    $update{$id} = $updateCodes{mostPep};
		}
	    }
	}
    }
    ##warn "!!! check f11678 -> FBgn0266452 !!!\n"; # checkdit mofo!

    my $preComments = "# identified baits in $in according to $baitIDFile and $flybaseFile\n";
    if (defined $baitKeyOut) {
	my @outCols = qw(search_id ms_inst_run_id prevBait newBait update);
	open my $OUT, ">", $baitKeyOut or die "Can't write to $baitKeyOut. $!";
	print $OUT $preComments;
	say $OUT join "\t", @outCols;

	for my $id (sort {$a <=> $b} keys %newBait) {
	    my $updateString;
	    if ($update{$id} eq $updateCodes{unknown} || 
		exists $rejectCodes{$update{$id}}) {
		$updateString = $update{$id};
	    } else {
		$updateString = join ",", keys %{ $updateFBgn{$sid2Bait{$id}} };
	    }
	    say $OUT (join "\t", $id, $msID{$id}, $sid2Bait{$id} // $dummyFBgn,
		      $newBait{$id}, $update{$id}, $updateString);
	}
	close $OUT;
    }
    
    say "update APMS";
    updateAPMS($apmsOut, $in, \%newBait, $preComments);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	baitid => '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.10-04-2016.tsv',
	flybase => '/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2015_04.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out apmsout < $defaultString -keyout ".
	"baitkey.log -nologp >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "baitid=s", "flybase=s", "keyout=s",
	"nologp");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{baitid});
    checkExist('f', $opts{flybase});

    return %opts;
}

# ret{oldfbgn} = {$fbgn1 => $symbol1, $fbgn2 => $symbol2,...}
sub makeUpdateMap {
    my ($ret, $flybaseFile) = @_;

    open my $IN, "<", $flybaseFile or die "Can't read from $flybaseFile. $!";
    while (my $line = <$IN>) {
	next if $line =~ /^#/;
	next if length($line) < 3;
	chomp $line;
	my ($symbol, $target, $prevString) = split /\t/, $line;
	say $line if ! defined $target;
	die "multiple maps for target '$target'" if exists $ret->{$target};
	# this shouldn't happen because $target is supposed to be the unique,
	#   most recent fbgn
	
	$ret->{$target}{$target} = $symbol;
	my @prev = split /,/, $prevString;
	for my $p (@prev) {
	    $ret->{$p}{$target} = $symbol;
	}
    }

    return;
}

# ret{search_id} = [{bait => 'given_bait', line => 'whole line'}]
sub readAPMS {
    my ($ret, $apmsFile, $msID) = @_;

    open my $IN, "<", $apmsFile or die "Can't read from $apmsFile. $!";

    my $reader = \&getLineDP4APMS;
    $reader = \&getLineDP4_1APMS unless exists $opts{nologp};
    my %row;
    while($reader->(\%row, $IN)) {
	next if $row{prey_ref} eq "TAG";
	#warn "multiple peptides ret->{$row{search_id}}{$row{prey_ref}}" 
	#    if exists $ret->{$row{search_id}}{$row{prey_ref}};
	$ret->{$row{search_id}}{$row{prey_ref}} += $row{total_peptides};
	$msID->{$row{search_id}} = $row{ms_inst_run_id};
    }

    return;
}

sub updateAPMS {
    my ($out, $in, $bait, $preComments) = @_;

    my @cols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date    ms_inst_run_id);
    
    my $reader = \&getLineDP4APMS;
    if (! exists $opts{nologp}) {
	$reader = \&getLineDP4_1APMS;
	push @cols, 'logp';
    }

    open my $OUT, ">", $out or die "can't write to $out. $!";
    print $OUT $preComments if defined $preComments;
    say $OUT join "\t", @cols;
    
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my %row;
    while($reader->(\%row, $IN, 'line')) {
	my $id = $row{search_id};
	$row{bait_ref} = $bait->{$id};
	say $OUT join "\t", map { $row{$_} // die Dumper(\%row) } @cols;
    }
    close $IN;
    close $OUT;

    return;
}

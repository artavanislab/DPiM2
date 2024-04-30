#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsHashRef);
use DpimLib qw(getLineDP4APMS getLineDP4_1APMS);

# identify the bait for this run programmatically 
# "programmatic" because I don't have a lookup table as of 03-30-2016,

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
    my $out = $opts{out};
    my $prevBaitFile = $opts{prevbait};
    my $flybaseFile = $opts{flybase};
    my $apmsOut = $opts{apmsout};
    my $flFile = $opts{fltab};
    
    say "reading previous baits";
    # prevBait{$ms_inst_run_id} = FBgnXXX
    my %prevBait;
    readColsHashRef(\%prevBait, $prevBaitFile, [qw(ms_inst_run_id oldfbgn)]);

    
    my %flMap;
    readColsHashRef(\%flMap, $flFile, [qw(ms_inst_run_id tap_id)]) 
	if defined $flFile;

    say "reading fbgn udate map";
    # updateFBgn{prev}{new} = new gene symbol
    my %updateFBgn;
    makeUpdateMap(\%updateFBgn, $flybaseFile);

    say "reading AP-MS data";
    my %apms;
    readAPMS(\%apms, $in);
    #die Dumper($apms{f05783}, $apms{f05783}{FBgn0265189});

    my %updateCodes = (
	unknown => 'unknown', # don't recognize this run id
	fl => 'tapFL', # this run has a tap ID starting "FL", i.e. it's dross
	nochange => 'nochange', # didn't update
	single => 'singleUpdate', # single possible update
	singlePep => 'singlePeptideMatch', # only one update matches peptides
	mostPep => 'mostPeptideMatches', # multiple peptide matches
	nonCG => 'nonCG', # the only symbol that's not CGXXXX 
	alphaCG => 'alphaGC', # alphanumerically first CG
	alphaNonCG => 'alphaNonCG', # alphanumerically first not CG
	);

    # newBait{$ms_inst_run_id} = FBgnXXX
    # *this is the work product of this script*
    my %newBait;
    my $dummyFBgn = 'FBgn0000000';
	
    say "identify baits";
    my %update; # update{ms_id} = $updateCode
    my %targets; # targets{ms_id} = number of targets
    for my $id (keys %apms) {
	my $prev = $prevBait{$id} // $dummyFBgn;
	my @targets = keys %{ $updateFBgn{$prev} };
	$targets{$id} = \@targets;
	if (1 == @targets) {
	    my $new = $targets[0];
	    $newBait{$id} = $new;
	    if ($new eq $prev) {
		$update{$id} = $updateCodes{nochange};
	    } else {
		$update{$id} = $updateCodes{single};
	    }
	} elsif ($prev = $dummyFBgn) {
	    my $new = $prev;
	    $newBait{$id} = $new;
	    if (exists $flMap{$id}) {
		$update{$id} = $updateCodes{fl};
	    } else {
		$update{$id} = $updateCodes{unknown};
	    }
	} else {
	    # seek peptides for all keys
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
	    #die Dumper($nonZero, $mostPep, $maxPep, \@targets, $apms{f05783}{FBgn0265189}) if ($id eq 'f05783');
	    if ($nonZero == 0) {
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
    #die Dumper(\%newBait);

    my $manCorrectWarning = 'WARNING: FBgn0262707 maps to multiple targets: FBgn0266452/CTPSyn for f06099 and f08108, and FBgn0266441/CG45071 for f11678.  Manually correct this to FBgn0266452/CTPSyn.';
    if ($newBait{f11678} ne 'FBgn0266441') {
	$manCorrectWarning='';
	#die "edge case disappeared?? newBait{f11678} ne 'FBgn0266441'";
    } else {
	warn $manCorrectWarning."\n";
	$newBait{f11678} = 'FBgn0266452';
	$update{f11678} = 'humanPreventsInconsistency_cf_f06099_and_f08108';
    }
    
    my @outCols = qw(ms_inst_run_id prevBait newBait update targets);
    my $preComments = "# identified baits in $in according to $prevBaitFile and $flybaseFile\n";
    $preComments.="# $manCorrectWarning\n" if $manCorrectWarning;
    
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    print $OUT $preComments;
    say $OUT join "\t", @outCols;

    for my $id (sort keys %newBait) {
	my $updateString;
	if ($update{$id} eq $updateCodes{unknown} || 
	    $update{$id} eq $updateCodes{fl}) {
	    $updateString = $update{$id};
	} else {
	    $updateString = join ",", keys %{ $updateFBgn{$prevBait{$id}} };
	}
	say $OUT (join "\t", $id, $prevBait{$id} // $dummyFBgn, $newBait{$id}, 
		  $update{$id}, $updateString);
    }
    close $OUT;

    if (defined $apmsOut) {
	say "update APMS";
	updateAPMS($apmsOut, $in, \%newBait, $preComments);
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	prevbait => '/home/glocke/DPiM/dpim4/withInstr/apmsData/baitMap.tsv',
	flybase => '/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2016_01.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $defaultString -apmsout ".
	"updated.apms -logp -fltab optionalFL.tab >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "prevbait=s", "flybase=s", "apmsout=s",
	"logp", "fltab=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{prevbait});
    checkExist('f', $opts{flybase});
    checkExist('f', $opts{fltab}) if exists $opts{fltab};

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

# ret{ms_inst_run_id} = [{bait => 'given_bait', line => 'whole line'}]
sub readAPMS {
    my ($ret, $apmsFile) = @_;

    open my $IN, "<", $apmsFile or die "Can't read from $apmsFile. $!";

    my $reader = \&getLineDP4APMS;
    $reader = \&getLineDP4_1APMS if exists $opts{logp};
    my %row;
    while($reader->(\%row, $IN, 'line')) {
	next if $row{prey_ref} eq "TAG";
	#warn "multiple peptides ret->{$row{ms_inst_run_id}}{$row{prey_ref}}" 
	#    if exists $ret->{$row{ms_inst_run_id}}{$row{prey_ref}};
	$ret->{$row{ms_inst_run_id}}{$row{prey_ref}} += $row{total_peptides};
    }

    return;
}

sub updateAPMS {
    my ($out, $in, $bait, $preComments) = @_;

    my @cols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date    ms_inst_run_id);
    
    my $reader = \&getLineDP4APMS;
    if (exists $opts{logp}) {
	$reader = \&getLineDP4_1APMS;
	push @cols, 'logp';
    }

    open my $OUT, ">", $out or die "can't write to $out. $!";
    print $OUT $preComments if defined $preComments;
    say $OUT join "\t", @cols;
    
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my %row;
    while($reader->(\%row, $IN, 'line')) {
	my $id = $row{ms_inst_run_id};
	$row{bait_ref} = $bait->{$id};
	say $OUT join "\t", map { $row{$_} } @cols;
    }
    close $IN;
    close $OUT;

    return;
}

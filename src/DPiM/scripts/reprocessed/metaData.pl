#!/usr/bin/env perl

##############################################################################80

use v5.10; 
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(min max sum);
use List::MoreUtils qw(each_array);
#use File::Temp qw(tempfile);
use DateTime;
use DateTime::Format::Strptime;

use HomeBrew::IO qw(checkExist readCol readColsRef readHeader readCols writeCols
                    readList readColsHashRef);
use DpimLib qw(getLineAPMS readHS getLineDP4APMS getLineRawAPMS  getLineHyperspecAPMS networkHashFromEdgeList getLineMayAPMS);

## find the meta-data associated with runs
## * Bait identity
## * Should it be scored/retained?

## Steps:
## 1. capture previously known baits
## 2. find baits corresponding to tap id's
## 3. handle disagreements between baits identified in steps #1, #2
##    - decisions were made by humans; this step applies those judgments.
## 4. remove experiments that are not DPIM1 experiments (marked "ignore" by Julian")
## 4.1. rescue some of these "ignore" runs 
## 5. figure out which of the "wonky" tap id's should be retained
## 6. publish final result

my $baseOut = $ARGV[0];
die "usage: $0 baseOutputFile";

## input
my $dataDir = '/home/glocke/DPiM/augRemap/apmsData';
my $annotationDir = '/home/glocke/DPiM/augRemap/apmsData/annotation';
my $prevBaitKeyFile = "$annotationDir/prevBaitKey.tsv";
my $TFRunFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/TF_Run_Numbers.txt';
my $wonkyDiscovered = "$annotationDir/RAO_rescueWonkyDiscovered.tsv";
my $dpim1Wonky = "/home/glocke/DPiM/augRemap/apmsData/dpim1/nodeydey.metadata.tsv";
my $haSummary = "$annotationDir/HA_summary_pipeline_slim.txt";

{
    my @dirs = ($dataDir, $annotationDir);
    checkExist('d', $_) for @dirs;
    my @files = ($prevBaitKeyFile, $TFRunFile, $wonkyDiscovered, $dpim1Wonky,
		 $haSummary);
    checkExist('f', $_) for @files;
}

## output
my $metaBase = "$annotationDir/metaDataTable6_10-04-2016";
my $meta1 = "$metaBase.tsv"; 
my $meta2 = "$metaBase.tap.tsv";
my $meta3 = "$metaBase.tap.ignore.tsv";
my $meta4 = "$metaBase.tap.ignore.byHand.tsv";
my $final = "$annotationDir/baitIdentityTable.10-04-2016.tsv";

if (! -e $prevBaitKeyFile) {
    ## collect previous bait identifications
    my $prevFile = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.r607.newFBgn';
    my %rid2Bait; ## $rid2Bait{ms_inst_run_id} = FBgn bait previously identified

    my %row;
    open my $IN, "<", $prevFile or die "Can't read from $prevFile: $!";
    while (getLineDP4APMS(\%row, $IN)) {
	if (exists $rid2Bait{$row{ms_inst_run_id}} && 
	    $rid2Bait{$row{ms_inst_run_id}} ne $row{bait_ref}) {
	    warn "rid2Bait{$row{ms_inst_run_id}} = $rid2Bait{$row{ms_inst_run_id}} and *$row{bait_ref}*";
	}
	$rid2Bait{$row{ms_inst_run_id}} = $row{bait_ref};
    }
    open my $OUT, ">", $prevBaitKeyFile or die "can't write to $prevBaitKeyFile: $!";
    say $OUT join "\t", qw(ms_inst_run_id bait_ref);
    say $OUT join "\t", $_, $rid2Bait{$_} for sort keys %rid2Bait;
    close $OUT;
}

if (! -e $meta1) {
    ## make a metadata table
    my @inFiles = qw(DPiM1_r607_protein_views.out_LDAmin7_PA_PSA01_160808
DPiM1_rejects_r607_protein_views.out_LDAmin7_PA_PSA01_160811
DPiM2_r607_protein_views.out_LDAmin7_PA_PSA01_160816
DPiM3_r607_protein_views.out_LDAmin7_PA_PSA01_160824);
    @inFiles = map { "$dataDir/$_" } @inFiles;
    my $dir = '/home/glocke/DPiM/augRemap/apmsData';
    my %prevBaitKey;
    readColsHashRef(\%prevBaitKey, $prevBaitKeyFile, 
		    [qw(ms_inst_run_id bait_ref)]);
    
    my @capture = qw(search_id ms_inst_run_id tap_id sample_date id_string);

    my %runs;
    for my $f (@inFiles) {
	say $f;
	open my $IN, "<", $f or die "Can't read from $f: $!";
	my %row;
	my $reject = 0;
	$reject = 1 if $f =~ /rejects/;
	while (getLineMayAPMS(\%row, $IN)) {
	    my $sid = $row{search_id};
	    my $rid = $row{ms_inst_run_id};
	    if (! exists $runs{$sid}) {
		$runs{$sid} = { map { $_ => $row{$_} } @capture };
		## convert date from MMDDyy to YYYY-MM-DD 
		$runs{$sid}{sample_date} =~ /^(\d\d?)(\d\d)(\d\d)$/ 
		    or die "can't parse date for $sid/$rid", Dumper(\%row);
		$runs{$sid}{sample_date} = "20$3-$1-$2";
		$runs{$sid}{bait} = $prevBaitKey{$rid} // 'FBgn0000000';
		$runs{$sid}{rejected} = $reject;
		$runs{$sid}{wonky} = 0;
		$runs{$sid}{wonky} = 1 if $row{tap_id} =~ /^hFH/;
		$runs{$sid}{wonky} = 1 if $row{tap_id} =~ /\D$/;
		$runs{$sid}{wonky} = 1 if $row{tap_id} =~ /^FL\d/;
	    }
	}
    }

    my @order = sort { $runs{$a}{ms_inst_run_id} cmp $runs{$b}{ms_inst_run_id} 
    } keys %runs;
    
    open my $OUT, ">", $meta1 or die "Can't write to $meta1. $!";
    my @cols = qw(search_id ms_inst_run_id tap_id bait sample_date id_string rejected wonky);
    say $OUT "# $0 make a metadata table";
    say $OUT join "\t", 'search_id', @cols;
    for my $id (@order) {
	say $OUT join "\t", $id, map { $runs{$id}{$_} // die Dumper($id, $runs{$id})} @cols;
    }
    close $OUT;
}

if (! -e $meta2) {
    ## update the metadata table with tap_id information
    my $tapFile = "$annotationDir/FH_Plate_Contents_for_JM_052715_short.updateFBgn.txt";

    my @cols = qw(search_id ms_inst_run_id bait tap_id sample_date id_string rejected wonky);
    my @meta;
    readColsRef(\@meta, $meta1, \@cols);
    my %tapHash;
    readColsHashRef(\%tapHash, $tapFile, [qw(TAP_Alias Flybase_ID)], "\t");

    my $unknown = 'FBgn0000000';

    @cols = qw(search_id ms_inst_run_id bait tap_id tap_bait sample_date id_string rejected wonky tap_test);

    open my $OUT, ">", $meta2 or die "Can't write to $meta2. $!";
    say $OUT "# $0 add tap_table info to metadata";
    say $OUT join "\t", @cols;

    my %eleven = map { $_ => 1 } qw(w33734 w33735 w33736 w33737 w33738 w33739 w33740 w33744 w33745 w33746 w33747); ## these are experiments which had previous records indicating doubled id's; they might look like disagreements, but GLL inspected them and the tap_id's agree for both duplicates 09-07-2016
    
    for my $expt (@meta) {
	my $tap = $expt->{tap_id};
	my $tapBait = $tapHash{$tap} // $unknown;
	$tapBait = $unknown if $tapBait eq 'N/A';
	if (! defined $tapHash{$tap} && $tap =~ /FH(\d+)/) {
	    my $tapN = sprintf("%04d", $1);
	    $tapBait = $tapHash{"FH$tapN"} // $unknown;
	}
	$expt->{tap_bait} = $tapBait;

	my $tapTest;
	if ($expt->{tap_bait} eq $unknown && 
	    $expt->{bait} eq $unknown) 
	{
	    $tapTest = "both_unknown";
	} 
	elsif ($expt->{tap_bait} eq $unknown 
		 && $expt->{bait} ne $unknown) 
	{
	    $tapTest = "known_but_not_in_table";
	} 
	elsif ($expt->{tap_bait} ne $unknown && 
		 $expt->{bait} eq $unknown) 
	{
	    $tapTest = "discovered_by_tap_table";
	} 
	elsif ($expt->{tap_bait} eq $expt->{bait}) {
	    $tapTest = "agree";
	} 
	elsif ($expt->{tap_bait} ne $expt->{bait}) {
	    if (exists $eleven{$expt->{ms_inst_run_id}}) {
		$tapTest = 'elevenPrevDupes';
	    } else {
		$tapTest = "DISAGREE";
	    }
	} else {
	    die "ehh??", Dumper($expt);
	}
	$expt->{tap_test} = $tapTest;

	say $OUT join "\t", map {$expt->{$_}} @cols;
    }
    close $OUT;
}

if (! -e $meta3) {
    my %ignoreRuns = readIgnoreRuns($haSummary);

    my @metaCols = qw(search_id ms_inst_run_id  bait    tap_id  tap_bait        sample_date     id_string       rejected        wonky     tap_test);
    my @meta;
    readColsRef(\@meta, $meta2, \@metaCols);

    open my $OUT, ">", $meta3 or die "can't write to $meta3: $!";
    say $OUT "# $0 correcting $meta2 for TF runs";
    say $OUT join "\t", @metaCols;
    for my $row (@meta) {
	my $rid = $row->{ms_inst_run_id};
	if (exists $ignoreRuns{$rid}) {
	    $row->{rejected}+= 0.5;
	    $row->{tap_test}.="_TF_run";
	}
	say $OUT join "\t", map {$row->{$_}} @metaCols;
    }
    close $OUT;
}

if (! -e $meta4) {
    my $cmd = "/home/glocke/DPiM/scripts/reprocessed/correctMetaData.pl -meta $meta3 -out $meta4";
    say $cmd;
    system($cmd);
}

if (! -e $final) {

    my @meta;
    readColsRef(\@meta, $meta4, [qw(search_id ms_inst_run_id tap_id newbait rejected wonky count tap_test)]);
    my %wonkyCorrects = readWonkyDiscovered($wonkyDiscovered);
    my %dpimBaits;
    {
	## A certain subset of "wonky" experiments were present in dpim1
	## some of these are good, some aren't
	## refer to email from Julian Mintseris Mon 9/12/2016 11:29 AM
	my %dontKeep = map {$_ => 1} 
	    qw(f14574 f14592 f14901 f14924 f07686 f07687 ltqg1755);
	
	## i have double-checked that none of the FBgn's have been updated
	my @read;
	readColsRef(\@read, $dpim1Wonky, [qw(ms_inst_run_id prevbait)]);
	for my $row (@read) {
	    next if exists $dontKeep{$row->{ms_inst_run_id}};
	    $dpimBaits{$row->{ms_inst_run_id}} = $row->{prevbait}
	}
    }
    
    for my $m (@meta) {
	my $rid = $m->{ms_inst_run_id};
	$m->{retain} = 'yes';

	if (exists $wonkyCorrects{$rid}) {
	    $m->{tap_test} = $wonkyCorrects{$rid};
	} elsif (exists $dpimBaits{$rid}) {
	    $m->{tap_test} = "wonky_dpim1";
	    $m->{newbait} = $dpimBaits{$rid};
	} else {
	    if ($m->{rejected}) {
		$m->{retain} = 'no';
		$m->{tap_test} .= "_rejected";
	    } elsif ($m->{wonky}) { 
		$m->{retain} = 'no';
		$m->{tap_test} .= "_wonky";
	    } elsif ($m->{tap_test} =~ /DISAGREE/ ||
		     $m->{tap_test} =~ /both_unknown/) 
	    {
		$m->{retain} = 'no';
		#$m->{tap_test} .= "_updateDisagreement";
		# removed this line because tap_test already records this info
	    }
	}
	$m->{notes} = $m->{tap_test};
	$m->{notes} =~ s/\s/_/g;
    }

    my @outCols = qw(search_id ms_inst_run_id bait_ref tap_id retain notes count);
    open my $OUT, ">", $final or die "can't write to $final: $!";
    say $OUT "# $0 computation assisted by humans to ";
    say $OUT join "\t", @outCols;
    $outCols[2] = 'newbait';
    for my $m (@meta) {
	say $OUT join "\t", map {$m->{$_} // die Dumper($m)} @outCols;
    }
}

sub readIgnoreRuns {
    my ($inFile) = @_;

    my @cols = ("Run number", "Julian's list");
    my %read;
    readColsHashRef(\%read, $inFile, \@cols, "\t");
    
    my %ret;
    for my $k (keys %read) {
	$ret{$k} = 1 if $read{$k} eq 'Ignore';
    }

    return %ret;
}

sub readWonkyDiscovered {
    my ($inFile) = @_;
    
    my @cols = qw(ms_inst_run_id	bait	tap_id	tap_bait	sample_date	id_string	rejected	wonky	tap_test	geneSymbol1	geneSymbol2	RAO_Determination	RAO_Confidence);
    my %colI = map { $cols[$_] => $_ } 0..$#cols;

    my %ret; ## ret{ms_inst_run_id} = RAO_determination iff this is to be retained
    open my $IN, "<", $inFile or die "Can't read from $inFile: $!";
    while (<$IN>) {
	chomp;
	my @spl = split /\t/;
	my $det = $spl[$colI{RAO_Determination}];
	if ($det =~ /retain/i) {
	    my $rid = $spl[$colI{ms_inst_run_id}];
	    $ret{$rid} = $det;
	}
    }

    return %ret;
}

#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef readColsHashRef readHeader);
use DpimLib qw(networkHashFromEdgeList getLineDP4APMS getLineHyperspecAPMS);

# for each new experiment, 
#   * how many proteins are seen which were never seen before?
#   * how many edges are seen which were never seen before?
# where "never before" means, no prior experiment has seen it

my %opts = getCommandLineOptions();

{
    my $apmsFile = $opts{apms};
    my $netFile = $opts{net};
    my $oldNetFile = $opts{oldnet};
    my $statsFile = $opts{stats};
    my $s2rFile = $opts{s2r};
    my $out = $opts{out};
    
    my %s2r;
    readColsHashRef(\%s2r, $s2rFile, [qw(Gene Protein)]);
    my %stats = readStats($statsFile);
    my %dataEra = readEra(map { $_ => $opts{$_} } qw(dp1 dp1rej dp2 dp3));
    my (%prey, %bait);
    readAPMS(\%prey, \%bait, $apmsFile);
    my %net;
    networkHashFromEdgeList(\%net, $netFile, undef);
    my %oldNet;
    networkHashFromEdgeList(\%oldNet, $oldNetFile, undef);

    my @sids = keys %prey;
    @sids = sort { $stats{$a}{rid} <=> $stats{$b}{rid} } @sids;
    @sids = sort { $dataEra{$a} <=> $dataEra{$b} } @sids;
    
    my @cols = readHeader($statsFile);
    push @cols, qw(newPrey newEdges newS2r oldEdges);

    

    # return the number of new proteins in this experiment
    my $seenPreyTest = nodeTest();
    my $s2rTest = nodeTest(\%s2r);
    my $seenEdgeTest = edgeTest(\%net);
    my $oldEdgeTest = edgeTest(\%oldNet);
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 found statistics on added proteins/edges ";
    say $OUT "# apms file $apmsFile";
    say $OUT "# network file $netFile";
    say $OUT "# previous network file $oldNetFile";
    say $OUT "# s2r proteome file $s2rFile";
    say $OUT "# stats file $statsFile";
    say $OUT join "\t", @cols;
    # main loop
    for my $sid (@sids) {
	my $newPrey = $seenPreyTest->($prey{$sid});
	my $newEdges = $seenEdgeTest->($prey{$sid});
	my $oldEdges = $oldEdgeTest->($prey{$sid});
	my $s2r = $s2rTest->($prey{$sid});
	chomp $stats{$sid}{line};
	say $OUT join "\t", $stats{$sid}{line}, $newPrey, $newEdges, $s2r
	    , $oldEdges;
    }
    
    close $OUT; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	apms => '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.rmDupes.pepFDR.sumIso.applyLC.untrans.tscCut.rebrand',
	net => '/home/glocke/DPiM/augRemap/nrBait_yRej_nBrand_11-16-2016/nrBait.net',
	oldnet => '/home/glocke/DPiM/prevDPIM/DPIM1_scores.r6.07.updateFBgn.tsv',
	stats => '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.rmDupes.pepFDR.sumIso.applyLC.untrans.tscCut.rebrand.statsBySID',
	dp1 => '/home/glocke/DPiM/augRemap/apmsData/DPiM1_r607_protein_views.out_LDAmin7_PA_PSA01_160808.dp4.statsBySID',
	dp1rej => '/home/glocke/DPiM/augRemap/apmsData/DPiM1_rejects_r607_protein_views.out_LDAmin7_PA_PSA01_160811.dp4.statsBySID',
	dp2 => '/home/glocke/DPiM/augRemap/apmsData/DPiM2_r607_protein_views.out_LDAmin7_PA_PSA01_160816.dp4.statsBySID',
	dp3 => '/home/glocke/DPiM/augRemap/apmsData/DPiM3_r607_protein_views.out_LDAmin7_PA_PSA01_160824.dp4.statsBySID',
	s2r => '/home/glocke/DPiM/cell_5871_mmc2_s2rProteome.r6.07.txt',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "apms=s", "net=s", "oldnet=s", "stats=s", 
	       "dp1=s", "dp1rej=s", "dp2=s", "dp3=s", "s2r=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{apms});
    checkExist('f', $opts{net});
    checkExist('f', $opts{stats});
    checkExist('f', $opts{dp1});
    checkExist('f', $opts{dp1rej});
    checkExist('f', $opts{dp2});
    checkExist('f', $opts{dp3});
    checkExist('f', $opts{s2r});

    return %opts;
}


sub readAPMS {
    my ($prey, $bait, $inFile) = @_;

    my $reader = \&getLineDP4APMS;
    if ($inFile =~ /nrBait$/) {
	$reader = \&getLineHyperspecAPMS;
    }
    
    #my @cols = qw(bait_ref prey_ref)
    my %row;
    open my $IN, "<", $inFile or die "can't read $inFile. $!";
    while ($reader->(\%row, $IN)) {
	my $sid = $row{search_id};
	$prey->{$sid}{$row{prey_ref}} = 1;
	$bait->{$sid} = $row{bait_ref};
    }	
    close $IN;
    
    return;
}

sub readStats {
    my ($inFile) = @_;

    my @read;
    readColsRef(\@read, $inFile, ['search_id', 'ms_inst_run_id'], 'line');
    my %ret;
    for my $row (@read) {
	my $line = $row->{line};
	chomp $line;
	$row->{ms_inst_run_id} =~ /([\d\.]+)/ or die "Can't parse run id ".
	    Dumper($row);
	$ret{$row->{search_id}} = {
	    rid => $1,
	    line => $line,
	};
    }

    return %ret;
}

sub readEra {
    my %inFiles = @_;
    my $inFile;

    my %eraVals = qw( dp1 1 dp1rej 2 dp2 3 dp3 4 );
    
    my %ret;
    
    for my $era (keys %eraVals) {
	my @read;
	readColsRef(\@read, $inFiles{$era}, ['search_id']);
	for my $row (@read) {
	    $ret{$row->{search_id}} = $eraVals{$era};
	}
    }

    return %ret;
}
    
# returns a closure telling you if you've seen a given node before
# if the optional argument $nodeList is define, only count nodes in $nodeList
sub nodeTest {
    my ($nodeList) = @_;

    my %seen;
    
    if (defined $nodeList) {
	return sub { 
	    my ($expt) = @_;
	    my $ret = 0;
	    for my $prey (keys %$expt) {
		if (! exists $seen{$prey} && exists $nodeList->{$prey}) {
		    $ret++;
		    $seen{$prey} = 1;
		}
	    }
	    return $ret;
	};
    } else {
	return sub {
	    my ($expt) = @_;
	    my $ret = 0;
	    for my $prey (keys %$expt) {
		if (! exists $seen{$prey}) {
		    $ret++;
		    $seen{$prey} = 1;
		}
	    }
	    return $ret;
	};
    }
}

# returns a closure telling you if you've seen a given edge before
# only count edges in $net
# keys in $net are assumed to be sorted; i.e. $net{$a}{$b} is only defined
#   if $a cmp $b
sub edgeTest {
    my ($net) = @_;

    my %seen;
    return sub { 
	my ($expt) = @_;

	my $ret = 0;
	
	my @prey = sort keys %$expt;
	for my $i (0..($#prey -1)) {
	    my $p1 = $prey[$i];
	    next unless exists $net->{$p1};
	    for my $j (($i+1)..$#prey) {
		my $p2 = $prey[$j];
		next unless exists $net->{$p1}{$p2};

		if (! exists $seen{$p1}{$p2}) {
		    $ret++;
		    $seen{$p1}{$p2} = 1;
		}
	    }
	}
	return $ret;
    };
}

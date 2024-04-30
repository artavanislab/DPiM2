#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(readGoDB networkHashFromEdgeList);

# shuffle annotations across the network
# compare the results to hyper-geometric test

my %opts = getCommandLineOptions();

{
    my $out = $opts{out};
    my $goDBFile = $opts{godb};
    my $goNameFile = $opts{goname};
    my $hgFile = $opts{hg};
    my $netFile1 = $opts{net1};
    my $netFile2 = $opts{net2};
    my $nTest = $opts{ntest};
    
    my %goDB;
    say "parsing $goDBFile";
    readGoDB(\%goDB, $goDBFile);



    say "retrieving goNameMap";
    my $goNameMap = makeGoNameMap($goNameFile);
    retrieve($goNameFile);    
    $goNameMap->{UNKNOWN} = { name => '"unknown"' };

    {
	# standardize go terms
	my %oldNew;
	for my $go (keys %goHist) {
	    if (exists $goNameMap->{$go}{standard_id}) {
		$oldNew{$go} = $goNameMap->{$go}{standard_id};
	    }
	}
	for my $oldGo (keys %oldNew) {
	    my $newGo = $oldNew{$oldGo};
	    $goHist{$newGo} += $goHist{$oldGo};
	    delete $goHist{$oldGo};
	}
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	godb => '/home/glocke/DPiM/flybase/gene_association.fb',
	goname => '/home/glocke/DPiM/flybase/goMap.storable',
	hg => 'edgesMissingFromDpim3.edge.GOHist.binomTest.tab',
	net1 => '/home/glocke/DPiM/dpim2/dpim2_nrtap.120123.nrBait.58.17.network',
	net2 => '/home/glocke/DPiM/analysis/nrDpim2/edgesMissingFromDpim3.tab',
	ntest => '10000',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "godb=s", "goname=s", "hg=s", "net1=s", 
	       "net2=s", "ntest=i");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

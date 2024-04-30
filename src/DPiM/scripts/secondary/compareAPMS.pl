#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4APMS);

# make a table comparing experiments from dpim1 vs dpim1a vs dpim4

my %opts = getCommandLineOptions();

{
    my $out = $opts{out};
    my $dp1 = $opts{dp1};
    my $dp1a = $opts{dp1a};
    my $dp4 = $opts{dp4};

    my @datasets = qw(dp1 dp1a dp4);
    my %parse; # parse{$ms_inst_run_id}{bait} = $FBgn
    #          # parse{$ms_inst_run_id}{date} = $sample_date
    #          # parse{               }{$preyFBgn}{dp1}  = TSC
    #          # parse{               }{$preyFBgn}{dp1a} = TSC
    #          # parse{               }{$preyFBgn}{dp4} = TSC
    my %runIDs = parse(\%parse, $dp1, 'dp1');
    parse(\%parse, $dp1a, 'dp1a');
    parse(\%parse, $dp4, 'dp4');

    open my $OUT, ">", $out or die "Can't write to $out: $!";
    say $OUT join "\t", qw(ms_inst_run_id date bait prey dp1TSC dp1aTSC dp4TSC);
    for my $rid (sort keys %runIDs) {
	my $bait = $parse{$rid}{bait};
	my $date = $parse{$rid}{date};
	my @prey = grep /FBgn/, keys %{ $parse{$rid} };
	for my $p (@prey) {
	    my $sum = 0;
	    for my $ds (@datasets) {
		$parse{$rid}{$p}{$ds} //= 0;
		die Dumper($rid, $parse{$rid}, \@prey) 
		    if $parse{$rid}{$p}{$ds} =~ /-/;
		
		$sum += $parse{$rid}{$p}{$ds};
	    }
	    $parse{$rid}{$p}{sum} = $sum;
	}
	#die Dumper(\@prey);
	@prey = sort { $parse{$rid}{$b}{sum} <=> $parse{$rid}{$a}{sum} }
	    @prey;
	for my $p (@prey) {
	    say $OUT join "\t", $rid, $date, $bait, $p
		, map {$parse{$rid}{$p}{$_}} @datasets;
	}
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	dp1 => '/home/glocke/DPiM/dpim1/dpim_all_lnSP_nrtap.101223.plusUniqMsInst.dp4.updateFBgn.sumIso',
	dp1a => '/home/glocke/DPiM/newDpim2/DPiM1_05-05-2016Map.dp4.newBait.logPCut.sumIso',
	dp4 => '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.newFBgn.rmDupes.sumIso',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "dp1=s", "dp1a=s", "dp4=s");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{dp1});
    checkExist('f', $opts{dp1a});
    checkExist('f', $opts{dp4});

    return %opts;
}


# parse{$ms_inst_run_id}{bait} = $FBgn
# parse{$ms_inst_run_id}{date} = $sample_date
# parse{               }{$preyFBgn}{dp1}  = TSC
# parse{               }{$preyFBgn}{dp1a} = TSC
# parse{               }{$preyFBgn}{dp4} = TSC
sub parse {
    my ($ret, $in, $dsKey, $runIDs) = @_;

    my %newRunIDs;
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN)) {
	my $rid = $row{ms_inst_run_id};
	next if defined $runIDs && ! exists $runIDs->{$rid};
	$newRunIDs{$rid} = 1;
	$ret->{$rid}{bait} = $row{bait_ref};
	$ret->{$rid}{date} = $row{sample_date};
	$ret->{$rid}{$row{prey_ref}}{$dsKey} = $row{total_peptides};
    }

    return %newRunIDs if ! defined $runIDs;
    return;
}

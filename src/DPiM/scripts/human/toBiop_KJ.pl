#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(each_array);
use HomeBrew::IO qw(checkExist readHeader readColsRef);
#use DpimLib qw(getLineDP4APMS);

# convert the bioplex to Biop format
# if -sepinst is set, write lines to different files depending on which 
#   instrument was used.  Further, separate any run ending in a letter

# output columns:
# search_id  bait_ref  prey_ref  total_peptides  sample_date  ms_inst_run_id  rep  plate
my %opts = getCommandLineOptions();

{
    my $outputFile = $opts{output}; # the main input file, but names...
    my $summaryFile = $opts{summary};
    my $out = $opts{out};

    # summary{bait}{A} = ms_inst_run_id, effectively
    my %summary;
    readSummary(\%summary, $summaryFile);

    my %add = (A => 1_000_000_000, B => 2_000_000_000 );

    my @inCols = qw(search_id Rep     bait_ref    prey_ref    total_peptides);
    my @outCols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date    ms_inst_run_id rep  plate);
    my $dummyDate = "1970-01-01"; # the epoch

    my $header = "# $0 munged $outputFile and $summaryFile\n".
	"# addA = $add{A}; addB = $add{B}\n".(join "\t", @outCols)."\n";
    
    open my $IN, "<", $outputFile or die "can't read $outputFile. $!";
    my %OUT;
    my $whichOut = 'a';
    if ($opts{sepinst}) {
	my $thisOut = "$out.outOfOrder";
	open $OUT{outOfOrder}, ">", $thisOut 
	    or die "can't write to $thisOut. $!";
	print { $OUT{outOfOrder} } $header;
    } else {
	open $OUT{$whichOut}, ">", $out or die "can't write to $out. $!";
	print { $OUT{$whichOut} } $header;	
    }
    while (<$IN>) {
	next if /^#/;
	next if /^APMS/;
	chomp;
	my @spl = split;
	my %row = map { $inCols[$_] => $spl[$_] } 0..$#spl;

	my $rep = $row{Rep};
	my $bait = $row{bait_ref};
	
	$row{search_id} += $add{$rep};
	$row{sample_date} = $dummyDate;
	my $msID = $summary{$bait}{$rep}{ms_inst_run_id} //
	    die "can't find summary{$bait}{$rep}";
	$row{ms_inst_run_id} = $msID;
	$row{rep} = $rep;
	$row{plate} = $summary{$bait}{$rep}{plate};
	if ($opts{sepinst}) {
	    if ($msID =~ /[a-z]+$/) {
		$whichOut = 'outOfOrder';
	    } else {
		$msID =~ /(^[a-z]+)/ // die "can't parse ".Dumper(\%row);
		$whichOut = $1;
		if (!exists $OUT{$whichOut}) {
		    my $thisOut = "$out.$whichOut";
		    open $OUT{$whichOut}, ">", $thisOut 
			or die "can't write to $thisOut. $!";
		    print { $OUT{$whichOut} } $header;
		}
	    }
	}
	say { $OUT{$whichOut} } join "\t", map { $row{$_} } @outCols;
    }
    close $IN; 
    close $_ for values %OUT; 
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	output => '/home/kli3/proj/Interactome/data/BioPlex2.0Nature/gpsRes4_output.tsv',
	summary => '/home/kli3/proj/Interactome/data/BioPlex2.0Nature/gpsRes4_runSummary.tsv'
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -out output < $defaultString -sepinst > \n";

    my %opts = ();
    GetOptions(\%opts, "out=s", "output=s", "summary=s", "sepinst");
    die $usage unless exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{output});
    checkExist('f', $opts{summary});

    return %opts;
}

# ret{bait}{A} = { ms_inst_run_id => $RAWFile, plate => $Plate }
# ret{bait}{B} = { ... }
sub readSummary {
    my ($ret, $in) = @_;
 
    my @cols = qw(BaitGeneID Plate Replicate RAWFile);
    my @read;
    readColsRef(\@read, $in, \@cols, undef, "\t");

    my $i=0;
    for my $row (@read) {
	$ret->{$row->{BaitGeneID}}{$row->{Replicate}} = {
	    ms_inst_run_id => $row->{RAWFile},
	    plate => $row->{Plate} 
	};
    }

    return;
}

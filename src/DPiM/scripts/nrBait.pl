#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#use Storable;
use DateTime;
#use DateTime::Format::Strptime;
use HomeBrew::IO qw(checkExist readColsHashRef readColsRef);
use DpimLib qw(getLineDP4APMS getLineBiopAPMS getLineHyperspecAPMS);

# for each bait, pick the replicate with the highest bait peptides

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $baitlessFile = $opts{baitless};
    my $curatedFile = $opts{curated};
    my $badYearStart = $opts{badyearstart};
    my $badYearEnd = $opts{badyearend};
    my $logFile = $opts{log};
   
    
    my %sort;
    # sort{bait}{search_id} = [ {prey_ref=>'fbgn', 'total_peptides'=>n},...]
    
    my @cols = qw(prey_ref total_peptides sample_date);
    
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my %row;
    say "parse";
    my $reader = \&getLineDP4APMS;
    if ($opts{mode} eq 'fourcols') {
	$reader = \&getLineHyperspecAPMS;
	@cols = qw(prey_ref total_peptides);
    } elsif ($opts{mode} eq 'human') {
	$reader = \&getLineBiopAPMS;
    }
    while ($reader->(\%row, $IN)) {
	my $bait = $row{bait_ref};
	my $id = $row{search_id};
	$sort{$bait}{$id} //= [];

	my $elem = {map { $_ => $row{$_} } @cols};
	push @{ $sort{$bait}{$id} }, $elem;
    }
    close $IN;

    if (defined $baitlessFile) {
	filterRebranded(\%sort, $baitlessFile);
    }
    %sort = removeTimeInterval(\%sort, $badYearStart, $badYearEnd)
	if exists $opts{badyear};
    #say 0+ keys %sort;
    
    my $LOG;
    if (defined $logFile) {
	open $LOG, ">", $logFile or die "can't write to $logFile. $!";
    }
    
    my %bestSID;
    for my $bait (keys %sort) {
	my %baitPep;
	for my $sid (keys %{ $sort{$bait} }) {
	    $baitPep{$sid} = 0;
	    for my $prey (@{ $sort{$bait}{$sid} }) {
		if ($prey->{prey_ref} eq $bait) {
		    $baitPep{$sid} = $prey->{total_peptides};
		    last;
		}
	    }
	}
	my @sort = sort {$baitPep{$a} <=> $baitPep{$b}} keys %baitPep;

	if (defined $logFile && @sort > 1) { 
	    say $bait;
	    for my $sid (@sort) { 
		say $LOG "\t$sid\t$baitPep{$sid}";
	    }
	    print $LOG "\n";
	}
	$bestSID{$bait} = $sort[-1]; 
    }
    close $LOG;

    if (-e $curatedFile) {
	curate(\%bestSID, $curatedFile);
    } else {
	warn "$0 can't find '$curatedFile', skipping\n";
    }
    
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    for my $bait (keys %sort) {
	my $sid = $bestSID{$bait} // die "no good bait for $bait?";
	for my $prey (@{ $sort{$bait}{$sid} }) {
	    say $OUT join("\t", $sid, $bait, $prey->{prey_ref},
			  $prey->{total_peptides});
	}
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    my @modes = qw(fly human fourcols);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    #my $start = DateTime->new(year=> 2009, month=> 4, day=> 1);
    #my $end   = DateTime->new(year=> 2010, month=> 6, day=> 1);

    my %defaults = (
	badyearstart => ("2009-04-01"),
	badyearend => ("2010-06-01"),
	curated => '/home/glocke/DPiM/augRemap/apmsData/annotation/RAO_Replicate_Winners_121916.short.short.txt',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString -baitless ".
	" apms.statsBySID $defaultString -badyear -log verbose.txt >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "baitless=s", 
	       "badyearstart=s", "badyearend=s", "curated=s", "badyear", 
	       "log=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};

    die "can't use -badyear with fourcol mode" if $opts{mode} eq 'fourcols' &&
	exists $opts{badyear};
    
    checkExist('f', $opts{in});
    checkExist('f', $opts{baitless}) if exists $opts{baitless};

    return %opts;
}

# prefer experiments where the bait was present to experiments where the bait
#   was rebranded 
sub filterRebranded {
    my ($sorted, $baitlessFile) = @_;

    my %pulldown;
    readColsHashRef(\%pulldown, $baitlessFile, [qw(search_id bait_peptides)]);
    for my $bait (keys %$sorted) {
	my %baitless;
	my ($foundBait, $foundBaitless) = (undef, undef);
	for my $sid (keys %{$sorted->{$bait}}) {
	    if ($pulldown{$sid} == 0) {
		$baitless{$sid} = 1;
		$foundBaitless = 1;
	    } else {
		$baitless{$sid} = undef;
		$foundBait = 1;
	    }
	}
	if ($foundBait && $foundBaitless) {
	    for my $sid (keys %baitless) {
		delete $sorted->{$bait}{$sid} if defined $baitless{$sid};
	    }
	}
    }

    return;
}

# remove experiments falling within a time interval
# unless $harsh is defined, only remove if additional experiments exist for 
#   that bait
sub removeTimeInterval {
    my ($sorted, $start, $end, $harsh) = @_;

    my $parseDate = DateTime::Format::Strptime->new(
	pattern   => '%Y-%m-%d');
    my $startDate = $parseDate->parse_datetime($start);
    my $endDate = $parseDate->parse_datetime($end);

    my $inInterval = sub { # return true if this date is in the interval
	my ($dateString) = @_;
	my $date = $parseDate->parse_datetime($dateString);
	die Dumper($dateString) if ! defined $date;
	    ##say Dumper($date, $startDate, $endDate);
	if (DateTime->compare($date, $startDate) >= 0 &&
	    DateTime->compare($date, $endDate) < 0) {
	    return 1;
	} 
	return undef;
    };

    my @drop; # drop[$i] = {bait=>$bait, sid=>$sid};
    for my $bait (keys %$sorted) {
	my @sids = keys %{$sorted->{$bait}};
	next unless 1 < @sids || $harsh;
	my @badDates = grep { 
	    $inInterval->($sorted->{$bait}{$_}[0]{sample_date});
	} @sids;
	next unless @badDates > 0 && (@badDates < @sids || $harsh);
	if (0) {
	    say "$bait";
	    for my $sid (@sids) {
		say "\t$sid\t$sorted->{$bait}{$sid}[0]{sample_date}";
	    }
	    exit;
	}
	push @drop, {bait => $bait, sid => $_} for @badDates;
	
    }

    for my $d (@drop) {
	delete $sorted->{$d->{bait}}{$d->{sid}};
	delete $sorted->{$d->{bait}} if 0 == keys %{ $sorted->{$d->{bait}}};
    }

    return %$sorted;
}


sub curate {
    my ($bestSID, $curatedFile) = @_;

    my @cols = qw(pickMe  search_id       bait_ref);
    my @curation;
    readColsRef(\@curation, $curatedFile, \@cols);
    for my $row (@curation) {
	next unless $row->{pickMe} == 1;
	$bestSID->{$row->{bait_ref}} = $row->{search_id};
    }
}
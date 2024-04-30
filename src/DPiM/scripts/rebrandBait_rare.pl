#!/usr/bin/env perl

use feature ':5.10';
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min);
use Data::Dumper;
use DateTime;
#use DateTime::Format::Strptime;
use HomeBrew::IO qw(checkExist readColsRef readColsHashRef);
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);

# for each experiment, find the prey meeting the following criteria
# * was available as a clone at sample_date time

# pick the lower fbgn in the case of a tie

my %opts = getCommandLineOptions();

{
    my $noBaitFile = $opts{nobait};
    my $out = $opts{out};
    my $freqFile = $opts{freq};
    my $cloneFile = $opts{clones};
    my $withBaitFile = $opts{withbait};

    my %freqs;
    readColsHashRef(\%freqs, $freqFile, [qw(fbgn Fraction)]);
    
    my %clones; # clones{fbgn} = earliest date
    readClones(\%clones, $cloneFile);

    
    my %apms;
    readAPMS(\%apms, $noBaitFile);
    # apms{search_id} = [ {all columns},...]

    rebrand(\%apms, \%freqs, \%clones);
	

    my $OUT;
    my @cols = qw( search_id bait_ref prey_ref total_peptides sample_date ms_inst_run_id );

    if (defined $withBaitFile) {
		open $OUT, ">", $out or die "can't write to $out. $!";
		say $OUT "# rebranded baits in $noBaitFile";
		close $OUT;
		my $cmd = "cat $withBaitFile >> $out";
		system($cmd);
		open $OUT, ">>", $out or die "can't write to $out. $!";
    } else {
		open $OUT, ">", $out or die "can't write to $out. $!";
		say $OUT "# rebranded baits in $noBaitFile";
		say $OUT join "\t", @cols;
    }

    for my $sid (sort {$a <=> $b} keys %apms) {
	for my $row (@{ $apms{$sid} }) {
	    say $OUT join "\t", map { $row->{$_} } @cols;
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
	freq => 'short.tscCut.protFreq',
	#clones => 'FH_Plate_Contents_By_Date.tidy.txt',
	clones => '/home/glocke/DPiM/augRemap/apmsData/annotation/FH_Plate_Contents_By_Date.tidy.updated.txt',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -nobait baitless.apms -out output < $defaultString".
	"-withbait apms > \n";

    my %opts = ();
    GetOptions(\%opts, "nobait=s", "out=s", "freq=s", "clones=s", "withbait=s");
    die $usage unless exists $opts{nobait} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{nobait});
    checkExist('f', $opts{freq});
    checkExist('f', $opts{clones});
    checkExist('f', $opts{withbait}) if exists $opts{withbait};

    return %opts;
}

sub readClones {
    my ($ret, $inFile) = @_;
    
    my @read;
    readColsRef(\@read, $inFile, [qw(Flybase_ID Earliest_Date_Received)], 
		undef, "\t");
    
    #my $parseDate = DateTime::Format::Strptime->new(
	#pattern   => '%m%d%y');
    for my $row (@read) {
	my $fbgn = $row->{Flybase_ID};
	my $dateString = $row->{Earliest_Date_Received};
	next if ! defined $fbgn && ! defined $dateString;
	#my $dateObj = $parseDate->parse_datetime($dateString);
	my $dateObj = $dateString;
	my @parts = split '', $dateObj;
	$dateObj = join( '', @parts[4,5,2,3,0,1]);
	$ret->{$fbgn} //= $dateObj;
	#$ret->{$fbgn} = $dateObj
	#    if 0 < DateTime->compare($ret->{$fbgn}, $dateObj);
	my $d =$ret->{$fbgn};
	@parts = split '', $d;
	$d = join( '', @parts[4,5,2,3,0,1]);
	$ret->{$fbgn} = $dateObj if $d < $dateObj;
	## if this date is earlier then use it
    }
    return;
}

sub readAPMS {
    my ($ret, $inFile) = @_;

    my @cols = qw( search_id bait_ref prey_ref total_peptides sample_date ms_inst_run_id );

    open my $IN, "<", $inFile or die "Can't read from $inFile. $!";
    my %row;
    say "parse";
    while (getLineDP4APMS(\%row, $IN)) {
	my $id = $row{search_id};
	$ret->{$id} //= [];

	push @{ $ret->{$id} }, { map { $_ => $row{$_} } @cols};
    }
    close $IN;
}


sub rebrand {
    my ($apms, $freq, $clones) = @_;

    #my $parseDate = DateTime::Format::Strptime->new(
	#pattern   => '%Y-%m-%d');

    for my $k (keys %$apms) {
		my $expt = $apms->{$k};
		#say "first bait = $expt->[0]{bait_ref}";
		my %maybes; # maybe{fbgn}=freq{fbgn} iff clones{fbgn} is defined
		#my $date = $parseDate->parse_datetime($expt->[0]{sample_date});
		my $date = $expt->[0]{sample_date};
		$date =~ s/-//g;
		#my %theseClones = map { $_ => 1 }
		#    grep { 1 > DateTime->compare($clones->{$_}, $date) } keys %$clones;
		my @v = values %$clones;
		my @nv = @v;
		my $count = 0;
		while ($count <= $#v) {
			my @parts = split '', $v[$count];
			push(@parts,2);
			push(@parts,0);
			$nv[$count] = join( '', @parts[6,7,0,1,4,5,2,3]);
			$count +=1;
			#print "$count\n";
		}
		my %theseClones = map { $_ => 1 }
			grep { $_ > $date } @nv;
		
		for my $row (@$expt) {
		    my $prey = $row->{prey_ref};
		    next unless exists $theseClones{$prey};
		    $maybes{$prey} = $freq->{$prey};
		}
		if (0 == keys %maybes) {
		    delete $apms->{$k};
		} else {
		    my $minFreq = min values %maybes;
		    my @newPrey = sort grep { $maybes{$_} == $minFreq } keys %maybes;
		    warn "experiment '$k' has more than one possible bait ( ".
			join(", ", @newPrey)." )\n" if @newPrey > 1;
		    $_->{bait_ref} = $newPrey[0] for @$expt;
		}
		#say "last bait = $expt->[0]{bait_ref}, maxTSC = $maxTSC; $maxPrey";
		#say "last bait = $expt->[0]{bait_ref}";
		#exit;
    }
}

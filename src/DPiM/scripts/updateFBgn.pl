#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readCols readColsHash readHeader);
use DpimLib qw(getLineDP4APMS);

# update old fbgn's to new ones

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $refMap = $opts{ref};
    my $mode = $opts{mode};
    my $baitMapFile = $opts{baitmap};
    
    if ($mode eq 'apms' || $mode eq 'dp1excel') {
	# repeat any prey with more than one match
	# select the "best" bait where necessary
	apmsUpdate($in, $out, $refMap, $baitMapFile);
    } elsif ($mode eq 'network') {
	# repeat any interaction with more than one match
	networkUpdate($in, $out, $refMap, $baitMapFile);
    } elsif ($mode eq 'rough') {
	# raw search/replace
	roughUpdate($in, $out, $refMap);
    } else {
	die "unknown mode '$mode'";
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    # mode:
    my @modes = qw(apms network rough dp1excel);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	ref=>'/home/glocke/DPiM/flybase/fbgn_annotation_ID_fb_2015_04.tsv',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString $defaultString ".
	"-findkeys -baitmap optinal.baitMap >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "ref=s", "findkeys", "baitmap=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{ref});
    checkExist('f', $opts{baitmap}) if exists $opts{baitmap};

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    return %opts;
}

# ret{old_fbgn} = new_fbgn
sub singleRefMap {
    my ($inFile, $baitMapFile) = @_;

    my %hashFormat;
    makeUpdateMap(\%hashFormat, $inFile);

    my %ret;
    for my $old (keys %hashFormat) {
	my @targets =  sort keys %{$hashFormat{$old}};
	if (1 == @targets) {
	    $ret{$old} = $targets[0];
	} else {
	    my @nonCG = grep { $hashFormat{$old}{$_} !~ /^CG\d+/ } @targets;
	    if (@nonCG > 0) {
		$ret{$old} = $nonCG[0];
	    } else {
		$ret{$old} = $targets[0];
	    }
	}
    }

    # if (defined $baitMapFile)... this has been taken care of by multiRefMap
    
    return %ret;
}

# ret{old_fbgn} = [all matches];
sub multiRefMap {
    my ($inFile, $baitMapFile) = @_;

    # this only slightly hacky.  makeUpdateMap was written for identifyBait.pl
    # since it works, but this script expects a different format, it's easy to
    # just reformat the output of makeUpdateMap
    my %hashFormat;
    makeUpdateMap(\%hashFormat, $inFile);
    
    my %ret;
    for my $old (keys %hashFormat) {
	$ret{$old} = [ keys %{ $hashFormat{$old} }];
    }
    
    if (defined $baitMapFile) {
	my %baitMap = readColsHash($baitMapFile, [qw(prevBait newBait)]);
	for my $old (keys %baitMap) {
	    $ret{$old} = [ $baitMap{$old} ];
	}
    }

    return %ret;
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


sub roughUpdate {
    my ($in, $out, $refMap) = @_;
    my %update = singleRefMap($refMap);
    say "\%update constructed";

    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";

    my $cnt = 0;
    while (my $line = <$IN>) {
	my $newLine = $line;
	while ($line =~ /(FBgn\d+)/g) {
	    my $old = $1;
	    my $new = $update{$old} // $old;
	    $newLine =~ s/$old/$new/;
	    $cnt++ if $new ne $old;
	}
	print $OUT $newLine;
    }
    close $IN;
    close $OUT;

    say "made $cnt substitutions";
}

sub apmsUpdate {
    my ($in, $out, $refMap) = @_;

    #say "peanuts!!";
    
    my %singleUpdate = singleRefMap($refMap);
    my %multiUpdate = multiRefMap($refMap); 
   
    my @cols = qw(search_id       bait_ref        prey_ref        total_peptides  sample_date     ms_inst_run_id);
    my $getLine = \&getLineDP4APMS;
    if ($opts{mode} eq 'dp1excel') {
	$getLine = \&getLineDP1ExcelAPMS;
	@cols = qw(search_id       bait_ref        prey_ref total_peptides)
    }
    
    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT "# sanitized $in according to $refMap";
    say $OUT join "\t", @cols;

    my %row;
    while ($getLine->(\%row, $IN)) {
	#say Dumper(\%row);
	$row{bait_ref} = $singleUpdate{$row{bait_ref}} // $row{bait_ref};
	my $cleanPrey = $row{prey_ref};
	$cleanPrey =~ s/reverse_//;
	if (exists $multiUpdate{$cleanPrey}) {
	    my $rev = '';
	    if ($row{prey_ref} =~ /reverse_/) {
		$rev = 'reverse_';
	    }
	    my $preys = $multiUpdate{$cleanPrey};
	    for my $p (@$preys) {
		$row{prey_ref} = $rev.$p;
		say $OUT join "\t", map { $row{$_} } @cols;
	    }
	} else {
	    say $OUT join "\t", map { $row{$_} } @cols;
	}
	die "found FBgn0003600" if $row{prey_ref} eq "FBgn0003600";
    }
    close $OUT;
}

sub networkUpdate {
    my ($in, $out, $refMap) = @_;

    my %multiUpdate = multiRefMap($refMap); 
   
    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";

    say $OUT "# sanitized $in according to $refMap";
    say $OUT join "\t", readHeader($in);

    my %row;
    while (my $line = <$IN>) {
	next unless $line =~ /FBgn/;
	#say Dumper(\%row);
	chomp $line;
	my ($old1, $old2, $score) = split /\t/, $line;
	if (exists $multiUpdate{$old1} && exists $multiUpdate{$old2}) {
	    my @fb1 = @{ $multiUpdate{$old1} };
	    my @fb2 = @{ $multiUpdate{$old2} };
	    for my $new1 (@fb1) {
		for my $new2 (@fb2) {
		    
		    say $OUT join "\t", ( sort ($new1, $new2) ), $score;
		}
	    }
	} elsif (exists $multiUpdate{$old1}) {
	    my @fb1 = @{ $multiUpdate{$old1} };
	    for my $new1 (@fb1) {
		say $OUT join "\t", (sort ($new1, $old2)), $score;
	    }	    
	} elsif (exists $multiUpdate{$old2}) {
	    my @fb2 = @{ $multiUpdate{$old2} };
	    for my $new2 (@fb2) {
		say $OUT join "\t", (sort ($old1, $new2)), $score;
	    }	    
	} else {
	    say $OUT $line;
	}

    }
    close $OUT;
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

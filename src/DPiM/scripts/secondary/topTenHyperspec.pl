#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use lib '/home/glocke/utilScripts/HomeBrew/lib';
use HomeBrew::IO qw(checkExist readList writeCols);

# look at the simulation output of hyperspec
# establish the frequency of each protein/interaction (node/edge)

my @modes = qw(edge protein);
my %modes = map {$_ => 1} @modes;

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $topN = $opts{n};


    my %rawTop; # rawTop{simNum} = [ {first, second, score}, {}, ...];
    my %nEdge;  # nEdge{simNum} = number of interactions 
    my %rawAll; # rawAll{protein_id} = # of appearances
    for my $f (readList($in)) {
	my @top;
	my $line;
	my $nRow=0;
	
	open my $IN, "<", $f or die "can't read from $f. $!";
	while (defined($line = readHS($IN))) {
	    $_ = $line;
	    my @spl = split;

	    if (@top < $topN) {
		push @top, {first => $spl[0], second => $spl[1], 
			    score => $spl[2]};
	    } elsif ($mode ne 'protein') {
		last;
	    }
	    $rawAll{$spl[0]}++;
	    $rawAll{$spl[1]}++;

	    $nRow++;
	}
	close $IN;

	my $simNum = getSimNum($f);
	say "finished $simNum";
	die "found $simNum again ($f)" if exists $rawTop{$simNum};
	$rawTop{$simNum} = \@top;
	$nEdge{$simNum} = $nRow;
    }
    
    my %topCounts = countUp(\%rawTop, $mode);
    writeResults($out, \%topCounts, \%nEdge);
    
    if ($mode eq 'protein') {
	my @proteins = sort keys %rawAll;
	
	my @data;
	push @data, \@proteins;
	push @data, [ map {$rawAll{$_}} @proteins ];
	my $header = "protein\tcount\n";
	my $format = "%s\t%d\n";
	writeCols("$out.raw", \@data, $header, $format);
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	n => 10,
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = join("/", @arr);

    my $usage = "usage: $0 -in sim.list -out report.tab < $modeString ".
	"$defaultString >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "n=i");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    return %opts;
}


sub readHS {
    my ($IN) = @_;

    return undef if eof($IN);

    my $goodLine = sub {
	return $_[0] =~ /^FBgn\d+\s+FBgn/;
    };
    
    my $line;
    do {
	$line = <$IN>;
    } while (!eof($IN) && ! $goodLine->($line));
    return undef if ! $goodLine->($line);

    $_ = $line;
    chomp;
    return $line;
}

sub getSimNum {
    my ($fileName) = @_;
    $fileName =~ /hyperspec.out.(\d+)/ or die "can't parse $fileName";
    return 0+ $1;
}

# ret{prot_id} = {count => number of appearnaces, sims=>"1,2,3"}
sub countUp {
    my ($rawCounts, $mode) = @_;

    my $adder = sub {
	my ($hash, $key, $sim) = @_;
	$hash->{$key} //= { count=> 0, sims=>'' };
	if (length($hash->{$key}->{sims})) {
	    $hash->{$key}->{sims}.=",$sim";
	} else {
	    $hash->{$key}->{sims} = $sim;
	}
	$hash->{$key}{count}++;
	return;
    };
    
    my $addProtein = sub {
	my ($hash, $prot1, $prot2, $sim) = @_;

	$adder->($hash, $prot1, $sim);
	$adder->($hash, $prot2, $sim);
	return;
    };
    my $addEdge = sub {
	my ($hash, $prot1, $prot2, $sim) = @_;

	my ($k1, $k2) = sort ($prot1, $prot2);
	$adder->($hash, $k1."_".$k2, $sim);
	return;
    };

    my $add;
    if ($mode eq 'protein') {
	$add = $addProtein;
    } elsif ($mode eq 'edge') {
	$add = $addEdge;
    } else {
	die "countUp doesn't recognize mode '$mode'";
    }
    
    my %ret = ();
    for my $sim (keys %$rawCounts) {
	for my $row (@{$rawCounts->{$sim}}) {
	    $add->(\%ret, $row->{first}, $row->{second}, $sim);
	}
    }
    
    for my $p (keys %ret) {
	my $simString = $ret{$p}{sims};
	my @spl = split /,/, $simString;
	my %uniq = map {$_=>1} @spl;
	$ret{$p}{nsims} = 0+ (keys %uniq);
    }

    return %ret;
}

sub writeResults {
    my ($out, $topCounts, $nEdge) = @_;

    open my $OUT, ">", $out or die "can't write to $out. $!";
    
    say $OUT join "\t", qw(protein count nsims sim);

    my $format = "%s\t%d\t%d\t%s\n";
    
    for my $prot (sort {$topCounts->{$a}{count} <=> $topCounts->{$b}{count}} 
		  keys %$topCounts)
    {
	
	printf $OUT $format, $prot, $topCounts->{$prot}{count},
	    $topCounts->{$prot}{nsims}, $topCounts->{$prot}{sims};
    }

    close $OUT;

    return;
}

#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef);
#use DpimLib qw(getLineAPMS);

# turn the output of hierarchy into a common cluster format
# example workflow:
#   $LOUVDIR/convert -i network -o -o network.bin -w network.weight
#   $LOUVDIR/community network.bin -w network.weight -l 1 > network.bin.louvain.out
#   $LOUVDIR/heirarchy network.bin.louvain.out # find max level
#   $LOUVDIR/heirarchy -l X network.bin.louvain.levelX


my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $keyFile = $opts{key};

    die "unimplemented mode $mode" unless $mode eq "mcl";

    my %clusters;
    open my $IN, "<", $in or die "Can't read from $in. $!";
    while (<$IN>) {
	chomp;
	next if /^0\s+0$/;
	my ($node, $cluster) = split;
	$clusters{$cluster} //= [];
	push @{$clusters{$cluster}}, $node;
    }

    my @order = sort {@{ $clusters{$b} } <=> @{ $clusters{$a} }} keys %clusters;

    my %nodeKey = readNodeKey($keyFile);
    
    open my $OUT, ">", $out or die "can't write to $out. $!";
    for my $cluster (@order) {
	my @nodes = @{ $clusters{$cluster} };
	for my $n (@nodes) {
	    next if $n == 0;
	    $n = $nodeKey{$n} // die "unrecognized node '$n'";
	}
	@nodes = sort @nodes;
	say $OUT join "\t", @nodes;
    }
    close $OUT;
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   


sub getCommandLineOptions {
    # mode:
    my @modes = qw(mcl cluster);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	key => '/home/glocke/DPiM/oldDpim/dpim3.1/louvain/dpim3.09-25-2015.nrBait.77.44.network.intkeys.log',
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;
    my $usage = "usage: $0 -in louvain.levelX -out output < $modeString ".
	"$defaultString >";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "key=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});
    checkExist('f', $opts{key});

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    return %opts;
}

sub readNodeKey {
    my $in = shift;
    
    my @cols = qw(nodeString nodeInt);
    my @read;
    readColsRef(\@read, $in, \@cols);
    return map { $_->{nodeInt} => $_->{nodeString} } @read;
}

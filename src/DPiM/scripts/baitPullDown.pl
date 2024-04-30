#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist writeCols);
use DpimLib qw(getLineAPMS getLineDP4APMS getLineDP4_1APMS getLineBiopAPMS getLineRawAPMS);

# for each $idCol, find the number of bait pulled down

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    my $mode = $opts{mode};
    my $idCol = $opts{idcol};
    
    my $reader = \&getLineDP4APMS;
    $reader = \&getLineDP4_1APMS if $mode eq 'logp';
    $reader = \&getLineBiopAPMS if $mode eq 'human';
    $reader = \&getLineRawAPMS if $mode eq 'raw';
    $reader = \&getLineAPMS if $mode eq 'clean';

    my %exp; # exp{sid}= {bait_ref=> $fbgn, bait_peptides=> $$, date}
    
    open my $IN, "<", $in or die "can't read from $in. $!";

    my %row;
    while ($reader->(\%row, $IN)) {
	my $id = $row{$idCol};
	$exp{$id} //= { bait_ref=>$row{bait_ref}, bait_peptides=>0, 
			date => $row{sample_date}, 
			ms_inst_run_id => $row{ms_inst_run_id} };
	$exp{$id}{bait_peptides}+= $row{total_peptides} 
	    if $row{bait_ref} eq $row{prey_ref};
	$exp{$id}{prey}{$row{prey_ref}} = 1 if $row{total_peptides} > 0;
	$exp{$id}{tsc}+=$row{total_peptides};
    }

    my @cols = qw(bait_ref bait_peptides tsc prey_count date ms_inst_run_id);
    unshift @cols, $idCol;
    my @sid;
    if ($idCol eq 'search_id') {
	@sid = sort {$a <=> $b} keys %exp;
    } else {
	@sid = sort keys %exp;
    }
    my @bait_ref = map { $exp{$_}{bait_ref} } @sid;
    my @bait_peptides = map { $exp{$_}{bait_peptides} } @sid;
    my @tsc = map { $exp{$_}{tsc} } @sid;
    my @prey_count = map { 0+ keys %{ $exp{$_}{prey} } } @sid;
    my @date = map { $exp{$_}{date} } @sid;
    my @rid = map { $exp{$_}{ms_inst_run_id} } @sid;

    my $header = join "\t", @cols;
    my $format = join "\t", qw(%s %s %d %d %d %s %s);
    writeCols($out, [\@sid, \@bait_ref, \@bait_peptides, \@tsc, \@prey_count,
		     \@date, \@rid], $header, $format);
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my @modes = qw(dp4 human logp raw clean);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	idcol => 'search_id',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output < $modeString $defaultString ".
	">\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "mode=s", "idcol=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    checkExist('f', $opts{in});

    return %opts;
}


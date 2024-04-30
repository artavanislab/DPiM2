#!/usr/bin/env perl

##############################################################################80

use v5.10; 
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Path qw(make_path);
use File::Spec;
use HomeBrew::IO qw(checkExist );
#use DpimLib qw(getLineAPMS getLineDP4APMS getLineRawAPMS );

my %opts = getCommandLineOptions();

{
    my $in = File::Spec->rel2abs($opts{in});
    my $baseDir = File::Spec->rel2abs($opts{dir});
    my $job = $opts{job};
    my $nSim = $opts{nsim};
    my $rt = $opts{rt};
    my $seed = $opts{seed};
    my $mode = $opts{mode};

    my $jdir = "$baseDir/json";
    my $qdir = "$baseDir/qdir";
    if (! -d $jdir) {
	make_path($jdir) or die "can't make $jdir";
    }
    if (! -d $qdir) {
	make_path($qdir) or die "can't make $qdir";
    }
    chdir($qdir);
    
    my $hyper = "module add boost; $ENV{DPCPP}/hyperspec --in $in";
    #my $hyper = "module add boost; ~/DPiM/cpp/hyperspec --in $in --protlen /home/glocke/DPiM/nsaf/dmel-all-translation-r6.09.aveLen.tsv";
    my $qsss = '/home/glocke/utilScripts/qsubWrap.pl';

    if ($mode eq 'human') {
	my $protLen = "/home/glocke/DPiM/human/nsaf/entrez2translatedLength.tsv";
	checkExist('f', $protLen);
	$hyper.= " --protlen $protLen";
    }
    
    if (! exists $opts{simonly}) {
	my $cmd .= "$hyper --out $jdir/$job.json";
	system("$qsss -cmd '$cmd' -job $job -rt '$rt'");
    }    

    exit if $nSim < 1;

    # run simulations
    $hyper.=" --sim 1";
    for my $n (1..$nSim) {
	my $thisJob = sprintf("$job\_sim%03d", $n);
	my $out = "$jdir/$thisJob.json";
	my $cmd = "module add boost; $hyper --out $out --seed $seed";
	system("$qsss -cmd '$cmd' -job $thisJob -rt '$rt'");
	sleep(1);
	$seed++;
    }
    exit;
}

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    my @modes = qw(fly human);
    my %modes = map {$_ => 1} @modes;
    my @arr = @modes;
    $arr[0] = "*".$arr[0]."*";
    my $modeString = "-mode ".join("/", @arr);

    my %defaults = (
	dir => "$ENV{PWD}/hyper",
	job => "hyper",
	nsim => 101,
	rt => '03:00:00',
	seed => time(),
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in apms.tsv <  $modeString $defaultString -simonly>\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "mode=s", "job=s", "dir=s", "nsim=i", "rt=s", 
	       "seed=i", "simonly");
    die $usage unless exists $opts{in};

    die "you must select one of the following modes: ", 
        join(", ",keys(%modes)), "\n" 
	if (exists $opts{mode} && ! exists $modes{$opts{mode}});
    $opts{mode} //= $modes[0];

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

#!/usr/bin/env perl

# wrap a command in a bash script and submit it via qsub

use strict;
use warnings;
use Getopt::Long;
use feature ':5.10'; 

my %opts = getCommandLineOptions();
my $cmd		= $opts{cmd};
my $job		= $opts{job};
my $rt		= $opts{rt};
my $modules	= $opts{modules};

{
    
    my $jobScript = $job.".qsub.bash";
    my $bashLoc = "/bin/bash";
    if ( -f $jobScript) {
	if (! defined $opts{rm}) {
	    print "script to submit, $jobScript, already exists (specify -rm to ignore this message). quitting\n";
	    exit -1;
	} else {
	    print "deleting job script ($jobScript) and making new one\n";
	}
    }

    open my $SCR, ">", $jobScript or die "can't write to $jobScript. $!\n";
    print $SCR "#!$bashLoc\n";
    print $SCR "#\$ -l h_rt=$rt";
    print $SCR '
#$ -V
#$ -cwd 
#$ -b y
';
    if (defined $modules) {
	my @spl = split /,/, $modules;
	say $SCR "module add all-collection";
	say $SCR "module add $_" for @spl;
    }
    
    print $SCR "$cmd\n";
    close $SCR;
    system("chmod 755 $jobScript");

    system("qsub -N $job $ENV{PWD}/$jobScript\n") unless (defined $opts{nosub});	
    
}
exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
    
    my %defaults = (
	rt => "00:59:59",
	);
    my $defaultString = join " ", map { "-$_ $defaults{$_}" } keys %defaults;
  
    my $usage = qq{usage: $0 -cmd "job.pl -opts" -job distance.dat < }.
	qq{$defaultString -modules mod1,mod2,... -rm -nosub>};

    my %opts;
    GetOptions(\%opts, "cmd=s", "job=s", "rt=s", "modules=s", "rm", "nosub");
    if (!exists $opts{cmd} || ! exists $opts{job}) {
	die "$usage\n";
    }

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    return %opts;
}

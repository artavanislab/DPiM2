#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readList);
use feature ':5.10'; 

# given a list of bash scripts, submit them to the queue, keeping a certain 
#   number in that queue at all times

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $nJobs = $opts{njob};
    my $logFile = $opts{log};
    my $addArgs = $opts{qsubargs};

    my %scripts;
    for my $s (readList($in)) {
	$s =~ /(.+)\.qsub\.bash$/ or die "can't find '.qsub.bash' in $s";
	my @spl = split /\//, $1;
	$scripts{pop(@spl)} = $s;
    }

    $nJobs = 10* (keys %scripts) if exists $opts{all};

    open my $LOG, ">", $logFile or die "can't write to $logFile. $!";

    while (keys %scripts) {
	while (queuedJobs(qstat()) < $nJobs && keys %scripts) {
	    my $job = (sort keys %scripts)[0];
	    my $cmd = "qsub -N $job $scripts{$job}";
	    $cmd.=" ".$addArgs if $addArgs;
	   
	    say $cmd;
	    say $LOG $cmd;
	    system($cmd);
	    delete $scripts{$job};
	    #TO DO: check if jobs completed properly, i.e. whether the output
	    # exists
	}

	sleep(20) while (queuedJobs(qstat()) >= $nJobs);
    }


}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	njob => 60,
	log => 'listQ.log',
	);

    my $usage = "usage: $0 -in bash.list < -njob $defaults{njob} -log ".
	"$defaults{log} -all -qsubargs \"-l etc etc\"> \n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "njob=i", "log=s", "all", "qsubargs=s");
    die $usage unless exists $opts{in};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# for now, just get job status and number for each job
# allow multiple jobs with same name to over-write each other cuz they're not
# mine
sub qstat {
    my $qstatTxt = `qstat -f`;
    my @jobs = split /Job Id: /, $qstatTxt;
    
    my %ret;
    for my $j (@jobs) {
	next if $j eq '';
	my @lines = split /\n/, $j;

	my @jobOwner = grep {/Job_Owner/} @lines;
	next unless $jobOwner[0] =~ /glocke/;

	$lines[0] =~ /(\d+)\.nikka/ 
	    or die "can't find job number in $lines[0]";
	my $id = $1;
	my @jobState = grep { /job_state/ } @lines;
	die "can't find job state in $id" if @jobState != 1;
	my @jobName =  grep { /Job_Name/ } @lines;
	die "can't find job name in $id" if @jobName != 1;
	$jobState[0] =~ /job_state = (.+)/ or die "can't parse $jobState[0]";
	my $state = $1;
	$jobName[0] =~ /Job_Name = (.+)/ or die "can't parse $jobName[0]";
	my $name = $1;
	#say "$jobName[0] -> $name";
	$ret{$name} = { id => $id, state => $state };
    }
    return %ret;
}

# how many of my jobs have status 'R' or 'Q' ?
sub queuedJobs {
    my %qstat = @_;

    my $ret = 0;

    for my $info (values %qstat) {
	$ret++ if $info->{state} eq 'R' || $info->{state} eq 'Q';
    }

    return $ret;
}

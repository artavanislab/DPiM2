#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use DateTime;
use DateTime::Format::Strptime;
use HomeBrew::IO qw(checkExist);
use DpimLib qw(getLineDP4APMS);

# make a matrix with prey per column and bait per row
# it will be large
# -perbait produces separate files for each bait

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};

    say "parse input";
    my (%allAPMS, %meta, %prey);
    parseAPMS(\%allAPMS, \%meta, \%prey, $in); 
    

    say "report";
    if (exists $opts{perbait}) {
	my %baits  = map { $meta{$_}{bait} => 1 } keys %meta;
	my @spl = split /\./, $out;
	my $last = pop @spl;
	my $first = $out;
	if (@spl) {
	    $first = join ".", @spl;
	} else {
	    undef $last;
	}
	for my $bait (keys %baits) {
	    my @expOrder = sort grep { ($meta{$_}{bait} // die "can't find meta{$_}{bait}") eq $bait } keys %meta;
	    if (@expOrder > 1) {
		@expOrder = sort { $meta{$a}{time} <=> $meta{$b}{time} } @expOrder;
	    }
	    my %theseAPMS = map { $_ => $allAPMS{$_} } @expOrder;
	    my %prey2;
	    for my $exp (@expOrder) {
		$prey2{$_} += $theseAPMS{$exp}{$_} 
		    for keys %{ $theseAPMS{$exp} };
	    }
	    #@die Dumper(\%prey2);
	    my @preyOrder = sort {$prey2{$b} <=> $prey2{$a}} keys %prey2;
	    my $thisOut = "$first.$bait";
	    $thisOut = "$first.$bait.$last" if $last;
	    writer($thisOut, \@expOrder, \@preyOrder, \%meta, \%theseAPMS);
	}
    } else {
	my @expOrder = sort keys %meta;
	@expOrder = sort { $meta{$a}{time} <=> $meta{$b}{time} } @expOrder;
	my @preyOrder = sort {$prey{$b} <=> $prey{$a}} keys %prey;
	writer($out, \@expOrder, \@preyOrder, \%meta, \%allAPMS);
    }
    
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in in.apms -out out.mat < -perbait >\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s", "perbait");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# $ret->{$search_id}{$prey}=$tsc
# $meta->{$search_id} = {ms_inst_run_id, bait, date, time}
# $prey->{$fbgn} = $tsc == total number of peptides found for this prey
#    defined only for prey that were actually found, naturally
sub parseAPMS {
    my ($ret, $meta, $prey, $in) = @_;

    my $parseDate = DateTime::Format::Strptime->new(
	pattern   => '%Y-%m-%d');
    
    
    open my $IN, "<", $in or die "can't read $in. $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN)) {
	my $sid = $row{search_id};
	$ret->{$sid}{$row{prey_ref}} = 
	    $row{total_peptides};
	$meta->{$sid} //= {
	    ms_inst_run_id => $row{ms_inst_run_id},
	    bait => $row{bait_ref}, 
	    date=>$row{sample_date}, 
	    time => $parseDate->parse_datetime($row{sample_date})
	};
	die Dumper(\%row) unless defined $meta->{$sid}{time};
	$prey->{$row{prey_ref}}+=$row{total_peptides};
    }
    close $IN;
    return;
}

sub writer {
    my ($out, $expOrder, $preyOrder, $meta, $allAPMS, $bait) = @_;
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 matrix showing TSC for each prey in each experiment";
    say $OUT "# showing only experiments with bait '$bait'" if defined $bait;
    say $OUT join "\t", "search_id", @$expOrder;
    say $OUT join "\t", "ms_inst_run_id", map { $meta->{$_}{ms_inst_run_id} } @$expOrder;
    say $OUT join "\t", "date", map { $meta->{$_}{date} } @$expOrder;
    say $OUT join "\t", "bait", map { $meta->{$_}{bait} } @$expOrder;
    for my $prey (@$preyOrder) {
	say $OUT join "\t", $prey, map { $allAPMS->{$_}{$prey} // 0 } @$expOrder;
    }
    close $OUT;
    

}

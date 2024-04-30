#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist readColsRef writeCols);
use Math::GSL::CDF qw(gsl_cdf_hypergeometric_P);
use Statistics::R;
use File::stat;
use Storable;

# do multiple hypergeometric tests to see if any GO terms are under or over
# represented

# unless -disjoint is set, this assumes that the smaller hist tests a subset
# of the larger hist.  (largeness is assessed automatically)

my %opts = getCommandLineOptions();

{
    my $in1 = $opts{in1};
    my $in2 = $opts{in2};
    my $out = $opts{out};
    my $out2 = $opts{out2};
    my $alpha = $opts{alpha};
    my $goNameFile = $opts{goname};
    
    my $total1 = countProteins($in1);
    my $total2 = countProteins($in2);
    #die ($total1, $total2);
    if ($total1 < $total2) {
	# in1 is the larger of the two
	($in1, $in2) = ($in2, $in1);
	($total1, $total2) = ($total2, $total1);
	say "$in1, $in2";
    }
    my %hist1 = readHist($in1);
    my %hist2 = readHist($in2);

    my @terms = grep { exists $hist1{$_} } keys %hist2;
    my $nTerms = 0+ @terms;
    
    my $goNameMap = retrieve($goNameFile);

    my @results;
    for my $t (@terms) {
	my $one = $hist1{$t};
	my $two = $hist2{$t};
	my $expect = $total2 * $one / $total1;

	my $name = $goNameMap->{$t}{name} // "UNKNOWN";

	my %row = (term => $t, one => $one, two => $two, expect => $expect, 
		   total1 => $total1, total2 => $total2, name => '"'.$name.'"');
	
	my ($found, $white, $black, $outOf) = ($two, $one, $total1 - $one,
					       $total2);
	if ($two < $expect) {
	    $row{test} = '-'; 
	} else {
	    $row{test} = '+';
	    ($white, $black) = ($black, $white);
	    $found = $outOf - $found;
	}
	#say "$t\t$name:\tgsl_cdf_hypergeometric_P($found, $white, $black, $outOf)";
	my $p = gsl_cdf_hypergeometric_P($found, $white, $black, $outOf);
	$row{p} = $p;
	$row{q} = $p * $nTerms;

	push @results, \%row;
    }
    #hyperGeometric();
    #my $R = Statistics::R->new();
    #$R->startR;
    

    
    @results = sort {$a->{p} <=> $b->{p}} @results;
    
    my @cols = qw( term p q test one two expect total1 total2 name);
    my @data;
    for my $c (@cols) {
	push @data, [map {$_->{$c}} @results];
    }
    
    my $header = join "\t", 'rank', @cols;
    my $format = join "\t", qw( %d %s %.3e %.3e %s %d %d %.2f %d %d %s );
	
    my $preComments = "# test = '-' is a test for extreme ".
	"underrepresentation\n";
    $preComments .= "# test = '+' is a test for extreme overrepresentation\n";
    $preComments .= "# q = p * number of tests\n";
    
    writeCols($out, \@data, $header, $format, $preComments, 1);

    if (defined $out2) {
	my @sigResults = grep { $_->{q} <= $alpha } @results;
	my @data;
	for my $c (@cols) {
	    push @data, [map {$_->{$c}} @sigResults];
	}
	$preComments .= "# selecting rows with q < $alpha";
	writeCols($out2, \@data, $header, $format, $preComments, 1);	
    }
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	alpha => 0.05,
	goname => '/home/glocke/DPiM/flybase/goMap.storable',
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in1 hist.tab -in2 hist.tab -out output < -disjoint".
	" -out2 qLTalpha.out $defaultString > \n";

    my %opts = ();
    GetOptions(\%opts, "in1=s", "in2=s", "out=s", "disjoint", "out2=s", 
	"alpha=f", "goname=s");
    die $usage unless exists $opts{in1} && exists $opts{in2} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in1});
    checkExist('f', $opts{in2});
    checkExist('f', $opts{goname});

    return %opts;
}

sub readHist {
    my ($in) = @_;
    my @cols = qw(term count);
    my @read;
    readColsRef(\@read, $in, \@cols);
    return map { $_->{term} => $_->{count} } @read;
}

sub countProteins {
    my ($in) = @_;

    open my $IN, "<", $in or die "Can't read from $in. $!";
    my $line;
    do {
	$line = <$IN> // die "can't count proteins in $in (eof)";
    } until ($line =~ /proteins found = (\d+)/ || 
	     $line =~ /total tests = (\d+)/);
    die "can't count proteins in $in (not eof)" unless defined $1 && $1 > 0;

    return $1;
}

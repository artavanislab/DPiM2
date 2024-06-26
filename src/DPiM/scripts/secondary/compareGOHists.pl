#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(each_array);
use HomeBrew::IO qw(checkExist readColsRef writeCols);
use Math::GSL::CDF qw(gsl_cdf_hypergeometric_P);
use Statistics::R;
use File::stat;
use Storable;

# do fisher tests to see if any GO terms are under- or overrepresented

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
	
	if ($two < $expect) {
	    $row{updown} = '-'; 
	} else {
	    $row{updown} = '+';
	}
	push @results, \%row;
	##die Dumper(\@results);
    }
    fisherTest(\@results);
    

    
    @results = sort {$a->{p} <=> $b->{p}} @results;
    
    my @cols = qw( term p q updown one two expect total1 total2 name);
    my @data;
    for my $c (@cols) {
	push @data, [map {$_->{$c}} @results];
    }
    
    my $header = join "\t", 'rank', @cols;
    my $format = join "\t", qw( %d %s %.3e %.3e %s %d %d %.2f %d %d %s );
	
    my $preComments = "# $0 compared $in1 (superset) to $in2 (subset)\n";
    $preComments = "# $0 compared $in1 to $in2 (disjoint sets)\n" 
	if exists $opts{disjoint};
    $preComments .= "# q = holm's correction\n";
    
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

# compute (two-tailed) fisher tests for every term
# then use holm's correction for "family-wise error rate" (multiple testing)
sub fisherTest {
    my ($results) = @_;
    
    my $R = Statistics::R->new();
    $R->startR;

    for my $term (@$results) {
	my ($white1, $white2, $black1, $black2);
	if (exists $opts{disjoint}) {
	    $white1 = $term->{one};
	    $white2 = $term->{two};
	    $black1 = $term->{total1} - $white1;
	    $black2 = $term->{total2} - $white2;
	} else {
	    ## 1 is superset, 2 is subset
	    $white2 = $term->{two};
	    $white1 = $term->{one} - $white2;
	    $black2 = $term->{total2} - $white2;
	    $black1 = $term->{total1} - $white1 - $black2;
	}
	$R->send(qq'tabl <- matrix(c($white1, $black1, $white2, $black2), nrow=2)');
	$R->send(qq'ft <- fisher.test(tabl)');
	$R->send(q'ft$p.value');
	my $pString = $R->read;
	$pString =~ s/\[1\] //;
	$term->{p} = $pString;
	#die Dumper($term);
    }
    my $vecString = join ",", map { $_->{p} } @$results;
    $R->send("p <- c($vecString)");
    $R->send(qq'q <- p.adjust(p, method="holm")');
    $R->send(qq'paste(q, collapse=",")');
    my $pString = $R->read;
    my $x = $pString;
    $pString =~ s/\[1\] //;
    $pString =~ s/"//g;
    my @q = split ",", $pString;

    my $it = each_array(@$results, @q);
    while (my ($term, $q) = $it->()) {
	$term->{q} = $q;
    }    
    return;
}

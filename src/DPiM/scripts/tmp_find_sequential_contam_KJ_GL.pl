#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Math::GSL::CDF qw(gsl_cdf_gaussian_Q);
use DpimLib qw(getLineDP4APMS);
use HomeBrew::IO qw(checkExist);

#my $SCRDIR = '/home/kli3/proj/Interactome/data/Fly_DPIM_data/DPIM3/newData';

my $outFile;
if (@ARGV != 4) {
    die "usage: $0 <fbgn_id2_2col.txt> <apms_input> <apms.commonContam> <out>\n";
}
$outFile = pop(@ARGV);
checkExist('f', $_) for @ARGV;

say "parse $ARGV[0]";
open my $IN, "<", $ARGV[0] or die "can't read from $ARGV[0]. $!";

my %fbgn2symbol;
while (my $line = <$IN>) {
    chomp($line);
    next if ($line =~ /^clone_name/); # skip header
    #($fbgn, $symbol, $cg_id, $id2, $clone_name, $avg_exp_tryptic) = 
    # split(/\s/, $line);
    my ($fbgn, $symbol) = split(/\s/, $line);
    $fbgn2symbol{$fbgn} = $symbol;
}
$fbgn2symbol{"FBgn0000000"} = "EmpVecCtrl";
$fbgn2symbol{"GFP"} = "GFPCtrl";
$fbgn2symbol{"turboGFP"} = "GFPCtrl";
$fbgn2symbol{"eGFP"} = "GFPCtrl";
$fbgn2symbol{"tGFP"} = "GFPCtrl";

close($IN);

my (%ms_total_by_inst, %ms_bait_by_inst, %ms_prey, %ms_prey_tsc_arr);

say "parse $ARGV[1]";
open $IN, "<", $ARGV[1] or die "Can't open $ARGV[1]. $!";
my %row;
while(getLineDP4APMS(\%row, $IN)){
    my $inst_id = $row{ms_inst_run_id};
    my $bait_ref = $row{bait_ref};
    my $prey_ref = $row{prey_ref};
    my $tot_pep = $row{total_peptides};
    $ms_total_by_inst{$inst_id}{$prey_ref} += $tot_pep;
    $ms_bait_by_inst{$inst_id} = $bait_ref;
    $ms_prey{$prey_ref} += $tot_pep;
}
close($IN);

# compute the mean and std of TSCs
foreach my $inst (keys %ms_total_by_inst) {
    foreach my $prey (keys %{ $ms_total_by_inst{$inst} }) {
	$ms_prey_tsc_arr{$prey} //= [];
	push @{ $ms_prey_tsc_arr{$prey} }, $ms_total_by_inst{$inst}{$prey};
    }
}

say "parse $ARGV[2]";
my (%common_contam, %common_contam_tsc);
#open($IN, "common_contaminants_nrtap_all.out");
open $IN, "<", $ARGV[2] or die "Can't open $ARGV[2]. $!";
while (my $line = <$IN>) {
    chomp($line);
    next if ($line =~ /Fraction/); # skip header
    my ($fbgn, $cnt, $frac, $mean_tsc) = split(/\t/, $line);
    $common_contam{$fbgn} = $frac;
    $common_contam_tsc{$fbgn} = $mean_tsc;
    #$common_contam_cnt{$fbgn} = $cnt;
}
close($IN);


my @sorted_inst = sort {$a <=> $b} keys %ms_total_by_inst;
my @sorted_prey = sort keys %ms_prey;

####print scalar(@sorted_inst)."\t".scalar(@sorted_prey)."\n"; exit(0);
####### CONSTRUCT SEQUENTIAL RUN X PREY MATRIX ########
say "Construct run x prey matrix";
my @inst_by_prey_ar;
for my $i (0..$#sorted_inst) {
    for my $j (0..$#sorted_prey) {
	$inst_by_prey_ar[$i][$j] = 
	    $ms_total_by_inst{ $sorted_inst[$i] }{ $sorted_prey[$j] } // 0;
    }
}

my $minTSC = 20;
######## LOOK FOR MONOTONICALLY DECREASING SEQUENTIAL PREY TSCS ########
#### Allow some leeway for washes/standards
#### Sequences have to start with at least $minTSC TSCs
say "Seek contaminant signuatures";
my %monotonic_seq_hash;
#for my $j (0..10) {
for my $j (0..$#sorted_prey) {
    ##say $j;
    my $prev_inst_tot = $inst_by_prey_ar[$#sorted_inst][$j];
    my $monotonic_seq = 0;
    my $monotonic_seq_last = $#sorted_inst;
    for (my $i = @sorted_inst-2; $i >= 0; $i--) {
	##say "\t$i" if $i % 100 == 0;
	####print $inst_by_prey_ar[$i][$j]."\t"; ####
	if ($inst_by_prey_ar[$i][$j] > $prev_inst_tot &&
	    ($sorted_inst[$i+1] - $sorted_inst[$i]) <= 4 &&
	    $ms_bait_by_inst{ $sorted_inst[$i] } ne 
	    $ms_bait_by_inst{ $sorted_inst[$i+1] }) 
	{
	    $monotonic_seq = 1;
	    if ($prev_inst_tot == 0) { 
		# avoid sequences ending with and including 0
		$monotonic_seq_last = $i;
	    }
	} elsif ($prev_inst_tot >= $minTSC && ($monotonic_seq_last > ($i+1))) {
	    my @inst_list = @sorted_inst[($i+1) .. $monotonic_seq_last];
	    my ($inst_line, $pep_sum) = 
		checkSeq(join(",", @inst_list), $sorted_prey[$j]);
	    @inst_list = split(/,/, $inst_line);
	    
	    ####print "\t".($i+1)." - ".$monotonic_seq_last."\t"; ####
	    ####print "$i|$monotonic_seq_last|\t|".join(",", @inst_list)."|\t|".$sorted_prey[$j]."|\n"; ####

	    if (scalar(@inst_list) > 1) {
		$monotonic_seq_hash{ $sorted_prey[$j]."\t".
					 join(",", @inst_list) } = $pep_sum;
	    }
	    $monotonic_seq_last = $i;
	} else {
	    ####print "zeroing\t".($i+1)." - ".$monotonic_seq_last."\t"; ####
	    for my $each (($i+1)..$monotonic_seq_last) {
		$inst_by_prey_ar[$each][$j] = 0;
	    }
	    $monotonic_seq = 0;
	    $monotonic_seq_last = $i;
	}
	$prev_inst_tot = $inst_by_prey_ar[$i][$j];
	####print "\n"; ####
    }
    ####exit(0); ####
}
    
my $DEBUG = undef;
if ($DEBUG) {
    for my $i (0..4) {
	print $sorted_inst[$i];
	for my $j (0..20) {
	    print "\t".$inst_by_prey_ar[$i][$j];
	}
	print "\n";
    }
}

open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
my $cnt = 0;
foreach my $one (sort {$monotonic_seq_hash{$a} <=> $monotonic_seq_hash{$b}} 
		 keys %monotonic_seq_hash) 
{
    $cnt++;
    my ($prey, $inst_line) = split(/\t/, $one);
    my @inst_list = split(/,/, $inst_line);

    printf $OUT "\t%0.2e\t%s\t%s\t%0.4f\t%0.2f\n"
	, $monotonic_seq_hash{$one}, $prey, $fbgn2symbol{$prey}//$prey
	, $common_contam{$prey}, $common_contam_tsc{$prey};
    for my $inst (sort {$a <=> $b} @inst_list) {
	#say "inst: ", isTrue($inst);
	#say "ms_bait_by_inst{$inst}: ", isTrue($inst);
	#say "fbgn2symbol{ $ms_bait_by_inst{$inst} }: ", isTrue($fbgn2symbol{ $ms_bait_by_inst{$inst} });
	#say "ms_total_by_inst{$inst}{$prey}: ", isTrue($ms_total_by_inst{$inst}{$prey});

	printf $OUT "%s\t%s\t%s\t%d\n", $inst, $ms_bait_by_inst{$inst}
	, $fbgn2symbol{ $ms_bait_by_inst{$inst} } // $ms_bait_by_inst{$inst}
	, $ms_total_by_inst{$inst}{$prey};
    }
    print $OUT "\n";
    my @bait_list = map{ $ms_bait_by_inst{$_} } @inst_list;
    print join("\t", @bait_list)."\n";
}
close($OUT);


sub checkSeq {
     my $instr_line = $_[0];
     my $prey = $_[1];
     my @instr_list = split(/,/, $instr_line);
     ####print "JMDEBUG: |$prey|$instr_line|\n"; ####

     my @sorted_instr_list = sort {$a <=> $b} @instr_list;

     my @tot_pep_ar = ();
     for my $k (0..$#sorted_instr_list) {
	 push @tot_pep_ar, $ms_total_by_inst{ $sorted_instr_list[$k] }{$prey};
     }
     ####print "JMDEBUG: ".join(",", @tot_pep_ar)."\n"; ####
     # compute probability of seeing this sum of TSCs by chance
     my @tsc_vec = @{ $ms_prey_tsc_arr{$prey} };
     my @ln_tsc_vec = &log_vec(@tsc_vec);
     my $mean_ln_tsc = &avg(@ln_tsc_vec);
     my $std_ln_tsc = &std(@ln_tsc_vec);
     my @ln_tot_pep_ar = &log_vec(@tot_pep_ar);
     my $expected_sum_mean = scalar(@tot_pep_ar)*$mean_ln_tsc;
     my $expected_sum_std = sqrt(scalar(@tot_pep_ar) * ($std_ln_tsc**2));

     my $cdf = gsl_cdf_gaussian_Q(sum(@ln_tot_pep_ar) - $expected_sum_mean, 
				  $expected_sum_std);

     # compute the probability of seeing a run this long given the avg frequency
     # discovered during December 2010 update that as the number of experiments got over 6K, calculations like 
     # run_prob 2 6538 0.7031 start returning nan instead of 1.  Unclear why but hack around it for now
     my $run_prob = "~/DPiM/cpp/run_prob";
     my $cmd = join ' ', $run_prob, 0+ @instr_list, 0+ keys %ms_total_by_inst
	 , $common_contam{$prey};
     my $run_p = `$cmd`;
     chomp($run_p);
     $run_p = 1 if ($run_p eq "nan"); ## hack around issue with run_prob for large number of experiments
     ####print "JMDEBUG: ".scalar(@instr_list)." ".scalar(keys %ms_total_by_inst)." ".$common_contam{$prey}." ".$run_p."\n"; ####
     ##$run_p = &run_prob(scalar(@instr_list), scalar(keys %ms_total_by_inst), $common_contam{$prey});
     
     my $score = $cdf * $run_p;
     if ($tot_pep_ar[0]/$tot_pep_ar[1] < 2) {
	 @instr_list = ();
     }

     ####printf "JMDEBUG: %s\t%0.4e\t%0.4e\t%0.4e\t%s\n", $prey, $cdf, $run_p, $score, join(",", @instr_list); ####
     return (join(",", @instr_list), $score);
}

sub sum {
    my @vec = @_;
    my $sum = 0;
    for my $each (@vec) {
	$sum += $each;
    }
    return $sum;
}

sub avg {
    my @vec = @_;
    return &sum(@vec)/scalar(@vec);
}

sub std {
    my @vec = @_;
    my $vecMean = &avg(@vec);
    my @mean_diff_squared = ();
    for my $i (0..$#vec) {
        push @mean_diff_squared, ($vec[$i] - $vecMean)**2;
    }
    my $SSE = &sum(@mean_diff_squared);
    if (scalar(@mean_diff_squared) > 1) {
        return sqrt($SSE/( scalar(@mean_diff_squared) - 1));
    } else {
        return sqrt($SSE/( scalar(@mean_diff_squared) ));
    }
}

sub log_vec {
    my @vec = @_;

    for my $i (0 .. $#vec) {
	$vec[$i] = log($vec[$i]);
    }
    return @vec;
}

# Returns the probability of a run of at least r consecutive successes in n
# independent trials where the probability of success at any one trial is p
sub run_prob {
    my $r = $_[0];
    my $n = $_[1];
    my $p = $_[2];
    my $Bnr = &calc_beta_ln($r, $n, $p);
    my $Bnrr = &calc_beta_ln($r, ($n-$r), $p);
    my $zn = $Bnr - ($p**$r)*$Bnrr;
    my $yn = 1 - $zn;
    $yn = 1 if ($yn > 1);
    ####print "JMDEBUG: sub run_prob: $r $n $p\t$yn\n"; ####
    return $yn;
}

# special sub-function of &run_prob
sub calc_beta_ln {
    my $r = $_[0];
    my $n = $_[1];
    my $p = $_[2];
    my $q = 1 - $p;
    my @Bnr_vec = ();
    for (my $l = 0; $l <= $n/($r+1); $l++) {
	my $tmp_ln = &ln_nchoosek($n-$l*$r, $l) + $l*log($q*$p**$r);
	my $tmp = exp($tmp_ln);
	push @Bnr_vec, ((-1)**$l) * $tmp;
    }
    return &sum(@Bnr_vec);
}

# Log count number of combinations C(n,k)
sub ln_nchoosek {
    my $n = $_[0];
    my $k = $_[1];
    my ($ln_c, $c, $i);
    if ($k > $n/2) {
	$k = $n - $k;
    }
    if ($k <= 1) {
	$c = $n**$k;
	$ln_c = log($c);
    } else {
	my @nums = ();
	for ($i = 0; $i < $k; $i++) {
	    push(@nums, ($n - $k + 1 + $i)/($i + 1));
	}
	$ln_c = &sum( &log_vec(@nums) );
    }
    return $ln_c;
}

sub product {
    my @vec = @_;
    my $p = 1;
    for my $each (@vec) {
	$p *= $each;
    }
    return $p;
}

sub isTrue {
    my $arg = shift;
    return "yes" if $arg;
    return "no";
}

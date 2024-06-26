#!/pkg/perlbrew/0.72/Linux/x86_64/perls/perl-5.20.1/bin/perl

use warnings;

my $SCRDIR = '/home/kli3/proj/Interactome/data/Fly_DPIM_data/DPIM3/newData';


##open(STDERR, ">&STDOUT");
#open(IN, "fbgn_id2.out");
open(IN, $ARGV[0]) or die "can't read from $ARGV[0]. $!";
while ($line = <IN>) {
    chomp($line);
    next if ($line =~ /^clone_name/); # skip header
    #($fbgn, $symbol, $cg_id, $id2, $clone_name, $avg_exp_tryptic) = split(/\s/, $line);
    ($fbgn, $symbol) = split(/\s/, $line);
    $fbgn2symbol{$fbgn} = $symbol;
}
$fbgn2symbol{"FBgn0000000"} = "EmpVecCtrl";
$fbgn2symbol{"GFP"} = "GFPCtrl";
$fbgn2symbol{"turboGFP"} = "GFPCtrl";
$fbgn2symbol{"eGFP"} = "GFPCtrl";
$fbgn2symbol{"tGFP"} = "GFPCtrl";

close(IN);

#open(IN, "dpim_tf_all.120126.plusUniqMsInst");
open(IN, $ARGV[1]);
while ($line = <IN>) {
    chomp($line);
	$line =~ s/\r//; #\KJ\ remove the \cM at the end of the line
    next if ($line =~ /tap_id/); # skip header
    my ($tap_id, $instr_run_id, $user, $search_id, $date, $tot_pep, $uniq_pep, $bait_ref, $prey_ref) = split(/\t/, $line);
    $instr_id = $instr_run_id;
    $instr_id =~ s/[fltqgw]//g;
    $instr_id2run{$instr_id} = $instr_run_id;
    $ms_total_by_search{$instr_id}{$prey_ref} += $tot_pep;
    $ms_bait_by_search{$instr_id} = $bait_ref;
    $ms_prey{$prey_ref} += $tot_pep;
}
close(IN);

# compute the mean and std of TSCs
foreach $instr (keys %ms_total_by_search) {
    foreach $prey (keys %{ $ms_total_by_search{$instr} }) {
	if (exists $ms_prey_tsc_arr{$prey}) {
	    push @{ $ms_prey_tsc_arr{$prey} }, $ms_total_by_search{$instr}{$prey};
	} else {
	    $ms_prey_tsc_arr{$prey} = [ $ms_total_by_search{$instr}{$prey} ];
	}
    }
}

#open(IN, "common_contaminants_nrtap_all.out");
open(IN, $ARGV[2]);
while ($line = <IN>) {
    chomp($line);
    next if ($line =~ /Fraction/); # skip header
    ($fbgn, $cnt, $frac, $mean_tsc) = split(/\t/, $line);
    $common_contam{$fbgn} = $frac;
    $common_contam_tsc{$fbgn} = $mean_tsc;
    #$common_contam_cnt{$fbgn} = $cnt;
}
close(IN);


@sorted_instr = sort {$a <=> $b} keys %ms_total_by_search;
@sorted_prey = sort keys %ms_prey;

####print scalar(@sorted_instr)."\t".scalar(@sorted_prey)."\n"; exit(0);
####### CONSTRUCT SEQUENTIAL RUN X PREY MATRIX ########
@instr_by_prey_ar = ();
for $i (0..$#sorted_instr) {
    for $j (0..$#sorted_prey) {
	if (exists $ms_total_by_search{ $sorted_instr[$i] }{ $sorted_prey[$j] }) {
	    $instr_by_prey_ar[$i][$j] = $ms_total_by_search{ $sorted_instr[$i] }{ $sorted_prey[$j] };
	} else {
	    $instr_by_prey_ar[$i][$j] = 0;
	}
    }
}

######## LOOK FOR MONOTONICALLY DECREASING SEQUENTIAL PREY TSCS ########
#### Allow some leeway for washes/standards
#### Sequences have to start with at least N TSCs (20)
#for $j (0..10) {
for $j (0..$#sorted_prey) {
    $prev_instr_tot = $instr_by_prey_ar[$#sorted_instr][$j];
    $monotonic_seq = 0;
    $monotonic_seq_last = scalar(@sorted_instr)-1;
    for ($i = scalar(@sorted_instr)-2; $i >= 0; $i--) {
	####print $instr_by_prey_ar[$i][$j]."\t"; ####
	if ($instr_by_prey_ar[$i][$j] > $prev_instr_tot &&
	    ($sorted_instr[$i+1] - $sorted_instr[$i]) <= 4 &&
	    $ms_bait_by_search{ $sorted_instr[$i] } ne $ms_bait_by_search{ $sorted_instr[$i+1] }) {
	    $monotonic_seq = 1;
	    if ($prev_instr_tot == 0) { # avoid sequences ending with and including 0
		$monotonic_seq_last = $i;
	    }
	} elsif ($prev_instr_tot >= 20 && ($monotonic_seq_last > ($i+1))) {
	    ####print "\t".($i+1)." - ".$monotonic_seq_last."\t"; ####
	    @instr_list = @sorted_instr[($i+1) .. $monotonic_seq_last];
	    ####print "$i|$monotonic_seq_last|\t|".join(",", @instr_list)."|\t|".$sorted_prey[$j]."|\n"; ####
	    ($instr_line, $pep_sum) = &checkSeq(join(",", @instr_list), $sorted_prey[$j]);
	    @instr_list = split(/,/, $instr_line);
	    if (scalar(@instr_list) > 1) {
		#$monotonic_seq_hash{ $sorted_prey[$j] }{ join(",", @instr_list) } = $pep_sum;
		$monotonic_seq_hash{ $sorted_prey[$j]."\t".join(",", @instr_list) } = $pep_sum;
	    }
	    $monotonic_seq_last = $i;
	} else {
	    ####print "zeroing\t".($i+1)." - ".$monotonic_seq_last."\t"; ####
	    for $each (($i+1)..$monotonic_seq_last) {
		$instr_by_prey_ar[$each][$j] = 0;
	    }
	    $monotonic_seq = 0;
	    $monotonic_seq_last = $i;
	}
	$prev_instr_tot = $instr_by_prey_ar[$i][$j];
	####print "\n"; ####
    }
    ####exit(0); ####
}
    

#for $i (0..$#sorted_instr) {
#    print $sorted_instr[$i];
#    for $j (0..20) {
#	print "\t".$instr_by_prey_ar[$i][$j];
#    }
#    print "\n";
#}

open(OUT, ">find_sequential_contam.out");
$cnt = 0;
foreach $one (sort {$monotonic_seq_hash{$a} <=> $monotonic_seq_hash{$b}} keys %monotonic_seq_hash) {
    $cnt ++;
    ($prey, $instr_line) = split(/\t/, $one);
    @instr_list = split(/,/, $instr_line);

    printf "\t%0.2e\t%s\t%s\t%0.4f\t%0.2f\n", $monotonic_seq_hash{$one}, $prey, $fbgn2symbol{$prey}, $common_contam{$prey}, $common_contam_tsc{$prey};
    for $instr (sort {$a <=> $b} @instr_list) {
	printf "%s\t%s\t%s\t%d\n", $instr_id2run{$instr}, $ms_bait_by_search{$instr}, $fbgn2symbol{ $ms_bait_by_search{$instr} }, $ms_total_by_search{$instr}{$prey};
    }
    print "\n";
    @bait_list = map{ $ms_bait_by_search{$_} } @instr_list;
    print OUT join("\t", @bait_list)."\n";
}
close(OUT);


sub checkSeq {
     my $instr_line = $_[0];
     my @instr_list = split(/,/, $instr_line);
     my $prey = $_[1];
     ####print "JMDEBUG: |$prey|$instr_line|\n"; ####

     my @sorted_instr_list = sort {$a <=> $b} @instr_list;

     my @tot_pep_ar = ();
     for my $k (0..$#sorted_instr_list) {
	 push @tot_pep_ar, $ms_total_by_search{ $sorted_instr_list[$k] }{$prey};
     }
     ####print "JMDEBUG: ".join(",", @tot_pep_ar)."\n"; ####
     # compute probability of seeing this sum of TSCs by chance
     @tsc_vec = @{ $ms_prey_tsc_arr{$prey} };
     @ln_tsc_vec = &log_vec(@tsc_vec);
     $mean_ln_tsc = &avg(@ln_tsc_vec);
     $std_ln_tsc = &std(@ln_tsc_vec);
     @ln_tot_pep_ar = &log_vec(@tot_pep_ar);
     $expected_sum_mean = scalar(@tot_pep_ar)*$mean_ln_tsc;
     $expected_sum_std = sqrt(scalar(@tot_pep_ar) * ($std_ln_tsc**2));

     
     open(CDF, "$SCRDIR/gsl_cdf_gaussian_Q ".(&sum(@ln_tot_pep_ar) - $expected_sum_mean)." ".$expected_sum_std." | ");
     my $cdf = <CDF>;
     close(CDF);
     chomp($cdf);

     # compute the probability of seeing a run this long given the avg frequency
     # discovered during December 2010 update that as the number of experiments got over 6K, calculations like 
     # run_prob 2 6538 0.7031 start returning nan instead of 1.  Unclear why but hack around it for now
     my $run_prob = "$SCRDIR/run_prob";
     open(RP, "$run_prob ".scalar(@instr_list)." ".scalar(keys %ms_total_by_search)." ".$common_contam{$prey}." | ");
     my $run_p = <RP>;
     close(RP);
     chomp($run_p);
     $run_p = 1 if ($run_p eq "nan"); ## hack around issue with run_prob for large number of experiments
     ####print "JMDEBUG: ".scalar(@instr_list)." ".scalar(keys %ms_total_by_search)." ".$common_contam{$prey}." ".$run_p."\n"; ####
     ##$run_p = &run_prob(scalar(@instr_list), scalar(keys %ms_total_by_search), $common_contam{$prey});
     
     $score = $cdf * $run_p;
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
    $SSE = &sum(@mean_diff_squared);
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
    $Bnr = &calc_beta_ln($r, $n, $p);
    $Bnrr = &calc_beta_ln($r, ($n-$r), $p);
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
    @Bnr_vec = ();
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
	$tmp = &sum( &log_vec(@nums) );
	$ln_c = $tmp;
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

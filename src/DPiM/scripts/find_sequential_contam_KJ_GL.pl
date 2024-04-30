#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Math::GSL::CDF qw(gsl_cdf_gaussian_Q);
use Data::Dumper;
use DpimLib qw(getLineDP4APMS getLineBiopAPMS);
use HomeBrew::IO qw(checkExist);

my %RUN_PROB_CACHE; # save time by remembering previous calculations

my %opts = getCommandLineOptions();

{
    my $apmsFile = $opts{apms};
    my $commFile = $opts{comm};
    my $symbFile = $opts{symb};
    my $outFile = $opts{out};
    my $logFile = "$outFile.log_seqContam";

    my %fbgn2symbol = mapFBgn2symb($symbFile);

    my (%tsc_by_runNum, %bait_by_runNum, %allPrey, %runNum2runID);
    {
	my ($tbr, $bbr, $ap, $rn2rid) = parseAPMS($apmsFile);
	%tsc_by_runNum = %$tbr;
	%bait_by_runNum = %$bbr;
	%allPrey = %$ap; 
	%runNum2runID = %$rn2rid;
    }

    my (%common_contam, %common_contam_tsc);
    {
	my ($cc, $cct) = parseCommContam($commFile);
	%common_contam = %$cc;
	%common_contam_tsc = %$cct;
    }

    ######## LOOK FOR MONOTONICALLY DECREASING SEQUENTIAL PREY TSCS ########
    #### Allow some leeway for washes/standards
    #### Sequences have to start with at least $minTSC TSCs
    my %monotonic_seq_hash = seekSeqs(\%tsc_by_runNum, \%bait_by_runNum, 
				      \%allPrey, \%common_contam);

    my @outputOrder = sort keys %monotonic_seq_hash;
    @outputOrder= sort {$monotonic_seq_hash{$a} <=> $monotonic_seq_hash{$b}} @outputOrder;
    
    open my $LOG, ">", $logFile
	or die "can't write to $logFile: $!";
    open my $OUT, ">", $outFile or die "can't write to $outFile: $!";
    foreach my $one (@outputOrder) {
	my ($prey, $instr_line) = split(/\t/, $one);
	my @instr_list = split(/,/, $instr_line);

	printf $OUT "\t%0.2e\t%s\t%s\t%0.4f\t%0.2f\n"
	    , $monotonic_seq_hash{$one}, $prey
	    , $fbgn2symbol{$prey} // "unknown_$prey"
	    , $common_contam{$prey}, $common_contam_tsc{$prey};
	
	for my $instr (sort {$a <=> $b} @instr_list) {
	    printf $OUT "%s\t%s\t%s\t%d\n"
		, $runNum2runID{$instr}, $bait_by_runNum{$instr}
	        , $fbgn2symbol{ $bait_by_runNum{$instr} } // "unknown_$prey"
		, $tsc_by_runNum{$instr}{$prey};
	}
	print $OUT "\n";
	my @bait_list = map{ $bait_by_runNum{$_} } @instr_list;
	print $LOG join("\t", @bait_list)."\n";
    }
    close($OUT);
    close($LOG);
}
exit;
   

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
	symb => '/home/glocke/DPiM/fbgn_id2_2col_03-17-2016.txt',
	## for human, use entrez2symb.concise.tsv
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -apms in.apms -comm in.commContam -out out.seqContam"
	." < $modeString $defaultString -decreasing >\n";

    my %opts = ();
    GetOptions(\%opts, "apms=s", "comm=s", "out=s", "mode=s", "symb=s", 
	"decreasing");
    die $usage unless exists $opts{apms} && exists $opts{comm} && 
	exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }
    $opts{mode} //= $modes[0];
    die "you must select one of the following modes: "
	, join(", ",keys(%modes)), "\n" if ! exists $modes{$opts{mode}};
    
    checkExist('f', $opts{apms});
    checkExist('f', $opts{comm});
    checkExist('f', $opts{symb});

    return %opts;
}

sub mapFBgn2symb {
    my ($symbFile) = @_;

    my %ret;
    open my $IN, "<", $symbFile or die "can't read from $symbFile. $!";
    while (my $line = <$IN>) {
	next unless $line =~ /^FBgn/ || $line =~ /^\d+/; 
	# latter condition reflects entrez id's, used in human data
	chomp($line);
	my ($fbgn, $symbol) = split(/\s/, $line);
	$ret{$fbgn} = $symbol;
    }
    $ret{"FBgn0000000"} = "EmpVecCtrl";
    $ret{"GFP"} = "GFPCtrl";
    $ret{"turboGFP"} = "GFPCtrl";
    $ret{"eGFP"} = "GFPCtrl";
    $ret{"tGFP"} = "GFPCtrl";

    close($IN);
    return %ret;
}

sub parseAPMS {
    my ($apmsFile) = @_;
    my $reader = \&getLineDP4APMS;
    $reader = \&getLineBiopAPMS if $opts{mode} eq 'human';

    my (%tsc_by_runNum, %bait_by_runNum, %allPrey, %runNum2runID);

    open my $IN, "<", $apmsFile or die "Can't open $apmsFile. $!";
    my %row;
    while($reader->(\%row, $IN)){
	my $runID = $row{ms_inst_run_id};
	my $bait_ref = $row{bait_ref};
	my $prey_ref = $row{prey_ref};
	my $tot_pep = $row{total_peptides};
	next if $tot_pep == 0;
	my $runNum = $runID;
	$runNum =~ s/[a-z]//g;
	$runNum2runID{$runNum} = $runID;
	$tsc_by_runNum{$runNum}{$prey_ref} += $tot_pep;
	$bait_by_runNum{$runNum} = $bait_ref;
	$allPrey{$prey_ref} = 1;
    }

    close($IN);
    return(\%tsc_by_runNum, \%bait_by_runNum, \%allPrey, \%runNum2runID);
}

sub parseCommContam {
    my ($commFile) = @_;

    my (%common_contam, %common_contam_tsc);    
	
    #open($IN, "common_contaminants_nrtap_all.out");
    open my $IN, "<", $commFile or die "Can't open $commFile. $!";
    while (my $line = <$IN>) {
	chomp($line);
	next if ($line =~ /Fraction/); # skip header
	my ($fbgn, $cnt, $frac, $mean_tsc) = split(/\t/, $line);
	$common_contam{$fbgn} = $frac;
	$common_contam_tsc{$fbgn} = $mean_tsc;
	#$common_contam_cnt{$fbgn} = $cnt;
    }
    close($IN);
    return(\%common_contam, \%common_contam_tsc);
}

######## LOOK FOR MONOTONICALLY DECREASING SEQUENTIAL PREY TSCS ########
#### Allow some leeway for washes/standards
#### Sequences have to start with at least $minTSC TSCs
sub seekSeqs {
    my ($tbr, $bbr, $ap, $cc) = @_;

    my %tsc_by_runNum = %$tbr;
    my %bait_by_runNum = %$bbr;
    my %allPrey = %$ap; 
    my %common_contam = %$cc;
	
    my %tsc_seq; # TSC for all appearances of each prey
    foreach my $instr (keys %tsc_by_runNum) {
	foreach my $prey (keys %{ $tsc_by_runNum{$instr} }) {
	    $tsc_seq{$prey} //= [];
	    push @{ $tsc_seq{$prey} }, $tsc_by_runNum{$instr}{$prey};
	}
    }
    my %tsc_stats = map { $_ => preyStats($tsc_seq{$_}) } keys %tsc_seq;

    my %monotonic_seq_hash; # return value
	
    my @sorted_instr = sort {$a <=> $b} keys %tsc_by_runNum;
    my @sorted_prey = sort keys %allPrey;

    ####### CONSTRUCT SEQUENTIAL RUN X PREY MATRIX ########
    ## GL note: this is a sparse data structure and ought to be stored sparsely
    my @instr_by_prey_ar;
    for my $i (0..$#sorted_instr) {
	for my $j (0..$#sorted_prey) {
	    $instr_by_prey_ar[$i][$j] = 
		$tsc_by_runNum{ $sorted_instr[$i] }{ $sorted_prey[$j] } // 0;
	}
    }

    my $minTSC = 20;
    my $cmp = sub {
	my ($tsc1, $tsc2) = @_;
	return $tsc1 >= $tsc2;
    };
    $cmp = sub {
	my ($tsc1, $tsc2) = @_;
	return $tsc1 > $tsc2;
    } if exists $opts{decreasing};
    
	
    my $nanCount = 0;
    #for $j (0..10) {
    for my $preyIndex (0..$#sorted_prey) {
	my $prey = $sorted_prey[$preyIndex];
	my $prev_instr_tot = $instr_by_prey_ar[$#sorted_instr][$preyIndex];
	##my $monotonic_seq = 0;
	my $monotonic_seq_last = $#sorted_instr;
	for (my $i = $#sorted_instr-1; $i >= 0; $i--) {
	    my $thisInstr = $sorted_instr[$i];
	    my $nextInstr = $sorted_instr[$i+1];
	    ####print $instr_by_prey_ar[$i][$preyIndex]."\t"; ####
	    if ($cmp->($instr_by_prey_ar[$i][$preyIndex], $prev_instr_tot) &&
		($nextInstr - $thisInstr) <= 4 &&
		$bait_by_runNum{ $thisInstr } ne $bait_by_runNum{ $nextInstr })
	    {
		##$monotonic_seq = 1;
		if ($prev_instr_tot == 0) { 
		    # avoid sequences ending with and including 0
		    $monotonic_seq_last = $i;
		}
	    } elsif ($prev_instr_tot >= $minTSC && ($monotonic_seq_last > ($i+1))) {
		my @instr_list = @sorted_instr[($i+1) .. $monotonic_seq_last];
		#my ($instr_line, $prey, $tbr, $ap, $cc) = @_;
		my ($instr, $seqPVal, $nanFlag) 
		    = checkSeq(\@instr_list, $sorted_prey[$preyIndex], 
			       \%tsc_by_runNum, \%common_contam, 
			       $tsc_stats{$prey}{mean}, $tsc_stats{$prey}{sd});
		$nanCount++ if $nanFlag;
		@instr_list = @$instr;
		
		if (scalar(@instr_list) > 1) {
		    $monotonic_seq_hash{ mshKey( $sorted_prey[$preyIndex], 
						 \@instr_list) } = $seqPVal;
		}
		$monotonic_seq_last = $i;
	    } else {
		for my $each (($i+1)..$monotonic_seq_last) {
		    $instr_by_prey_ar[$each][$preyIndex] = 0;
		}
		##$monotonic_seq = 0;
		$monotonic_seq_last = $i;
	    }
	    $prev_instr_tot = $instr_by_prey_ar[$i][$preyIndex];
	}
    }
    
    warn "nanCount = $nanCount" if $nanCount > 0;
    return %monotonic_seq_hash;
}

# take log of every tsc for every appearance of this prey
# find mean and stddev for this vector
sub preyStats {
    my ($tsc) = @_;

    my @tsc_vec = @$tsc;
    my @ln_tsc_vec = log_vec(@tsc_vec);
    my $mean = avg(@ln_tsc_vec);
    my $sd = std(@ln_tsc_vec);
    return { mean=>$mean, sd=>$sd };
}

sub mshKey {
    my ($sorted_prey, $instr_list) = @_;
    return $sorted_prey."\t".join(",", @$instr_list);
}

# combine two statistical tests
# test1: is the sum of spectral counts in this sequence unusually high?
# test2: is it unusual to have so many consecutive appearances of this protein?
# test1 is a z-test, while test2 is a more complicated formula
# each gives you a p-value, and the product of these is returned
sub checkSeq {
    my ($instr, $prey, $tbr, $cc, $mean_ln_tsc, $std_ln_tsc) = @_;

    my %tsc_by_runNum = %$tbr;
    my %common_contam = %$cc;

    my @instr_list = @$instr;
    ####print "JMDEBUG: |$prey|$instr_line|\n"; ####

    my @sorted_instr_list = sort {$a <=> $b} @instr_list;

    # compute probability of seeing this sum of TSCs by chance
    my @tot_pep_ar = ();
    for my $k (0..$#sorted_instr_list) {
	push @tot_pep_ar, $tsc_by_runNum{ $sorted_instr_list[$k] }{$prey};
    }
    my @ln_tot_pep_ar = &log_vec(@tot_pep_ar);
    my $expected_sum_mean = scalar(@tot_pep_ar)*$mean_ln_tsc;
    my $expected_sum_std = sqrt(scalar(@tot_pep_ar) * ($std_ln_tsc**2));

    my $cdf = gsl_cdf_gaussian_Q(sum(@ln_tot_pep_ar) - $expected_sum_mean, 
				 $expected_sum_std);

    # compute the probability of seeing a run this long given the avg frequency
    ## discovered during December 2010 update that as the number of experiments 
    ##   got over 6K, calculations like `run_prob 2 6538 0.7031` start returning
    ##   nan instead of 1.  Unclear why but hack around it for now 
    my $run_p;
    my $nanFlag = undef;
    if ($common_contam{$prey} == 1) {
	$run_p = 1; 
    } else {
	my $run_prob = "$ENV{DPSCR}/run_prob";
	my $consec = 0+ @instr_list;
	my $trials = 0+ keys %tsc_by_runNum;
	my $prob = $common_contam{$prey};
	if (exists $RUN_PROB_CACHE{$consec}{$trials}{$prob}) {
	    $run_p = $RUN_PROB_CACHE{$consec}{$trials}{$prob};
	} else {
	    my $cmd = "$run_prob $consec $trials $prob";
	    $run_p = `$cmd`;
	    chomp($run_p);
	    $RUN_PROB_CACHE{$consec}{$trials}{$prob} = $run_p;
	}
	if ($run_p =~ /nan/i) {
	    ## hack around issue with run_prob for large number of experiments
	    $run_p = 1;
	    $nanFlag=1;
	    #warn "nan: $cmd\n";
	}
    }
    
    my $score = $cdf * $run_p;
    if ($tot_pep_ar[0]/$tot_pep_ar[1] < 2) {
	@instr_list = ();
    }

    return (\@instr_list, $score, $nanFlag);
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



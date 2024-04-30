package HomeBrew::Bio;

use strict;
use warnings;
use List::Util qw(sum);

require Exporter;

# makeOcc(\@occ,\@prob,$objSize) first arg is return value.  must be empty.
# makeOcc(\@prob,$objSize) # as prev, but returns array-ref instead of void
## \_> do not normalize!
# reverseComplement($seq); # alias rc
# makeWeightMatrix(\@output, \@trainingSeq, \%background)
## $output[$i]{$word} = ln( freqOfWordAtPosnI / bgFreqOfWord )
# applyWeightMatrix($seq, \@matrix) = score # $seq may be string or array ref
# seqToDinucs($seq, $mode) $seq is pbr.  current modes: (AT GC AATT)
# wordCount($seq, $wordLength, $counts) count all words of length N in seq
## optionally pass an existing counts to add to
# makeWords(@alphabet, $length); # make all words of $length from @alpha
# hamming($word1, $word2); # hamming distance between two words **OF SAME LEN**

our @ISA = qw(Exporter);

#our @EXPORT = qw(checkExist);
our @EXPORT_OK = qw(makeOcc reverseComplement rc makeWeightMatrix 
    applyWeightMatrix seqToDinucs monoNucToWord wordCount makeWords hamming);

our $VERSION = '0.01';

our $argIllDefined = -2; # if arg is not the correct reference type or 

sub makeOcc {
    my $occNotProvided = @_ < 3;
    my ($occ, $prob, $objSize);
    if ($occNotProvided) {
	$occ = [];
	($prob, $objSize) = @_;
    } else {
	($occ, $prob, $objSize) = @_;
    }
    
    # verify arguments
    unless ($objSize>0) {
	warn "makeOcc: objSize too small ($objSize)\n";
	return $argIllDefined;
    }
    unless (ref($occ) eq "ARRAY" && @$occ == 0) {
	warn "makeOcc: occ arg must be an empty array ref\n";
	return $argIllDefined;
    }    
    unless (ref($prob) eq "ARRAY" && @$prob >= $objSize ) {
	warn "makeOcc: prob arg must be a ref to array at least as large as objSize ($objSize)\n";
	return $argIllDefined;
    }
    
    my $i;
    my $o=0;
    for $i (0..$objSize-1) {
	$o+= $prob->[$i];
	push(@$occ, $o);
    } 
    for my $i ($objSize..$#$prob) {
	$o -= $prob->[$i-$objSize];
	$o += $prob->[$i];
	push(@$occ, $o);
    }

    return $occ if $occNotProvided; # let user de-reference rather than spamming stack
    return;
}

sub reverseComplement {
    return &rc;
}

sub rc {
    my ($seq) = @_;
    my %pairs = ();
    $pairs{'A'} = 'T';
    $pairs{'T'} = 'A';
    $pairs{'C'} = 'G';
    $pairs{'G'} = 'C';
    $pairs{'a'} = 'T';
    $pairs{'t'} = 'A';
    $pairs{'c'} = 'G';
    $pairs{'g'} = 'C';
    my @bits_r = reverse split('',$seq);
    my @bits_rc = ();
    foreach my $bit (@bits_r) {
	if (defined $pairs{$bit}) {
	    push @bits_rc,$pairs{$bit};
	}else{
	    die "Unknown base in rc: '$bit'!\n";
	}
    }
    return join '',@bits_rc;
}

# score = sum_over_positions_i[log(freq of base at i/ base freq of base)]
sub makeWeightMatrix {
    my ($ret, $trainingSeq, $background) = @_;

    # verify input
    die "Bio::makeWeightMatrix expects first arg to be empty array ref\n"
	unless ref($ret) eq "ARRAY" && @$ret==0;
    die "Bio::makeWeightMatrix expects second arg to be a non-empty array of strings\n"
	unless ref($trainingSeq) eq "ARRAY" && @$trainingSeq && 
	length($trainingSeq->[0]);
    die "Bio::makeWeightMatrix expects third arg to be non-empty hash ref\n" 
	unless ref($background) eq "HASH" && scalar(keys %$background);
    my @bases = keys(%$background);
    for (@bases) {
	die "Bio::makeWeightMatrix background must be frequencies (0<=x<=1)\n"
	    unless $background->{$_} >=0 && $background->{$_} <= 1;
    }
    

    my @freq; # freq[pos]{key} = # of counts / number of training sequences
    # keys are expected to match background (ACGT or amino acid)
    my $expectedLength = length($trainingSeq->[0]);
    for my $seq (@$trainingSeq) {
	die "Bio::makeWeightMatrix found a training sequence of non-uniform length.\n\tExpected length $expectedLength.  Found $seq (length ", length($seq), ")\n"
	    if length($seq) != $expectedLength;

	my $i=0;
	for my $c ($seq =~ /./g) { # loop over each letter in 
	    $freq[$i]{$c}++;
	    $i++;
	}
    }
    #use Data::Dumper;
    #die Dumper(\@freq);

    # check for unrecognized keys
    my $pos;
    for $pos (@freq) {
	for (keys %$pos) {
	    die "can't find key $_ in background\n" 
		unless exists $background->{$_};
	}
    }

    # as of now, @freq is really counts, not frequencies.
    # fix this, and look out for zeros
    my $i=0;
    for $pos (@freq) {
	for (@bases) {
	    unless (exists $pos->{$_} && $pos->{$_}) {
		warn "found no $_ at position $i\n";
		$pos->{$_} = 0.1;
	    }
	    $pos->{$_} /= @$trainingSeq;
	}
	$i++;
    }
    
    
    $i=0;
    for $pos (@freq) {
	for (@bases) {
	    $ret->[$i]{$_} = log($pos->{$_} / $background->{$_});
	}
	$i++;
    }
    return;
}

sub applyWeightMatrix {
    
    my ($seq, $weightMatrix) = @_;
    
    die "Bio::applyWeightMatrix: seq may be a string or an array ref\n" 
	if (ref($seq) && ref($seq) ne 'ARRAY');

    my $score = 0;
    if (ref($seq)) {
	for my $i (0..$#$seq) {
	    if (exists $weightMatrix->[$i]{$seq->[$i]}) {
		$score += $weightMatrix->[$i]{$seq->[$i]};
	    } else {
		warn "Bio::applyWeightMatrix: can't find base ", $seq->[$i], 
		    " at position $i in weight matrix\n";
		die $argIllDefined;
	    }
	}
    } else {
	my $i=0;
	for my $c ($seq =~ /./g) { # loop over each letter in 
	    if (exists $weightMatrix->[$i]{$c}) {
		$score += $weightMatrix->[$i]{$c};
	    } else {
		warn "Bio::applyWeightMatrix: can't find base $c", 
		    " at position $i in weight matrix\n";
		die $argIllDefined;
	    }
	    $i++;
	}
    }
    return $score;
}

# return a string of zeroes and ones where ones match according to mode
sub seqToDinucs {
    my ($seq, $mode) = @_;
    die "Bio::seqToDinucs first arg must be scalar or scalar ref\n" 
	if ref($seq) && ref($seq) ne 'SCALAR';
    $seq = \$_[0] unless ref($seq);
    
    my %modes = qw(AT 1 GC 1 AATT 1);
    die "Bio::seqToDinucs doesn't recognize mode = $mode.\n".
	"\tAcceptable modes are ", join(', ', keys(%modes)), "\n" 
	unless exists $modes{$mode};

    if ($mode eq 'AT') {
	# turn AA AT TA TT into 1.  all others are 0.
	$$seq =~ s/(A|T|a|t)(?=(A|T|a|t))/1/g;
	# regex matches A or T followed by A or T
    } elsif ($mode eq 'GC') {
	# turn GG GC CG CC into 1.  all others are 0.
	$$seq =~ s/(G|C|g|c)(?=(G|C|g|c))/1/g;
	# regex matches G or C followed by G or C
    } elsif ($mode eq 'AATT') { 
	# turn AA and TT to 1, others to 0
	$$seq =~ s/(A|a)(?=(A|a))/1/g;
	$$seq =~ s/(T|t)(?=(T|t))/1/g;

	# can also do with single regex, but keep it simple for maintenance
	# s/((A|a)(?=(A|a))|(T|t)(?=(T|t)))/1/g;
    }
    $$seq =~ s/[^1]/0/g; # nullify everything that ain't a match

    chop($$seq); # number of dinucs = length of sequence - 1
}

# get an array of words from an array of mononucs
sub monoNucToWord {
    my ($ret, $seq, $N) = @_;
    die "Bio::monoNucToWord second arg (ret) must be array ref" 
	if ref($seq) ne 'ARRAY';
    die "Bio::monoNucToWord second arg (seq) must be array ref" 
	if ref($seq) ne 'ARRAY';
    die "Bio::monoNucToWord third arg (word length) must be greater than one" 
	if $N < 2;

    @$ret = @$seq;
    for my $n (1..($N-1)) {
	pop @$ret;
	for my $i ($n..$#$seq) {
	    $ret->[$i-$n].=$seq->[$i];
	}
    }
    return;
}

sub wordCount {
    my ($seq, $wordLength, $counts) = @_;
    
    die "Bio::wordCount optional third arg must be hash-ref\n" 
	if defined $counts && ref($counts) ne 'HASH';
    die "Bio::wordCount word length is longer than sequence\n" 
	if $wordLength > length($seq);

    my $cnt = $counts // {};

    my $w;
    for my $i (0 .. length($seq) - $wordLength) {
	$w = substr($seq, $i, $wordLength);
	$cnt->{$w}++;
    }
    
    return %$cnt unless defined $counts;
    return;
}

sub makeWords {
    my ($alphabet, $n) = @_;
    
    return @$alphabet if ($n==1);
    
    my (@words,$w);
    for my $nuc (@$alphabet) {
        for $w (makeWords($alphabet, $n-1)) {
            push(@words, $nuc.$w);
        }
    }

    return @words;
}

sub hamming {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

"no text after this line"


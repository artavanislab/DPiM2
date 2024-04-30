package libBindSter;

use strict;
use warnings;
use Exporter;
#use Bio::SeqIO; # it's not hard to work around this, but it's convenient

# IO related
# checkExist(mode,path) dies on fail
# readCol($file,$col) read named column from file, return @values
# readColRef(\@values, $file,$col) as above, but pass empty \@values
# /--
# readCols(file, @cols, $wholeLine) = @ret; $ret[i] = hash ref keyed by @cols
# readColsRef(\@ret, file, @cols, $wholeLine); as above.
# readCols2(file, @cols, $wholeLine) = %ret; $ret[$col] = array ref
# readCols2Ref(\%values, file, @cols, $wholeLine) as above but pass empty \%ret
# \_.> if whole line $wholeLine is defined, it names a column that holds
####   the whole (chomped) line of the file

# bio related
# makeWords(@alphabet, $length); # make all words of $length from @alpha
# rc($string) = reverse complement.  will die if $string =~ /![ACGTacgt]/

# these subroutines help in handling fullProfile output
# readFullProfile($ret, $file, $seqIDList) <- read output of fullProfile
# sumProfile($ret, $profile) <- for each probe, sum all binding modes
# findObjSize($singleProbe), given a fullProfile output, find the footprint
# \_.> input would be $rfp->{seqID} where $rfp is set by readFullProfile

our @ISA = qw(Exporter);

#our @EXPORT = qw(checkExist);
our @EXPORT_OK = qw(checkExist readList readListRef readList2 readCol readColRef readCols readColsRef 
    readCols2 readCols2Ref readFasta readFastaRef readHeader makeWords rc
readFullProfile sumProfile findObjSize);

our $pathNotFound = -3;
our $emptyFile = -4;

# die unless path isn't found or if file is empty ('f' only)
sub checkExist {
  my ($type, $path) = @_;
  die "checkExist: path variable not defined\n" unless defined $path;
  if ($type eq 'd') {
    if (! -d $path) {
      warn("checkExist: dir not found: $path\n");
      die $pathNotFound;
    }
  }
  elsif ($type eq 'f') {
    if (! -f $path) {
      warn("checkExist: file not found: $path\n");
      die $pathNotFound;
    }
    elsif (! -s $path) {
      warn("checkExist: empty file: $path\n");
      die $emptyFile;
    }
  } else {
      warn "checkExist: unrecognized mode '$type'\n";
  }
  return;
}


# This function reads in a list of items from a file (one item per line, rest of line ignored;
# '#' comment symbols recognized).
sub readList {
    my @ret;
    readListRef(\@ret, @_);
    return @ret;
}

sub readListRef {
  my ($ret, $file, $verbose) = @_;
  die "readList: undefined file\n" unless defined $file;
  die "readListRef: expect first arg to be empty array ref\n" 
      unless ref($ret) eq 'ARRAY' && @$ret==0;

  my $cnt = 0;
  open my $IN, "<", $file or warn "readList can't read from $file: $!\n" and 
      return;
  while (<$IN>) {
    next if /^#/;
    chomp;
    my @tt = split;
    next if (@tt==0);
    my $token = $tt[0];
    push @$ret, $token;
    $cnt++;
  }
  close $IN;
  print "Read $cnt token(s) from file: $file ..\n" if $verbose;
}

# This function reads in a list of items from a file 
# ret[line][entry] one entry per space separated item on line
# '#' comment symbols recognized if they appear first in a line
sub readList2 {
  my ($file, $verbose) = @_;
  my @tokens = ();
  my $cnt = 0;
  open my $IN, "<", "$file" or die "Can't open file: $!\n";
  while (<$IN>) {
    next if /^#/;
    chomp;
    my @tt = split;
    next if (@tt==0);
    push @tokens, \@tt;
    $cnt+=@tt;
  }
  close $IN;
  print "Read $cnt token(s) from file: $file ..\n" if $verbose;
  return @tokens;
}

sub readCols {
    my @ret;
    readColsRef(\@ret, @_);
    return @ret;
}

# third arg is an array, but its elements may be array references.
#  cols = ([c11,c12],c2,[c31,c32...],...)
# the return value is an array of hash references:
#  ret[i] = { c11 -> val1, c2 -> val2, c31-> val3,...}
#   note that the FIRST element in the array ref [cn1,cn2] is always the key for valn
# the purpose of these array references is that sometimes a certain column
## in the source file may have been named inconsistently (e.g. filter was 
## briefly written as counts).  $cols = ["seqN",["filter","counts"]]
## gets you the column named 'filter' or 'counts', preferring the former.
## Note that the return value will be named 'filter' in either case.
sub readColsRef {
    my ($ret, $f, $cols, $wholeLine) = @_;

    die "IO::readColsRef expects first arg to be empty array ref\n" 
	unless ref($ret) eq 'ARRAY' && @$ret==0;
    die "IO::readColsRef expects third arg to be non-empty array ref\n" 
	unless ref($cols) eq 'ARRAY' && @$cols;
    die "IO::readColsRef input filename arg not defined\n"
	if ! defined $f || ! length($f);

    my $IN;
    if ( (split /\./, $f)[-1] eq 'gz') {
	die "can't open zipped files yet!\n";
	open $IN, "<", $f or die "readCol couldn't read from '$f'. $!\n";
    } else {
	open $IN, "<", $f or die "readCol couldn't read from '$f'. $!\n";
    }
    
    my (%loc,$c,$cc,$line);
    for $c (@$cols) {
	if (ref($c) eq 'ARRAY') {
	    $loc{$c->[0]} = -1;
	} else {
	    $loc{$c} = -1;
	}
    }

    # find header
    do {
	$line = <$IN>;
    } while (substr($line,0,1) eq '#');
    chomp $line;

    $_ = $line;
    my @spl = split;
    
    # parse header to find each target column
    my $flag=0;
    for my $i (0..$#spl) {
	for $c (@$cols) {
	    if (ref($c) eq 'ARRAY') {
		$flag = 0;
		for $cc (@$c) {
		    if ($spl[$i] eq $cc) {		    
			$loc{$c->[0]} = $i;
			$flag=1;
			last;
		    }
		}
		last if $flag;
	    } else {
		if ($spl[$i] eq $c) {
		    $loc{$c} = $i;
		    last;
		}
	    }
	}
    }
    
    for $c (@$cols) {
	if (ref($c) eq 'ARRAY') {
	    for $cc (@$c) {
		die "readColsRef: can't find '$cc' column in header of '$f'\n" if $loc{$c->[0]} < 0;
	    } 
	} else {
	    die "readColsRef: can't find '$c' column in header of '$f'\n" if $loc{$c} < 0;
	}
    }

    while ($line = <$IN>) {
	my $ent = {};
	next if substr($line,0,1) eq '#';
	chomp $line;
	$_ = $line;
	@spl = split;
	for $c (@$cols) {
	    if (ref($c) eq 'ARRAY') {
		$ent->{$c->[0]} = $spl[$loc{$c->[0]}];	    
	    } else {
		$ent->{$c} = $spl[$loc{$c}];	    
	    }
	}
	$ent->{$wholeLine} = $line if defined $wholeLine;

	push(@$ret,$ent);
    }
    close IN;
    return;
}

sub readCols2 {
    my %ret;
    readCols2Ref(\%ret, @_);
    return %ret;
}

# as readCols EXCEPT it returns a hash of array refs instead of an array of 
# hash refs
sub readCols2Ref {
    my ($ret, $f,$colsA, $wholeLine) = @_;
    # verify input
    die "readCols2Ref: first arg (return value) must be empty hash ref\n" 
	unless ref($ret) eq "HASH" && scalar(keys(%$ret)) == 0;
    die "readCols2Ref: third arg (column list) must be non-empty array ref\n"
	unless ref($colsA) eq "ARRAY" && @$colsA>0;

    my @cols = @$colsA;

    open(IN,$f) or die "readCols couldn't read from '$f'. $!\n";
    
    my (%loc,$c,$cc,$line);    
    for $c (@cols) {
	if (ref($c) eq 'ARRAY') {
	    $loc{$c->[0]} = -1;
	    $ret->{$c->[0]} = [];
	} else {
	    $loc{$c} = -1;
	    $ret->{$c} = [];
	}
    }

    # find header
    do {
	$line = <IN>;
    } while (substr($line,0,1) eq '#');
    chomp $line;

    $_ = $line;
    my @spl = split;
    
    # parse header to find each target column
    my $flag=0;
    for my $i (0..$#spl) {
	for $c (@cols) {
	    if (ref($c) eq 'ARRAY') {
		$flag = 0;
		for $cc (@$c) {
		    if ($spl[$i] eq $cc) {		    
			$loc{$c->[0]} = $i;
			$flag=1;
			last;
		    }
		}
		last if $flag;
	    } else {
		if ($spl[$i] eq $c) {
		    $loc{$c} = $i;
		    last;
		}
	    }
	}
    }
    
    # verify that you have found all columns
    for $c (@cols) {
	if (ref($c) eq 'ARRAY') {
	    for $cc (@$c) {
		die "readCols2Ref: can't find '$cc' column in header of '$f'\n"
		    if $loc{$c->[0]} < 0;
	    } 
	} else {
	    die "readCols2Ref: can't find '$c' column in header of '$f'\n"
		if $loc{$c} < 0;
	}
    }


    while ($line = <IN>) {
	my $ent = {};
	next if substr($line,0,1) eq '#';
	chomp $line;
	$_ = $line;
	@spl = split;
	for $c (keys %loc) {
	    push(@{$ret->{$c}}, $spl[$loc{$c}])
	}
	push(@{$ret->{$wholeLine}}, $line) if defined $wholeLine;
    }
    close IN;
    return;
}

sub readCol {
    my @col;
    readColRef(\@col,@_);
    return @col;
}

sub readColRef {
    my ($ret, $file, $col) = @_;
    die "readColRef requires first argument to be an empty array ref" 
	if ref($ret) ne "ARRAY" || @$ret != 0;

    open(IN,$file);

    my $line;
    do {
	$line = <IN>;
    } while (substr($line,0,1) eq '#');
    chomp ($line);

    $_ = $line;
    my @spl = split;
    my $loc = -1;
    for my $i (0..$#spl) {
	$loc = $i if $spl[$i] eq $col;
    }
    die "readColRef: can't find $col in $file\n" if $loc < 0;

    while ($line = <IN>) {
	next if substr($line,0,1) eq '#';
	chomp($line);
	$_ = $line;
	@spl = split;
	die "wrong number of cols in $line\n" unless defined $spl[$loc];
	push(@$ret,$spl[$loc]);
    }

    close IN;
    return;
}


# return an array of sequences
# second arg is a PBR array-ref that we fill with names of sequences
sub readFastaRef {
    my ($seqs, $fastaFile, $chr, $nocheck, $mode) = @_;

    die "can't read fasta files!!";
    
    die "readFastaRef: first and optional third arg must be empty array ref\n" 
	if ref($seqs) ne 'ARRAY' ||  @$seqs > 0 || 
	@_>=3 && (ref($chr) ne 'ARRAY' || @$chr > 0);

    my $seqio = undef;
    #my $seqio = Bio::SeqIO->new( -file   => $fastaFile, -format => "fasta" );
    my $seqobj;
    while ($seqobj = $seqio->next_seq()) {
	push(@$seqs, uc($seqobj->seq()));
	if (!defined $nocheck && $seqs->[-1] =~ m/([^ACGT])/) {
	    die "found illegal nucleotide '$1' in ", $seqobj->display_id(),"\n";
	}
	    
	push(@$chr, $seqobj->display_id()) if defined $chr;
    }

    if (defined $mode) {
	if ($mode eq "checkWidth") { # do all seqs have the same length?
	    for my $i (1 .. $#$seqs) {
		die "readFasta: Seqs read from $fastaFile are not of the same length!\n"
		    if ((length($seqs->[$i]) - length($seqs->[$i-1]) ) != 0);
	    }
	} else {
	    warn "readFasta: unrecognized mode '$mode'\n";
	}
    }

    return;
}
sub readFasta { 
    my @seqs;
    readFastaRef(\@seqs, @_); 
    return @seqs;
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

# these subroutines help in handling fullProfile output

# read the output of fullProfile
# ret{seqID} = [ probs for all species ]
# note that prob{seqID} [ 2i ] is the revcom of [2i + 1]
sub readFullProfile {
    my ($profile, $file, $seqIDList) = @_;
    
    die "seqIDList must be an empty array ref" if defined $seqIDList &&
	( ref($seqIDList) ne 'ARRAY' || @$seqIDList != 0 );
    
    open my $IN, "<", $file or die "can't read $file. $!";
    
    my $seqID;
    my $profiles; # binding profiles of each mode on this probe
    while (my $line = <$IN>) {
	next if $line =~ /^#/;
	chomp $line;
	if ($line =~ /^>/) {
	    if ($seqID) { # unless this is the first probe in the file
		$profile->{$seqID} = $profiles;
		push @$seqIDList, $seqID if defined $seqIDList;
	    }
	    $seqID = substr($line, 1);
	    $profiles = [];
	} else {
	    push @$profiles, [ split /,/, $line ];
	}
    }
    $profile->{$seqID} = $profiles; # last one
    push @$seqIDList, $seqID if defined $seqIDList;
    return;
}

# for each probe, sum all binding modes, retaining position-information
# return $ret->{seqID} = [sum at 0, sum at 1, sum at 2,.... ];
sub sumProfile {
    my ($ret, $profile) = @_;

    for my $seqID (keys %$profile) {
	my @profiles = @{ $profile->{$seqID} };
	my $nProfiles = @profiles;
	my @sum;
	for my $i (0..$#{$profiles[0]}) {
	    my $s = 0;
	    for my $j (0..($nProfiles-1)) {
		$s += $profiles[$j][$i] // 
		    die "can't find profile{$seqID}[$j][$i]";
	    }
	    push @sum, $s;
	}
	$ret->{$seqID} = \@sum;
    }
    return;
}

# What is the size of the object?
# Input is the binding profile of a single probe
# assume that the zeroes at the end of the prob profile is objSize-1
# assume the number of zeroes at the end of each sequence doesn't change
sub findObjSize {
    my ($profile) = @_;

    my $i = -1;
    while($profile->[$i] == 0) { $i--; }

    return -$i;
}

"no text after this line"

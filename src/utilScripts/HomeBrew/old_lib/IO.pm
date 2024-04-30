package HomeBrew::IO;

use strict;
use warnings;
use List::Util qw(min);

require Exporter;

# checkExist(mode,path) dies on fail
# checkExist2(mode,path) returns undef on success, error on fail
# readList(file) = @firstElemInEachLine # at start of line recognized
# readList2(file) = @arr = [each item in each line]
# readList3(file) = @arr = [list of array refs, one for each column]
# readCol($file,$col) read named column from file, return @values
# readColRef(\@values, $file,$col) as above, but pass empty \@values 
# readCols(file, @cols) = @ret; $ret[i] = hash ref keyed by @cols
# readCols2(file, @cols) = %ret; $ret[$col] = array ref, as readCol($col)
# readFastaGlocke(file,$chr,$nocheck,$mode) second arg is PBR array ref.
# readFasta alias for above
# writeCols(file,\@cols, $header, $colPrintf, $preComments, $seqN);
## cols = [ [col1], [col2]...]
## header = sprintf "%12s%12s...", "col1", "col2"...
## colPrintf = "%12d%12.4e..." must have as many fields as there are cols
## precomments = "# these are some comments\n#that preced the data\n"
## seqN: if defined, add a seqN column 
### \- Note that you must account for seqN in colPrintf!
# writeFasta(file,seq,nameOfSeq,length=60) seq may be string or array ref
# writeHist(file, hist, binMin, header) # use with HomeBrew::Stats::histogram

# for legacy code: underscored names for subroutines:
# all are mere aliases for their counterparts
# read_list(file) = @firstElemInEachLine # at start of line recognized
# read_list_2(file) = @arr = [each item in each line]
# read_fasta_glocke(file,$chr,$mode) second arg is PBR array ref.

###
### TODO: 
### for each readColX add readColXRef which will take a ref argument
#### make the originals mere aliases for this version.


our @ISA = qw(Exporter);

#our @EXPORT = qw(checkExist);
our @EXPORT_OK = qw(checkExist checkExist2 readList readList2 readList3 readCol
    readColRef readCols readCols2 readCols2Ref readFastaGlocke readFasta 
    read_list read_list_2 read_fasta_glocke writeCols writeFasta writeHist);
    

our $VERSION = '0.01';

our $pathNotFound = -3;
our $emptyFile = -4;

# die unless path isn't found or if file is empty ('f' only)
sub checkExist {
  my ($type, $path) = @_;
  if ($type eq 'd') {
    if (! -d $path) {
      warn("dir not found: $path\n");
      die $pathNotFound;
    }
  }
  elsif ($type eq 'f') {
    if (! -f $path) {
      warn("file not found: $path\n");
      die $pathNotFound;
    }
    elsif (! -s $path) {
      warn("empty file: $path\n");
      die $emptyFile;
    }
  } else {
      warn "checkExist: unrecognized mode '$type'\n";
  }
  return;
}

# safer. return error value on fail.
sub checkExist2 {
  my ($type, $path) = @_;
  if ($type eq 'd') {
    if (! -d $path) {
      warn("dir not found: $path\n");
      return $pathNotFound;
    }
  }
  elsif ($type eq 'f') {
    if (! -f $path) {
      warn("file not found: $path\n");
      return $pathNotFound;
    }
    elsif (! -s $path) {
      warn("empty file: $path\n");
      return $emptyFile;
    }
  } else {
      warn "checkExist2: unrecognized mode '$type'\n";
  }
  return;
}


# This function reads in a list of items from a file (one item per line, rest of line ignored;
# '#' comment symbols recognized).
sub readList {
  my ($file) = @_;
  die "readList: undefined file\n" unless defined $file;
  my @tokens = ();
  my $cnt = 0;
  open(I,"$file") or (warn "readList can't read from $file: $!\n" && return 0); 
  while (<I>) {
    chomp;
    my @tt = split ' ';
    next if (@tt==0);
    next if (substr($tt[0],0,1) eq '#');
    my $token = $tt[0];
    push @tokens,$token;
    $cnt++;
  }
  close I;
  print "Read $cnt token(s) from file: $file ..\n";
  return @tokens;
}
sub read_list { return &readList; }

# This function reads in a list of items from a file 
# ret[line][entry] one entry per space separated item on line
# '#' comment symbols recognized if they appear first in a line
sub readList2 {
  my ($file) = @_;
  my @tokens = ();
  my $cnt = 0;
  open(I,"$file") or die "Can't open file: $!\n";
  while (<I>) {
    chomp;
    my @tt = split(/ +/,$_);
    next if (@tt==0);
    next if (substr($tt[0],0,1) eq '#');
    push @tokens,\@tt;
    $cnt++;
  }
  close I;
  print "Read $cnt token(s) from file: $file ..\n";
  return @tokens;
}
sub read_list_2 { return &readList2; }

# This function reads in a list of items from a file 
# ret[$column][$line] the entry on $line from $column
# '#' comment symbols recognized iff they appear first in a line
sub readList3 {
    my ($file) = @_;
    my @tokens = ();
    my $cnt = 0; # number of lines actually parsed
    open(I,"$file") or die "Can't open file: $!\n";
    my $i;
    while (<I>) {
	next if m/^#/; # skip commented lines
	chomp;
	my @tt = split(/ +/,$_);
	next if (@tt==0);
	for $i (0..$#tt) {
	    $tokens[$i] = [] if !defined $tokens[$i];
	    $tokens[$i][$cnt] = $tt[$i];
	}
	$cnt++;
    }
    close I;
    print "Read $cnt token(s) from file: $file ..\n";
    return @tokens;
}


# second arg is an array, but its elements may be array references.
#  cols = ([c11,c12],c2,[c31,c32...],...)
# the return value is an array of hash references:
#  ret[i] = { c11 -> val1, c2 -> val2, c31-> val3,...}
#   note that the FIRST element in the array ref [cn1,cn2] is always the key for valn
# the purpose of these array references is that sometimes a certain column
## in the source file may have been named inconsistently (e.g. filter was 
## briefly written as counts).  colsA = ["seqN",["filter","counts"]]
## yields a return value of 
sub readCols {
    my ($f,$colsA) = @_;
    my @cols = @$colsA;

    open(IN,$f) or die "readCols couldn't read from '$f'. $!\n";
    
    my (%loc,$c,$cc,$line);
    for $c (@cols) {
	if (ref($c) eq 'ARRAY') {
	    $loc{$c->[0]} = -1;
	} else {
	    $loc{$c} = -1;
	}
    }

    # find header
    do {
	$line = <IN>;
    } while (substr($line,0,1) eq '#');
    chomp $line;

    my @spl = split(/ +/," ".$line);
    
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
    
    for $c (@cols) {
	if (ref($c) eq 'ARRAY') {
	    for $cc (@$c) {
		die "can't find '$cc' column in header of '$f'\n" if $loc{$c->[0]} < 0;
	    } 
	} else {
	    die "can't find '$c' column in header of '$f'\n" if $loc{$c} < 0;
	}
    }

    my @ret = ();
    while ($line = <IN>) {
	my $ent = {};
	next if substr($line,0,1) eq '#';
	chomp $line;
	@spl = split(/ +/," ".$line);
	for $c (@cols) {
	    if (ref($c) eq 'ARRAY') {
		$ent->{$c->[0]} = $spl[$loc{$c->[0]}];	    
	    } else {
		$ent->{$c} = $spl[$loc{$c}];	    
	    }
	}
	# $ent->{line} = $line; unsafe, what if one of the columns is "line"?

	push(@ret,$ent);
    }
    close IN;
    return @ret;
}

sub readCols2 {
    my %ret;
    readCols2Ref(\%ret, @_);
    return %ret;
}

# as readCols EXCEPT it returns a hash of array refs instead of an array of 
# hash refs
sub readCols2Ref {
    my ($ret, $f,$colsA) = @_;
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

    my @spl = split(/ +/," ".$line);
    
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
		die "can't find '$cc' column in header of '$f'\n" if $loc{$c->[0]} < 0;
	    } 
	} else {
	    die "can't find '$c' column in header of '$f'\n" if $loc{$c} < 0;
	}
    }


    while ($line = <IN>) {
	my $ent = {};
	next if substr($line,0,1) eq '#';
	chomp $line;
	@spl = split(/ +/," ".$line);
	for $c (keys %loc) {
	    push(@{$ret->{$c}}, $spl[$loc{$c}])
	}
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

    my @spl = split(/ +/," ".$line);
    my $loc = -1;
    for my $i (0..$#spl) {
	$loc = $i if $spl[$i] eq $col;
    }
    die "can't filed $col in $file\n" if $loc < 0;

    while ($line = <IN>) {
	next if substr($line,0,1) eq '#';
	chomp($line);
	@spl = split(/ +/," ".$line);
	die "wrong number of cols in $line\n" unless defined $spl[$loc];
	push(@$ret,$spl[$loc]);
    }

    close IN;
    return;
}

sub readFasta {
    return &readFastaGlocke;
}

#This sub loads fasta seqs into an array WITHOUT BioPerl
#functions. It can optionally perform various checks on
#the seqs it read.
#WARNING: this sub is NOT thoroughly tested!
sub readFastaGlocke {
    my ($fasta_file,$chr,$nocheck,$mode) = @_;
    die "optional second arg to read_fasta_glocke must be array ref\n" if @_>1 && ref($chr) ne 'ARRAY';
    my @seqs = ();
    my $cnt = -1;
    #my $nlines = 0;
    open(II,"$fasta_file") or die "Can't open file: $!\n";
    while (<II>) {
	chomp;
	my @tt = split ' ';
	next if (@tt==0);
	next if substr($tt[0],0,1) eq '#';
	if (substr($tt[0],0,1) eq '>') {
	    if (defined $chr) {
		push(@$chr,substr($tt[0],1))
	    }
		
	    $cnt++;
	    next;
	}
	if (!defined $seqs[$cnt]) {
	    $seqs[$cnt] = uc($tt[0]); # make everything uppercase!!!
	}else{
	    $seqs[$cnt] .= uc($tt[0]); # make everything uppercase!!!
	}
	#$nlines++;
	#print "ahoy nlines = $nlines\n" unless $nlines % 1000;
	die "found bad base $& in $fasta_file'\n" 
	    if !(defined $nocheck && $nocheck) && $seqs[$cnt] =~ m/[^ACGT]/;
    }

    close II;
    if (defined $mode) {
	if ($mode eq "check_width") { # do all seqs have the same length?
	    for my $i (1 .. $#seqs) {
		die "Seqs read from $fasta_file are not of the same length!\n"
		    if ((length($seqs[$i]) - length($seqs[$i-1]) ) != 0);
	    }
	} else {
	    warn "called readFastaGlocke with unrecognized mode '$mode'\n";
	}
    }

    return @seqs;
}
sub read_fasta_glocke { return &readFastaGlocke; }

sub writeCols {
    my ($file, $cols, $header, $colPrintf, $preComments, $seqN) = @_;
    # verify input
    die "writeCols: cannot write to $file. $!\n" if -f $file && ! -w $file;
    die "writeCols: cols variable must be non-empty array ref\n" 
	unless ref($cols) eq "ARRAY" && @$cols>0;
    my $size;
    use Data::Dumper;
    for my $c (@$cols) {
	die "writeCols: cols variable must contain array refs\n" 
	    unless ref($c) eq "ARRAY";
	die "writeCols: cols must all be same size\n" 
	    if defined $size && @$c != $size;
	$size = @$c;
    }
    
    my $nFields = 0;
    $nFields++ while $colPrintf =~ m/% ?\+?-?0?#?(0x|0X|0b|0B)?\d*(\.\d*)?(c|s|d|u|o|x|e|f|g|X|E|G|b|B|p|n|i|D|U|O|F)/g;
    die "writeCols: colPrintf arg must have same number of fields as cols arg (+1 if seqN)\nfound $nFields, expected ", @$cols + ((defined $seqN)?1:0), "\n"
	unless $nFields==(@$cols + ((defined $seqN)?1:0));
    
    $header.="\n" unless $header =~ m/\n$/;
    $colPrintf.="\n" unless $colPrintf =~ m/\n$/;
    my @preComments;
    @preComments = split(/\n/, $preComments) 
	if defined $preComments && $preComments;
    for (@preComments) {
	$_.="\n"   unless /\n$/;
	$_ = "# $_" unless /^#/;
    }

    open (my $OUT, ">", $file) or 
	die "writeCols: couldn't write to $file. $!\n";
    
    print {$OUT} @preComments if @preComments;
    print {$OUT} $header;
    my ($i, @data);
    for $i (0..$#{$cols->[0]}) {
	@data = ();
	push(@data, $i+1) if defined $seqN;
	push(@data, $cols->[$_]->[$i]) for (0..$#$cols);
	printf {$OUT} $colPrintf, @data;
    }
    close $OUT;

    return;
}

sub writeFasta {
    my ($file, $seq, $nameOfSeq, $lineLength) = @_;
    # verify input
    die "cannot write to $file. $!\n" if -f $file && ! -w $file;
    die "empty sequence\n" unless ref($seq) || $seq; # if seq is empty string
    die "illegal sequence\n" if ref($seq) && (ref($seq) != "ARRAY" || @$seq<1);
    die "nameOfSeq cannot be blank\n" unless defined $nameOfSeq && $nameOfSeq;
    $lineLength //= 60;
    die "lineLength must be positive\n" if $lineLength<1;

    open( my $OUT, ">", $file) or die "couldn't write to $file.  $!\n";
    print {$OUT} ">$nameOfSeq\n";
    
    my $i = 0;
    my $m;
    if (ref($seq)) { # $seq is array ref
	while ($i < $#$seq) {
	    $m = min($i+$lineLength, $#$seq);
	    print {$OUT} join("", @$seq[$i..$m]), "\n";
	    $i+=$lineLength;
	}
    } else { # $seq is string
	# a more perly version that surely suffers for the hefty join
	#my $format = "(A$lineLength)*";
	#my @seq = unpack($format, $seq);
	#print {$OUT} join("\n", @seq);
	while ($i < length($seq)) {
	    $m = min($lineLength, length($seq) - $i);
	    print {$OUT} substr($seq,$i,$m), "\n";
	    $i+=$lineLength;
	}
    }
    close $OUT;

    return;
}

# args: 
#  $out: filename
#  $hist: histogram (array-ref)
#  $binMin: minimum to be in hist bin (array-ref)
#  $header: optional header to put in print-out
sub writeHist {
    my ($out,$hist,$binMin,$header) = @_;
    die "cannot write to $out\n" if -f $out && ! -w $out;
    die "printHist: expected array refs for 2nd and 3rd args\n" 
	unless ref($hist) eq "ARRAY" && ref($binMin) eq "ARRAY";
    die "printHist: size of binMin doesn't match size of hist\n" 
	unless $#$binMin == $#$hist;
    die "printHist: histogram has ", scalar(@$hist), " bin(s)\n" if @$hist<2;
    

    my $bin; # counts in bin
    my $binSize = $binMin->[1] - $binMin->[0];

    open (OUT,">$out") or die "couldn't write to $out.  $!\n";
    print OUT $header if defined $header;
    printf OUT "%10s%10s%10s\n", "binMin","binMid","count";    
    for my $i (0..$#$hist) {
	printf OUT ("%10.3f%10.3f%10d\n", $binMin->[$i], 
		    $binMin->[$i]+$binSize/2, $hist->[$i]);
    }
    close OUT;
    return;
}

"no text after this line"

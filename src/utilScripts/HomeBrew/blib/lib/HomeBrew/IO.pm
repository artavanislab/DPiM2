package HomeBrew::IO;

use strict;
use warnings;
use List::Util qw(min);
use Text::ParseWords;

require Exporter;

# checkExist(mode,path) dies on fail
# checkExist2(mode,path) returns undef on success, error on fail
# readList(file) = @firstElemInEachLine # at start of line recognized
# readList2(file) = @arr = [each item in each line]
# readList3(file) = @arr = [list of array refs, one for each column]
# readCol($file,$col) read named column from file, return @values
# readColRef(\@values, $file,$col) as above, but pass empty \@values
# /--
# readCols(file, @cols, $wholeLine) = @ret; $ret[i] = hash ref keyed by @cols
# readColsRef(\@ret, file, @cols, $wholeLine); as above.
# readCols2(file, @cols, $wholeLine) = %ret; $ret[$col] = array ref
# readCols2Ref(\%values, file, @cols, $wholeLine) as above but pass empty \%ret
# \_.> if whole line $wholeLine is defined, it names a column that holds
####   the whole (chomped) line of the file
# readFasta(file,$chr,$nocheck,$mode) second arg is PBR array ref.
# readFastaRef as above but first arg is \@seqs;
# readSAGE(file, countCol, libSize) like readCols(file, qw(Tag $countCol source pos. strand gene locus desc.))
## second arg is optional, reports normalization
# readSAGERef(\@sage, file, libSize) as readSAGE but ret is passed by ref
# writeCols(file,\@cols, $header, $colPrintf, $preComments, $seqN);
## cols = [ [col1], [col2]...]
## header = sprintf "%12s%12s...", "col1", "col2"...
## colPrintf = "%12d%12.4e..." must have as many fields as there are cols
## precomments = "# these are some comments\n#that preced the data\n"
## seqN: optionally (ifdef) add a seqN column starting at 1 
##   \_.> Note that you must account for seqN in colPrintf!
# writeFasta(file,seq,nameOfSeq,length=60) seq may be string or array ref
# writeHist(file, hist, binMin, preComments) # use with HomeBrew::Stats::histogram
# read2DArray(ret, file); # more notes at sub definition
# write2DArray(file, struct, keyFormat, dataFormat, (optional: preComments)) 
#  struct must be array of hash or array of array
#  all elements of struct must be formatted like one another
# readHeader($file, \$raw) return an array of all the columns listed in the header.  last arg optionally gets set to the raw string header

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
our @EXPORT_OK = qw(checkExist checkExist2 readList readListRef readList2 
    readList3 readCol readColRef readCols readColsRef readCols2 readCols2Ref 
    readFasta readFastaRef readSAGE readSAGERef 
    writeCols writeFasta writeMultiFasta writeHist
    read2DArray write2DArray readHeader);
    

our $VERSION = '0.01';

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

# safer. return error value on fail.
sub checkExist2 {
  my ($type, $path) = @_;
  die "checkExis2t: path variable not defined\n" unless defined $path;
  if ($type eq 'd') {
    if (! -d $path) {
      warn("checkExist2: dir not found: $path\n");
      return $pathNotFound;
    }
  }
  elsif ($type eq 'f') {
    if (! -f $path) {
      warn("checkExist2: file not found: $path\n");
      return $pathNotFound;
    }
    elsif (! -s $path) {
      warn("checkExist2: empty file: $path\n");
      return $emptyFile;
    }
  } else {
      warn "checkExist2: unrecognized mode '$type'\n";
  }
  return;
}

# any input routine that splits lines will call this
# the purpose is essentially to wrap the quote checking rather than repeat
#   the same code
# it returns a subroutine so that callers will only check whether it needs to
#   check for quotes on every line
sub lineSplitter {
    my ($checkQuotes) = @_;
    print "checkQuotes = ";
    if (defined $checkQuotes) {
	print "$checkQuotes\n";
    } else {
	print "NOT DEFINED\n";
    }
    my $ret;
    if ($checkQuotes) {
	print "\tquotewords style\n";
	$ret = sub {
	    my $line = shift;
	    chomp $line;
	    my @ret = quotewords('\s+', 'keep', $line);
	    shift @ret if length($ret[0]) < 1;
	    return @ret;
	}
    } else {
	$ret = sub {
	    my $line = shift;
	    $_ = $line;
	    chomp;
	    return split;
	};
    }
    return $ret;
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
  open(I,"$file") or die "Can't open file: $!\n";
  while (<I>) {
    next if /^#/;
    chomp;
    my @tt = split;
    next if (@tt==0);
    push @tokens, \@tt;
    $cnt+=@tt;
  }
  close I;
  print "Read $cnt token(s) from file: $file ..\n" if $verbose;
  return @tokens;
}

# This function reads in a list of items from a file 
# ret[$column][$line] the entry on $line from $column
# '#' comment symbols recognized iff they appear first in a line
sub readList3 {
    my ($file, $verbose) = @_;
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
## Note that the result will be named 'filter' in either case.
sub readColsRef {
    my ($ret, $f, $cols, $wholeLine, $checkQuotes) = @_;

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

    my $splitter = lineSplitter($checkQuotes);
    
    # find header
    do {
	$line = <$IN>;
    } while (substr($line,0,1) eq '#');
    chomp $line;

    my @spl = $splitter->($line, $checkQuotes);
    sanitizeHeader(\@spl);
    
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
	next if length($line) < 2;
	@spl = $splitter->($line, $checkQuotes);
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
    my ($ret, $f,$colsA, $wholeLine, $checkQuotes) = @_;
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

    my $splitter = lineSplitter($checkQuotes);

    # find header
    do {
	$line = <IN>;
    } while (substr($line,0,1) eq '#');
    chomp $line;

    my @spl = $splitter->($line, $checkQuotes);
    sanitizeHeader(\@spl);
    
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
	next if length($line) < 2;

	@spl = $splitter->($line, $checkQuotes);
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
    my ($ret, $file, $col, $checkQuotes) = @_;
    die "readColRef requires first argument to be an empty array ref" 
	if ref($ret) ne "ARRAY" || @$ret != 0;

    my $splitter = lineSplitter($checkQuotes);
    
    open(IN,$file);

    my $line;
    do {
	$line = <IN>;
    } while (substr($line,0,1) eq '#');

    my @spl = $splitter->($line);
    #die join "\t", @spl;
    
    sanitizeHeader(\@spl);
    my $loc = -1;
    for my $i (0..$#spl) {
	$loc = $i if $spl[$i] eq $col;
    }
    die "readColRef: can't find $col in $file\n" if $loc < 0;

    while ($line = <IN>) {
	next if substr($line,0,1) eq '#';
	next if length($line) < 2;
	@spl = $splitter->($line);
	die "wrong number of cols in $line\n" unless defined $spl[$loc];
	push(@$ret, $spl[$loc]);
    }

    close IN;
    return;
}

# return an array of sequences
# second arg is a PBR array-ref that we fill with names of sequences
sub readFastaRef {
    #use Bio::SeqIO;
    die "HomeBrew::IO::readFasta disabled due to Bio::SeqIO";
    my ($seqs, $fastaFile, $chr, $nocheck, $mode) = @_;

    die "readFastaRef: first and optional third arg must be empty array ref\n" 
	if ref($seqs) ne 'ARRAY' ||  @$seqs > 0 || 
	@_>=3 && (ref($chr) ne 'ARRAY' || @$chr > 0);

    #my $seqio = Bio::SeqIO->new( -file   => $fastaFile, -format => "fasta" );
    my $seqio;
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

sub readSAGE {
    my @sage;
    readSAGERef(\@sage, @_);
    return @sage;
}

sub readSAGERef {
    my ($ret, $file, $countCol, $libSize) = @_;
    
    die "IO::readSAGERef first arg must be ref to array\n" 
	if ref($ret) ne 'ARRAY';

    use List::Util qw(min);

    open(my $IN, "<", $file) 
	or die "IO::readSAGERef couldn't read from $file. $!\n";
    
    my $line;
    # find header
    do {
	$line = <$IN>;
    } while ($line =~ /^#/ || $line =~ /^ +\*/ || length($line) <= 1 
	     || $line =~ /^ +$/ );
    chomp $line;
    
    # verify format
    # name of second column varies from file to file depending on assay
    # you may specify it to assure yourself that you've got the right assay
    # or you may leave it blank. 
    # in either case, the output will name this column 'counts'
    my @cols = ('Tag', $countCol, qw(source pos. strand gene locus desc.));

    $_ = $line;
    my @spl = split;
    
    #use Data::Dumper;
    #die "IO::readSAGERef header does not appear to have the desired format\n"
    #, "line = '$line'\n", Dumper(\@spl) if @spl!=@cols;

    for my $i (0..$#spl) {
	die "IO::readSAGERef can't find $cols[$i] at $i in $file\n" 
	    if defined $cols[$i] && $spl[$i] ne $cols[$i] && $cols[$i];
	# cols[i] will be defined and evaluate true unless
	#   i==1 && $countCol is not defined or empty (or ill defined i guess)
    }
    $cols[1] = 'counts';

    $line = <$IN>;
    chomp $line;
    die "IO::readSAGERef line after header does not begin 'Lib. size...' (found '$line')\n"
	if $line !~ /^Lib\. size( +| *\t *)\d+/;
    if (ref($libSize) eq 'SCALAR') {
	$_ = $line;
	@spl = split;
	$$libSize = $spl[-1];
    }

    # now parse data
    while (<$IN>) {
	next if /^#/;
	next unless length($_);
	@spl = split;
	my $ent = {};
	# if @spl < @cols, not all cols are filled in
	# if @spl > @cols, this is because of a block-text description
	for my $i (0..min($#spl, $#cols)) {
	    $ent->{$cols[$i]} = $spl[$i];
	}
	if (@spl>@cols) { # description entry contains spaces
	    $ent->{$cols[-1]} .= ' '.join(" ", @spl[@cols..$#spl]);
	}
	push(@$ret, $ent);
    }
    
    return;
}

sub writeCols {
    my ($file, $cols, $header, $colPrintf, $preComments, $seqN) = @_;
    # verify input
    die "writeCols: cannot write to $file. $!\n" if -f $file && ! -w $file;
    die "writeCols: cols variable must be non-empty array ref\n" 
	unless ref($cols) eq "ARRAY" && @$cols>0;
    my $size;

    for my $c (@$cols) {
	die "writeCols: cols variable must contain array refs\n" 
	    unless ref($c) eq "ARRAY";
	die "writeCols: cols must all be same size\n" 
	    if defined $size && @$c != $size;
	$size = @$c;
    }
    
    my $nFields = 0;
    $nFields++ while $colPrintf =~ m/% ?\+?-?0?#?(0x|0X|0b|0B)?\d*(\.\d*)?(c|s|d|u|o|x|e|f|g|X|E|G|b|B|p|n|i|D|U|O|F)/g;
    die "writeCols: colPrintf arg must have same number of fields as cols arg (+1 if seqN)\n\tfound $nFields, expected ", @$cols + ((defined $seqN)?1:0), "\n"
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

    my $OUT;
    if (grep /^$file$/, qw(stdout STDOUT)) {
	$OUT = \*STDOUT;
    } else {
	open ($OUT, ">", $file) or 
	    die "writeCols: couldn't write to $file. $!\n";
    }
    
    print $OUT @preComments if @preComments;
    print $OUT $header;
    my ($i, @data);
    for $i (0..$#{$cols->[0]}) {
	@data = ();
	push(@data, $i+1) if defined $seqN;
	push(@data, $cols->[$_]->[$i]) for (0..$#$cols);
	printf $OUT $colPrintf, @data;
    }
    close $OUT;

    return;
}

sub writeFasta {
    my ($file, $seq, $nameOfSeq, $lineLength) = @_;
    # verify input
    die "cannot write to $file. $!\n" if -f $file && ! -w $file;
    die "empty sequence\n" unless ref($seq) || $seq; # if seq is empty string
    die "illegal sequence\n" if ref($seq) && (ref($seq) ne "ARRAY" || @$seq<1);
    die "nameOfSeq cannot be blank\n" unless defined $nameOfSeq && $nameOfSeq;
    if (ref($nameOfSeq) eq 'ARRAY') {
	die "if nameOfSeq is an array, it must be the same size as seq!\n"
	    if ref($seq) ne 'ARRAY' || @$seq != @$nameOfSeq;
    }
    $lineLength //= 60;
    die "lineLength must be positive\n" if $lineLength<1;

    my (@seqs, @names);
    if (ref($nameOfSeq) eq 'ARRAY') {
	@seqs = @$seq;
	@names = @$nameOfSeq;
    } else {
	@seqs = ($seq);
	@names = ($nameOfSeq);
    }

    open my $OUT, ">", $file or die "couldn't write to $file.  $!\n";

    for my $i (0..$#seqs) {
	my $seq = $seqs[$i];
	my $nameOfSeq = $names[$i];
	print $OUT ">$nameOfSeq\n";
	
	my $i = 0;
	my $m;
	if (ref($seq)) { # $seq is array ref
	    while ($i < $#$seq) {
		$m = min($i+$lineLength-1, $#$seq);
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
    }
    close $OUT;

    return;
}

sub writeMultiFasta {
    my ($file, $seqs, $names, $lineLength) = @_;
    # verify input
    die "cannot write to $file. $!\n" if -f $file && ! -w $file;
    die "empty sequence\n" unless ref($seqs) eq 'ARRAY' && @$seqs > 0;
    die "expect seqs to be a ref to array of strings" if ref($seqs->[0]);
    die "names cannot be blank\n" unless ref($names) eq 'ARRAY' && @$names > 0;
    die "names and seqs must be same size" if @$names != @$seqs;

    $lineLength //= 60;
    die "lineLength must be positive\n" if $lineLength<1;

    open my $OUT, ">", $file or die "couldn't write to $file.  $!\n";

    for my $i (0..$#$seqs) {
	my $seq = $seqs->[$i];
	my $nameOfSeq = $names->[$i];
	print $OUT ">$nameOfSeq\n";
	
	my $i = 0;
	my $m;
	if (ref($seq)) { # $seq is array ref
	    while ($i < $#$seq) {
		$m = min($i+$lineLength-1, $#$seq);
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
    }
    close $OUT;

    return;
}

# args: 
#  $out: filename
#  $hist: histogram (array-ref)
#  $binMin: minimum to be in hist bin (array-ref)
#  $preComments: optional preComments to put in print-out
sub writeHist {
    my ($out,$hist,$binMin,$preComments) = @_;
    die "cannot write to $out\n" if -f $out && ! -w $out;
    die "printHist: expected array refs for 2nd and 3rd args\n" 
	unless ref($hist) eq "ARRAY" && ref($binMin) eq "ARRAY";
    die "printHist: size of binMin doesn't match size of hist\n" 
	unless $#$binMin == $#$hist;
    die "printHist: histogram has ", scalar(@$hist), " bin(s)\n" if @$hist<2;
    

    my $bin; # counts in bin
    my $binSize = $binMin->[1] - $binMin->[0];

    $preComments.="\n"   unless $preComments=~ m/\n$/;
    $preComments = "# $preComments" unless $preComments=~ m/^#/;

    open (OUT,">$out") or die "couldn't write to $out.  $!\n";
    print OUT $preComments if defined $preComments;
    printf OUT "%15s%15s%10s\n", "binMin","binMid","count";    
    for my $i (0..$#$hist) {
	printf OUT ("%15.3f%15.3f%10d\n", $binMin->[$i], 
		    $binMin->[$i]+$binSize/2, $hist->[$i]);
    }
    close OUT;
    return;
}

# return format is ret[column] = {rowName1 => X1, rowName2 => X2,....}
# colNames are the same for every column
# input file must be formatted as follows
#       0   1   2 ...
# colA a0  a1  a2 ...
# colB b0  b1  b2 ...
# ...
sub read2DArray {
    my ($ret, $file, $checkQuotes) = @_;
    die "IO::read2DArray expected first arg to be empty array ref\n"
	unless ref($ret) eq 'ARRAY' && @$ret==0;
    
    open my $IN, "<", $file 
	or die "IO::read2DArray couldn't read from $file. $!\n";

    my $splitter = lineSplitter($checkQuotes);

    #my $line;
    # find header
    do {
	$_ = <$IN>;
    } while (/^#/);
    my @spl = $splitter->($_);
    my $nElem = @spl;
    #push(@$ret, {}) for 1..$nElem;

    while (<$IN>) {
	next if /^#/;
	next if length($_) == 0;
	@spl = $splitter->($_);
	die "IO::read2DArray found wrong number of elements in line $_\n" 
	    unless @spl == $nElem+1;
	my $k = shift(@spl);	
	for my $i (0..$#spl) {
	    $ret->[$i]{$k} = $spl[$i];
	}
    }
    return;
}

sub write2DArray {
    my ($out, $data, $keyFormat, $dataFormat, $preComments) = @_;
    die "IO::write2DHash second arg must be array ref\n" 
	unless ref($data) eq "ARRAY";
    die "IO::write2DHash formats must be valid printf codes '%12d/s/e'"
	if defined $keyFormat && $keyFormat !~ m/^%\d+/ 
	&& $keyFormat !~ m/^%\d+/;

    # verify that formats have the same length 
    my $nspaces;
    {
	($nspaces) = ($keyFormat =~ m/(\d+)/);
	my ($nspaces2) = ($dataFormat =~ m/(\d+)/);
	die "key format and data format must be the same length" 
	    unless $nspaces == $nspaces2;
    }

    

    # verify that elements of $data are all the same format
    if (ref($data->[0]) eq "ARRAY") {
	my $size = @{$data->[0]};
	my $i=0;
	for (@$data) {
	    die "IO::write2DHash $i'th element of data arg has wrong length.",
	        "  got ", scalar(@$_), ", expected $size.\n" 
		    unless @$_==$size;
	    $i++;
	}
    } elsif (ref($data->[0]) eq "HASH") {
	my @k = keys %{$data->[0]};
	my %k = map {$_ => 1} @k;
	my $i=0;
	for my $h (@$data) {
	    die "IO::write2DHash $i'th hash in data arg has wrong number of ".
	        " keys\n" if scalar(keys %$h)!=@k;
	    for (keys %$h) {
		die "IO::write2DHash $i'th hash in data arg has wrong key $_\n" 
		    unless $k{$_};
	    }	    
	    $i++;
	}
    } else {
	die "IO::write2DHash second arg must be 2D array (array->array/hash)" 
    }

    my $spaces = " " x $nspaces;
    
    open my $OUT, ">", $out or die "couldn't write to $out. $!\n";

    if (defined $preComments) {
	print $OUT $preComments;
	print $OUT "\n" unless $preComments =~ /\n$/;
    }
    print $OUT $spaces;
    printf $OUT $keyFormat, $_ for (0..$#$data);
    print $OUT "\n";
    if (ref($data->[0]) eq "ARRAY") {
	for my $j (0..$#{$data->[0]}) {
	    printf $OUT $keyFormat, $j;
	    for my $i (0..$#$data) {
		printf $OUT $dataFormat, $data->[$i][$j];
	    }
	    print $OUT "\n";
	}
    } elsif (ref($data->[0]) eq "HASH") {
	for my $k (sort keys %{$data->[0]}) {
	    printf $OUT $keyFormat, $k;
	    for my $i (0..$#$data) {
		printf $OUT $dataFormat, $data->[$i]{$k};
	    }
	    print $OUT "\n";
	}
    } 
    close $OUT;
    return;
}

sub readHeader {
    my ($file, $raw) = @_;
    
    open my $IN, "<", $file 
	or die "readHeader: couldn't open file '$file'. $!\n";
    do {
	$_ = <$IN>;
    } while (/^#/);
    chomp;
    my @spl = split;
    
    $$raw = $_ if ref($raw) eq 'SCALAR';
    return @spl;
}

# remove quotes from header columns
# accounts for R's write.table behavior
sub sanitizeHeader {
    my ($header) = @_;
    $_ =~ s/"//g for @$header;
    return;
}

"no text after this line"

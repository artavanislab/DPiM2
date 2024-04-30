package HomeBrew::Stats;

use strict;
use warnings;
use List::Util qw(sum min max);

require Exporter;

# mean(list) return mean
# sumsquare(list) return sum of squares
# sd(list) return standard deviation
# ssd(list) return sample standard deviation

# histogram( data => \@data, hist => \@ret, 
#  (optional: max=> $max, min => $min, 
#    ( binSize => $binSize // 1 or nBins => $nBins), 
#    binMin => \@retBinMin, mode=> int (no other modes yet)
#    bottomOverflow => \$bottomOverflow, topOverflow => \$topOverflow
#  )
# )     
# getHistTails(\@counts, $fraction) # find bottom and top fraction of hist.
##  return (bottomIndex, topIndex);
####  bottomIndex is the largest index within the bottom $fraction    
# histMean(\@counts (optional \@binMid)) returns index (or real value)
# histSD(\@counts, (optional \@binMid)) 

# make a 2D histogram from an array, array[X] = Y
# hist2D_array( data => \@data, hist => \@ret 
#  (optional: filter=> \@, intY => 1 // 0, minX => $, maxX => $, minY => $, 
#   maxY $, binMinX => \@, binMinY => \@, (binSizeX => $ or nBinsX => $), 
#   (binSizeY => $ or nBinsY => $), bottomOverflowX => \$, topOverflowX => \$, 
#   bottomOverflowY => \$, topOverflowY => \$, 
#  ) 
# )

# make a 2D histogram from a hash of array refs
# hist2D_hash( data => {x => [], 'y' => []}, hist => \@ret
#  (optional: , x=>'x', y=>'y', filter=> \@, intX => 1 // 0, intY => 1 // 0, 
#   minX => $, maxX => $, minY => $, maxY $, binMinX => \@, binMinY => \@, 
#   (binSizeX => $ or nBinsX => $), (binSizeY => $ or nBinsY => $), 
#   bottomOverflowX => \$, topOverflowX => \$, 
#   bottomOverflowY => \$, topOverflowY => \$, 
#  ) 
# )

our @ISA = qw(Exporter);

our @EXPORT_OK = qw( mean sd ssd sumsquare histogram  getHistTails 
                     histMean histSD hist2D_array hist2D_hash);

our $VERSION = 0.01;

our $argNotDefined = -1; # required arg has not been passed
our $argIllDefined = -2; # if arg is not the correct reference type or 
#                        # is/isn't empty when it should(n't) be

sub mean {
    if (ref($_[0]) eq "ARRAY") {
	my $list = shift;
	die "not enough elements in list to get mean\n"
	    if  @$list < 1;
	return sum(@$list)/@$list;
    } else {
	die "not enough elements in list to get mean\n"
	    if  @_ < 1;
	return sum(@_)/@_;
    }
}

sub sumSquare {
    my $sum=0;
    if (ref($_[0]) eq "ARRAY") {
	my $list = shift;
	$sum += $_*$_ for (@$list);
    } else { 
	$sum += $_*$_ for (@_);
    }
    return $sum;
}

sub sd {
    my @args = @_;
    my $mean = &mean; # messes @_ up.
    my $sumSquare = sumSquare(@args); 

    my $size;
    if (ref($args[0]) eq "ARRAY") {
	$size = scalar(@{$args[0]});
    } else {
	$size = scalar(@_);
    }

    return sqrt($sumSquare/$size - $mean*$mean);
}

sub ssd {
    my @args = @_;
    my $sumSquare = &sumSquare;

    my ($sum,$size);
    if (ref($args[0]) eq "ARRAY") {
	$sum = sum(@{$args[0]});
	$size = scalar(@{$args[0]});
    } else {
	$sum = sum(@args);
	$size = scalar(@args);
    }
    die "not enough elements in list to get sample standard deviation\n"
	if $size < 2;

    return sqrt( ($sumSquare - $sum*$sum/$size)/($size-1) );
}

# histogram( data => \@data, hist => \@ret, 
#  (optional: max=> $max, min => $min, 
#    ( binSize => $binSize // 1 or nBins => $nBins), 
#    binMin => \@retBinMin, mode=> int (no other modes yet)
#    bottomOverflow => \$bottomOverflow, topOverflow => \$topOverflow
#  )
# )     
sub histogram {
    my %args =  @_;

    # first check for unrecognized args
    my @keys = qw( data hist max min binSize nBins binMin mode
                   bottomOverflow topOverflow filter);
    my %keys = map {$_ => 1} @keys;
    for (keys(%args)) {
	die "unrecognized arg '$_'\n" unless $keys{$_};
    }
    
    # verify mandatory args:
    if (!exists $args{data}) {
	warn "histogram: called histogram without defining the data (input) argument!\n";
	return $argNotDefined;
    }
    if (ref($args{data}) ne "ARRAY") {
	warn "histogram: data argument is not an array ref\n";
	return $argIllDefined;
    }
    if (@{$args{data}} < 2) {
	warn "histogram: data argument has less than 2 elems\n";
	return $argIllDefined;
    }

    if (!exists $args{hist}) {
	warn "histogram: called histogram without defining the hist (output) argument!\n";
	return -1;
    }
    if (ref($args{hist}) ne "ARRAY") {
	warn "histogram: hist argument is not an array ref\n";
	return $argIllDefined;
    }
    if (@{$args{hist}} > 0) {
	warn "histogram: hist argument is not empty\n";
	return $argIllDefined;
    }

    # verify optional args:
    if (exists $args{binMin}) {
	if (ref($args{binMin}) ne "ARRAY") {
	    warn "histogram: binMin argument is not an array ref\n";
	    return $argIllDefined;
	}
	if (@{$args{binMin}} > 0) {
	    warn "histogram: binMin argument is not empty\n";
	    return $argIllDefined;
	}
    }

    if (exists $args{binSize} && $args{binSize} <= 0) { # defined by default
	warn "histogram: illegal binSize<=0 ", $args{binSize}, "\n";
	return $argIllDefined;
    }
    if (exists $args{nBins}) {
	if ($args{nBins}<2 || $args{nBins} != int($args{nBins})) {
	    warn "histogram: illegal value of nBins: ", $args{nBins}, "\n";
	    return $argIllDefined;
	}
	if (exists $args{binSize}) {
	    warn "histogram: WARNING: binSize and nBins both defined.  nBins overrides!\n";
	}
    }
    if (exists $args{bottomOverflow} && 
	ref($args{bottomOverflow}) ne "SCALAR") 
    {
	warn "histogram: bottomOverflow must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{topOverflow} && 
	ref($args{topOverflow}) ne "SCALAR") 
    {
	warn "histogram: topOverflow must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{min} && exists $args{max} && $args{min}>=$args{max}) {
	warn "histogram: min must be less than max\n";
	return $argIllDefined;
    }
    # done verifying.  phew!

    my $data = $args{data};
    my $hist = $args{hist};

    # initialize optional arguments
    my $binMin = $args{binMin} // [];
    my $bottomOverflow = 0;
    my $topOverflow = 0;
    my $min = $args{min} // min(@$data);
    my $max = $args{max} // max(@$data);
    if (exists $args{mode} && $args{mode} eq "int") { 
	$min-=0.5;
	$max+=0.5;
    }
    my $binSize = $args{binSize} // 1;
    if (exists $args{nBins}) { 
	$binSize = ($max-$min)/$args{nBins};
    }

    # note loop condition prevents floating point errors
    for (my $thisBin = $min; $thisBin < $max - $binSize*0.01;
	 $thisBin+=$binSize) {
	push(@$binMin, $thisBin);
    }

    die "Stats::histogram failed to create the right number of bins.  "
	, "nBins arg = ", $args{nBins}, ", size of binMin = ", @$binMin+0, "\n" 
	if exists $args{nBins} && $args{nBins} != @$binMin;

    my $bin;
    for my $d (@$data) {

	if ($d<$min) {
	    $bottomOverflow++;
	    next;
	} elsif ($d>$max) {
	    $topOverflow++;
	    next;
	}
	$bin = int(($d-$min)/$binSize);
	$bin-- if $bin>$#$binMin; # in case you forgot to use int mode
	$hist->[$bin]++;
    }
    
    $hist->[$_] //= 0 for 0..$#$binMin; # set uninitialized entries to zero
    
    warn "lingering bug?  binMin and hist are not the same size\n" 
	unless @$binMin == @$hist;
    
    ${$args{bottomOverflow}} = $bottomOverflow;
    ${$args{topOverflow}} = $topOverflow;
    
    return;
}

# getHistTails(\@counts, $fraction) # find bottom and top fraction of hist.
sub getHistTails {
    my ($counts, $fraction) = @_;
    die "Stats::getHistTails expects first arg to be array ref\n" 
	unless defined $counts && ref($counts) eq 'ARRAY';
    die "Stats::getHistTails expects second arg to be scalar between zero and one half\n"
	if ref($fraction) ne '' || $fraction < 0 || $fraction > 0.5;

    my $total = sum(@$counts);
    my $i = 0;
    my $sum = 0;
    while ($sum < $fraction * $total) {
	$sum += $counts->[$i];
	$i++;
    }
    my $bottomIndex = $i-1;
    $i = $#$counts;
    $sum = 0;
    while ($sum < $fraction * $total) {
	$sum += $counts->[$i];
	$i--;
    }
    my $topIndex = $i+1;
    return ($bottomIndex, $topIndex);
}

sub histMean {
    my ($counts, $binMid) = @_;
    die "Stats::histMean expects first arg to be array ref\n" 
	unless defined $counts && ref($counts) eq 'ARRAY';
    die "Stats::histMean expects optional second arg to be array ref\n"
	if defined $binMid && ref($binMid) ne 'ARRAY';
    
    my $binz;
    if (defined $binMid) {
	$binz = $binMid;
    } else {
	$binz = [0..$#$counts];
    }

    my $sum=0;
    for my $i (0..$#$counts) {
	$sum += $counts->[$i] * $binz->[$i];
    }
    return $sum / sum(@$counts);
}

sub histSD {
    my ($counts, $binMid) = @_;
    die "Stats::histSD expects first arg to be array ref\n" 
	unless defined $counts && ref($counts) eq 'ARRAY';
    die "Stats::histSD expects optional second arg to be array ref\n"
	if defined $binMid && ref($binMid) ne 'ARRAY';
    
    my $binz;
    if (defined $binMid) {
	$binz = $binMid;
    } else {
	$binz = [0..$#$counts];
    }

    my $sum1=0;
    my $sum2=0;
    for my $i (0..$#$counts) {
	$sum1 += $counts->[$i] * $binz->[$i];
	$sum2 += $counts->[$i] * $binz->[$i] * $binz->[$i];
    }
    my $sum3 = sum(@$counts);
    $sum1 /= $sum3;
    $sum2 /= $sum3;
    return sqrt($sum2 - $sum1*$sum1);
}

# make a 2D histogram from an array, array[X] = Y
# hist2D_array( data => \@data, hist => \@ret 
#  (optional: filter=> \@, intY => 1 // 0, minX => $, maxX => $, minY => $, 
#   maxY $, binMinX => \@, binMinY => \@, (binSizeX => $ or nBinsX => $), 
#   (binSizeY => $ or nBinsY => $), bottomOverflowX => \$, topOverflowX => \$, 
#   bottomOverflowY => \$, topOverflowY => \$, 
#  ) 
# )
sub hist2D_array {
    my %args = @_;

    # first check for unrecognized args
    my @keys = qw( data hist filter intY maxX minX binSizeX nBinsX binMinX
                   bottomOverflowX topOverflowX
                   maxY minY binSizeY nBinsY binMinY 
                   bottomOverflowY topOverflowY
                   );
    my %keys = map {$_ => 1} @keys;
    for (keys(%args)) {
	die "hist2D_array: unrecognized arg '$_'\n" unless $keys{$_};
    }
    
    # verify mandatory args:
    if (!exists $args{data}) {
	warn "hist2D_array: data (input) argument undefined\n";
	return $argNotDefined;
    }
    if (ref($args{data}) ne "ARRAY") {
	warn "hist2D_array: data argument is not an array ref\n";
	return $argIllDefined;
    }
    if (@{$args{data}} < 2) {
	warn "hist2D_array: data argument has less than 2 elems\n";
	return $argIllDefined;
    }

    if (!exists $args{hist}) {
	warn "hist2D_array: hist argument undefined\n";
	return $argNotDefined;
    }
    if (ref($args{hist}) ne "ARRAY") {
	warn "hist2D_array: hist argument is not an array ref\n";
	return $argIllDefined;
    }
    if (@{$args{hist}} > 0) {
	warn "hist2D_array: hist argument is not empty\n";
	return $argIllDefined;
    }

    # verify optional args:
    if (exists $args{filter}) {
	if (ref($args{filter}) ne 'ARRAY') {
	    warn "hist2D_array: filter arg is not array ref\n";
	    return $argIllDefined;
	}
	my $fSize = @{$args{filter}} +0;
	my $dSize = @{$args{data}} +0;
	if ($fSize != $dSize) {
	    warn "hist2D_array: size of filter ( $fSize ) does not match size "
		, "of data $dSize\n";
	    return $argIllDefined;
	}
    }

    if (exists $args{binMinX}) {
	if (ref($args{binMinX}) ne "ARRAY") {
	    warn "hist2D_array: binMinX argument is not an array ref\n";
	    return $argIllDefined;
	}
	if (@{$args{binMinX}} > 0) {
	    warn "hist2D_array: binMinX argument is not empty\n";
	    return $argIllDefined;
	}
    }

    if (exists $args{binSizeX} && $args{binSizeX} <= 0) { # defined by default
	warn "hist2D_array: illegal binSizeX<=0: '", $args{binSizeX}, "'\n";
	return $argIllDefined;
    }
    if (exists $args{nBinsX}) {
	if ($args{nBinsX}<2 || $args{nBinsX} != int($args{nBinsX})) {
	    warn "hist2D_array: illegal value of nBinsX: ", $args{nBinsX}, "\n";
	    return $argIllDefined;
	}
	if ($args{nBinsX} > @{ $args{data} }) {
	    warn "hist2D_array: nBinsX is large than input array\n";
	    return $argIllDefined;
	}
	if (exists $args{binSizeX}) {
	    warn "hist2D_array: WARNING: binSizeX and nBinsX both defined.  nBinsX overrides!\n";
	}
    }
    if (exists $args{bottomOverflowX} && 
	ref($args{bottomOverflowX}) ne "SCALAR") 
    {
	warn "hist2D_array: bottomOverflowX must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{topOverflowX} && 
	ref($args{topOverflowX}) ne "SCALAR") 
    {
	warn "hist2D_array: topOverflowX must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{minX} && exists $args{maxX}) {
	if ($args{minX}>=$args{maxX}) {
	    warn "hist2D_array: minX must be less than maxX\n";
	    return $argIllDefined;
	} 
	if (exists $args{nBins} && $args{maxX} - $args{minX} >= $args{nBinsX}) {
	    warn "hist2D_array: maxX - minX >= nBinsX\n";
	}
    }

    if (exists $args{binMinY}) {
	if (ref($args{binMinY}) ne "ARRAY") {
	    warn "hist2D_array: binMinY argument is not an array ref\n";
	    return $argIllDefined;
	}
	if (@{$args{binMinY}} > 0) {
	    warn "hist2D_array: binMinY argument is not empty\n";
	    return $argIllDefined;
	}
    }

    if (exists $args{binSizeY} && $args{binSizeY} <= 0) { # defined by default
	warn "hist2D_array: illegal binSizeY<=0: '", $args{binSizeY}, "'\n";
	return $argIllDefined;
    }
    if (exists $args{nBinsY}) {
	if ($args{nBinsY}<2 || $args{nBinsY} != int($args{nBinsY})) {
	    warn "hist2D_array: illegal value of nBinsY: ", $args{nBinsY}, "\n";
	    return $argIllDefined;
	}
	if (exists $args{binSizeY}) {
	    warn "hist2D_array: WARNING: binSizeY and nBinsY both defined.  nBinsY overrides!\n";
	}
    }
    if (exists $args{bottomOverflowY} && 
	ref($args{bottomOverflowY}) ne "SCALAR") 
    {
	warn "hist2D_array: bottomOverflowY must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{topOverflowY} && 
	ref($args{topOverflowY}) ne "SCALAR") 
    {
	warn "hist2D_array: topOverflowY must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{minY} && exists $args{maxY}) {
	if ($args{minY}>=$args{maxY}) {
	    warn "hist2D_array: minY must be less than maxY\n";
	    return $argIllDefined;
	} 
    }
    # done verifying.  phew!

    my $data = $args{data};
    my $hist = $args{hist};

    # initialize optional arguments
    my $filter = $args{filter};
    my $binMinX = $args{binMinX} // [];
    my $bottomOverflowX = 0;
    my $topOverflowX = 0;
    my $minX = $args{minX} // 0;
    my $maxX = $args{maxX} // @$data +0;

    my $binMinY = $args{binMinY} // [];
    my $bottomOverflowY = 0;
    my $topOverflowY = 0;
    my $minY = $args{minY} // min(@$data);
    my $maxY = $args{maxY} // max(@$data);

    if (exists $args{intY}) { 
	$minY-=0.5;
	$maxY+=0.5;
    }

    my $binSizeX = $args{binSizeX} // 1;
    if (defined $args{nBinsX}) { 
	$binSizeX = ($maxX-$minX)/$args{nBinsX};
    }

    my $binSizeY = $args{binSizeY} // 1;
    if (defined $args{nBinsY}) { 
	$binSizeY = ($maxY-$minY)/$args{nBinsY};
    }

    # note loop condition prevents floating point errors
    for (my $thisBin = $minX; $thisBin < $maxX - $binSizeX*0.01;
	 $thisBin+=$binSizeX) {
	push(@$binMinX, $thisBin);
    }

    for (my $thisBin = $minY; $thisBin < $maxY - $binSizeY*0.01;
	 $thisBin+=$binSizeY) {
	push(@$binMinY, $thisBin);
    }

    #use Data::Dumper;
    #print Dumper($binMinX);
    #print Dumper($binMinY);

    die "Stats::hist2D_array failed to create the right number of bins in X.  "
	, "nBins arg = ", $args{nBinsX}, ", size of binMinX = ", @$binMinX+0
	, "\n" if defined $args{nBinsX} && $args{nBinsX} != @$binMinX;

    my ($binX, $binY);
    for my $i (0..$#$data) {
	next if defined $filter && ! $filter->[$i];
	# X overflow supercedes Y overflow.  no double counting.
	if ($i < $minX) {
	    $bottomOverflowX++;	    
	    next;
	} elsif ($i > $maxX) {
	    $topOverflowX++;
	    next;
	}
	my $d = $data->[$i];
	if ($d<$minY) {
	    $bottomOverflowY++;
	    next;
	} elsif ($d>$maxY) {
	    $topOverflowY++;
	    next;
	}
	
	$binX = int( ($i-$minX) / $binSizeX);
	$binX-- if $binX>$#$binMinX; # in case you forgot to use int mode
	
	$binY = int( ($d-$minY) / $binSizeY);
	$binY-- if $binY>$#$binMinY; # in case you forgot to use int mode

	$hist->[$binX][$binY]++;
    }
    
    for my $x (0..$#$binMinX) {
	for my $y (0..$#$binMinY) {
	    $hist->[$x][$y] //= 0;
	}
    }
    
    warn "lingering bug?  binMin and hist are not the same size\n" 
	unless @$binMinX == @$hist;
    
    ${$args{bottomOverflowX}} = $bottomOverflowX;
    ${$args{topOverflowX}} = $topOverflowX;
    ${$args{bottomOverflowY}} = $bottomOverflowY;
    ${$args{topOverflowY}} = $topOverflowY;
    
    return;
}

# make a 2D histogram from a hash of array refs
# hist2D_hash( data => {x => [], 'y' => []}, hist => \@ret
#  (optional: , x=>'x', y=>'y', filter=> \@, intX => 1 // 0, intY => 1 // 0, 
#   minX => $, maxX => $, minY => $, maxY $, binMinX => \@, binMinY => \@, 
#   (binSizeX => $ or nBinsX => $), (binSizeY => $ or nBinsY => $), 
#   bottomOverflowX => \$, topOverflowX => \$, 
#   bottomOverflowY => \$, topOverflowY => \$, 
#  ) 
# )
sub hist2D_hash {
    my %args = @_;

    # first check for unrecognized args
    my @keys = qw( data hist x y filter intX maxX minX binSizeX nBinsX
                   binMinX bottomOverflowX topOverflowX
                   intY maxY minY binSizeY nBinsY binMinY 
                   bottomOverflowY topOverflowY
                 );
    my %keys = map {$_ => 1} @keys;
    for (keys(%args)) {
	die "hist2D_hash: unrecognized arg '$_'\n" unless $keys{$_};
    }
    
    # verify mandatory args:
    if (!exists $args{data}) {
	warn "hist2D_hash: data (input) argument undefined\n";
	return $argNotDefined;
    }
    if (ref($args{data}) ne "HASH") {
	warn "hist2D_hash: data argument is not hash ref\n";
	return $argIllDefined;
    }

    if (!exists $args{hist}) {
	warn "hist2D_hash: hist argument undefined\n";
	return $argNotDefined;
    }
    if (ref($args{hist}) ne "ARRAY") {
	warn "hist2D_hash: hist argument is not an array ref\n";
	return $argIllDefined;
    }
    if (@{$args{hist}} > 0) {
	warn "hist2D_hash: hist argument is not empty\n";
	return $argIllDefined;
    }

    # verify optional args:
    my $x =  $args{x} // 'x';
    my $y =  $args{y} // 'y';
    if (!exists $args{data}{$x}) {
	warn "hist2D_hash: can't find data{$x}\n";
	return $argIllDefined;
    }
    if (ref($args{data}{$x}) ne 'ARRAY') {
	warn "hist2D_hash: data{$x} not an array ref\n";
	return $argIllDefined;
    }
    if (!exists $args{data}{$y}) {
	warn "hist2D_hash: can't find data{$y}\n";
	return $argIllDefined;
    }
    if (ref($args{data}{$y}) ne 'ARRAY') {
	warn "hist2D_hash: data{$y} not an array ref\n";
	return $argIllDefined;
    }
    if (@{ $args{data}{$x} } != @{ $args{data}{$y} } ) {
	warn "hist2D_hash: sizes of data{$x} and data{$y} don't match\n";
	return $argIllDefined;
    }

    if (exists $args{filter}) {
	die "hist2D_hash: filter arg is not implemented\n";
	if (ref($args{filter}) ne 'ARRAY') {
	    warn "hist2D_hash: filter arg is not array ref\n";
	    return $argIllDefined;
	}
	my $fSize = @{$args{filter}} +0;
	my $dSize = @{$args{data}} +0;
	if ($fSize != $dSize) {
	    warn "hist2D_hash: size of filter ( $fSize ) does not match size "
		, "of data $dSize\n";
	    return $argIllDefined;
	}
    }

    if (exists $args{binMinX}) {
	if (ref($args{binMinX}) ne "ARRAY") {
	    warn "hist2D_hash: binMinX argument is not an array ref\n";
	    return $argIllDefined;
	}
	if (@{$args{binMinX}} > 0) {
	    warn "hist2D_hash: binMinX argument is not empty\n";
	    return $argIllDefined;
	}
    }

    if (exists $args{binSizeX} && $args{binSizeX} <= 0) { # defined by default
	warn "hist2D_hash: illegal binSizeX<=0: '", $args{binSizeX}, "'\n";
	return $argIllDefined;
    }
    if (exists $args{nBinsX}) {
	if ($args{nBinsX}<2 || $args{nBinsX} != int($args{nBinsX})) {
	    warn "hist2D_hash: illegal value of nBinsX: ", $args{nBinsX}, "\n";
	    return $argIllDefined;
	}
	if (exists $args{binSizeX}) {
	    warn "hist2D_hash: WARNING: binSizeX and nBinsX both defined.  nBinsX overrides!\n";
	}
    }
    if (exists $args{bottomOverflowX} && 
	ref($args{bottomOverflowX}) ne "SCALAR") 
    {
	warn "hist2D_hash: bottomOverflowX must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{topOverflowX} && 
	ref($args{topOverflowX}) ne "SCALAR") 
    {
	warn "hist2D_hash: topOverflowX must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{minX} && exists $args{maxX}) {
	if ($args{minX}>=$args{maxX}) {
	    warn "hist2D_hash: minX must be less than maxX\n";
	    return $argIllDefined;
	} 
	if (exists $args{nBins} && $args{maxX} - $args{minX} >= $args{nBinsX}) {
	    warn "hist2D_hash: maxX - minX >= nBinsX\n";
	}
    }

    if (exists $args{binMinY}) {
	if (ref($args{binMinY}) ne "ARRAY") {
	    warn "hist2D_hash: binMinY argument is not an array ref\n";
	    return $argIllDefined;
	}
	if (@{$args{binMinY}} > 0) {
	    warn "hist2D_hash: binMinY argument is not empty\n";
	    return $argIllDefined;
	}
    }

    if (exists $args{binSizeY} && $args{binSizeY} <= 0) { # defined by default
	warn "hist2D_hash: illegal binSizeY<=0: '", $args{binSizeY}, "'\n";
	return $argIllDefined;
    }
    if (exists $args{nBinsY}) {
	if ($args{nBinsY}<2 || $args{nBinsY} != int($args{nBinsY})) {
	    warn "hist2D_hash: illegal value of nBinsY: ", $args{nBinsY}, "\n";
	    return $argIllDefined;
	}
	if (exists $args{binSizeY}) {
	    warn "hist2D_hash: WARNING: binSizeY and nBinsY both defined.  nBinsY overrides!\n";
	}
    }
    if (exists $args{bottomOverflowY} && 
	ref($args{bottomOverflowY}) ne "SCALAR") 
    {
	warn "hist2D_hash: bottomOverflowY must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{topOverflowY} && 
	ref($args{topOverflowY}) ne "SCALAR") 
    {
	warn "hist2D_hash: topOverflowY must be a scalar reference\n";
	return $argIllDefined;
    }
    if (exists $args{minY} && exists $args{maxY}) {
	if ($args{minY}>=$args{maxY}) {
	    warn "hist2D_hash: minY must be less than maxY\n";
	    return $argIllDefined;
	} 
    }
    # done verifying.  phew!
    
    my $data = $args{data};
    my $hist = $args{hist};

    # initialize optional arguments
    #my $filter = $args{filter};
    my $binMinX = $args{binMinX} // [];
    my $bottomOverflowX = 0;
    my $topOverflowX = 0;
    my $minX = $args{minX} // min( @{ $data->{$x} } );
    my $maxX = $args{maxX} // max( @{ $data->{$x} } );

    my $binMinY = $args{binMinY} // [];
    my $bottomOverflowY = 0;
    my $topOverflowY = 0;
    my $minY = $args{minY} // min( @{ $data->{$y} } );
    my $maxY = $args{maxY} // max( @{ $data->{$y} } );

    if (exists $args{intX}) { 
	$minX-=0.5;
	$maxX+=0.5;
    }
    if (exists $args{intY}) { 
	$minY-=0.5;
	$maxY+=0.5;
    }

    my $binSizeX = $args{binSizeX} // 1;
    if (defined $args{nBinsX}) { 
	$binSizeX = ($maxX-$minX)/$args{nBinsX};
    }

    my $binSizeY = $args{binSizeY} // 1;
    if (defined $args{nBinsY}) { 
	$binSizeY = ($maxY-$minY)/$args{nBinsY};
    }

    # note loop condition prevents floating point errors
    for (my $thisBin = $minX; $thisBin < $maxX - $binSizeX*0.01;
	 $thisBin+=$binSizeX) {
	push(@$binMinX, $thisBin);
    }
    for (my $thisBin = $minY; $thisBin < $maxY - $binSizeY*0.01;
	 $thisBin+=$binSizeY) {
	push(@$binMinY, $thisBin);
    }

    #use Data::Dumper;
    #print Dumper($binMinX);
    #print Dumper($binMinY);

    die "Stats::hist2D_hash failed to create the right number of bins in X.  "
	, "nBins arg = ", $args{nBinsX}, ", size of binMinX = ", @$binMinX+0
	, "\n" if defined $args{nBinsX} && $args{nBinsX} != @$binMinX;

    my ($binX, $binY);
    for my $i (0..$#{ $data->{$x} }) {
	#next if defined $filter && ! $filter->[$i];
	# X overflow supercedes Y overflow.  no double counting.
	my $xVal = $data->{$x}[$i];
	my $yVal = $data->{$y}[$i];
	if ($xVal < $minX) {
	    $bottomOverflowX++;	    
	    next;
	} elsif ($xVal > $maxX) {
	    $topOverflowX++;
	    next;
	}
	if ($yVal<$minY) {
	    $bottomOverflowY++;
	    next;
	} elsif ($yVal>$maxY) {
	    $topOverflowY++;
	    next;
	}
	
	$binX = int( ($xVal-$minX) / $binSizeX);
	$binX-- if $binX>$#$binMinX; # in case you forgot to use int mode
	
	$binY = int( ($yVal-$minY) / $binSizeY);
	$binY-- if $binY>$#$binMinY; # in case you forgot to use int mode

	$hist->[$binX][$binY]++;
    }
    
    for my $x (0..$#$binMinX) {
	for my $y (0..$#$binMinY) {
	    $hist->[$x][$y] //= 0;
	}
    }
    
    warn "lingering bug?  binMin and hist are not the same size\n" 
	unless @$binMinX == @$hist;
    
    ${$args{bottomOverflowX}} = $bottomOverflowX;
    ${$args{topOverflowX}} = $topOverflowX;
    ${$args{bottomOverflowY}} = $bottomOverflowY;
    ${$args{topOverflowY}} = $topOverflowY;
    
    return;
}

"no text below this line"

package HomeBrew::StructManip;

use strict;
use warnings;

require Exporter;

# sortArrayOfHashes(\@array, $col, (optional: $sortFunction))
# sortHashOfArrays(\%hash, $col, (optional: $sortFunction))
#  must convert %hash{cols} to @array[{cols}] then call sAoH
# convertHash2Array($hash, (optional $ret))
# convertArray2Hash($hash, (optional $ret))

our @ISA = qw(Exporter);

our @EXPORT_OK = qw(sortArrayOfHashes sortHashOfArrays convertHash2Array
    convertArray2Hash);

our $VERSION = '0.01';

sub sortArrayOfHashes {
    my ($array, $col, $sortFunction) = @_;
    my $fRef = \$_[-1];
    $sortFunction = $$fRef;
    die "StructManip::sortArrayOfHashes first array must be nonempty array of".
	"hashes\n" unless ref($array) eq 'ARRAY' && @$array>0 
	&& ref($array->[0]) eq 'HASH';
    die "StructManip::sortArrayOfHashes can't find $col in array\n"
	unless exists $array->[0]{$col};

    my %sortable; # sortable{colVal} = [ all elems of $array with colVal ]
    for my $h (@$array) {
	$sortable{$h->{$col}} //= [];
	push @{ $sortable{ $h->{$col} } }, $h;
    }

    @$array = (); # clear array so that you can refill it in order
    my @sortedK;
    if (defined $sortFunction) {
	@sortedK = sort $sortFunction keys %sortable;
    } else {
	@sortedK = sort keys %sortable;
    }

    for my $k (@sortedK) {
	for my $h (@{$sortable{$k}}) {
	    push @$array, $h;
	}
    }
    
    return;
}

sub sortHashOfArrays {
    my ($hash, $col, $sortFunction) = @_;
    die "StructManip::sortHashOfArrays first arg must be hash ref\n" 
	unless ref($hash) eq 'HASH';
    die "StructManip::sortHashOfArrays can't find hash{$col} or it's not an ".
	"aref\n" unless exists $hash->{$col} && ref($hash->{$col}) eq 'ARRAY';
    my $size = @{$hash->{$col}};
    for (my($k, $v) = each %$hash) {
	die "StructManip::sortHashOfArrays hash{$k} is not an aref\n" 
	    unless ref($v) eq 'ARRAY';

	die "StructManip::sortHashOfArrays hash{$k} has wrong size.  "
	    , "got ", @$v+0, ", expected $size.\n" 
	    unless @$v == $size;
    }

    my $ar = [];
    convertHash2Array($hash, $ar);
    sortArrayOfHashes($ar, $col, $sortFunction);
    convertArray2Hash($ar, $hash);
}

sub convertHash2Array {
    my ($hash, $ar) = @_; # second arg is optional
    die "StructManip::convertHash2Array first arg must be hash ref\n" 
	unless ref($hash) eq 'HASH';
    die "StructManip::convertHash2Array optional second arg must be array ref\n"
	if defined $ar && ref($ar) ne 'ARRAY';

    my $ret = [];
    $ret = $ar if defined $ar; # second arg is optional
    @$ret = ();

    my $size;
    {
	my @keys = keys %$hash;
	$size = @{$hash->{$keys[0]}};
    }
    push(@$ret, {}) for 1..$size;
    while (my ($k, $v) = each %$hash) {
	die "StructManip::covertHash2Array hash{$k} is not an aref\n" 
	    unless ref($v) eq 'ARRAY';

	die "StructManip::covertHash2Array hash{$k} has wrong size.  "
	    , "got ", @$v+0, ", expected $size.\n" 
	    unless @$v == $size;

	$ret->[$_]{$k} = $v->[$_] for 0..$#$v;
    }

    return @$ret unless defined $ar;
    return;
}

sub convertArray2Hash {
    my ($ar, $hash) = @_; # second arg is optional
    die "StructManip::convertArray2Hash optional second arg must be array ref\n"
	unless ref($ar) eq 'ARRAY';
    die "StructManip::convertArray2Hash first arg must be hash ref\n" 
	if defined $hash && ref($hash) ne 'HASH';

    my $ret = {};
    $ret = $hash if defined $hash; # second arg is optional
    %$ret = ();

    my @cols = keys %{$ar->[0]};

    $ret->{$_} = [] for @cols;
    for my $h (@$ar) {
	die "StructManip::covertArray2Hash arr element '$h' is not an href\n" 
	    unless ref($h) eq 'HASH';

	for (@cols) {
	    die "StructManip::covertArray2Hash cannot find key $_ in array ".
		"element 'h'.\n" unless exists $h->{$_};

	    push @{$ret->{$_}}, $h->{$_};
	}
    }

    return %$ret unless defined $hash;
    return;
}


"no text below this line"

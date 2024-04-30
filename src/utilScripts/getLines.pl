#!/usr/bin/env perl

unless(@ARGV==3 || @ARGV==4) {    
    die "usage:\t<script> <in> <out> <start> <stop> \n\tOR\n\t<script> <in> <out> <nlines>\n";
}

$in = $ARGV[0];
$out = $ARGV[1];
$start = $ARGV[2];
$stop = $ARGV[3];
if (@ARGV==3) {    
    $stop = $start;
    $start = 1;
}

open(IN,$in) or die "can't open infile '$in'.  $!\n";
open(OUT,">$out") or die "can't open outfile '$out'.  $!\n";

$count = 1;

while ($line = <IN>) {
    print OUT $line if ($count>=$start && $count<=$stop);
    last if ($count>=$stop);
    $count++;
}
close IN;
close OUT;

die "file ended before line $stop\n" if ($count<$stop);

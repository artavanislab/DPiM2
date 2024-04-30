#!/usr/bin/env perl

$out = "tmp.sh";
unless (@ARGV == 2 || @ARGV == 3) {
    print "usage: <script> first last <out = $out>\n";
    exit;
}

$first = $ARGV[0];
$last = $ARGV[1];
$out = $ARGV[2] if defined $ARGV[2];

@arr = ($first..$last);
open(OUT,">$out");
print OUT "#!/bin/bash\n";
for $i (@arr) {
    print OUT "qdel $i\n";
}
close OUT;
system("chmod 755 $out");
exit;

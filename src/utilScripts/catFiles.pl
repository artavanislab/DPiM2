#!/usr/bin/env perl

use strict;
use Data::Dumper;

# take a list of files and concatenate them.


# getCommandLineOptions
# checkExist(option,path)
# read_list($file) read the first string on each line of file, return array
# readOcc(file)
# writeOcc(\%occ,$file);

my %opts = &getCommandLineOptions();
my $fileList = $opts{l}; # input filelist to concatenate
my $out = $opts{out};

if (defined $opts{dum}) {
    dum();
} else {
    main();
}
exit;

sub main {

    my @files = ();
    @files = read_list($fileList);
    if (@files<2) {
	die "need at least two files!\n";
    }
    foreach my $f (@files) {
	checkExist('f',$f);
    }

    my $header;
    my $file = $files[0];
    my %occ = readTab($file,\$header);

    # check headers first
    my @arr = @files;
    shift(@arr);
    while ($file = shift(@arr)) {
	checkHeader($file,$header);
    }

    
    @arr = @files;
    shift(@arr);
    my $prevMax = scalar(keys(%occ));
    my $i;
    while ($file = shift(@arr)) {
	my %occ2 = readTab($file,\$header);
	#print Dumper(\%occ2);
	for $i (1..keys(%occ2)) {
	    $occ{$i+$prevMax} = $occ2{$i};
	}
	$prevMax = scalar(keys(%occ));
    }

    open(OFILE,">$out");
    print OFILE "# concatenating files: ( @files )\n";
    print OFILE "$header\n";
    for $i (1..keys(%occ)) {
	printf OFILE "%8d   ".$occ{$i}."\n",$i;
    }
    close OFILE;

    return;
}


# as cat > $out execpt you omit headers and comments
sub dum {
    my @files = ();
    @files = read_list($fileList);
    if (@files<2) {
	die "need at least two files!\n";
    }
    foreach my $f (@files) {
	checkExist('f',$f);
    }

    open(OUT,">$out");
    my ($line,$header);
    foreach my $f (@files) {
	open(IN,$f) or die "couldn't read from $f\n";
	print "reading $f...\n";
	# find header
	do {
	    $line = <IN>;
	} while (substr($line,0,1) eq '#');

	if (!defined $header) {
	    $header = $line;
	    print OUT $header;
	}

	while ($line = <IN>) {
	    print OUT $line; # assume no more comments in file
	}
	close IN;
    }
    close OUT;
    
    return;
}

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   


sub getCommandLineOptions {
  use Getopt::Long;
  my $usage = qq{usage: $0 -l occ.list -out concatenated.occ < -dum >\n};
  # Get args
  &GetOptions(\%opts, "l=s", "out=s","dum");
  if (!defined $opts{l} || !defined $opts{out}) {
    print STDERR "$usage\n";
    exit -1;
  }
  
  &checkExist('f', $opts{l});

  return %opts;
}

# checkExist()
#
sub checkExist {
  my ($type, $path) = @_;
  if ($type eq 'd') {
    if (! -d $path) {
      warn("dir not found: $path\n");
      exit -3;
    }
  }
  elsif ($type eq 'f') {
    if (! -f $path) {
      warn("file not found: $path\n");
      exit -3;
    }
    elsif (! -s $path) {
      warn("empty file: $path\n");
      exit -3;
    }
  }
}


# This function reads in a list of items from a file (one item per line, rest of line ignored;
# '#' comment symbols recognized).
sub read_list {
  my ($file) = @_;
  my @tokens = ();
  my $cnt = 0;
  open(I,"$file") or die "Can't open file: $!\n";
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


# read occupancy file (name of file is arg)
# returns an hash of hash $occ{seqN}{'Prob1'} and $occ{seqN}{'Occ1'}
sub readOcc {
    my ($occFile) = @_;
    open (OFILE,$occFile) or die "can't open occupancy file $!\n";
    
    my $line1;

    #remove comments
    do {
	$line1 = <OFILE>;
    } while (substr($line1,0,1) eq '#');
    chomp ($line1);

    my @header = split(/ +/,$line1);
    if ($header[0] eq '') {
	shift(@header);
    }
    #print Dumper(\@header);

    
    my @line = ();
    my %occ = ();
    my $seqN=0;
    my $j=1;
    my $ent;
    my $i=0;
    while (<OFILE>) {
	chomp;
	if (length($_)==0) {
	    next;
	}
	$i++;
	@line = split(/ +/,$_);
	if (@line==0) {
	    next;
	}

	$seqN = shift(@line);
	if (!length($seqN)) { 
	    $seqN = shift(@line);
	}
	foreach $ent (@line) {
	    #print_r($ent);
	    $occ{$seqN}{$header[$j]} = $ent;
	    $j++;
	}
	$j=1;
	
    }
    
    close OFILE;

    return %occ;
}

sub writeOcc {
    my ($occA,$outFile) = @_;
    my %occ = %$occA;

    open (OUTF,">$outFile") or die "can't open output file $outFile\n";

    my $header = "    seqN      Prob1       Occ1\n";
    
    print OUTF $header;
    
    my $line;
    for my $i (1..keys(%occ)) {	
	$line = sprintf("%8d   %12.8e   %12.8e\n",($i),$occ{$i}{Prob1},$occ{$i}{Occ1});
	print OUTF $line;
    }
    
    close OUTF;
    return;
}

# read any file formatted like this
#  "# do not parse coments"
#  "X junk"
#  "Y junk2"
# return a hash tab{X} = junk, tab{Y} = junk2
sub readTab {
    my ($inFile,$header) = @_;
    open (OFILE,$inFile) or die "can't open occupancy file $!\n";
    
    my $line1;

    #remove comments
    do {
	$line1 = <OFILE>;
    } while (substr($line1,0,1) eq '#');

    chomp ($line1);    
    my $myHeader = $line1;
    if (defined $$header) {
	die "should have already caught this error\n" if $$header ne $myHeader;
    } else {
	$$header = $myHeader;
    }
    
    my @line = ();
    my %tab;
    my @spl;
    my $seqN=0;                                        # '      -4 23 45 23 -4 3 2'...
    my $len1; # number of blank chars preceding the index      ^ 6 
    my $len2; # end position of index string                      ^ 8
    my $i=0;
    my $first=1;
    my $E;
    my $lump;
    while ($line1 = <OFILE>) {
	chomp $line1;
	if (length($line1)==0) {
	    next;
	}
	$len1=0;
	while (substr($line1,$len1,1) eq ' ') {
	    $len1++;
	}
	$len2 = $len1;
	do {
	    $len2++;
	} while (substr($line1,$len2,1) ne ' ');	    
	#$len++;
	$i =  substr($line1,$len1,($len2-$len1));
	if ($i =~ /\D/)  {
	    die "found a non-digit in key: '$i'\nQuitting\n";
	}
	$lump = substr($line1,$len2);
	$tab{$i} = $lump;

    }

    close OFILE;
    
    my @arr = sort numerically (keys(%tab));
    my $prev = shift(@arr);
    while ($i = shift @arr) {
	if ($i-$prev!=1) {
	    die "nonconsecutive keys: $prev , $i \n";
	}
	$prev = $i;
    }

    return %tab;
}

# find first non-comment line in file, 
#  compare this line with $header, 
#   die if different
#    Tuez la difference!
sub checkHeader {
    my ($inFile,$header) = @_;
    open (OFILE,$inFile) or die "can't open occupancy file $!\n";
    
    my $line1;

    #remove comments
    do {
	$line1 = <OFILE>;
    } while (substr($line1,0,1) eq '#');

    chomp ($line1);    
    my $myHeader = $line1;
    die "old header '$header' doesn't match new header '$myHeader' in file $inFile\n" if $header ne $myHeader;
    close OFILE;
    return;
}
sub numerically { $a <=> $b; }

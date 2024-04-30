#!/usr/bin/env perl

use strict;
use Data::Dumper;

# find subdirectories
# find .XX files in them
# produce a list of all .XX files with full path.  

# use -local for local path
# use -rec to search subdirectories

my $DEF_EXT = "nz";

my %opts = &getCommandLineOptions();
my $out  = $opts{out};
my $home   = $opts{dir};
my $ext  = $opts{ext};
my $grep = $opts{grep};
my $except = $opts{except};
$out = $ext.".list" if defined $ext && ! defined $out;
$out = $grep.".list" if ! defined $ext && ! defined $out;

main();
exit;

sub main {

    if (! defined $home) {
	$home = `pwd`;
	chomp $home;
    }
    
    my @dirs = ();
    unless (defined $opts{rec}) {
	@dirs = ( $home );
    } else {
	@dirs = split(/\n/, `find $home -type d`);   
    }

    my @matches;
    for my $dir (@dirs) {
	next if defined $except && $dir =~ /$except/;
	opendir(my $DIR, $dir);
	my @files;
	if (defined $ext) {
	    @files = grep(/\.$ext$/, readdir($DIR));
	} else {
	    @files = grep(/$grep/, readdir($DIR));
	}
	closedir($DIR);
	@files = grep(/$grep/, @files) 
	    if defined $grep && defined $ext && $dir !~ /$grep/;
	@files = grep(!/$except/, @files) if defined $except;
	@files = grep(!/$out/, @files);
	foreach my $file (@files) {
	    my $path = (defined $opts{local})?$file:$dir.'/'.$file;
	    push(@matches,$path);
	}
    }
    print "found ".@matches." files\n";
    return if @matches == 0;

    my $OUT;
    if ($out eq 'stdout' || $out eq 'STDOUT') {
	$OUT = *STDOUT;
    } else {
	open ($OUT, ">", $out) or die "couldn't open outfile $out ( $! )";
    }
    print $OUT "# finding files in $home";    
    print $OUT " ending in .$ext files starting at $home" if defined $ext;
    print $OUT " requiring '$grep' " if defined $grep;
    print $OUT " and excluding '$except'" if defined $except;
    print $OUT "\n";
    print $OUT join("\n",sort @matches), "\n";
    close $OUT;
    return;
}

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {
  use Getopt::Long;
  my $usage = qq{usage: $0 (-ext $DEF_EXT and/or -grep "regex") <-dir (default .) -out $DEF_EXT.list -rec -local -except "regex">};
  # Get args
  my %opts = ();
  &GetOptions(\%opts, "out=s", "dir=s", "ext=s", "rec", "local", "except=s",
	      "grep=s");
  die $usage, "\n" unless (exists $opts{ext} || exists $opts{grep});

  return %opts;
}

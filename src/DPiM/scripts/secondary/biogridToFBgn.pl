#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use HomeBrew::IO qw(checkExist);
#use DpimLib qw(getLineAPMS);

# convert a BioGRID tab2 file into FBgns

my %opts = getCommandLineOptions();

{
    my $in = $opts{in};
    my $out = $opts{out};
    
    my @keep = (
	'BioGRID Interaction ID',
	'Systematic Name Interactor A',
	'Systematic Name Interactor B',
	'Official Symbol Interactor A',
	'Official Symbol Interactor B',
	'Experimental System',
	'Experimental System Type',
	'Author',
	'Pubmed ID',
	'Score',
	'Modification',
	'Phenotypes',
	'Qualifications',
	'Tags',
	)
	
}

exit;

   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

sub getCommandLineOptions {

    my %defaults = (
	);
    my $defaultString = 
	join " ", map { "-$_ $defaults{$_}" } sort keys %defaults;

    my $usage = "usage: $0 -in input -out output\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# sets the $ret hash keyed as specified
sub getLineBioGRID {
    my ($ret, $IN, $wholeLine) = @_;

    return undef if eof($IN);

    my @cols = (
	'BioGRID Interaction ID',
	'Entrez Gene Interactor A',
	'Entrez Gene Interactor B',
	'BioGRID ID Interactor A',
	'BioGRID ID Interactor B',
	'Systematic Name Interactor A',
	'Systematic Name Interactor B',
	'Official Symbol Interactor A',
	'Official Symbol Interactor B',
	'Synonyms Interactor A',
	'Synonyms Interactor B',
	'Experimental System',
	'Experimental System Type',
	'Author',
	'Pubmed ID',
	'Organism Interactor A',
	'Organism Interactor B',
	'Throughput',
	'Score',
	'Modification',
	'Phenotypes',
	'Qualifications',
	'Tags',
	'Source Database'
        );

                  prey_ref);
    my @spl;
    my $line;
    do {{ # perl's next statement doesn't work inside a do block.
	$line = <$IN>;
	next if $line =~ /^#/;
	next if $line =~ /tap_id/;
	next if length($line) < 4;
	$_ = $line;
	chomp;
	@spl = split;
	}
    } while (!eof($IN) && @spl != @cols);
    return undef if @spl != @cols; # checks for eof

    $ret->{$cols[$_]} = $spl[$_] for 0..$#cols ;
    $ret->{$wholeLine} = $line if defined $wholeLine;
    return 1;
}

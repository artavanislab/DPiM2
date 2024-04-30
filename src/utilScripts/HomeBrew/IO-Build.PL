#!/usr/bin/env perl

use warnings;
use strict;
use Module::Build;
use lib '/home/glocke/perl5/lib/perl5';

my $build = Module::Build->new
    (
     module_name => "HomeBrew::IO",
     dist_author => "George Locke",
     dist_abstract => "Various utilities for reading files",
     build_requires => {
	 #'Bio::SeqIO' => '0', # to read fasta sequences
	 'Test::More' => '0.10',
	 'Test::Exception' => '0', # to test that checkExist dies
	 'List::Util' => '0', # for min in printFasta
	 'Text::ParseWords' => '0', # for quotewords
     },
    );
$build->create_build_script();

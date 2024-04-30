#!/usr/bin/env perl

use warnings;
use strict;
use Module::Build;

my $build = Module::Build->new
    (
     module_name => "HomeBrew::Bio",
     dist_author => "George Locke",
     dist_abstract => "bioinformatical utilities",
     build_requires => {
	 'Test::More' => '0.10',
	 'Test::Exception' => '0', 
	 'List::Util' => '0',
     },
    );
$build->create_build_script();

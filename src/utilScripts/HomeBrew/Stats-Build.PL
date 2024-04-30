#!/usr/bin/env perl

use warnings;
use strict;
use Module::Build;
my $build = Module::Build->new
    (
     module_name => "HomeBrew::Stats",
     dist_author => "George Locke",
     dist_abstract => "a few statistics type functions",
     build_requires => {
        'Test::More' => '0.10',
	'List::Util' => '0', # for sum
     },
    );
$build->create_build_script();

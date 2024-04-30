#!/usr/bin/env perl

use warnings;
use strict;
use Module::Build;

my $build = Module::Build->new
    (
     module_name => "HomeBrew::StructManip",
     dist_author => "George Locke",
     dist_abstract => "",
     build_requires => {
	 'Test::More' => '0.10',
     },
    );
$build->create_build_script();

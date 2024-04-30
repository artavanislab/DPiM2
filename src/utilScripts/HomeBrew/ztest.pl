#!/usr/bin/env perl

use warnings;
use strict;
use Test::More qw(no_plan);

# Verify module can be included via "use" pragma
BEGIN { use_ok('HomeBrew::IO') };

# Verify module can be included via "require" pragma
require_ok('HomeBrew::IO');

my $retVal;
print "checkExist('f', $0)\n";
$retVal = HomeBrew::IO::checkExist('f', $0);
ok(!$retVal, "checkExist('f') OK test");

print "\n\n";

$retVal = HomeBrew::IO::checkExist('d', $ENV{PWD});
ok(!$retVal, "checkExist('d') OK test");


print "retVal == undef? ", ($retVal eq undef)?"yes":"no", "\n";
print "OK!\n";




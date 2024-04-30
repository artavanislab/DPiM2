#!/bin/bash 
module add perl
module add gcc
module add gsl
module add boost
module add CPAN
module add CPAN/BioPerl
module add CPAN/DateTime
module add CPAN/Math-GSL-CDF

DPSCR="/home/glocke/DPiM/scripts/"; export DPSCR;
DPCPP="/home/glocke/DPiM/cpp/"; export DPCPP;
UTILSCR="/home/glocke/DPiM/scripts/"; export UTILSCR;
PERL5LIB="/home/glocke/DPiM/scripts/lib:/home/glocke/utilScripts/HomeBrew/lib:/home/glocke/perl5/lib/perl5${PERL5LIB+:}${PERL5LIB}"; export PERL5LIB;
PATH="/home/glocke/perl5/bin${PATH+:}${PATH}"; export PATH;
PERL_LOCAL_LIB_ROOT="/home/glocke/perl5${PERL_LOCAL_LIB_ROOT+:}${PERL_LOCAL_LIB_ROOT}"; export PERL_LOCAL_LIB_ROOT;
PERL_MB_OPT="--install_base \"/home/glocke/perl5\""; export PERL_MB_OPT;
PERL_MM_OPT="INSTALL_BASE=/home/glocke/perl5"; export PERL_MM_OPT;

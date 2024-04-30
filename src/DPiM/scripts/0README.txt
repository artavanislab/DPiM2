##                                      ##
 ## Biogen Interactome Software Pipeline ##
  ##                                      ##
   ##     Notes on installation and use    ##
    ##                                      ##

This document written by George Locke August 2016

Software written by George Locke, Julian Mintseris, and Kejie Li

##
 ## Installation
  ##

The software is located in https://bitbucket.org/biogenidec/compbiobatk/branch/interactome-feature

## Environment

I have created a script called 0setup.sh that implements the following changes
to the environment.

You will need to define/edit the following environmental variables:
DPSCR - points to the directory containing this readme, preFilterMaster.pl, etc
UTILSCR - points to the directory containing qsubWrap.pl
PERL5LIB - add pointers to HomeBrew and DpimLib

Here are the lines I've added to my .bashrc:
DPSCR="/home/glocke/DPiM/scripts/"; export DPSCR;
DPCPP="/home/glocke/DPiM/cpp/"; export DPCPP;
UTILSCR="/home/glocke/DPiM/scripts/"; export UTILSCR;
PERL5LIB="/home/glocke/DPiM/scripts/lib:/home/glocke/utilScripts/HomeBrew/lib:/home/glocke/perl5/lib/perl5${PERL5LIB+:}${PERL5LIB}"; export PERL5LIB;

If for some reason you don't want to mess with your environmental variables,
you could add "use lib '/path/to/lib';" lines to the perl, and hard-code 
locations where necessary.

If you're working on Biogen's camhpc cluster, you will need the following
modules added:
module add perl
module add gcc
module add gsl
module add boost
module add CPAN
module add CPAN/BioPerl
module add CPAN/DateTime
module add CPAN/Math-GSL-CDF

Another note for biogen, is that the adding the boost module did not enable g++
was to find the boost dll's, so there is a line in the makefile giving a full
path to where the .so files are.  Whenever the sysadmins update boost, this
address will need to be updated as well.

If you're not working there, you may want to change the way the scripts interact
with qsub.  Relevant scripts include utilScripts/qsubWrap.pl, runHyper.pl,
and preFilterMaster.pl.

I'll note here that the readFasta routine in HomeBrew::IO uses Bio::SeqIO.
Installing BioPerl can be a pain, and you may find that you wish to alter that
routine so that it doesn't require BioPerl.  While there are several routines
in that module that aren't called, readFasta is indeed called to establish
whether there are proteins with unknown length (length correction is part of
hyperspec).


## External libraries

The perl code includes a variety of stuff from CPAN.  

The C++ code requires boost/program_options, gsl/gsl_histogram, gsl/gsl_cdf,
and rapidjson.  (The json output is sort of a verbose log, so, if you were so
inclined, you could simply delete all the code that has anything json related
and be left with functioning code.)

We also use MCL for clustering.  http://micans.org/mcl/

## External data

* for fly:
GLL note 2016-12-19 these links are incorrect
use r6.07 data

//ftp://ftp.flybase.net/releases/FB2016_04/precomputed_files/genes/fbgn_annotation_ID_fb_2016_04.tsv.gz
//ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.12_FB2016_04/fasta/dmel-all-translation-r6.12.fasta.gz

* for human:
all gpff files from
ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/

##
 ## Use
  ##

## general note on software
executing any script without any arguments will produce a usage message
informing you of the names of the arguments that looks something like this:
"script.pl -in required -out required < -mode *default*/alt1/alt2 -opt defaultVal -optflag >"
the arguments inside the brackets are optional; the script assumes default
values if the user doesn't specify the argument.

The C++ executable, which you may never access directly, has a usage message
that is more verbose and hopefully self-explanatory.

## input format

The perl expects apms data in tsv format with the following columns: 
search_id	 bait_ref  prey_ref	 total_peptides	    sample_date	ms_inst_run_id

This is "dp4 format".  There may be comments in the file (first character in
line being '#'), and the first non-commented line must name the columns.

The first column is a unique numeric id for each experiment; aside from
uniqueness, its value has no importance.  bait_ref and prey_ref are id's
(FBgn expected for fly, entrez id for human). total_peptides is the abundance of
this prey.  sample_date is the date when this experiment took place.
ms_inst_run_id is a run identifier - the LC carryover procedure requires that
these id's accurately represent the order in which experiments took place (the
code interprets f1234 as preceding f1235).

## before you can run the pipeline

Heaps of work went in to curating the data properly.  We boiled it all down to
a single flat file presently located at /home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_112116.txt

You need to create a file that matches each protein id to its length (number of
amino acids in protein). for flies, this is performed by
secondary/fastaToProteinLength.pl; for human, use
human/secondary/parseRefseqGFF.pl

You also need to create a file that maps from protein id to gene symbol.  For
FBgn, this is from secondary/FBgn2SymbTab.pl .   For human I recommend the
org.Hs.eg.db R library; I have not yet implemented a script for it.

## main pipeline

Here is the order of operations:

0. reprocessed/reorderMultipleIDs.pl
1. reprocessed/preFilterMaster.pl
2. runHyper.pl
3. compute_fdr_GL.pl

First, you prepare the data to be scored, then you score it, then you interpret
the scores and finalize the network.

# 0. reprocessed/reorderMultipleIDs.pl

Some run id's refer to multiple experiments.  This script invents new run id's
for those experiments, thus disambiguating them.  Key fact for it is 
disambiguates experiments with run id's that refer 

# 1. preFilterMaster.pl

This script prepares the data to be scored.  Its arguments are:
in - the apms data formatted as specified above
out - names the output files (see below)
id2symb - file mapping from protein id to gene symbol
lenfile - file mapping from protein id to translated length
scrdir - the directory where all scripts are found
fbgnstandard - points to fbgn_annotation_ID_fb_2016_04.tsv
mode - fly or human

Note on output files: all output files will be have names starting with the
value specified by the -out argument, and this process produces many files.  The
main thing to know is that **the final output of the process is named
$out.toscore.nrBait** (where $out is what you specified).

# 2. runHyper.pl

This script scores every interaction and produces a set of null (simulated)
networks for comparison.  Its arguments are

in - the filtered input (note only four columns instead of 6)
mode - fly or human
dir - output goes here
job - the name attached to the output files and the jobs in the queue
nsim - the number of null networks to produce
rt - the maximum run time for each job
seed - random seed (automatically initialized to system time)
simonly - do not compute the true network, just the simulations (use if you want more simulated networks, i.e. probably never)

the output goes to subdirectories within the specified -dir (the script will
create the directories if they do not already exist).  json files go into one
directory, "json", and files related to the queue go into another, "qdir".
You will probably only be interested in qdir.  

NB: You will probably never touch the json files.

# 3. compute_fdr_GL.pl

create a list pointing the the files with all simulation output, e.g.
ls $PWD/*sim*.o* > sim.list


This interprets the true scores in light of the null scores and outputs a
network of "high confidence interactions".  Its arguments are:

real - the output log of the job producing the real scores
sim - file list of output logs for the null (simulated network) jobs
out - final network goes here
onesim - set this flag if your -sim argument is one file rather than a list
       of files
mode - fly or human

before executing this script, go into the qdir from the previous step and make
the list of simulation files, e.g. "ls $PWD/*sim*.o* > sim.list"

Now you've got your very own network.  Bob's your uncle!



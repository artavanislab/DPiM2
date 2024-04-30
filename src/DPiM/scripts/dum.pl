#!/usr/bin/env perl

##############################################################################80

use v5.10; 
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(min max sum);
use List::MoreUtils qw(each_array);
use File::Path qw(make_path);
#use File::Temp qw(tempfile);
use Text::ParseWords;
use Storable qw(retrieve);
use DateTime;
use DateTime::Format::Strptime;
use Math::GSL::CDF qw(gsl_cdf_hypergeometric_P);
use Statistics::R;

use HomeBrew::IO qw(checkExist readCol readColsRef readHeader readCols writeCols
                    readList readList2 readColsHashRef);
use DpimLib qw(getLineAPMS readHS getLineDP4APMS getLineDP4_1APMS getLineRawAPMS
  getLineHyperspecAPMS networkHashFromEdgeList getLineMayAPMS);

{
    my $in = "/home/glocke/allFiles.txt";
    open my $IN, "<", $in or die "can't read from $in. $!";
    while (<$IN>) {
	chomp;
	next if /^\.+$/;
	my $cmd = "chmod +r $_";
	#say $cmd;
	#last;
	system($cmd);
    }

    exit;
}

{
    my $in = '/home/glocke/DPiM/augRemap/apmsData/annotation/RAO_Replicate_Winners_121916.txt';
    my @cols = qw(ms_inst_run_id	search_id	bait_ref	sample_date	numProteins	baitPeptides	nrBaitWinner	Bobs_Choice	Rebranded_from	Rebranded_to);

    my $badyearstart = "2009-04-01";
    my $badyearend = "2010-06-01";

    my %choose;
    {
	my @read;
	readColsRef(\@read, $in, \@cols, 'line', "\t");
	for my $row (@read) {
	    $row->{sample_date} =~ m:(\d+)/(\d+)/(\d\d\d\d):
		or die "can't parse date ".Dumper($row);
	    $row->{sample_date} = sprintf("%04d-%02d-%02d", $3, $1, $2);
	    $choose{$row->{bait_ref}}{$row->{search_id}} = $row;
	}
    }
    #die Dumper($choose{FBgn0035753});

    # return true if this date is inside the badyear
    my $inBadYear = sub {
	my ($checkDate) = @_;

	return 1 if $checkDate eq $badyearstart;
	#return 1 if $checkDate eq $badyearend;

	return ($badyearstart lt $checkDate && $checkDate lt $badyearend);
    };
    
    
    my @allCols = qw(ms_inst_run_id	search_id	bait_ref	sample_date	numProteins	baitPeptides	nrBaitWinner	Bobs_Choice	Rebranded_from	Rebranded_to	Rebranded_Annotation	Rebranded_Name	Rebranded_Symbol	Final_Updated_Bait_ID	tap_id	tap_bait	rejected	wonky	geneSymbol1	search_id	sample_date	TAP_ID	CG_ID	CG-X_ID	Flybase_ID	Unique_Peptides	Total_Peptides	File_Name	Data_Receipt_Date	Comments	Problems	Contamination_Status	num_contaminating_peptides);
    say "# $0 selected only those corrections not due to rebranding or badyear from $in";
    ##say join "\t", @allCols;
    my @outCols = qw(search_id bait_ref baitPeptides);
    say join "\t", "pickMe", @outCols;
    
    for my $bait (sort keys %choose) {
	my @sids = keys %{ $choose{$bait} };
	die "should be two keys\n" if 2 != @sids;
	@sids = sort {$choose{$bait}{$a}{Bobs_Choice} cmp 
			  $choose{$bait}{$b}{Bobs_Choice}} @sids;

	# ignore failure to account for rebranded baits
	next if length($choose{$bait}{$sids[1]}{Rebranded_to}) && 
	    ! length($choose{$bait}{$sids[0]}{Rebranded_to});

	# ignore the bad year if that's the only problem
	next if $inBadYear->($choose{$bait}{$sids[0]}{sample_date}) && 
	    ! $inBadYear->($choose{$bait}{$sids[1]}{sample_date}) && 	    
	    $choose{$bait}{$sids[1]}{baitPeptides} < 
	    $choose{$bait}{$sids[0]}{baitPeptides};

	
	say join "\t", 1, map {$choose{$bait}{$sids[0]}{$_}} @outCols;
	say join "\t", 0, map {$choose{$bait}{$sids[1]}{$_}} @outCols;
	##print $choose{$bait}{$sids[0]}{line};
	##print $choose{$bait}{$sids[1]}{line};
	
    }
    
    ##my $parseDate = DateTime::Format::Strptime->new(
    ## pattern   => '%Y-%m-%d');

    
    exit;
}


{
    my $in = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.rmDupes.pepFDR.sumIso.applyLC.untrans.tscCut.rebrand.working.noBait';

    my %seen = (); # have i seen this sid?
    my %topPrey = (); # how many times was this prey the top prey?
    open my $IN, "<", $in or die "can't read from $in. $!";

    my %row;
    while(getLineDP4APMS(\%row, $IN, 'line')){
	my $sid = $row{search_id};
	if (! exists $seen{$sid}) {
	    $seen{$sid} = 1;
	    $topPrey{$row{prey_ref}}++;
	    warn $row{line};
	}
    }

    my @ranked = sort { $topPrey{$b} <=> $topPrey{$a} } keys %topPrey;
    say "# found ".(0+ @ranked)." different top prey";
    say "# found ".(sum(values %topPrey))." total experiments";
    say join "\t", qw(rebrandTo times);
    say join "\t", $_, $topPrey{$_} for @ranked;
    exit;
}


{
    # collect DroID interactions into a single file
    my $listFile = '/home/glocke/DPiM/droid/d2u.list';
    my @f = readList($listFile);
    my %colNames = qw(
/camhpc/home/glocke/DPiM/droid/curagen_yth.txt.d2u curagen_yth
/camhpc/home/glocke/DPiM/droid/finley_yth.txt.d2u finley_yth
/camhpc/home/glocke/DPiM/droid/fly_other_physical.txt.d2u multidb
/camhpc/home/glocke/DPiM/droid/flybase_noDPIM1_ppi.txt.d2u flybase
/camhpc/home/glocke/DPiM/droid/human_interologs.txt.d2u hs_interologs
/camhpc/home/glocke/DPiM/droid/hybrigenics_yth.txt.d2u hybrigenics_yth
/camhpc/home/glocke/DPiM/droid/perrimon_coapcomplex.txt.d2u perrimon_apms
/camhpc/home/glocke/DPiM/droid/worm_interologs.txt.d2u ce_interologs
/camhpc/home/glocke/DPiM/droid/yeast_interologs.txt.d2u sc_interologs
); ## qw()
    my @cols = qw(multidb flybase curagen_yth finley_yth hybrigenics_yth
perrimon_apms sc_interologs hs_interologs ce_interologs);

    $colNames{$_} // die "can't find colname for '$_'" for @f;
    
    my %net;
    for my $f (@f) {
	say $f;
	my $thisCol = $colNames{$f};
	my %thisNet;
	open my $IN, "<", $f or die "can't read from $f. $!";
	while (<$IN>) {
	    next unless /^FBgn/;
	    my @spl = split;
	    my ($fb1, $fb2) = sort @spl[0..1];
	    die "'$fb1' not FBgn ($f, '$_')" unless $fb1 =~ /^FBgn\d+$/;
	    die "'$fb2' not FBgn ($f, '$_')" unless $fb2 =~ /^FBgn\d+$/;
	    $thisNet{$fb1}{$fb2} = 1;
	}
	close $IN;
	for my $fb1 (keys %thisNet) {
	    for my $fb2 (keys %{ $thisNet{$fb1} }) {
		$net{$fb1}{$fb2}{$thisCol} = 1;
	    }
	}
    }

    my $out = '/home/glocke/DPiM/droid/comprehensive.net';
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# comprehensive droid network made by dum.pl";
    say $OUT join "\t", qw(protein1 protein2), @cols, 'all';
    for my $fb1 (sort keys %net) {
	for my $fb2 (sort keys %{$net{$fb1}}) {
	    my @row = map { $net{$fb1}{$fb2}{$_} // 0 } @cols;
	    my $sum = sum(@row);
	    say $OUT join "\t", $fb1, $fb2, @row, $sum;
	}
    }
    close $OUT;
    exit;
}

{
    my @a;
    for (my $i=1.5; $i<=3.01; $i+=0.1) { 
	push @a, sprintf("%.1f", $i);
	
    }
    say join ",", @a;
    exit;
}

{
    ## 
    my $in = '/home/glocke/DPiM/prevDPIM/dpim1Net/publishedClusters/cell_5871_mmc4.sanitize.csv';
    my $out = '/home/glocke/DPiM/prevDPIM/dpim1Net/publishedClusters/cell_5871_mmc4.sanitize.mcl';
    open my $OUT, ">", $out or die "Can't write to $out: $!";
    open my $IN, "<", $in or die "Can't read from $in. $!";
    <$IN>; # skip first line
    while (my $line = <$IN>) {
	my @spl = split /,/, $line;
	my $id = shift @spl;;
	my @fbgns = grep /FBgn/, @spl;
	die "can't find fbgns in '$line'" unless 1 == @fbgns;
	$fbgns[0] =~ s/"//g;
	my @members = split /;/, $fbgns[0];
	say $OUT join "\t", $id, @members;
    }

    exit;
}

{
    my $parseAPMS = sub {
	my ($in, $reader) = @_;

	open my $IN, "<", $in or die "Can't read from $in. $!";
	my $x = <$IN>;
	
	my (%row, %ret);
	while ($reader->(\%row, $IN, 'line')) {
	    $ret{$row{search_id}}{$row{prey_ref}} = $row{total_peptides};
	}
	return %ret;
    };

    my $in = '/home/glocke/DPiM/prevDPIM/dpim1apms/dpim1_excel2tab.r6.07.updateFBgn.sumIso';
    my $out = '/home/glocke/DPiM/prevDPIM/dpim1apms/dpim1_excel2tab.r6.07.updateFBgn.sumIso.protFreq';
    my %apms = $parseAPMS->($in, \&getLineHyperspecAPMS);

    my %appear;
    my %tsc;
    for my $sid (keys %apms) {
	for my $prey (keys %{ $apms{$sid} }) {
	    $appear{$prey}++;
	    $tsc{$prey}+= $apms{$sid}{$prey};
	}
    }
    my @allPrey =  sort {$appear{$b} <=> $appear{$a}} keys %appear;

    my $nExpts = 0+ keys %apms;
    my %freq = map { $_ => $appear{$_}/$nExpts } @allPrey;
    

    
    my @cols = qw(fbgn    Count   Fraction        avg_tot_pep);
    open my $OUT, ">", $out or die "Can't write to $out: $!";
    say $OUT "# $0 found protein frequencies in $in";
    say $OUT "# totalExpts = $nExpts";
    say $OUT join "\t", @cols;
    for my $prey (@allPrey) {
	#die "can't find tsc{$prey}" unless exists $tsc{$prey};
	#die "can't find appear{$prey}" unless exists $appear{$prey};
	#die "can't find freq{$prey}" unless exists $freq{$prey};
	say $OUT join("\t", $prey, $appear{$prey}, $freq{$prey}, $tsc{$prey});
    }
    
    exit;
}


{
    my %h = ( a=> 1, b=>2, c=>3);
    my @k = qw(a b);
    my %h2 = %h{@k};
    die Dumper(\%h2);
}

{

    my $readClusters = sub {
	my ($ret, $clFile, $spr) = @_;

	my @clusters = readList2($clFile);
	#my @x = map { sprintf($spr, $_+1) } 0..$#clusters;
	my @newIDs;
	for my $i (0..$#clusters) {
	    my $k = sprintf($spr, $i+1);
	    push @newIDs, $k;
	    $ret->{$k} = $clusters[$i];
	}

	return @newIDs;
    };
    my $f = '/tmp/8niZA1snSy';

    my %cl;
    my @newIDs = $readClusters->(\%cl, $f, '%04d');
    #die Dumper(\%cl, \@newIDs);
    my $cl2 = $cl{'0002'};
    my $grepEdge = '~/DPiM/scripts/secondary/findEdgesWithNodes.pl';
    my $cmd = '';
    exit;
}

{
    my @a = ([1..10], [2..10], [3..4]);
    shift @$_ for @a;
    die Dumper(\@a);
}

{
    my $in  = '/home/glocke/DPiM/augRemap/nrBait_yRej_nBrand_11-16-2016/mclx/mcl.clusters.i2.txt';
    open my $IN, "<", $in or die "can't read from $in $!";
    while (<$IN>) {
	say $_ if /\d\d\d\d\.\d\d\d/;
    }
    exit;
}

{
    my @a = qw(1 2 3 4);
    my %h = (1=>1, 2=>2);
    @a = grep {! exists $h{$_} } @a;
    die Dumper(\@a);
}

{
    my %x = ( 1=> [0..10], 2=> [1..5] );
    my $one = $x{1};
    delete $x{1};
    die Dumper($one, \%x);
    exit;
}

{
    # parse the pfam dat file
    # their own stupid hmmfetch program doesn't work!!
    my $pfamFile = '/home/glocke/DPiM/pfam/Pfam-A.hmm.dat';
    my $outFile =  '/home/glocke/DPiM/pfam/Pfam-A.hmm.dat.my.tsv';

    open my $OUT, ">", $outFile or die "can't write to $outFile: $!";
    say $OUT "# $0 parsed $pfamFile into column format";
    say $OUT join "\t", qw(term name uniqName);
    open my $IN, "<", $pfamFile or die "Can't read from $pfamFile : $!";
    my $term;
    while (<$IN>) {
	chomp;
	if (/^#=GF AC/) {
	    my $acc = (split)[-1];
	    $acc =~ /(PF\d+)\.\d+/ or die "can't parse '$acc' '$_'";
	    $term = $1;
	    print $OUT $term, "\t";
	} elsif (/^#=GF DE/) {
	    s/#=GF DE   // or die "Can't can't clean '$_'";
	    say $OUT qq("$_"\t"$_ $term");
	}
    }

    exit;
}


{
    ## annotate nrBait and rebranded status on the final metadata table
    my %passNrb; # passNrb{$sid} = 1 iff this sid passes nr-bait
    {
	my $nrbFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.rmDupes.pepFDR.sumIso.applyLC.untrans.tscCut.rebrand.mindTheBranding.nrBait';
	open my $IN, "<", $nrbFile or die "Can't read from $nrbFile : $!";
	while (<$IN>) {
	    my $sid = (split)[0];
	    $passNrb{$sid} = 1;
	}
	close $IN;
    }

    my %rebrand;
    {
	my $rebrandFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.rmDupes.pepFDR.sumIso.applyLC.untrans.tscCut.rebrand.working.noBait.dum.rebrand.statsBySID';
	readColsHashRef(\%rebrand, $rebrandFile, [qw(search_id bait_ref)]);
    }

    my $annFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_110816.txt';
    my $outFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_112116.txt';
    open my $IN, "<", $annFile or die "Can't read from $annFile : $!";
    open my $OUT, ">", $outFile or die "can't write to $outFile: $!";
    while (<$IN>) {
	chomp;
	if (/^id_string/) {
	    say $OUT join "\t", $_, qw(nrBaitWinner rebrandedBait);
	} else {
	    my $sid = (split /\t/)[-2];
	    die Dumper($sid, $_) if $sid < 100000;
	    say $OUT join "\t", $_, $passNrb{$sid}//0, $rebrand{$sid}//'NA';
	}
    }
    close $OUT;
    close $IN;
    exit;
}

{
    my $s = "<seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <clan>";
    my @spl = split /> </, $s;
    $_ =~ s/ /_/ for @spl;
    $spl[0] =~ s/<//;
    $spl[-1] =~ s/>//;
    $_ = "'$_'" for @spl;
    say "c(", (join ", ", @spl), ")";

    exit;
}

{
    my $f = '/home/glocke/DPiM/prevDPIM/dpim1Net/DPIM1_scores.r6.07.updateFBgn.i4.annotated.preyPrey';
    my @cols = readHeader($f);
    my @read = readCols($f, \@cols, undef, "\t");
    say "ncols = ", 0+ @read;

    my %baitSup;
    $baitSup{$_->{baitSupport}}++ for @read;
    die Dumper(\%baitSup);
    exit;
}

{
    my $f = '/home/glocke/DPiM/prevDPIM/dpim1Net/DPIM1_scores.r6.07.updateFBgn.i4.annotated.preyPrey';
    my $s = `grep 'reticulum","cytoplasm","endoplasmic' $f`;
    my @spl = split /\t/, $s;
    die Dumper(\@spl);
}

{
    ## select the key columns
    my $metaFile  = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_110816.txt';
    my $outFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.11-14-2016_reject.tsv';
    my @inCols = ("search_id", "Final_Updated_Bait_ID", "RAO_Determination", 
		  "rejected");
    my @read;
    readColsRef(\@read, $metaFile, \@inCols, undef, "\t");
    pop @inCols;
    my @outCols = qw(search_id bait retain);
    open my $OUT, ">", $outFile or die "can't write to $outFile: $!";
    say $OUT "# $0 selected key columns from $metaFile";
    say $OUT join "\t", @outCols;
    for my $row (@read) {
	if ('1' eq $row->{rejected} // die Dumper($row)) {
	    $row->{RAO_Determination} = 'Remove'
	}
	say $OUT join "\t", map { $row->{$_} // die Dumper($row) } @inCols;
    }

    exit;
}


{
    my $s = "'lipid particle','alcohol metabolic process','protein homodimerization activity','acetaldehyde metabolic process','alcohol dehydrogenase (NAD) activity','ethanol oxidation','protein complex','acetaldehyde dehydrogenase (acetylating) activity','alcohol catabolic process','microtubule associated complex','ethanol metabolic process','behavioral response to ethanol','cytosol'";
    my @spl = split /'/, $s;
    my @values = grep { length($_) > 1 } @spl;
    die Dumper(\@values);
}

{
    my $clFile = '/home/glocke/DPiM/augRemap/nrBait11-08-2016/mcl/mcl.clusters/mcl.clusters.nrBait.ppi.abc-format.i4.txt.R_notEASE.GOTest';
    my @read;
    readColsRef(\@read, $clFile, ["minP", "minQ", "term", "sig", "size", 
				  "bestCnt", "allCnt", "nWayTie", "name", "i"],
		undef, "\t", "checkQuotes");

    my @x = @read[0..2];
    die Dumper(\@x);
    exit;
}

{
    ## revert bob's fbgn updates
    my $bobFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/bobsUpdate_plus.tsv';
    my %baitSwitch;
    readColsHashRef(\%baitSwitch, $bobFile, [qw(CurrentID SubmittedID)]);

    my $metaFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_110116.txt';
    my $outFile  = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_110816.txt';

    open my $IN, "<", $metaFile or die "can't read from $metaFile: $!";
    open my $OUT, ">", $outFile or die "can't write to $outFile: $!";
    while (<$IN>) {
	chomp;
	my @spl = split /\t/;
	if (@spl > 1 && exists $baitSwitch{$spl[3]}) {
	    say "changing '$spl[3]' to '$baitSwitch{$spl[3]}'";
	    $spl[3] = $baitSwitch{$spl[3]};
	    say $OUT join "\t", @spl;
	    
	} else {
	    say $OUT $_;
	}
    }
    close $IN;
    close $OUT;
    exit;
}

{
    my @a = qw(
FBgn0003663
FBgn0033466
FBgn0045949
); 
    my $aveLen = "REVdmel-all-translation-r6.07_TAGS_sorted_trEl_vir.aveLen.tsv"; 
    for (@a) { 
	my $grep = `grep $_ $aveLen`; 
	say "$_ \"$grep\"" if length($grep);
    }
    exit;
}

{
    my %test = ( a => 1==0, b => 1==0, c => 1==1, d => 1==1);
    my @winners = grep { $test{$_} } keys %test;
    say join ", ", @winners;
    exit;
}

{
    my %h = (a => { a => 1, b => 1, c => 1 }, x => { x=> 1, y => 1, z=> 1});
    $h{b} = $h{a};
    $h{c} = $h{a};
    $h{c}{d} = 1;
    $h{d} = $h{c};
    say " h{d} == h{c} ?", isTrue($h{d} == $h{c});
    say " h{d} == h{a} ?", isTrue($h{d} == $h{a});
    say " h{d} == h{x} ?", isTrue($h{d} == $h{x});
    
    for my $v (values %h) {
	say Dumper("huuu", $v);
    }
    die Dumper(\%h);
}


{
    ## map id_string to search_id
    
    my @files = readList('/home/glocke/DPiM/augRemap/apmsData/raw.list');
    my %ids;
    my %row;
    for my $f (@files) {
	open my $IN, "<", $f or die "can't read from $f: $!";
	while (getLineMayAPMS(\%row, $IN)) {
	    $ids{$row{id_string}}{$row{search_id}} = $row{sample_date};
	}
    }

    my @multiKeys = grep { 1 < keys %{ $ids{$_} }} keys %ids;
    say "# $0 compiled a map from id string to search_id";
    say "# $_" for @multiKeys;
    say "# number of multi keys = ", 0+ @multiKeys;
    say join "\t", qw(search_id id_string sample_date);
    for my $id (sort keys %ids) {
	my @sid = sort {$a <=> $b} keys %{ $ids{$id} };
	for my $sid (@sid) {
	    my $date = $ids{$id}{$sid};
	    my ($mm, $dd, $yy);
	    if (length($date) == 5) {
		$date =~ /^(\d)(\d\d)(\d\d)$/ 
		    or die "can't parse ids{$id}{$sid} (length 5)\n", 
		    Dumper($ids{$id});
		($mm, $dd, $yy) = ($1, $2, $3);
	    } else {
		$date =~ /^(\d\d)(\d\d)(\d\d)$/ 
		    or die "can't parse ids{$id}{$sid} (length 6)\n", 
		    Dumper($ids{$id});
		($mm, $dd, $yy) = ($1, $2, $3);
	    }
	    $date = sprintf("2%03d-%02d-%02d", $yy, $mm, $dd);
	    say join "\t", $sid, $id, $date
	}
    }
    
    exit;
}

{
    my @missing = readList('/home/glocke/DPiM/augRemap/apmsData/annotation/tmp.newBaitIDs.txt');

    my %symbols;
    {
	my @inCols = qw(Final_Updated_Bait_ID geneSymbol1     geneSymbol2 id_string);
	my @symbols;
	my $annFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_102216_noComments.txt';
	readColsRef(\@symbols, $annFile, \@inCols, undef, "\t");
	
	for my $row (@symbols) {
	    $symbols{$row->{Final_Updated_Bait_ID}} = 
		[$row->{geneSymbol1}, $row->{geneSymbol2}, $row->{id_string}];
	}
    }
    say join "\t", qw(fbgn symbol1 symbol2 id_string);
    for my $fbgn (@missing) {
	say join "\t", $fbgn, @{ $symbols{$fbgn} // die "can't find $fbgn" };
    }

}

{
    ## find any bait identities that aren't in the aveLen File
    my $baitFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.10-31-2016.tsv';
    my $lenFile = '/home/glocke/DPiM/nsaf/dmel-all-translation-r6.09.aveLen.tsv';
    my %len;
    readColsHashRef(\%len, $lenFile, [qw( fbgn fbgn )]);
    my @baits;
    readColsRef(\@baits, $baitFile, [qw(bait retain)]);
    my %notFound;
    for my $row (@baits) {
	next unless $row->{retain} eq 'Retain';
	my $b = $row->{bait};
	$notFound{$b} = 1 if ! exists $len{$b};
    }

    say $_ for sort keys %notFound;
    exit;
}



{
    my $f = '/home/glocke/DPiM/augRemap/apmsData/tmp.log';
    open my $IN, "<", $f or die "can't read from $f: $!";
    my %ids;
    my %row;
    while (getLineMayAPMS(\%row, $IN)) {
	$ids{$row{id_string}}{$row{search_id}} = 1;
    }
    die Dumper(\%ids);
    exit;
}

{
    my %id;
    my $in = '/home/glocke/DPiM/augRemap/apmsData/DPiM3_r607_protein_views.out_LDAmin7_PA_PSA01_160824.reorderDupes.dp4';
    open my $IN, "<", $in or die "Can't open $in: $!";
    my %row;
    while (getLineDP4_1APMS(\%row, $IN)) {
	$id{$row{ms_inst_run_id}}{$row{search_id}}=1;
    }

    for my $rid (keys %id) {
	my @k = keys %{ $id{$rid} };
	if (1 < @k) {
	    say "$rid $_" for @k;
	}
    }

    exit;
}

{
    my %h1 = qw( 1 a 2 b);
    my %h2 = qw( 3 c 4 d);
    %h1 = (%h1, %h2);
    die Dumper(\%h1);
}


{
    my @a = 0..10;
    splice(@a,2,0,'hi there!');
    say join ",", @a;
    exit;
}


{
    my %x = map { $_ => $_ } 1..10;

    #my @y;
    my %y = %x{3..5};
    say join ", ", keys %y;
    exit;
}


{
    my $ignoreFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/tmp.ignoreRef.txt';
    my %ignore = map {$_ => 1} readList($ignoreFile);

    my $header = join "\t", qw'search_id       bait_ref        prey_ref        total_peptides  sample_date     ms_inst_run_id';
    say $header;
    my $dataFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.pepFDR.sumIso';
    open my $IN, "<", $dataFile or die "can't read from $dataFile: $!";
    my %row;
    my %found;
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	next unless exists $ignore{$row{ms_inst_run_id}};
	next if $row{bait_ref} eq 'FBgn0000000';
	print $row{line};
	$found{$row{ms_inst_run_id}} = 1;
    }

    
    {
	my $metaFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.09-28-2016.tsv';
	my $out =  '/home/glocke/DPiM/augRemap/apmsData/ignore.butKept.metaData.txt';
	open my $IN, "<", $metaFile or die "can't read from $metaFile: $!";

	open my $OUT, ">", $out or die "can't write to $out: $!";
	while (<$IN>) {
	    my @spl = split;
	    next unless $spl[0] eq 'search_id' || exists $found{$spl[1]};
	    print $OUT $_;
	}
	close;
    }
    
    exit;
}

{
    my %check = map {$_ => 1} qw(f16937 f17567 f18511 f18512 f18613 f18514 f18516 f18517 f18581);
    my $metaFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable4.tap.correct1.tsv';
    open my $IN, "<", $metaFile or die "can't read from $metaFile: $!";
    while (<$IN>) {
	my @spl = split;
	next unless $spl[0] eq 'search_id' || exists $check{$spl[1]};
	print;
    }

    exit;
}

{
    my %check = map {$_ => 1} qw(f16937 f17567 f18511 f18512 f18613 f18514 f18516 f18517 f18581);
    my $header = 'search_id       bait_ref        prey_ref        total_peptides  sample_date     ms_inst_run_id';
    say $header;
    my $dataFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.pepFDR.sumIso.applyLC';
    open my $IN, "<", $dataFile or die "can't read from $dataFile: $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	print $row{line} if exists $check{$row{ms_inst_run_id}};
    }

    exit;
}

{
    ## get the APMS data for the runs that seem to be adding WAY too many edges
    ## unknown previous baits, in 2011
    my $findUnknown = sub {
	my ($metaFile, $baitFile) = @_;
	
	my @cols = qw(search_id prevbait sample_date);
	my @read;
	readColsRef(\@read, $metaFile, \@cols);
	@cols = qw(search_id retain);
	my %retain;
	readColsHashRef(\%retain, $baitFile, \@cols);

	my $nullBait = 'FBgn0000000';
	my %unknown; # unkonwn{sid} = 1 if prevbait == FBgn0000000
	for my $row (@read) {
	    my $sid = $row->{search_id};
	    next unless $row->{prevbait} eq $nullBait;
	    next unless $retain{$row->{search_id}} eq 'yes';
	    next unless $row->{sample_date} =~ /2011/;
	    $unknown{$sid} = 1;
	}
	    
	return %unknown;
    };

    my $metaFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable4.tap.correct1.tsv';
    my $baitFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.09-14-2016.tsv';
    my %sids = $findUnknown->($metaFile, $baitFile);
    
    my $dataFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.pepFDR.sumIso.applyLC.trans.0_05bottomTSC.rebranded';
    open my $IN, "<", $dataFile or die "can't read from $dataFile: $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	print $row{line} if exists $sids{$row{search_id}};
    }
    exit;
}


{
    my $metaFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable4.tap.correct1.doubles.txt';
    open my $IN, "<", $metaFile or die "can't read from $metaFile. $!";
    my %sids;
    while (<$IN>) {
	my @spl = split;
	my $sid = shift @spl;
	$sids{$sid} = 1;
    }
    close $IN;

    my $dataFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.pepFDR.sumIso';
    open $IN, "<", $dataFile or die "can't read from $dataFile: $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	print $row{line} if exists $sids{$row{search_id}};
    }
    exit;
}

{
    # run modularity.pl on all mcl cluster sets
    my @nets = qw( DPIM1 DPIM1_2 DPIM2 DPIM2_2 nrBait dropNet );
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20);
    my $outDir = '/home/glocke/DPiM/augRemap/apmsData/nrBait09-14-2016';
    my $date = '09-15-2016';

    my %dirs = (
	DPIM1 => '/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters',
	DPIM1_2 => '/home/glocke/DPiM/newDpim1_2/nrBait/mcl/mcl.clusters',
	DPIM2 => '/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters',
	DPIM2_2 => '/home/glocke/DPiM/newDpim2/nrBait_08-23-2016/mcl/mcl.clusters',
	nrBait => '/home/glocke/DPiM/dpim4/withInstr/nrBait/mcl/mcl.clusters',
	dropNet => '/home/glocke/DPiM/augRemap/apmsData/nrBait09-14-2016/mcl/mcl.clusters',
	);
    my %nets = (
	DPIM1 => '/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv',
	DPIM1_2 => '/home/glocke/DPiM/newDpim2/nrBait/newDpim2.nrBait.net',
	DPIM2 => '/home/glocke/DPiM/prevDPIM/dpim2_nrtap.120123.nrBait.58.17.network',
	DPIM2_2 => '/home/glocke/DPiM/newDpim2/nrBait_08-23-2016/newDpim2_08-2016.net',
	nrBait => '/home/glocke/DPiM/dpim4/withInstr/nrBait/nrBait.net',
	dropNet => '/home/glocke/DPiM/augRemap/apmsData/nrBait09-14-2016/nrBait.net',
	);

    my $meanMod = sub {
	my ($in, $connected) = @_;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my $ret = -1;
	while (<$IN>) {
	    if (($connected && /^# mean connected modularity ([\d\.-]+)/) ||
		(! $connected && /^# mean modularity ([\d\.-]+)/)) 
	    {
		$ret = $1;
		last;
	    }
	}
	die "can't find mean mod in $in." if $ret < -0.5;
	return $ret;
    };
    
    
    my $connected = undef; # seek mean modularity only among connected clusters
    my %modularity;
    my $scr = '~/DPiM/scripts/secondary/modularity.pl -mode mcl -mincluster 3';
    my $ext = 'mod';
    for my $n (@nets) {
	my $dir = $dirs{$n};
	my $net = $nets{$n};
	my %mods;
	for my $i (@i) {
	    my @f = glob("$dir/*i$i.txt");
	    die Dumper($n, $i, \@f) if 1 != @f;
	    my $in = $f[0];
	    my $out = "$in.$ext";
	    if (-e $out) {
		$mods{$i} = $meanMod->($out, $connected);
	    } else {
		my $cmd = "$scr -net $net -mod $in -out $out";
		$cmd .= " -human" if $n eq 'human' || $n eq 'CompPASS';
		say $cmd;
		$mods{$i} = `$cmd`;
		if ($connected) {
		    $mods{$i} = $meanMod->($out, $connected);
		}
	    }
	    $mods{$i} = sprintf("%.4f", $mods{$i});
	}
	$modularity{$n} = \%mods;
    }


    {
	my $outFile = "$outDir/modularity.$date.tsv";
	say $outFile;
	$outFile =~ s/modularity\./modularity.connected./ if $connected;
	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	say $OUT join "\t", 'i', @nets;
	for my $i (@i) {
	    my @row = map { $modularity{$_}{$i} } @nets;
	    say $OUT join "\t", $i, @row;
	}
    }
    exit;
}

{
    # collect clusterGOTest.pl results for all mcl cluster sets
    my @nets = qw( DPIM1 DPIM1_2 DPIM2 DPIM2_2 nrBait dropNet );
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20);
    my $outDir = '/home/glocke/DPiM/augRemap/apmsData/nrBait09-14-2016/';
    my $date = 'max5.09-15-2016';
    
    my %dirs = (
	DPIM1 => '/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters',
	DPIM1_2 => '/home/glocke/DPiM/newDpim1_2/nrBait/mcl/mcl.clusters',
	DPIM2 => '/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters',
	DPIM2_2 => '/home/glocke/DPiM/newDpim2/nrBait_08-23-2016/mcl/mcl.clusters',
	nrBait => '/home/glocke/DPiM/dpim4/withInstr/nrBait/mcl/mcl.clusters',
	dropNet => '/home/glocke/DPiM/augRemap/apmsData/nrBait09-14-2016/mcl/mcl.clusters',
	);
    
    my $getPercent = sub {
	my $in = shift;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my $ret = -1;
	while (<$IN>) {
	    if (/^# ([\d\.]+)% of clusters/) {
		$ret = $1;
		last;
	    }
	}
	die "can't find percent in $in." if $ret < 0;
	return $ret;
    };
    
    my $getNSig = sub {
	my $in = shift;
	my @sig = readCol($in, 'sig');
	return sum(@sig);
    };
    
    my %percents;
    my %nSigs;
    for my $n (@nets) {
	say $n;
	my $dir = $dirs{$n};
	
	my %percs;
	my %nSig;
	for my $i (@i) {
	    my @f = glob("$dir/*i$i.txt.max5.GOTest");
	    die Dumper($n, $i, 'max5', \@f) if 1 != @f;
	    $percs{$i} = $getPercent->($f[0]);
	    $nSig{$i} = $getNSig->($f[0]);
	}
	$percents{$n} = \%percs;
	$nSigs{$n} = \%nSig;
    }
    
    {
	my $outFile = "$outDir/GOTest.$date.tsv";
	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	say $OUT "# $0 grep'd and joined";
	say $OUT join "\t", 'i', @nets;
	for my $i (@i) {
	    my @row = map { $percents{$_}{$i} } @nets;
	    say $OUT join "\t", $i, @row;
	}
    }
    
    {
	my $outFile = "$outDir/GOCount.$date.tsv";
	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	say $OUT "# $0 grep'd and joined";
	say $OUT join "\t", 'i', @nets;
	for my $i (@i) {
	    my @row = map { $nSigs{$_}{$i} } @nets;
	    say $OUT join "\t", $i, @row;
	}
    }
    exit;
}


{

    my $readSeqContam = sub {
	my ($experiments, $seqFile) = @_;

	my $totalToRemove=0;
	open(INPUT, "<$seqFile") or die "Cannot open $seqFile\n";
	my $prey;
	my $firstInSeries = 1;
	while(my $buf = <INPUT>){
	    chomp $buf;
	    next if ($buf eq '');

	    my @tokens = split("\t", $buf);
	    if ($buf =~ "^\t") {
		$prey = $tokens[2];
		$firstInSeries = 1;
	    } else {
		if ($firstInSeries) {
		    $firstInSeries = 0;
		} else {
		    my $instr_run_id = $tokens[0];
		    my $bait = $tokens[1];
		    if ($bait ne $prey) {
			$totalToRemove++;
			$experiments->{$instr_run_id}{$prey} = 1;
		    }
		}
	    }
	}
	return $totalToRemove;
    };
    my $scEqual = sub {
	my ($source, $target) = @_;
	my $equal = 1;
	
	for my $rid (keys %$source) {
	    if (! exists $target->{$rid}) {
		warn "can't find $rid in target";
		$equal = undef;
		next;
	    }
	    for my $prey (keys %{$source->{$rid}}) {
		if (! exists $target->{$rid}{$prey}) {
		    warn "can't find target->{$rid}{$prey}";
		    $equal = undef;
		}
	    }
	}
	return $equal;
    };
    
    my $prevSeq = '/home/glocke/DPiM/newDpim1_2/DPiM1_05-05-2016Map.dp4.newBait.logPCut.sumIso.trans.seqContam';
    my $newSeq = '/home/glocke/DPiM/newDpim1_2/repeat09-14-2016/DPiM1_05-05-2016Map.dp4.newBait.logPCut.sumIso.trans.seqContam';
    my (%prev, %new);
    my $prevTotal = $readSeqContam->(\%prev, $prevSeq);
    my $newTotal  = $readSeqContam->(\%new, $newSeq);

    say "prevTotal = newTotal? ", isTrue($prevTotal == $newTotal);
    say "scEqual(prev, new)? ", isTrue($scEqual->(\%prev, \%new));
    say "scEqual(new, prev)? ", isTrue($scEqual->(\%new, \%prev));
    exit;
}

{
    ## find supporting evidence for edges in augRemap/dpim1/oldNrBait
    ## not in ../nrBait
    ## this is a comparison between versions of the LC carryover detection
    my $missing = '/home/glocke/DPiM/augRemap/apmsData/dpim1/oldNrBait09-12-2016/oldNrBait.net.missingFromNrBait.net';
    my $data0 = '/home/glocke/DPiM/augRemap/apmsData/dpim1/augDpim1.dp4.newBait.pepFDR.sumIso.old.applyLC.trans.0_05tscFilter.nrBait';
    my $data1 = '/home/glocke/DPiM/augRemap/apmsData/dpim1/augDpim1.dp4.newBait.pepFDR.sumIso.applyLC.trans.0_05tscFilter.nrBait';
    my $outDir = '/home/glocke/DPiM/augRemap/apmsData/dpim1/oldNrBait09-12-2016/missing';

    my $coA = '~/DPiM/scripts/secondary/coAppearance.pl';
    ## -in apms -prot FBgn1,FBgn2  -out output < -mode *4cols*/6cols -sidout print.search_ids >

    open my $IN, "<", $missing or die "Can't read from $missing. $!";
    while (my $line = readHS($IN)) {
	my @spl = split /\s+/, $line;
	my $prot = join ",", @spl[0..1];
	my $protName = join "-", @spl[0..1];
	my $out0 = "$outDir/apms0$protName.tab";
	my $sidout0 = "$outDir/apms0$protName.json";
	my $cmd = "$coA -in $data0 -prot $prot -out $out0 -sidout $sidout0";
	say $cmd;
	system $cmd;
	my $out1 = "$outDir/apms1$protName.tab";
	my $sidout1 = "$outDir/apms1$protName.json";
	$cmd = "$coA -in $data1 -prot $prot -out $out1 -sidout $sidout1";
	say $cmd;
	system $cmd;
    }

    exit;
}

    

{
    ## get metadata for baits that were rejected from augRemap dpim1 data
    my $missing = '/home/glocke/DPiM/augRemap/apmsData/dpim1/augDpim1.dp4.newBait.pepFDR.sumIso.nodeydey';
    my %missingIDs;
    {
	open my $IN, "<", $missing or die "can't read from $missing: $!";
	my %row;
	while(getLineDP4APMS(\%row, $IN)) {
	    $missingIDs{$row{search_id}}=1;
	}
    }
    my $prev = '/home/glocke/DPiM/newDpim1_2/DPiM1_05-05-2016Map.dp4.newBait.logPCut.sumIso';
    my %prevIDs;
    {
	open my $IN, "<", $prev or die "can't read from $prev: $!";
	my %row;
	while(getLineDP4APMS(\%row, $IN)) {
	    $prevIDs{$row{ms_inst_run_id}}=1;
	}
    }
    
    
    my $metaFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable3.tap.correct1.tsv';
    my @cols = qw(search_id       ms_inst_run_id  prevbait        tap_id  newbait sample_date     id_string       rejected     wonky   tap_test        count);
    my @meta;
    readColsRef(\@meta, $metaFile, \@cols, 'line');

    say join "\t", @cols;
    for my $row (@meta) {
	print $row->{line} if exists $missingIDs{$row->{search_id}} && 
	    exists $prevIDs{$row->{ms_inst_run_id}};
    }
    exit;
}


{
    ## update the metadata table with tap_id information
    my $dir = '/home/glocke/DPiM/augRemap/apmsData/annotation';
    my $metaFile = "$dir/metaDataTable2.txt";
    my $tapFile = "$dir/FH_Plate_Contents_for_JM_052715_short.updateFBgn.txt";
    my $outFile = "$dir/metaDataTable2.tap.txt";

    my @cols = qw(ms_inst_run_id  bait    tap_id  sample_date     id_string       rejected        wonky);
    my @meta;
    readColsRef(\@meta, $metaFile, \@cols);
    my %tapHash;
    readColsHashRef(\%tapHash, $tapFile, [qw(TAP_Alias Flybase_ID)], "\t");

    my $unknown = 'FBgn0000000';

    @cols = qw(ms_inst_run_id  bait    tap_id tap_bait  sample_date     id_string       rejected        wonky tap_test);

    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    say $OUT "# $0 add tap_table info to metadata";
    say $OUT join "\t", @cols;

    for my $expt (@meta) {
	my $tap = $expt->{tap_id};
	my $tapBait = $tapHash{$tap} // $unknown;
	if (! defined $tapHash{$tap} && $tap =~ /FH(\d+)/) {
	    my $tapN = sprintf("%04d", $1);
	    $tapBait = $tapHash{"FH$tapN"} // $unknown;
	}
	$expt->{tap_bait} = $tapBait;

	my $tapTest;
	if ($expt->{tap_bait} eq $unknown && 
	    $expt->{bait} eq $unknown) 
	{
	    $tapTest = "both_unknown";
	} 
	elsif ($expt->{tap_bait} eq $unknown 
		 && $expt->{bait} ne $unknown) 
	{
	    $tapTest = "known_but_not_in_table";
	} 
	elsif ($expt->{tap_bait} ne $unknown && 
		 $expt->{bait} eq $unknown) 
	{
	    $tapTest = "discovered_by_tap_table";
	} 
	elsif ($expt->{tap_bait} eq $expt->{bait}) {
	    $tapTest = "agree";
	} 
	elsif ($expt->{tap_bait} ne $expt->{bait}) {
	    $tapTest = "DISAGREE";
	} else {
	    die "ehh??", Dumper($expt);
	}
	$expt->{tap_test} = $tapTest;

	say $OUT join "\t", map {$expt->{$_}} @cols;
    }
    close $OUT;
    exit;
}


{
    ## make a metadata table
    my @inFiles = qw(DPiM1_r607_protein_views.out_LDAmin7_PA_PSA01_160808
DPiM1_rejects_r607_protein_views.out_LDAmin7_PA_PSA01_160811
DPiM2_r607_protein_views.out_LDAmin7_PA_PSA01_160816
DPiM3_r607_protein_views.out_LDAmin7_PA_PSA01_160824);
    my $baitKeyFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.baitKey';
    my $dir = '/home/glocke/DPiM/augRemap/apmsData';
    my %baitKey;
    readColsHashRef(\%baitKey, $baitKeyFile, [qw(ms_inst_run_id newBait)]);
    
    my @capture = qw(ms_inst_run_id tap_id sample_date id_string);

    my %runs;
    for my $f (@inFiles) {
	say $f;
	open my $IN, "<", $f or die "Can't read from $f: $!";
	my %row;
	my $reject = 0;
	$reject = 1 if $f =~ /rejects/;
	while (getLineMayAPMS(\%row, $IN)) {
	    my $run = $row{search_id};
	    if (! exists $runs{$run}) {
		$runs{$run} = { map { $_ => $row{$_} } @capture };
		## convert date from MMDDyy to YYYY-MM-DD 
		$runs{$run}{sample_date} =~ /^(\d\d?)(\d\d)(\d\d)$/ 
		    or die "can't parse date for $run", Dumper(\%row);
		$runs{$run}{sample_date} = "20$3-$1-$2";
		$runs{$run}{bait} = $baitKey{$run} // 'FBgn0000000';
		$runs{$run}{rejected} = $reject;
		$runs{$run}{wonky} = 0;
		$runs{$run}{wonky} = 1 if $row{tap_id} =~ /\D$/;
		$runs{$run}{wonky} = 1 if $row{tap_id} =~ /[^FHL\d]/;
	    }
	}
    }

    my $outFile = "$dir/metaDataTable2.txt";
    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    my @cols = ('bait', @capture, 'rejected', 'wonky');
    say $OUT "# $0 make a metadata table";
    say $OUT join "\t", 'search_id', @cols;
    for my $id (sort keys %runs) {
	say $OUT join "\t", $id, map { $runs{$id}{$_} } @cols;
    }
    close $OUT;
    exit;
}

{
    ## count the number of baits in the nrBait/hyperspec input
    #my $in = 'augRemap.dp4.newBait.logP.sumIso.applyLC.trans.rebranded.0_05tscFilter.nrBait';
    my $in = 'augRemap.dp4.newBait.logP.sumIso.applyLC.trans';
    open my $IN, "<", $in or die "Can't read from $in. $!";
    
    
    my %baits = ();
    my %expts = ();
    my %row;
    ##while (getLineHyperspecAPMS(\%row, $IN)) {
    while (getLineDP4APMS(\%row, $IN)) {
	$baits{$row{bait_ref}} = 1;
	$expts{$row{search_id}} = 1;
    }

    say "number of baits = ", (0+ keys %baits);
    say "number of experiments = ", (0+ keys %expts);
    exit;
}


{
    ## check if Bob Obar's corrected bait id's are found in r6.07
    my @correct = qw(
FBgn0036117
FBgn0265192
FBgn0010292
FBgn0016974
FBgn0038577
FBgn0033968
FBgn0263121
FBgn0010905
FBgn0259721
FBgn0031939
FBgn0031016
FBgn0003721
FBgn0031038
FBgn0263121
FBgn0263121
FBgn0263121
FBgn0263121
FBgn0262516
FBgn0263121
FBgn0000046
FBgn0037177
FBgn0010359
FBgn0283462
FBgn0263750
FBgn0000043
FBgn0053502
FBgn0051274
FBgn0040233
FBgn0033439
FBgn0039944
FBgn0050015
FBgn0024754
FBgn0013278
FBgn0032321
FBgn0028864
FBgn0004103
FBgn0250869
FBgn0036606
FBgn0083967
FBgn0000046
FBgn0263750
FBgn0001091
FBgn0004103
FBgn0003889
FBgn0016974
FBgn0053293
FBgn0037912
FBgn0037912
FBgn0259195
FBgn0000557
FBgn0031728
FBgn0040233
FBgn0003683
FBgn0029825
FBgn0051835
FBgn0086611
FBgn0000078
FBgn0261549
FBgn0028956
FBgn0036791
FBgn0036872
FBgn0004103
FBgn0259721
FBgn0263121
FBgn0010905
FBgn0003721
FBgn0003889
FBgn0000043
FBgn0004403
FBgn0001091
FBgn0261353
FBgn0037185
FBgn0064126
FBgn0029676
FBgn0029525
FBgn0040233
FBgn0024754
FBgn0013278
FBgn0036606
FBgn0038462
FBgn0029676
FBgn0259916
FBgn0260388
FBgn0259711
FBgn0083967
FBgn0050015
FBgn0261673
FBgn0263121
FBgn0035488
FBgn0259211
FBgn0069354
FBgn0259139
FBgn0264442
FBgn0033842
FBgn0026403
FBgn0033733
FBgn0038058
FBgn0037037
FBgn0039637
FBgn0032250
FBgn0283510
FBgn0036691
FBgn0050392
FBgn0053057
FBgn0038903
FBgn0040890
FBgn0034129
FBgn0010395
FBgn0038903
FBgn0029535
FBgn0266451
FBgn0030809
FBgn0069913
FBgn0036014
FBgn0037086
FBgn0032296
FBgn0028487
FBgn0036180
FBgn0260461
FBgn0051161
FBgn0039141
FBgn0036217
FBgn0031861
FBgn0035014
FBgn0034849
FBgn0050278
FBgn0029596
FBgn0034129
FBgn0031451
FBgn0028942
FBgn0028848
FBgn0030319
FBgn0261714
FBgn0053116
FBgn0003371
FBgn0036708
FBgn0037713
FBgn0035593
FBgn0053293
FBgn0032340
FBgn0042175
FBgn0033733
FBgn0061359
FBgn0039114
FBgn0283472
FBgn0086757
FBgn0037146
FBgn0029154
FBgn0032753
FBgn0052243
FBgn0264357
FBgn0010762
FBgn0263198
FBgn0000153
FBgn0259678
FBgn0052104
FBgn0024558
FBgn0037573
FBgn0030240
FBgn0029958
FBgn0010434
FBgn0050395
FBgn0001987
FBgn0036534
FBgn0013772
FBgn0004117
FBgn0052675
FBgn0031395
FBgn0039068
FBgn0035265
FBgn0033481
FBgn0037543
FBgn0263933
FBgn0038135
FBgn0063493
FBgn0261353
FBgn0030196
FBgn0039358
FBgn0053293
FBgn0040066
FBgn0038126
FBgn0033836
FBgn0264742
FBgn0000644
FBgn0053105
FBgn0030374
FBgn0259168
FBgn0040305
FBgn0026620
FBgn0011476
FBgn0085248
FBgn0264747
FBgn0045073
FBgn0259785
);
    my $nsafFile = '/home/glocke/DPiM/nsaf/REVdmel-all-translation-r6.07_TAGS_sorted_trEl_vir.aveLen.tsv';
    my %inDat;
    readColsHashRef(\%inDat, $nsafFile, [qw(fbgn fbgn)]);

    for my $fb (@correct) {
	say $fb if ! exists $inDat{$fb};
    }

    exit;
}

{
    ## correct baits incorrectly labeled as coming from non melanogaster species
    ## see "Re: Correction to Previous Corrections", from Bob Obar,
    ##  ## Thu 9/1/2016 10:29 AM
    ## from same thread, from me, George Locke, esq., Thu 9/1/2016 6:16 PM
    
    my $manualFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/FH_Plate_Contents_for_JM_052715_short_manual.txt';
    my $prevFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/FH_Plate_Contents_for_JM_052715_short.txt';
    my @corrected = `grep -i manual $manualFile`;
    my %correct;
    for my $row (@corrected) {
	my @spl = split /\s+/, $row;
	$correct{$spl[0]} = $spl[1];
    }

    my %prev;
    readColsHashRef(\%prev, $prevFile, [qw(TAP_Alias Flybase_ID)], "\t");

    my %incorrect = map { $_ => $prev{$_} } keys %correct;
    die Dumper(\%correct, \%incorrect);
    exit;
}

{
    my $p = 'MGNKCCSKRQDQELALAYPTGGYKKSDYTFGQTHINSSGGGNMGGVLGQKHNNGGSLDSRYTPDPNHRGPLKIGGKGGVD
IIRPRTTPTGVPGVVLKRVVVALYDYKSRDESDLSFMKGDRMEVIDDTESDWWRVVNLTTRQEGLIPLNFVAEERSVNSE
DWFFENVLRKEADKLLLAEENPRGTFLVRPSEHNPNGYSLSVKDWEDGRGYHVKHYRIKPLDNGGYYIATNQTFPSLQAL
VMAYSKNALGLCHILSRPCPKPQPQMWDLGPELRDKYEIPRSEIQLLRKLGRGNFGEVFYGKWRNSIDVAVKTLREGTMS
TAAFLQEAAIMKKFRHNRLVALYAVCSQEEPIYIVQEYMSKGSLLDFLREGDGRYLHFEDLIYIATQVASGMEYLESKQL
IHRDLAARNVLIGENNVAKICDFGLARVIADDEYCPKQGSRFPVKWTAPEAIIYGKFSIKSDVWSYGILLMELFTYGQVP
YPGMHSREVIENIERGFRMPKPTNHYFPDNIYQLLLQCWDAVPEKRPTFEFLNHYFESFSVTSEVPYREVQD';
    $p =~ s/\s//g;

    
    
    exit;
}

{
    ## test if gene symbol assignments have changed from r6.07 to r6.09
    my $r607File = '/home/glocke/DPiM/fbgn_id2_2col_08-30-2016.txt';
    my $r609File = '/home/glocke/DPiM/fbgn_id2_2col_03-17-2016.txt';
    my (%r607, %r609);
    readColsHashRef(\%r607, $r607File, [qw(fbgn symbol)]);
    readColsHashRef(\%r609, $r609File, [qw(fbgn symbol)]);

    my $mismatch = 0;
    for my $fb (sort keys %r607) {
	next unless ((! exists $r609{$fb}) || $r607{$fb} ne $r609{$fb});
	print "r607{$fb} = '".($r607{$fb} // 'not_there')."'";
	say "\tr609{$fb} = '".($r609{$fb} // 'not_there')."'";
	$mismatch++;
    }

    say "$mismatch mismatches out of ", (0+ keys %r607);
    
    exit;
}

{
    ## check to see which of the 207 disagreements were not present in the 
    ## r6.09 version
    my $dir = '/home/glocke/DPiM/augRemap/apmsData/annotation';
    my $r607File = "$dir/metaDataTable.tap.DISAGREE.txt";
    my $r609File = "$dir/r6.09/metaDataTable.tap.DISAGREE.txt";

    my $reader = sub {
	my ($in) = @_;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my %ret;
	while (<$IN>) {
	    my @spl = split;
	    $ret{$spl[0]} = $spl[3];
	}
	return %ret;
    };

    my %r607 = $reader->($r607File);
    my %r609 = $reader->($r609File);

    my %x = map { $_ => 1 } (keys(%r607), keys(%r609));
    my @keys = keys %x;
    for my $k (sort @keys) {
	next unless ((! exists $r607{$k}) || (! exists $r609{$k})) ||
	    $r607{$k} ne $r609{$k};
	print "r607{$k} = '".($r607{$k} // 'not_there')."'";
	say "\tr609{$k} = '".($r609{$k} // 'not_there')."'";
    }

    exit;
}

{
    ## update the metadata table with tap_id information
    my $dir = '/home/glocke/DPiM/augRemap/apmsData/annotation';
    my $metaFile = "$dir/metaDataTable.txt";
    my $tapFile = "$dir/FH_Plate_Contents_for_JM_052715_short.updateFBgn.txt";
    my $outFile = "$dir/metaDataTable.tap.txt";

    my @cols = qw(ms_inst_run_id  bait    tap_id  sample_date     id_string       rejected        wonky);
    my @meta;
    readColsRef(\@meta, $metaFile, \@cols);
    my %tapHash;
    readColsHashRef(\%tapHash, $tapFile, [qw(TAP_Alias Flybase_ID)], "\t");

    my $unknown = 'FBgn0000000';

    @cols = qw(ms_inst_run_id  bait    tap_id tap_bait  sample_date     id_string       rejected        wonky tap_test);

    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    say $OUT "# $0 add tap_table info to metadata";
    say $OUT join "\t", @cols;

    for my $expt (@meta) {
	my $tap = $expt->{tap_id};
	my $tapBait = $tapHash{$tap} // $unknown;
	if (! defined $tapHash{$tap} && $tap =~ /FH(\d+)/) {
	    my $tapN = sprintf("%04d", $1);
	    $tapBait = $tapHash{"FH$tapN"} // $unknown;
	}
	$expt->{tap_bait} = $tapBait;

	my $tapTest;
	if ($expt->{tap_bait} eq $unknown && 
	    $expt->{bait} eq $unknown) 
	{
	    $tapTest = "both_unknown";
	} 
	elsif ($expt->{tap_bait} eq $unknown 
		 && $expt->{bait} ne $unknown) 
	{
	    $tapTest = "known_but_not_in_table";
	} 
	elsif ($expt->{tap_bait} ne $unknown && 
		 $expt->{bait} eq $unknown) 
	{
	    $tapTest = "discovered_by_tap_table";
	} 
	elsif ($expt->{tap_bait} eq $expt->{bait}) {
	    $tapTest = "agree";
	} 
	elsif ($expt->{tap_bait} ne $expt->{bait}) {
	    $tapTest = "DISAGREE";
	} else {
	    die "ehh??", Dumper($expt);
	}
	$expt->{tap_test} = $tapTest;

	say $OUT join "\t", map {$expt->{$_}} @cols;
    }
    close $OUT;
    exit;
}


{
    ## count the number of baits in the nrBait/hyperspec input
    my $in = '/home/glocke/DPiM/newDpim2/newDpim2.dp4.newBait.logP.sumIso.updateFBgn.applyLC.trans.0_05.tscFilter.fixDate2.nrBait';
    open my $IN, "<", $in or die "Can't read from $in. $!";
    
    
    my %baits = ();
    my %row;
    while (getLineHyperspecAPMS(\%row, $IN)) {
	$baits{$row{bait_ref}} = 1;
    }

    say "number of baits = ", (0+ keys %baits);
    exit;
}

{
    my @excludedTerms = qw( GO:0003674 GO:0005515 GO:0005575 GO:0005634 
GO:0005737 GO:0005739 GO:0005829 GO:0008150 GO:0005829 GO:0008150 );
    my $s = join ", ", map { "'$_'" } @excludedTerms;
    say "c($s)";
    exit;
}


{
    my $dir = '/home/glocke/DPiM/human/biopStyleData/nrBait/qdir';
    my @read = glob("$dir/*sim*.o*");

    my @nonZero = grep { -s $_ } @read;

    say for @nonZero;
    
    
    exit;
}

{
    # look for ms_inst_run_id that are non-sequential but look like they are
    # found none!
    my $in = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.newFBgn.rmDupes.sumIso';

    open my $IN, "<", $in or die "cna't read from $in. $!";
    my %row;
    my %runNum; 
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	my $rid = $row{ms_inst_run_id};
	$rid =~ /(\d+)$/ or die "can't parse $row{line}";
	$runNum{$1+0}{$rid} = 1;
    }

    for my $n (sort {$a <=> $b} keys %runNum) {
	my @k = keys %{ $runNum{$n} };
	say "$n - ", join ", ", @k;
    }
    exit;
}

{
    # collect DroID interactions into a single file
    my $listFile = '/home/glocke/DPiM/droid/d2u.list';
    my @f = readList($listFile);

    my %net;
    for my $f (@f) {
	say $f;
	my %thisNet;
	open my $IN, "<", $f or die "can't read from $f. $!";
	while (<$IN>) {
	    next unless /^FBgn/;
	    my @spl = split;
	    my ($fb1, $fb2) = sort @spl[0..1];
	    die "'$fb1' not FBgn ($f, '$_')" unless $fb1 =~ /^FBgn\d+$/;
	    die "'$fb2' not FBgn ($f, '$_')" unless $fb2 =~ /^FBgn\d+$/;
	    $thisNet{$fb1}{$fb2} = 1;
	}
	close $IN;
	for my $fb1 (keys %thisNet) {
	    for my $fb2 (keys %{ $thisNet{$fb1} }) {
		$net{$fb1}{$fb2}++;
	    }
	}
    }

    my $out = '/home/glocke/DPiM/droid/support.net';
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# compiled droid support net using dum.pl";
    say $OUT join "\t", qw(protein1 protein2 support);
    for my $fb1 (sort keys %net) {
	for my $fb2 (sort keys %{$net{$fb1}}) {
	    say $OUT join "\t", $fb1, $fb2, $net{$fb1}{$fb2};
	}
    }
    close $OUT;
    exit;
}

{
    # cuz dos2unix never works
    # try mac2unix next time!
    # aka dos2unix -c Mac
    my $dir = '/home/glocke/DPiM/droid';
    my @files = `ls *txt`;
    #@files = ('/home/glocke/DPiM/human/CompPASS/BioPlex.v5884_HCIs.txt');
    chomp @files;
    for my $f (@files) {
	my $f2 = $f.".d2u";
	open my $IN, "<", $f or die "Can't read from $f: $!";
	open my $OUT, ">", $f2 or die "Can't write to $f2: $!";
	while (<$IN>) {
	    s/
/\n/;
	    say $OUT $_;
	}
	close $IN;
	close $OUT;
    }
    exit;
}


{
    my $list = '/home/glocke/DPiM/dpim4/withInstr/apmsData/sepBait/txt.list';
    my $single = '/home/glocke/DPiM/dpim4/withInstr/apmsData/sepBait/single';
    my $multi = '/home/glocke/DPiM/dpim4/withInstr/apmsData/sepBait/multi';
    for my $in (readList($list)) {
	open my $IN, "<", $in or die "Can't read from $in: $!";
	my %row;
	my %sid;
	while (getLineDP4APMS(\%row, $IN, 'line')) {
	    $sid{$row{search_id}} = 1;
	    last if 1 < keys %sid;
	}
	close $IN;
	if (1 < keys %sid) {
	    system("mv $in $multi/.")
	} else {
	    system("mv $in $single/.")
	} 
    }

    exit;
}

exit;
{
    my $in = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.newFBgn.rmDupes.sumIso';
    my $outDir = '/home/glocke/DPiM/dpim4/withInstr/apmsData/sepBait';
    
    open my $IN, "<", $in or die "Can't read from $in: $!";
    my %row;
    my ($OUT, $prevBait);
    $prevBait='';
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	my $bait = $row{bait_ref};
	if ($bait ne $prevBait) {
	    close $OUT if defined $OUT;
	    my $out = "$outDir/perBait.$bait.txt";
	    open $OUT, ">>", $out or die "Can't write to $out: $!";
	}
	print $OUT $row{line};
	$prevBait = $bait;
    }

    exit;
}

{
    my $pathFile = '/home/glocke/sysbio/IPA/kli3IPA/Ingenuity_PathwayNodes_03.28.2016.txt';
    my $out = '/home/glocke/sysbio/IPA/Ingenuity_PathwayNodes_03.28.2016.Human.txt';
    my $entrezCol = 7;
    open my $IN, "<", $pathFile or die "Can't read from $pathFile. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    <$IN>;
    $_ = <$IN>;
    print $OUT $_; # print header

    my $i = 0;
    while (my $line = <$IN>) {
	next unless $line =~ /^ING/;
	my @spl = split /\t/, $line;
	#say join ", ", @spl;
	next unless $#spl >= $entrezCol && $spl[$entrezCol] =~ /^\d+$/;
	print $OUT $line;
    }
    exit;
}



{
    say "hi there!";
    exit;
}

{
    my ($k, $W, $B, $N) = (0, 5, 10, 5);
    #gsl_cdf_hypergeometric_P($found, $white, $black, $outOf);
    my $p = gsl_cdf_hypergeometric_P($k, $W, $B, $N);

    say "p = $p";
    exit;
}

{
    # make a many-to-one entrez to symbol map 

    my $startFile = '/home/glocke/DPiM/human/nsaf/entrez2symbol.NP.tsv';
    my $secondFile = '/home/glocke/DPiM/human/nsaf/entrez2symbol.XP.YP.tsv';
    my $thirdFile = '/home/glocke/DPiM/human/nsaf/tmp/symbol2entrez.biomart.tsv';
    my $fourthFile = '/home/glocke/DPiM/human/CompPASS/symb2entrez.8467.manual.txt';
    my $out = '/home/glocke/DPiM/human/nsaf/entrez2symbol.concise.tsv';

    my %map;

    my $preferRegex = sub {
	my ($s1, $s2, $regex, $not) = @_;
	my $ret = $s1;
	if ($s1 =~ /$regex/ && $s2 !~ /$regex/) {
	    #say "$regex $s1 but not $s2";
	    return $s2 unless $not;
	    return $s1;
	} elsif ($s2 =~ /$regex/ && $s1 !~ /$regex/) {
	    #say "$regex $s2 but not $s1";
	    return $s1 unless $not;
	    return $s2;
	} else {
	    #say "$regex neither or both";
	    return undef;
	}
    };

    # pick the shorter symbol; if same length, pick alphanumerically
    sub picker {
	my ($s1, $s2) = @_;

	my %regex = qw( LRG_ 0 [A-Z] 1 orf 0);
	for my $reg (keys %regex) {
	    if (my $got = $preferRegex->($s1, $s2, $reg, $regex{$reg})) {
		return $got;
	    }
	}
	
	my $ret = $s1;
	if (length($s2) < length($s1)) {
	    $ret = $s2;
	} elsif (length($s1) == length($s2)) {
	    $ret = (sort($s1, $s2))[0];
	}
	return $ret;
    }
    
    {
	my @read;
	readColsRef(\@read, $startFile, [qw(entrez  symbol)]);
	for my $row (@read) {
	    my $en = $row->{entrez};
	    my $sy = $row->{symbol};
	    next if exists $map{$en} && $map{$en} eq $sy;
	    my $pick = $sy;
	    if (exists $map{$en}) {
		$pick = picker($map{$en}, $sy);
		warn "1: map{$en} = $map{$en} and $sy (pick $pick)\n";
	    }
	    $map{$en} = $pick;
	}
    }
    
    
    {
	my %map2;
	my @read;
	readColsRef(\@read, $secondFile, [qw(entrez  symbol)]);
	for my $row (@read) {
	    my $en = $row->{entrez};
	    next if exists $map{$en};
	    my $sy = $row->{symbol};
	    my $pick = $sy;
	    if (exists $map2{$en} && $map2{$en} ne $sy) {
		$pick = picker($map2{$en}, $sy);
		warn "2: map2{$en} = $map2{$en} and $sy (pick $pick)\n";
	    }
	    $map2{$en} = $pick;
	}
	for my $k (keys %map2) {
	    $map{$k} //= $map2{$k};
	}
    }
    
    {
	my %map2;
	my @read;
	readColsRef(\@read, $thirdFile, [qw(entrezgene  external_gene_name)]);
	for my $row (@read) {
	    my $en = $row->{entrezgene};
	    next if exists $map{$en};
	    my $sy = $row->{external_gene_name};
	    my $pick = $sy;
	    if (exists $map2{$en} && $map2{$en} ne $sy) {
		$pick = picker($map2{$en}, $sy);
		warn "3: map2{$en} = $map2{$en} and $sy (pick $pick)\n" 
	    }
	    $map2{$en} = $pick;
	}
	for my $k (keys %map2) {
	    $map{$k} //= $map2{$k};
	}
    }

    {
	my %map2;
	my @read;
	readColsRef(\@read, $fourthFile, [qw(entrez  symbol)]);
	for my $row (@read) {
	    my $en = $row->{entrez};
	    next if exists $map{$en};
	    my $sy = $row->{symbol};
	    my $pick = $sy;
	    if (exists $map2{$en} && $map2{$en} ne $sy) {
		$pick = picker($map2{$en}, $sy);
		warn "4: map2{$en} = $map2{$en} and $sy (pick $pick)\n";
	    }
	    $map2{$en} = $pick;
	}
	for my $k (keys %map2) {
	    $map{$k} //= $map2{$k};
	}
    }

    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# $0 made map from $startFile , $secondFile , and $thirdFile (in decreasing order of precedence)";
    say $OUT join "\t", qw(entrez symbol);
    for my $en (sort {$a <=> $b} keys %map) {
	if (! defined $map{$en}) {
	    warn "map{$en} undefined\n";
	    next;
	}
	say $OUT join "\t", $en, $map{$en};
    }
    close $OUT;
    exit;
}

{
    my $s1 = 'TLA123';
    my @s2 = qw(LRG_heyo tla456 C9orf72 TLA0123);

    # if $s1 matches and $s2 doesn't, prefer $s2, and vice versa
    # if both match or neither match, prefer neither
    my $preferRegex = sub {
	my ($s1, $s2, $regex, $not) = @_;
	my $ret = $s1;
	if ($s1 =~ /$regex/ && $s2 !~ /$regex/) {
	    #say "$regex $s1 but not $s2";
	    return $s2 unless $not;
	    return $s1;
	} elsif ($s2 =~ /$regex/ && $s1 !~ /$regex/) {
	    #say "$regex $s2 but not $s1";
	    return $s1 unless $not;
	    return $s2;
	} else {
	    #say "$regex neither or both";
	    return undef;
	}
    };
    
    my %regex = qw( LRG_ 0 [A-Z] 1 orf 0);
    #die Dumper(\%regex);
    #my @regex = qw( handi wipes);
    for my $s2 (@s2) {
	my $ret = $s1;
	for my $reg (keys %regex) {
	    #say "\t\t$s1 =~ /$reg/? ", isTrue( $s1 =~ /$reg/);
	    #say "\t\t$s2 =~ /$reg/? ", isTrue( $s2 =~ /$reg/);
	    if (my $got = $preferRegex->($s1, $s2, $reg, $regex{$reg})) {
		#say "\t$s2 catching on $reg";
		$ret = $got;
		last;
	    }
	}
	say "prefer $ret among $s1, $s2";
    }
    exit;
}

{
    # run modularity.pl on all mcl cluster sets
    my @nets = qw( nrBait meanBait DPIM1 DPIM1_2 human CompPASS );
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20);

    my $meanMod = sub {
	my ($in, $connected) = @_;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my $ret = -1;
	while (<$IN>) {
	    if (($connected && /^# mean connected modularity ([\d\.-]+)/) ||
		(! $connected && /^# mean modularity ([\d\.-]+)/)) 
	    {
		$ret = $1;
		last;
	    }
	}
	die "can't find mean mod in $in." if $ret < -0.5;
	return $ret;
    };
    
    my %dirs = (
	nrBait => '/home/glocke/DPiM/dpim4/withInstr/nrBait/mcl/mcl.clusters',
	meanBait => '/home/glocke/DPiM/dpim4/withInstr/meanBait/mcl/mcl.clusters',
	DPIM1 => '/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters',
	DPIM1_2 => '/home/glocke/DPiM/newDpim2/nrBait/mcl/mcl.clusters',
	human => '/home/glocke/DPiM/human/biopStyleData/meanBait/mcl/mcl.clusters',
	CompPASS => '/home/glocke/DPiM/human/CompPASS/mcl/mcl.clusters',
	);
    my %nets = (
	nrBait => '/home/glocke/DPiM/dpim4/withInstr/nrBait/nrBait.net',
	meanBait => '/home/glocke/DPiM/dpim4/withInstr/meanBait/meanBait.net',
	DPIM1 => '/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv',
	DPIM1_2 => '/home/glocke/DPiM/newDpim2/nrBait/newDpim2.nrBait.net',
	human => '/home/glocke/DPiM/human/biopStyleData/meanBait/human.06-23-2016.meanBait.net',
	CompPASS => '/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.known.net'
	);
    
    my $connected = 1; # seek mean modularity only among connected clusters
    my %modularity;
    my $scr = '~/DPiM/scripts/secondary/modularity.pl -mode mcl -mincluster 3';
    my $ext = 'mod';
    for my $n (@nets) {
	my $dir = $dirs{$n};
	my $net = $nets{$n};
	my %mods;
	for my $i (@i) {
	    my @f = glob("$dir/*i$i.txt");
	    die Dumper($n, $i, \@f) if 1 != @f;
	    my $in = $f[0];
	    my $out = "$in.$ext";
	    if (-e $out) {
		$mods{$i} = $meanMod->($out, $connected);
	    } else {
		my $cmd = "$scr -net $net -mod $in -out $out";
		$cmd .= " -human" if $n eq 'human' || $n eq 'CompPASS';
		say $cmd;
		$mods{$i} = `$cmd`;
		if ($connected) {
		    $mods{$i} = $meanMod->($out, $connected);
		}
	    }
	    $mods{$i} = sprintf("%.4f", $mods{$i});
	}
	$modularity{$n} = \%mods;
    }


    {
	my $outFile = "/home/glocke/DPiM/dpim4/withInstr/modularity.min3.06-27-2016.tsv";
	$outFile =~ s/modularity\./modularity.connected./ if $connected;
	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	say $OUT join "\t", 'i', @nets;
	for my $i (@i) {
	    my @row = map { $modularity{$_}{$i} } @nets;
	    say $OUT join "\t", $i, @row;
	}
    }
    exit;
}


{
    # count experiments and baits in dpim4 # prev dpim1
    my %ids;
    my %baits;
    #my $in = '/home/glocke/DPiM/prevDPIM/dpim1apms/dpim1_excel2tab.tsv';
    #my $in = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.newFBgn';
    my $in = '/home/glocke/DPiM/newDpim2/DPiM1_05-05-2016Map.dp4.newBait';
    open my $IN, "<", $in or die "can't read from $in. $!";
    while (<$IN>) {
	next unless /^\d/; 
	my @a=split;
	die "can't understand '$_'" if @a < 2;
	my ($id, $bait) = ($a[0], $a[1]); 
	$ids{$id} = 1;
	$baits{$bait} = 1;
    }

    say "n ids = ", 0+ keys %ids;
    say "n baits = ", 0+ keys %baits;
    exit;
}

{
    # collect clusterGOTest.pl results for all mcl cluster sets
    my @nets = qw( nrBait meanBait DPIM1 DPIM1_2 human CompPASS );
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20);

    my %dirs = (
	nrBait => '/home/glocke/DPiM/dpim4/withInstr/nrBait/mcl/mcl.clusters',
	meanBait => '/home/glocke/DPiM/dpim4/withInstr/meanBait/mcl/mcl.clusters',
	DPIM1 => '/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters',
	DPIM1_2 => '/home/glocke/DPiM/newDpim2/nrBait/mcl/mcl.clusters',
	human => '/home/glocke/DPiM/human/biopStyleData/meanBait/mcl/mcl.clusters',
	CompPASS => '/home/glocke/DPiM/human/CompPASS/mcl/mcl.clusters',
	);
    
    my $getPercent = sub {
	my $in = shift;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my $ret = -1;
	while (<$IN>) {
	    if (/^# ([\d\.]+)% of clusters/) {
		$ret = $1;
		last;
	    }
	}
	die "can't find percent in $in." if $ret < 0;
	return $ret;
    };
    
    my $getNSig = sub {
	my $in = shift;
	my @sig = readCol($in, 'sig');
	return sum(@sig);
    };
    
    my %percents;
    my %nSigs;
    for my $n (@nets) {
	my $dir = $dirs{$n};
	
	my %percs;
	my %nSig;
	for my $i (@i) {
	    my @f = glob("$dir/*i$i.txt.max5.GOTest");
	    die Dumper($n, $i, 'max5', \@f) if 1 != @f;
	    $percs{$i} = $getPercent->($f[0]);
	    $nSig{$i} = $getNSig->($f[0]);
	}
	$percents{$n} = \%percs;
	$nSigs{$n} = \%nSig;
    }
    
    {
	my $outFile = "/home/glocke/DPiM/dpim4/withInstr/GOTest.max5.06-27-2016.tsv";
	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	say $OUT join "\t", 'i', @nets;
	for my $i (@i) {
	    my @row = map { $percents{$_}{$i} } @nets;
	    say $OUT join "\t", $i, @row;
	}
    }
    
    {
	my $outFile = "/home/glocke/DPiM/dpim4/withInstr/GOCount.max5.06-27-2016.tsv";
	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	say $OUT join "\t", 'i', @nets;
	for my $i (@i) {
	    my @row = map { $nSigs{$_}{$i} } @nets;
	    say $OUT join "\t", $i, @row;
	}
    }
    exit;
}



{
    my $log = '/home/glocke/DPiM/human/biopStyleData/meanBait/GO0016021.log';
    my $netFile = '/home/glocke/DPiM/human/biopStyleData/meanBait/human.06-23-2016.meanBait.net';
    my %net;
    networkHashFromEdgeList(\%net, $netFile, undef, 'symmetric', undef, 'human');
    my %goNodes;
    {
	open my $IN, "<", $log or die "Can't read from $log. $!";
	while (<$IN>) {
	    my @a = split;
	    next unless @a > 1;
	    $goNodes{$a[0]}=1;
	}
    }
    my $cnt = 0;
    for my $k (sort {$a <=> $b} keys %goNodes) {
	if (exists $net{$k}) {
	    say "I found $k" ;
	}
    }
    exit;
}
    

{
    my $transFile = '/home/glocke/DPiM/human/GO/uniprot2entrez.merge.txt';
    my %trans;
    readColsHashRef(\%trans, $transFile, [qw(uniprot entrez)]);

    my $goFile = '/home/glocke/DPiM/human/GO/uniprot.go.tsv';
    my @go;
    readColsRef(\@go, $goFile, [qw(uniprot go)]);

    my $outFile = '/home/glocke/DPiM/human/GO/entrez.go.tsv';
    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    say $OUT "# $0 converted $goFile using $transFile";
    say $OUT join "\t", qw(entrez term);
    my %cantFind;
    for my $row (@go) {
	my ($uni, $term) = ($row->{uniprot}, $row->{go});
	if (! exists $trans{$uni} && ! exists $cantFind{$uni}) {
	    $cantFind{$uni} = 1;
	    warn ((join "\t", $uni, $term)."\n");
	    next;
	} elsif (! exists $trans{$uni}) {
	    next;
	}
	my $ent = $trans{$uni};
	say $OUT join "\t", $ent, $term;
    }
    close $OUT;
    exit;
    
}

{
    ## merge uniprot and biomart translations from uniprot to entrez

    my $picker = sub {
	my ($s1, $s2) = @_;
	return $s2 if ! defined $s1;
	my $ret = $s1;
	if (length($s2) < length($s1)) {
	    $ret = $s2;
	} elsif (length($s1) == length($s2)) {
	    $ret = (sort($s1, $s2))[0];
	}
	return $ret;
    };
    my $uniprotSite = '/home/glocke/DPiM/human/GO/allUni2entrez.txt';
    my $biomart = '/home/glocke/DPiM/human/GO/biomart.uniprot2entrez.txt';
    my $outFile = '/home/glocke/DPiM/human/GO/uniprot2entrez.merge.txt';
    my %trans;
    {
	open my $IN, "<", $uniprotSite or die "Can't read from $uniprotSite. $!";
	<$IN>; ## skip header
	while (<$IN>) {
	    chomp;
	    my @spl = split;
	    my ($uni, $ent) = @spl;
	    $trans{$uni} = $picker->($trans{$uni}, $ent)
	}
	close $IN;
    }

    {
	open my $IN, "<", $biomart or die "Can't read from $biomart. $!";
	<$IN>; ## skip header
	my %trans2;
	while (<$IN>) {
	    chomp;
	    my @spl = split;
	    my ($uni, $ent) = @spl;
	    next if exists $trans{$uni};
	    $trans2{$uni} = $picker->($trans2{$uni}, $ent)
	}
	close $IN;
	$trans{$_} = $trans2{$_} for keys %trans2;
    }


    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    say $OUT "# $0 merged $uniprotSite and $biomart conversions";
    say $OUT join "\t", qw(uniprot entrez);
    say $OUT join "\t", $_, $trans{$_} for sort keys %trans;
    close $OUT;
    exit;
}



{
    my $in = '/home/glocke/DPiM/dpim4/withInstr/meanBait/meanBait.net.node.GOHist';
    my @read;
    readColsRef(\@read, $in, [qw(term count name)]);
    my @cnt = map {$_ ->{count}} @read;
    say for @cnt;
    exit;
}



{
    # see if there are gene symbols with more than one entrez id in my file
    
    my $f = '/home/glocke/DPiM/human/nsaf/entrez2symbol.tsv';
    
    my @read;
    readColsRef(\@read, $f, [qw(entrez  symbol)]);

    my %map;
    for my $row (@read) {
	my ($en, $sy) = ($row->{entrez}, $row->{symbol});
	if (exists $map{$sy} && $map{$sy} != $en) {
	    say "map{$sy} = $map{$sy} and $en";
	} else {
	    $map{$sy} = $en;
	}
    }

    my $out = '/home/glocke/DPiM/human/nsaf/symbol2entrez.tsv';
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# $0 made map from $f";
    say $OUT join "\t", qw(entrez symbol);
    for my $sy (sort keys %map) {
	say $OUT join "\t", $map{$sy}, $sy;
    }
    close $OUT;
    exit;
}

{
    # check if there are entrez id's in the CompPASS network that aren't in my
    # refseq set
    my $compFile = '/home/glocke/DPiM/human/CompPASS/symb2entrez.8467.manual.txt';
    my $refFile = '/home/glocke/DPiM/human/nsaf/symbol2entrez.tsv';
    my %comp;
    #readColsHashRef(\%comp, $compFile, [qw(entrez symbol)]);
    my %ref;
    readColsHashRef(\%ref, $refFile, [qw(entrez symbol)]);
    exit;

    my %not;
    for my $e (keys %comp) {
	$not{$e} = $comp{$e} if ! exists $ref{$e};
    }

    my $out = '/home/glocke/DPiM/human/CompPASS/translate/notInRefSeq.txt';
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# $0 found ".(0+ keys %not)." entrez id's in $compFile not in $refFile";
    say join "\t", qw(entrez symbol);
    for my $k (sort {$a <=> $b} keys %not) {
	say $OUT join "\t", $k, $not{$k};
    }
    close $OUT;
    exit;
}    


{
    # cuz dos2unix never works
    # try mac2unix next time!
    my $dir = '/home/glocke/DPiM/droid';
    my @files = `ls *txt`;
    #@files = ('/home/glocke/DPiM/human/CompPASS/BioPlex.v5884_HCIs.txt');
    chomp @files;
    for my $f (@files) {
	my $f2 = $f.".d2u";
	open my $IN, "<", $f or die "Can't read from $f: $!";
	open my $OUT, ">", $f2 or die "Can't write to $f2: $!";
	while (<$IN>) {
	    s/
/\n/;
	    say $OUT $_;
	}
	close $IN;
	close $OUT;
    }
    exit;
}



{
    ## move failed simulations to a special place where they will never be seen
    ## again
    my $listFile = '/home/glocke/DPiM/human/biopStyleData/meanBait/qdir/sim.list';
    my @files = readList($listFile);
    my $rt = "10:00:00";
    my $newDir = '/home/glocke/DPiM/human/biopStyleData/meanBait/qdir/failedSims';
    
    for my $f (@files) {
	my $wc = `wc -l $f`;
	my $lines = (split /\s+/, $wc)[0];
	if ($lines < 100) {
	    my $cmd = "mv $f $newDir/.";
	    say $cmd;
	    system($cmd);
	    $f =~ s/(sim\d+\.)o/$1e/;
	    $cmd = "mv $f $newDir/.";
	    say $cmd;
	    system($cmd);
	}
    }


    exit;
}

{
    ## re-run failed simulations
    my $listFile = '/home/glocke/DPiM/human/biopStyleData/nrBait/qdir/sim.list';
    my @files = readList($listFile);
    my $rt = "10:00:00";
    
    for my $f (@files) {
	my $wc = `wc -l $f`;
	my $lines = (split /\s+/, $wc)[0];
	if ($lines < 100) {
	    $f =~ /qdir\/(.+_sim\d+)/ or die "can't parse $f 1";
	    my $job = $1;
	    my $scr = $f;
	    $scr =~ s/\.o\d+/.qsub.bash/ or die "can't parse $f";
	    my $cmd = "qsub -N $job $scr -l rt='$rt'";
	    say $cmd;
	    ###system($cmd);
	}
    }


    exit;
}

{
    my $notFile = 'notAFile';
    $notFile .= 'x' while -e $notFile;
    my $ret = system("cp $notFile whoWouldNameAFileLikeThis");
    say "got '$ret'";

    exit;
}

{
    my %out;
    my $outA = 'outA'; my $outB = 'outB';
    open $out{A}, ">", $outA or die "Can't write to $outA: $!";
    open $out{B}, ">", $outB or die "Can't write to $outB: $!";

    my @x = qw(a b c d);
    for my $x (@x) {
	say { $out{A} } $x;
	say { $out{B} } $x;
    }
    close $_ for values %out;
    exit;
}

{
    ## compare koko to Vti1b
    my $vtFile = '/home/glocke/DPiM/newDpim2/nrBait/Vti1b.edges';
    my $koFile = '/home/glocke/DPiM/newDpim2/nrBait/comparisons/koko.edges';
    my $vtID = 'FBgn0264751';
    my $koID = 'FBgn0264816';

    sub edgesf {
	my ($inFile, $fb) = @_;
	
	open my $IN, "<", $inFile or die "Can't read from $inFile: $!";
	my %ret;
	while (<$IN>) {
	    my @spl = split;
	    my ($p1, $p2) = @spl[0..1];
	    if ($p1 eq $fb) {
		$ret{$p2} = 1;
	    } elsif ($p2 eq $fb) {
		$ret{$p1} = 1;
	    }
	}
	return %ret;
    }

    my %vtEdges = edgesf($vtFile, $vtID);
    my %koEdges = edgesf($koFile, $koID);
    my $vtCnt = 0;
    for my $k (sort keys %vtEdges) {
	$vtCnt++ if exists $koEdges{$k};
	say $k if ! exists $koEdges{$k};
    }
    say "$vtCnt out of ".(0+ keys %vtEdges)." Vti1b edges are koko edges";
    my $koCnt = 0;
    for my $k (sort keys %koEdges) {
	$koCnt++ if exists $vtEdges{$k};
	say $k if ! exists $vtEdges{$k};
    }
    say "$koCnt out of ".(0+ keys %koEdges)." ko edges are Vti1b edges";

    exit;
}

{
    # look for injection numbers that end in 'a' or 'b' where there is *not* a
    #   run with the same number lacking the a/b
    my $summFile = '/home/kli3/Interactome/data/BioPlex_8400_03282016/gpsRes6_runSummary.tsv';

    my @read;
    readColsRef(\@read, $summFile, [qw(RAWFile InjectionNumber)]);
    my %runID;
    my %letterEnd;
    for my $row (@read) {
	$row->{InjectionNumber} =~ /(\d+)/ 
	    or die "can't find number in $row->{RAWFile}";
	my $n = $1;
	$runID{$n} //= [];
	push @{ $runID{$n} }, $row->{RAWFile};
	if ($row->{InjectionNumber} =~ /[a-z]$/) {
	    $letterEnd{$n} = 1;
	}
    }

    for my $n (sort {$a <=> $b} keys %letterEnd) {
	say Dumper($n, $runID{$n-1}, $runID{$n}, $runID{$n+1})
    }
    exit;
}

{
    # find no baits with more than two replicates (none found)
    my $summFile = '/home/kli3/Interactome/data/BioPlex_8400_03282016/gpsRes6_runSummary.tsv';
    my @baits = readCol($summFile, 'BaitGeneID');
    my %cnt = ();
    for my $b (@baits) {
	$cnt{$b}++;
    }

    for my $b (sort {$a <=> $b} keys %cnt) {
	say "$b\t$cnt{$b}" if $cnt{$b} > 2;
    }
    exit;
}

{
    # look for injection numbers that are non-sequential but look like they are
    my $summFile = '/home/kli3/Interactome/data/BioPlex_8400_03282016/gpsRes6_runSummary.tsv';
    my @read;
    readColsRef(\@read, $summFile, [qw(RAWFile InjectionNumber)]);
    my %runID;
    my %letterEnd;
    for my $row (@read) {
	$row->{InjectionNumber} =~ /(\d+)/ 
	    or die "can't find number in $row->{RAWFile}";
	my $n = $1;
	$runID{$n} //= [];
	push @{ $runID{$n} }, $row->{RAWFile};
	if ($row->{InjectionNumber} =~ /[a-z]$/) {
	    $letterEnd{$n} = 1;
	}
    }

    for my $n (sort {$a <=> $b} keys %letterEnd) {
	say Dumper($n, $runID{$n-1}, $runID{$n}, $runID{$n+1})
    }
    exit;
}

{
    # look for injection numbers that are non-sequential but look like they are
    my $summFile = '/home/kli3/Interactome/data/BioPlex_8400_03282016/gpsRes6_runSummary.tsv';
    my @read;
    readColsRef(\@read, $summFile, [qw(RAWFile InjectionNumber)]);
    my %runID;
    for my $row (@read) {
	$row->{InjectionNumber} =~ /(\d+)/ 
	    or die "can't find number in $row->{RAWFile}";
	$runID{$row->{RAWFile}} = $1;
    }

    my @k = sort { $runID{$a} <=> $runID{$b} } keys %runID;
    my @history;
    for my $id (@k) {
	$id =~ /^([a-z]+)/ or die "can't get alphanumerics from $id";
	push @history, $1;
	if (@history == 3) {
	    warn "$id: history = ".( join ", ", @history)."\n"
		if $history[0] ne $history[1] && $history[1] ne $history[2];
	    shift @history;
	} elsif (@history > 3) {
	    die "programmer error";
	}
    }

    exit;
}

{
    # look for ms_inst_run_id's that are non-sequential but look like they are
    
    #my $apms = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016';
    my $apms = '/home/glocke/DPiM/human/dpimStyleData/gpsRes6_dp4_05-13-2016.tsv';
    my %runID;
    {
	my @runID = readCol($apms, 'ms_inst_run_id');
	for my $id (@runID) {
	    $id =~ /(\d+)/ or die "can't parse $id";
	    $runID{$id} = $1;
	}
    }
    my @k = sort { $runID{$a} <=> $runID{$b} } keys %runID;
    my @history;
    for my $id (@k) {
	$id =~ /^([a-z]+)/ or die "can't get alphanumerics from $id";
	push @history, $1;
	if (@history == 3) {
	    warn "$id: history = ".( join ", ", @history)."\n"
		if $history[0] ne $history[1] && $history[1] ne $history[2];
	    shift @history;
	} elsif (@history > 3) {
	    die "programmer error";
	}
    }

    exit;
}


{
    my $x = 'RpLP0,Aprt,Argk,arm,Arr1,BicD,Cam,CkIIbeta,Pka-C1,chic,comt,cos,crn,CycA,CycB,D,Dip-B,dsh,ecd,Ef1alpha48D,Ef1alpha100E,EF2,E
no,exu,Fas2,dup,fu,fzy,Gapdh1,Gapdh2,Gdh,DCTN1-p150,Gpdh,gsb,GstD1,hook,Hsp26,Hsp27,Hsp83,ImpL3,Jra,kay,Khc,CIAPIN1,l(2)37Cc,Aats-asp,l(2)gl,dco,hyd,RpS5a,R
pS3,Rcc1,mago,mod,Atpalpha,Pgi,polo,r,ref(2)P,sbr,RanGAP,sgg,Sod,sqh,Su(dx),tub,alphaTub84B,betaTub56D,twi,Pp1-87B,boca,mts,unk,Cyp1,Rop,drk,N,sim,Su(H),Hrb
27C,hop,enc,tws,14-3-3zeta,Top1,Aats-glupro,GstD2,GstD4,GstD5,GstD7,gammaTub37C,hdc,GstS1,Myo61F,Ccs,dock,Aats-trp,Snap25,Uch-L5,cpb,Moe,rpr,twin,Hem,Gem3,A
ldh,l(2)dtl,mtd,FK506-bp2,MAPk-Ak2,Pdi,Rel,Rho1,Ugt,fax,Cctgamma,CkIalpha,Rpt2,RpS21,alpha-Est7,ktub,Rab10,Rab11,Nlp,Nurf-38,Rab1,dos,Rga,Dat,Drice,14-3-3ep
silon,Rpt6,ade5,Trxr-1,l(2)k01209,dnk,Pka-R2,Df31,qkr54B,Pp4-19C,Taldo,Ocrl,vig,san,IKKbeta,chico,WASp,MED6,Flo1,Clc,NTPase,Syx4,Pli,cry,scf,Inos,Nsun2,CDC45L,mei-P26,gus,msk,antdh,Dronc,Hus1-like,Hsc70Cb,Apc2,Eb1,Aats-tyr,Aats-ile,Aats-gly,Aats-gln,Aats-arg,Kap-alpha3,Kank,Akap200,l(1)G0255,Nxt1,gammaSnap1,Rpt5,Rpt4,Rpt3,Rpt1,Dak1,tant,Uba2,Ef1gamma,Atg5,alpha-PheRS,CG7033,CG15739,Upf1,CG1673,BthD,CG4872,Nprl2,CG13001,Pvf1,CG17259,Elp3,mxt,CG3008,CG31917,Oscillin,baf,Klp31E,CG17124,Csl4,CG12264,Ski6,Prosalpha6T,loqs,bsf,CG10417,Not3,Pngl,Strica,p47,tsu,trsn,Git,Prp8,RpS11,CG8525,CG6145,Rpn13,Rrp42,clu,Atg9,CG7997,Mapmodulin,MetRS,EndoB,TBCB,CG10527,Rrp4,sigmar,Pym,DCP1,yki,Cct2,CG11486,PHGPx,Sc2,Fit1,mad2,CG8209,GstO3,GstO2,GstO1,CG6282,CG11811,CG11652,viaf,Pop2,CG10688,Tdrd3,CG17027,nxf2,Smn,CG5577,CG5567,Mtr3,RhoGDI,CG9391,CG7611,CG7414,CG9804,GstZ1,bocks,FBXO11,Art4,Rrp46,CG6465,Tctp,Elp1,GstD9,Hexim,CG7265,Sra-1,Trax,CG12320,Smu1,CG4813,CG6726,atl,CG5913,CG6330,CG4849,Nph,krz,bor,grsm,Slbp,Atx2,key,SCAR,GstD10,Fsn,S-Lap3,btz,GstT1,Sk2,4EHP,Rtnl1,CG33156,GstE8,GstE7,GstE6,GstE5,GstE3,GstE10,Not1,mib2,Aats-asn,cbs,egg,psidin,Vrp1,alphaSnap,dUTPase,Prosalpha6,Pgk,Nedd4,Pp2A-29B,Diap1,Rrp40,Bet3,pic,TER94,Prosalpha3,hpo,capt,lic,larp,RanBPM,AGO1,Fs(2)Ket,DnaJ-1,Acsl,vret,Acn,Hk,Cdk7,nwk,BubR1,tgo,Flo2,CkIIalpha,Lst8,mgr,Cdep,sud1,CTPsyn,Nsf2,Gp210,Hsc70-4,gammaSnap2,Ptp61F,Pen,HnRNP-K';
    my $cnt=0;
    $cnt++ while $x =~ /,/g;
    say "found $cnt commas";
    exit;
}

{
    my @a = split //, 'MSRGSSAGFDRHITIFSPEGRLYQVEYAFKAIAQENITTVALKSGDCAVVATQKKVTEKNIVPETVTHLFRITKDIGCAMTGRIADSRSQVQKARYEAANFRYKYGYEMPVDVLCRRIADINQVYTQNAEMRPLGCSMVLIAYDNEIGPSVYKTDPAGYFSGFKACSVGAKTLEANSYLEKKYKPNLSEEKAIQLAISCLSSVLAIDFKPNGIEIGVVSKSDPTFRILDEREIEEHLTKIAEKD';
    my @b = split //, 'MSRGSSAGFDRHITIFSPEGRLYQVEYAFKAIAQENITTVALKSGDCAVVATQKKVTEKNIVPETVTHLFRITKDIGCAMTGRIADSRSQVQKARYEAANFRYKYGYEMPVDVLCRRIADINQVYTQNAEMRPLGCSMVLIAYDNEIGPSVYKTDPAGYFSGFKACSVGAKTLEANSYLEKKYKPNLSEEKAIQLAISCLSSVLAIDFKPNGIEIGVVSKSDPTFRILDEREIEEHLTKIAEKD';
    for my $i (0..$#a) {
	say "mismatch $i ( $a[$i], $b[$i] )" if $a[$i] ne $b[$i];
    }
    exit;
}

{
    my $new = '/home/glocke/DPiM/newDpim2/tmp.log';
    my $old = '/home/glocke/DPiM/prevDPIM/dpim1apms/tmp.log';
    my @prey;
    {
	my $minTSC = 3;
	
	my $in = $new;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	while (<$IN>) {
	    my @a = split;
	    next unless $a[3] >= $minTSC;
	    push @prey, $a[2];
	}
	close $IN;
    }

    my $in = $old;
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my @apms = <$IN>;
    close $IN;
    for my $p (@prey) {
	my @f = grep { /$p/ } @apms;
	#die Dumper(@f);
	say $p if 0==@f;
    }

    exit;
}

{
    my $missFile = '/home/glocke/DPiM/newDpim2/nrBait/dpim1_edges_missing.net';
    my %missing;
    networkHashFromEdgeList(\%missing, $missFile, undef, 'keepScore');
    my $prosA1 = 'FBgn0263121';
    my @pairs;
    my $minScore = 400;
    for my $p1 (keys %missing) {
	next if $p1 eq $prosA1;
	for my $p2 (keys %{ $missing{$p1} }) {
	    next if $p2 eq $prosA1;
	    push @pairs, [$p1, $p2] if $missing{$p1}{$p2} >= $minScore;
	}
    }


    my $scr = "~/DPiM/scripts/secondary/coAppearance.pl";
    my $outDir = "/home/glocke/DPiM/newDpim2/nrBait/tests";
    my $prev = '/home/glocke/DPiM/prevDPIM/dpim1apms/dpim1_excel2tab.updateFBgn.nrBait.tsv';
    my $new = '/home/glocke/DPiM/newDpim2/DPiM1_05-05-2016Map.dp4.newBait.logPCut.sumIso.trans.applyLC.nrBait';
    for my $pair (@pairs) {
	my $k = join "_", @$pair;
	my $p = join ",", @$pair;
	my $cmd = "$scr -in $prev -out $outDir/old.$k.tsv -prot $p -sidout $outDir/old.$k.sid.tsv";
	say $cmd;
	system($cmd);
	$cmd = "$scr -in $new -out $outDir/new.$k.tsv -prot $p -sidout $outDir/new.$k.sid.tsv";
	say $cmd;
	system($cmd);
    }
    exit;
}

{
    my $base = 'new2_hyper';
    my $dir = '/home/glocke/DPiM/newDpim2/hyper/qdir';
    my $cmd = "qsub -N $base $dir/$base.qsub.bash";
    say $cmd;
    system($cmd);
    for my $i (1..11) {
	my $job = sprintf("$base\_sim%03d", $i);
	my $cmd = "qsub -N $job $dir/$job.qsub.bash";
	say $cmd;
	system($cmd);
    }
    exit;
}

{
    my $mapFile = '/home/glocke/DPiM/human/nsaf/entrez2translatedLength.tsv';
    my $protFile = '/home/glocke/DPiM/human/dpimStyleData/gpsRes6_dp4_05-13-2016.inNCBI.countPrey';
    my @prot = readList($protFile);
    #say "read proteins";
    shift(@prot); # remove header
    my %symMap;
    #say "read map";
    readColsHashRef(\%symMap, $mapFile, [qw(entrez entrez)]);

    for my $p (@prot) {
	say $p unless exists $symMap{$p};
    }
    exit;
}

{
    ## download human fasta files
    my $url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot';
    sub hgFaa {
	my $chr = shift;
	return "human.$chr.protein.faa.gz";
    }
    my $outDir = "/home/glocke/DPiM/human/nsaf";
    for my $i (1..26) {
	my $fa = hgFaa($i);
	my $cmd = "wget $url/$fa &"; 
	say $cmd;
	system($cmd);
    }
    exit;
}

{
    # pull out Smn and TBPH experiments for harsha
    my $idFile = '/home/glocke/DPiM/dpim4/withInstr/apmsData/SmnTBPH.experiments';
    my %ids = map { $_=>1 } readList($idFile);
    #die Dumper(\%ids);
    
    my $in = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.newFBgn.rmDupes.sumIso';
    my $out = '/home/glocke/tmp/DPIM4_Smn_TBPH.txt';
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN, 'line')) {
	print $OUT $row{line} if exists $ids{$row{search_id}};
    }
    close $IN;
    close $OUT;
    exit;
}

{
    # run clusterGOTest.pl on all mcl cluster sets
    my @supp = qw( support0001
support0500
support1000
support1500
support1800
support1900
support1980
support2000 
DPIM1);
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20 30);

    sub clustFile {
	my ($supp, $i) = @_;
	if ($supp eq 'DPIM1') {
	    return "/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters/mcl.clusters.DPIM1_scores.updateFBgn.ppi.abc-format.i$i.txt";	    
	}
	return "/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/$supp/mcl2/mcl.clusters/mcl.clusters.$supp.ppi.abc-format.i$i.txt";
    }
    sub goHistFile {
	my $supp = shift;
	if ($supp eq 'DPIM1') {
	    return '/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.node.GOHist.tsv';
	}
	return "/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/$supp/$supp.net.node.GOHist.tsv";
    }
    my $scr = '~/DPiM/scripts/secondary/clusterGOTest.pl -godb /home/glocke/DPiM/flybase/gene_association.noComponent.fb';
    # usage: /home/glocke/DPiM/scripts/secondary/clusterGOTest.pl -in cluster.list -gohist wholeNet.nodewise.gohist -out output < -mode *mcl*/mode2>


    my %percents; 
    for my $s (@supp) {
	my $hist = goHistFile($s);
	#say $hist;
	my %percs;
	for my $i (@i) {
	    my $cl = clustFile($s, $i);
	    my $cmd = "$scr -in $cl -out $cl.GOtest -gohist $hist";
	    say $cmd;
	    my $perc = `$cmd`;
	    chomp $perc;
	    $percs{$i} = $perc;
	}
	$percents{$s} = \%percs;
    }

    my $outFile = "/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/allSupports.GOtest.tsv";
    open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
    say $OUT join "\t", 'i', @supp;
    for my $i (@i) {
	my @row = map { $percents{$_}{$i} } @supp;
	say $OUT join "\t", $i, @row;
    }

    exit;
}

{
    # apply MCL to all levels of support
    my $mcl = '/home/glocke/DPiM/MCL_network-enrichment-PPI/runAllFromPPI_GL.sh';
    my @supp = qw( support0001
support0500
support1000
support1500
support1800
support1900
support1980
support2000 
);

    my $baseDir = '/home/glocke/DPiM/dpim4/withInstr/consensus5_nets';
    for my $s (@supp) {
	my $net = "$baseDir/$s/$s.net";
	my $outDir = "$baseDir/$s/mcl2";
	my $job = $s.'_mcl';
	my $cmd = "$mcl $net $outDir $job";
	say $cmd;
	system($cmd);
    }
    exit;
}

{
    my $in = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.nrBait';
    my $outDir = '/home/glocke/DPiM/dpim4/withInstr/nrBait';
    my $hyper = "~/DPiM/cpp/hyperspec --in $in --sim 1";
    my $qsss = '/home/glocke/utilScripts/qsubWrap.pl';
    my $seed = time();
    for my $n (1..11) {
	my $job = sprintf("nrBaitSim%02d_05-04-2016", $n);
	my $out = "$outDir/$job.json";
	my $cmd = "$hyper --out $out --seed $seed";
	system("$qsss -cmd '$cmd' -job $job -rt '00:10:00'");
	$seed++;
    }
    exit;
}

{
    # find modularity against corum clusters
    my $modScript = '~/DPiM/scripts/secondary/modularity.pl';
    my $quoScript = '~/utilScripts/quoteCol.pl -col name';
    my @nets = readList('/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/net.list');
    push @nets, '/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv';
    my @supp = qw( support0001
support0500
support1000
support1500
support1800
support1900
support1980
support2000 
DPIM1);


    say "# got modularity using dum 'find modularity against corum clusters'";
    say "net\tmodularity";
    my $it = each_array(@nets, @supp);
    while (my ($netFile, $name) = $it->()) {
	my $modFile = "$netFile.corum.allMin3.clust";
	my $out1 = "$modFile.modularity.noquotes";
	my $out2 = "$modFile.modularity";
	if ($netFile =~ /DPIM1/) {
	    $modFile = $netFile;
	    $modFile =~ s/tsv/corum.allMin3.clust/;
	}
	my $cmd = "$modScript -net $netFile -mod $modFile -out $out1";
	my $score = `$cmd`;
	$cmd = "$quoScript -in $out1 -out $out2";
	system($cmd);
	
	chomp $score;
	say join "\t", $name, sprintf("%.4f", $score);
    }

    exit;
}

{
    my @min = qw( 1 500 1000 1500 1800 1900 1980 2000 );
    my $scr = "$ENV{DPSCR}/supportCutoff.pl";
    my $baseNet = '/home/glocke/DPiM/dpim4/withInstr/consensus5_05-04-2016/nets/support.network';
    my $baseOut = 'support';
    my $baseDir = '/home/glocke/DPiM/dpim4/withInstr/consensus5_nets';
    for my $m (@min) {
	my $mString = sprintf("%04d", $m);
	my $dir = "$baseDir/support$mString/";
	make_path($dir);
	my $cmd = "$scr -in $baseNet -out $dir/support$mString.net -minsupport $m";
	say $cmd;
	system($cmd);
    }
    exit;
}
{
    # test for duplicate CG's in the file
    my $in = '/home/glocke/DPiM/nsaf/Compiled_Dmel_Protein_Sequences_for_George.txt';
    my @read = readCol($in, 'CG_ID', "\t");

    my %cnt;
    for my $iso (@read) {
	$iso =~ /(CG\d+)/ or die "Can't parse $iso";
	$cnt{$1}++;
	# last if 10 < keys %cnt;
    }

    #my @dupes = sort grep {$cnt{$_} > 1} keys %cnt;
    #say "cnt{$_} = $cnt{$_}" for @dupes;
    say for sort keys %cnt;
    exit;
}

{
    sub parseHeeb {
	my ($in) = @_;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my %ret;
	while (<$IN>) {
	    chomp;
	    my @spl = split;
	    $ret{$spl[2]} = $spl[3];
	}
	return %ret;
    }
    
    # doublecheck findPairInAPMS.pl
    my $dp1File = '/home/glocke/DPiM/prevDPIM/dpim1apms/tmp.tsv';
    my $dp4File = '/home/glocke/DPiM/dpim4/withInstr/apmsData/tmp.tsv';

    my %dp1 = parseHeeb($dp1File); # dp1{bait} = tsc
    
}

{
    my $updated = '/home/glocke/DPiM/prevDPIM/dpim1apms/dpim1_excel2tab.updateFBgn.match_dpim2.tsv';
    my $orig = '/home/glocke/DPiM/prevDPIM/dpim1apms/dpim1_excel2tab.tsv';
    my $out = '/home/glocke/DPiM/prevDPIM/dpim1apms/missingBaits.txt';

    my @cols1 = qw(bait dp1ID dp4ID geneSymbol);
    my @cols2 = qw(search_id bait);
    my @matchTab = readCols($updated, \@cols1);
    my @orig = readCols($orig, \@cols2);
    my @missing = grep { $_->{dp4ID} == -1 } @matchTab;

    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT join "\t", qw(dp1ID updatedBait origBait updatedGeneSymbol);
    for my $miss (@missing) {
	my @sameID = grep { $_->{search_id} == $miss->{dp1ID} } @orig;
	my $origBait = $sameID[0]{bait};
	say $OUT join "\t", $miss->{dp1ID}, $miss->{bait}, $origBait
	    , $miss->{geneSymbol}
    }
    exit;
}

{
    # http://flybase.org/reports/FBgg0000189.html
    # VACUOLAR ATPASE V0 DOMAIN SUBUNITS
    my @nodes = qw(FBgn0028664
FBgn0039058
FBgn0032294
FBgn0028669
FBgn0262513
FBgn0262514
FBgn0035521
FBgn0028670
FBgn0028663
FBgn0028667
FBgn0032373
FBgn0028665
FBgn0028662
FBgn0038613
FBgn0038458
FBgn0028671
FBgn0028668
FBgn0262736);
    say join ",", @nodes;
    exit;
}

{
    #annotate notchome.csv
    my $in1 = "notchome.csv"; 
    my $in2 = "notchome.geneSymbol.tsv"; 
    my $out = "notchome.tsv"; 
    my @prefix;
    {
	open my $IN, "<", $in2 or die "can't read from $in2. $!";
	while (<$IN>) {
	    next if /^id/;
	    my @spl = split /\t/;
	    push @prefix, join "\t", @spl[0..1];
	}
    }

    open my $IN, "<", $in1 or die "can't read from $in1. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT join "\t", qw(id name proteins);
    my $i = 0;
    while (<$IN>) {
	print $OUT join "\t", $prefix[$i++], $_;
    }
    close $IN;
    close $OUT;
    exit;
}

{
    # given Harsha's notchome cluster list, write out each cluster in a csv file
    my $in = '/home/glocke/DPiM/notchome/cluster.membership.tsv';
    my $out = '/home/glocke/DPiM/notchome/notchome.csv';
    my %clusters;
    open my $IN, "<", $in or die "can't read from $in. $!";

    while (<$IN>) {
	chomp;
	my @a = split;
	$clusters{$a[1]}{$a[0]} = 1;
    }
    close $IN;

    open my $OUT, ">", $out or die "Can't write to $out. $!";
    for my $n (sort {$a <=> $b} keys %clusters) {
	say $OUT join ",", sort keys %{ $clusters{$n} };
    }
    close $OUT;
    exit;
}

{
    ## doublecheck that there is no search_id associated with multiple and ms_inst_run
    my $in = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.01-30-2016.newFBgn.rmDupes.sumIso';

    my %sidHash;
    my %iidHash;
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN)) {
	my $sid = $row{search_id};
	my $iid = $row{ms_inst_run_id};
	if (exists $sidHash{$sid}) {
	    warn "found $sid with $sidHash{$sid} and $iid" unless
		$sidHash{$sid} eq $iid;
	} else {
	    $sidHash{$sid} = $iid;
	}

	$iidHash{$iid}{$sid} = { bait => $row{bait_ref}, 
				 date => $row{sample_date} };
    }

    for my $iid (sort keys %iidHash) {
	my @sids = sort {$a <=> $b} keys %{ $iidHash{$iid} };
	next unless @sids > 1;
	say "found $iid with ";
	for my $sid (@sids) {
	    say "\t $sid bait $iidHash{$iid}{$sid}{bait}, date $iidHash{$iid}{$sid}{date}";
	}
    }
    exit;
}



{
    # run minScore on subsets of the consensus dataset
    my $outDir = '/home/glocke/DPiM/dpim4/withInstr/consensus4_04-04-2016/qdir/minScore';
    my @fullList = readList('/home/glocke/DPiM/dpim4/withInstr/consensus4_04-04-2016/qdir/real.list');
    my $smallSize = @fullList / 20;

    my $minScore = '~/DPiM/scripts/consensus/minScore.pl';
    my $qsss = '/home/glocke/utilScripts/qsubWrap.pl';

    my $n = 0;
    while (@fullList) {
	$n++;
	my @thisList;
	my $job = sprintf("minScore%02d", $n);
	my $thisOut = "$outDir/$job.net";
	my $listFile = "$outDir/$job.list";
	push @thisList, shift @fullList while @thisList < $smallSize;
	if (0 < @fullList && @fullList < $smallSize) {
	    push @thisList, @fullList;
	    @fullList = ();
	}
	open my $OUT, ">", $listFile or die "Can't write to $listFile. $!";
	say $OUT $_ for @thisList;
	close $OUT;
	my $cmd = "$minScore -in $listFile -out $thisOut";
	say $cmd;
	system("$qsss -job $job -cmd '$cmd' -rt '01:20:00'")
    }

    exit;
}
    #my $out = '/tmp/oogly_boogly';


{
    my @supp = qw( support0001
support0500
support1000
support1500
support1800
support1900
support1980
support2000 
DPIM1);
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20 30);
    
    say "c( ", (join ", ", map {qq("$_")} @supp), ")";
    say "c( ", (join ", ", @i), ")";
    exit;
}
{
    # list Clusters
    my @supp = qw( support0001
support0500
support1000
support1500
support1800
support1900
support1980
support2000 
DPIM1);
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20 30);
    

    sub clustFileY {
	my ($supp, $i) = @_;
	if ($supp eq 'DPIM1') {
	    return "/home/glocke/DPiM/prevDPIM/DPIM1_MCL/mcl.clusters/mcl.clusters.DPIM1_scores.updateFBgn.ppi.abc-format.i$i.txt";	    
	}
	return "/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/$supp/mcl2/mcl.clusters/mcl.clusters.$supp.ppi.abc-format.i$i.txt";
    }

    for my $s (@supp) {
	my $out = "/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/clust.$s.list";
	open my $OUT, ">", $out or die "can't write to $out. $!";
	for my $i (@i) {
	    say $OUT clustFileY($s, $i);
	}
	close $OUT;
    }
    exit;
}


{
    # double-check compareNetToDroID.pl
    my @sup3;
    {
	my $droid = '/home/glocke/DPiM/droid/support.updateFBgn.net';
	my %net;
	networkHashFromEdgeList(\%net, $droid, undef, 'score');
	for my $n1 (keys %net) {
	    for my $n2 (keys %{ $net{$n1} }) {
		push @sup3, [$n1, $n2] if $net{$n1}{$n2} == 3;
	    }
	}
    }
    
    my $dp1F = '/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv';
    my $dp4F = '/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/support0001/support0001.net';
    my (%dp1, %dp4);
    networkHashFromEdgeList(\%dp1, $dp1F, 'symm');
    networkHashFromEdgeList(\%dp4, $dp4F, 'symm');
    my ($cnt1, $cnt4) = (0,0);
    for my $edge (@sup3) {
	my ($n1, $n2) = @$edge;
	$cnt1++ if exists $dp1{$n1}{$n2};
	$cnt4++ if exists $dp4{$n1}{$n2};
    }

    say "cnt1, cnt4 = $cnt1, $cnt4";
    exit;
}

{
    # double-check compareNetToDroID.pl
    my @sup4 = map { [split] } split /\n/, 'FBgn0000083     FBgn0032140 
FBgn0000083     FBgn0035112 
FBgn0000426     FBgn0261068 
FBgn0003510     FBgn0034362 
FBgn0011335     FBgn0036248 
FBgn0020249     FBgn0028546 
FBgn0026088     FBgn0034461 
FBgn0026141     FBgn0033081 
FBgn0031106     FBgn0261049 
FBgn0031437     FBgn0031781 
FBgn0032004     FBgn0036919 
FBgn0032329     FBgn0033457 
FBgn0033467     FBgn0034842 
FBgn0034051     FBgn0263106 
FBgn0034461     FBgn0039413 
FBgn0039404     FBgn0261597 
FBgn0053193     FBgn0261456';
    
    my $dp1F = '/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv';
    my $dp4F = '/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/support0001/support0001.net';
    my (%dp1, %dp4);
    networkHashFromEdgeList(\%dp1, $dp1F, 'symm');
    networkHashFromEdgeList(\%dp4, $dp4F, 'symm');
    my ($cnt1, $cnt4) = (0,0);
    for my $edge (@sup4) {
	my ($n1, $n2) = @$edge;
	$cnt1++ if exists $dp1{$n1}{$n2};
	$cnt4++ if exists $dp4{$n1}{$n2};
    }

    say "cnt1, cnt4 = $cnt1, $cnt4";
    exit;
}


{
    # collect DroID interactions into a single file
    my @f = qw(
/home/glocke/DPiM/droid/curagen_yth.txt
/home/glocke/DPiM/droid/finley_yth.txt
/home/glocke/DPiM/droid/hybrigenics_yth.txt
/home/glocke/DPiM/droid/perrimon_coapcomplex.txt
);

    my %net;
    for my $f (@f) {
	say $f;
	open my $IN, "<", $f or die "can't read from $f. $!";
	while (<$IN>) {
	    next unless /^FBgn/;
	    my @spl = split;
	    my ($fb1, $fb2) = sort @spl[0..1];
	    die "'$fb1' not FBgn ($f, '$_')" unless $fb1 =~ /^FBgn\d+$/;
	    die "'$fb2' not FBgn ($f, '$_')" unless $fb2 =~ /^FBgn\d+$/;
	    $net{$fb1}{$fb2}++;
	}
	close $IN;
    }

    my $out = '/home/glocke/DPiM/droid/support.net';
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# compiled using dum.pl";
    say $OUT join "\t", qw(protein1 protein2 support);
    for my $fb1 (sort keys %net) {
	for my $fb2 (sort keys %{$net{$fb1}}) {
	    say $OUT join "\t", $fb1, $fb2, $net{$fb1}{$fb2};
	}
    }
    close $OUT;
    exit;
}

{
    # run clusterGOTest.pl on all mcl cluster sets
    my @supp = qw( support0001
support0500
support1000
support1500
support1800
support1900
support1980
support2000 );
    my @i = qw(1.2 1.4 1.6 1.8 2 3 4 6 10 15 20 30);

    sub clustFileX {
	my ($supp, $i) = @_;
	return "/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/$supp/mcl2/mcl.clusters/mcl.clusters.$supp.ppi.abc-format.i$i.txt";
    }
    sub goHistFileX {
	my $supp = shift;
	return "/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/$supp/$supp.net.GOHist.tsv";
    }
    my $scr = '~/DPiM/scripts/secondary/clusterGOTest.pl';
    # usage: /home/glocke/DPiM/scripts/secondary/clusterGOTest.pl -in cluster.list -gohist wholeNet.nodewise.gohist -out output < -mode *mcl*/mode2>

    my $header = "i\tpercent";
    for my $s (@supp) {
	my $hist = goHistFileX($s);
	my $outFile = "/home/glocke/DPiM/dpim4/withInstr/consensus4_nets/clust.$s.GOTest.tsv";
	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	say $OUT $header;
	for my $i (@i) {
	    my $cl = clustFileX($s, $i);
	    my $cmd = "$scr -in $cl -out $cl.GOtest -gohist $hist";
	    say $cmd;
	    my $perc = `$cmd`;
	    chomp $perc;
	    say $OUT "$i\t$perc";
	}
	close $OUT;
    }
    exit;
    
}

exit;
{
    sub count {
	my $s = shift;
	my %ret;
	while ($s =~ /(reverse_FBgn\d+)/g) {
	    $ret{$1}++;
	}
	$s =~ s/reverse_FBgn\d+//g;
	while ($s =~ /(FBgn\d++)/g) {
	    $ret{$1}++;
	}

	return \%ret;
    }
    my $s1 = 'reverse_FBgn123 FBgn456 reverse_FBgn789 ';
    my $s2 = 'FBgn123 FBgn456 FBgn789 ';

    say "count s1";
    say Dumper(count($s1));
    say "count s2";
    say Dumper(count($s2));
    exit;
}



{
    my $newStorable = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.01-30-2016.storable';
    # apms{bait}{search_id} = [ {prey_ref=>'fbgn', 'total_peptides'=>n},...]
    my $apms = retrieve($newStorable);

    my %uniqFBgn = ();
    for my $bait (keys %$apms) {
	$uniqFBgn{$bait}=1;
	for my $exp (values %{$apms->{$bait}}) {
	    for my $prey (@$exp) {
		$uniqFBgn{$prey->{prey_ref}}=1;
	    }
	}
    }

    die "found ".(0+ keys %uniqFBgn)." keys in $newStorable";
}

{
    # check if ambig.baitMap multiples appear correctly assigned
    my $ambig = '/home/glocke/DPiM/dpim4/withInstr/apmsData/ambig.baitMap.tsv';
    my $single = '/home/glocke/DPiM/dpim4/withInstr/apmsData/baitMap.tsv';
    my %ambigRuns = map {$_ => 1} readCol($ambig, 'ms_inst_run_id');

    for my $id (sort keys %ambigRuns) {
	say "**single**";
	system("grep $id $single");
	say "**ambig**";
	system("grep $id $ambig");
	say "";
    }

    exit;
}

{
    my $untidy = '/home/glocke/DPiM/newDpim1/oldMissing.txt';
    open my $IN, "<", $untidy or die "can't read from $untidy. $!";

    my @parsed;
    while (<$IN>) {
	chomp;
	push @parsed, [split];
	$parsed[-1][1] =~ s/T00:00:00//;
    }
    @parsed = sort {$a->[0] cmp $b->[0]} @parsed;
    @parsed = sort {$a->[1] cmp $b->[1]} @parsed;
    say join "\t", @$_ for @parsed;
    exit;
}

{
    # seek ms_inst_run_id that are not in the file julian sent on 03-16-2016
    my $old = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.01-30-2016';
    my $new = '/home/glocke/DPiM/newDpim1/DPiM1_r607_protein_views.out';
    die "parseDate not working";
    my $parseDate = DateTime::Format::Strptime->new(
       pattern   => '%Y-%m-%d');
    my %row;
    open my $IN, "<", $old or die "Can't read from $old. $!";
    my %date;
    warn "parse old";
    while (getLineDP4APMS(\%row, $IN)) {
	$date{$row{ms_inst_run_id}} = 
	    $parseDate->parse_datetime($row{sample_date});
    }

    # seek runs earlier than june 26, 2008
    my $earliest = $parseDate->parse_datetime("2011-11-20");
    my @k = grep {$date{$_} < $earliest} keys %date;

    warn "seek ids in new";
    my @notFound;   
    for my $id (@k) {
	my @g = `grep $id $new`;
	push @notFound, $id if @g == 0;
    }
    say "$_\t$date{$_}" for @notFound;
    exit;
}

{
    # make a table reporting the number of replicates for each bait
    sub reportNRep {
	my ($outFile, $repHash, $comments) = @_;

	my @keys = sort {$repHash->{$b} <=> $repHash->{$a}} keys %$repHash;

	open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";

	say $OUT $comments if $comments;
	say $OUT join "\t", qw(bait replicates);
	for my $k (@keys) {
	    say $OUT join "\t", $k, $repHash->{$k};
	}
	close $OUT;
    }
    
    my $newStorable = '/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.01-30-2016.storable';
    my $prevStorable = '/home/glocke/DPiM/dpim4/withInstr/forgotBait/all.combined.01-30-2016.rmDupes.sumIso.0_05.tscFilter.applyLC.storable';
    
    for my $f ($newStorable, $prevStorable) {
	
	my $apms = retrieve($f);
    
	my %nRep = map { $_ => 0+ keys %{ $apms->{$_} } } keys %$apms;
	my @oneRepl = grep {$nRep{$_} == 1} keys %nRep;
	my @twoRepl = grep {$nRep{$_} == 2} keys %nRep;
	my @multiRepl = grep {$nRep{$_} > 2} keys %nRep;

	my $outF = 'new.replicateTable.tsv';
	$outF = 'old.replicateTable.tsv' if $f =~ /forgot/;
	reportNRep($outF, \%nRep, 
		   "# ".(0+ @oneRepl)." single-replicate baits; ".
		   (0+ @twoRepl)." double-replicate baits; ".
		   (0+ @multiRepl)." multi-replicate baits; ");
    }
    exit;
}

{
    my @modFiles = readList('/home/glocke/DPiM/dpim4/withInstr/nets/allOrthos.modularity.list');
    say "support\tmodularity";
    for my $f (@modFiles) {
	$f =~ m/support\.(\d+)\./ or die "Can't parse $f";
	my $support = 0+ $1;
	open my $IN, "<", $f or die "can't read from $f";
	my $mod = 0;
	while (<$IN>) {
	    if (/mean modularity ([-\d\.]+)/) {
		$mod = $1;
		last;
	    }
	}
	say join "\t", $support, $mod;
    }
    exit;
}

{

    my @files = qw'
/home/glocke/DPiM/dpim4/withInstr/all.combined.01-30-2016
/home/glocke/DPiM/dpim4/withInstr/forgotBait/all.combined.01-30-2016.rmDupes.sumIso.0_05.tscFilter.applyLC
/home/glocke/DPiM/dpim4/withInstr/all.combined.01-30-2016.rmDupes.sumIso.withBait.0_05.tscFilter
/home/glocke/DPiM/dpim4/withInstr/forgotBait/cons1_01-30-2016.tsv
/home/glocke/DPiM/dpim4/withInstr/forgotBait/bestNet_02-01-2016.tsv
';

    sub countFBgns {
	my $in = shift;

	my %fbgn = ();
	
	open my $IN, "<", $in or die "Can't read from $in. $!";
	while (my $line = <$IN>) {
	    while ($line =~ /(FBgn\d+)/g) {
		$fbgn{$1} = 1;
	    }
	}

	return 0+ keys %fbgn
    }

    for my $f (@files) {
	say join "\t", $f, countFBgns($f);
    }

    exit;
}

{
    my $in = '/home/glocke/DPiM/dpim4/witInstr/all.combined.01-30-2016';
    #my $in = '/home/glocke/DPiM/dpim4/witInstr/all.combined.01-30-2016.rmDupes.sumIso.0_05.tscFilter.applyLC';

    my (%baits, %preys, %sids);

    open my $IN, "<", $in or die "Can't read from $in. $!";
    my %row;
    while (getLineDP4APMS(\%row, $IN)) {
	$baits{$row{bait_ref}} = 1;
	$preys{$row{prey_ref}} = 1;
	$sids{$row{search_id}} = 1;
    }

    say 0+ keys %baits;
    say 0+ keys %preys;
    say 0+ keys %sids;
    exit;
}

{
    sub min4 {
	my $in = shift;
	my $out = "$in.min4";

	my $header;
	readHeader($in, \$header);
	open my $OUT, ">", $out or die "Can't write to $out. $!";
	say $OUT "# selected entries with at least four proteins from $in";
	say $OUT $header;

	open my $IN, "<", $in or die "Can't read from $in. $!";
	while (<$IN>) {
	    my $line = $_;
	    next if /^#/;
	    next if /^term/;
	    my @spl = split /\t/;
	    my $prot = $spl[2];
	    my $cnt=0;
	    $cnt++ while $prot =~ /,/g;
	    print $OUT $line if $cnt >= 3;
	}
	close $IN;
	close $OUT;
    }

    # select clusters that have at least four members
    my @in = readList('/home/glocke/DPiM/dpim4/witInstr/GOcomplex.list');
    min4($_) for @in;

    exit;
}


{
    # parse csv files to extract a map from search_id to date and instr_id

    my @csv = qw'
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_011615.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_040215.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_042015.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_050815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_052215.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_060815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_061815.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_071015.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_072415.csv
/home/glocke/DPiM/dpim4/excelCsv/DPiM3_Raw_data_first_analysis_2015_073115.csv
'; #

    die "parseDate not working";
    my $parseDate;
    #my $parseDate = DateTime::Format::Strptime->new(
    #    pattern   => '%m%d%y');

    my %dates; # dates{$dateTime}{$search_id}=1
    my %dateTime;
    for my $f (@csv) {
	$f =~ /_(\d+)\.csv/ or die "can't parse $f.";
	my $date = $parseDate->parse_datetime($1);
	#say "$f $date";
	$dateTime{$date} = $date;

	open my $IN, "<", $f or die "Can't read from $f. $!";
	my $header = <$IN>; # strip header
	my @spl = split /,/, $header;
	my $sidCol = (grep { $spl[$_] =~ /search_id/i } 0..$#spl)[0];
	my $instrCol = (grep { $spl[$_] =~ /name/i } 0..$#spl)[0];
	die "can't find search_id in $f" unless defined $sidCol;
	$sidCol += 2 if $f =~ /072415.csv/; ## grumble grumble
	## these excel files are unique snowflakes!

	while (<$IN>) {
	    my @spl = split /,/;
	    my $sid = $spl[$sidCol] or die "f = $f, spl = ", join "\t", @spl;
	    my $instr = $spl[$instrCol];
	    $dates{$date}{$sid}=$instr;
	}
    }
    
    die "writeDate not working";
    my $writeDate;
    #my $writeDate = DateTime::Format::Strptime->new(
     #   pattern   => '%Y-%m-%d');
    
    say "search_id\tdate\tms_inst_run_id";
    for my $d ( sort {$dateTime{$a} <=> $dateTime{$b}} keys %dates) {
	my $date = $dateTime{$d};
	for my $sid (sort {$a <=> $b} keys %{$dates{$d}} ) {
	    say join "\t", $sid, $writeDate->format_datetime($date)
		, $dates{$d}{$sid};
	}
    }
    

    exit;
}

{
    # find rows removed after find_sequ that are different between new and old
    sub parseAPMS {
	my ($in, $reader) = @_;

	open my $IN, "<", $in or die "Can't read from $in. $!";
	my $x = <$IN>;
	
	my (%row, %ret);
	while ($reader->(\%row, $IN, 'line')) {
	    $ret{$row{search_id}}{$row{prey_ref}} = $row{line};
	}
	return %ret;
    }
    my $prevFile = '/home/glocke/DPiM/dpim4/lcTest/prev.apply_lc.removed.out';
    my $newFile = '/home/glocke/DPiM/dpim4/lcTest/my.apply_lc.removed.out';
    
    my %prev = parseAPMS($prevFile, \&getLineRawAPMS);
    my %new = parseAPMS($newFile, \&getLineDP4APMS);
    my @prevMiss;
    my @newMiss;
    my %k1 = map {$_ => 1} (keys %prev, keys %new);
    for my $k1 (keys %k1) {
	if (! exists $prev{$k1}) {
	    warn "can't find $k1 in prev" ;
	    push @prevMiss, join "", values %{$new{$k1}};
	    next;
	}
	if (! exists $new{$k1}) {
	    warn "can't find $k1 in new" ;
	    push @newMiss, join "", values %{$prev{$k1}};
	    next;
	}
	my %k2 = map {$_ => 1} (keys %{$prev{$k1}}, keys %{$new{$k1}});
	for my $k2 (keys %k2) {
	    if (! exists $prev{$k1}) {
		warn "can't find $k1->$k2 in prev" ;
		push @prevMiss, join "",  $new{$k1}{$k2};
		next;
	    }
	    if (! exists $new{$k1}) {
		warn "can't find $k1->$k2 in new" ;
		push @newMiss, join "", $prev{$k1}{$k2};
		next;
	    }
	}
    }

    say "unique to new";
    print for @prevMiss;
    print "\n\n\n\n";
    say "unique to prev";
    print for @newMiss;
    print "\n\n\n\n";
    say "newMiss = ", (0+ @newMiss), "\tprevMiss = ", (0+ @prevMiss);
    exit;
}

{
    # validate the output of find_sequential_contaminants against previous results
    my $prevFile = '/home/glocke/DPiM/dpim4/lcTest/test.seqContam';
    my $newFile = '/home/glocke/DPiM/dpim4/lcTest/my.seqContam';

    sub parseSeqContam {
	my ($in) = @_;

	open my $IN, "<", $in or die "can't read from $in. $!";
	my %ret; # $ret{$bait_$n} = { first, second, third, 
	my %count;
	my $contam; # the potential contaminant
	my $key;
	while (my $line = <$IN>) {
	    next if length($line) < 3;
	    chomp $line;
	    $_ = $line;
	    my @spl = split;
	    if ($line =~ /^\t/) {
		$contam = $spl[1];
		$count{$contam}++;
		if (0) {
		    $key = "$contam"."_$count{$contam}";
		    $ret{$contam} = {
			first => $spl[0],
			second => $spl[3],
			third => $spl[4],
		    };
		}
	    } else {
		my $run = $spl[0];
		my $bait = $spl[1];
		my $count = $spl[-1];
		$ret{$contam}{$bait}{$run} = $count;
	    }
	}
	return %ret;
    }

    my %prev = parseSeqContam($prevFile);
    my %new = parseSeqContam($newFile);

    my %k1 = map {$_ => 1} (keys %prev, keys %new);
    for my $k1 (keys %k1) {
	die "can't find $k1 in prev" if ! exists $prev{$k1};
	die "can't find $k1 in new" if ! exists $new{$k1};
	my %k2 = map {$_ => 1} (keys %{$prev{$k1}}, keys %{$new{$k1}});
	for my $k2 (keys %k2) {
	    die "can't find $k1->$k2 in prev" if ! exists $prev{$k1}{$k2};
	    die "can't find $k1->$k2 in new" if ! exists $new{$k1}{$k2};
	    my %k3 = map {$_ => 1} (keys %{$prev{$k1}{$k2}}, 
				    keys %{$new{$k1}{$k2}});  
	    for my $k3 (keys %k3) {
		die "can't find $k1->$k2->$k3 in prev" 
		    if ! exists $prev{$k1}{$k2}{$k3};
		die "can't find $k1->$k2->$k3 in new" 
		    if ! exists $new{$k1}{$k2}{$k3};
		die "mismatch at $k1->$k2->$k3" 
		    if $new{$k1}{$k2}{$k3} != $prev{$k1}{$k2}{$k3};
	    }
	}
    }
    say "Success!";
    exit;
}


{
    
    my $prevFile = '/home/glocke/DPiM/dpim4/lcTest/prev.commonContam';
    my $newFile = '/home/glocke/DPiM/dpim4/lcTest/my.commonContam';
    my @cols = qw(fbgn Count Fraction avg_tot_pep);

    my (@prev, @new);
    readColsRef(\@prev, $prevFile, \@cols, 'line');
    readColsRef(\@new, $newFile, \@cols, 'line');
    my %prev = map { $_->{fbgn} => $_ } @prev;
    my %new = map { $_->{fbgn} => $_ } @new;
    shift(@cols);
    my %k = map {$_ => 1} (keys %prev, keys %new);
    for my $k (keys %k) {
	die "can't find $k in prev" if ! exists $prev{$k};
	die "can't find $k in new" if ! exists $new{$k};
	for my $col (@cols) {
	    die "failure to match $k $col\n$prev{$k}{line}\n$new{$k}{line}"
		if $prev{$k}{$col} != $new{$k}{$col};
	}
    }
    say "Success!";
    exit;
}

{
    my $dir = '/home/glocke/DPiM/dpim4/consTest2_01-28-2016/qdir';
    my @files = `ls $dir/*.o*`;
    chomp @files;
    for my $f (@files) {
	open my $IN, "<", $f or die "Can't read from $f. $!";
	my $line;
	my $found = undef;
	while ($line = <$IN>) {
	    if ($line =~ /FBgn/) {
		$found = 1;
		last;
	    }
	}
	die "read eof in $f" unless $found;
	chomp $line;
	$_ = $line;
	my @spl = split;
	if ($spl[-1] == 745) {
	    #say "$f is real ($line)";
	} else {
	    #say "$f is simulated ($line)";
	    @spl = split /\./, $f;
	    my $last = pop(@spl);
	    push @spl, "sim001";
	    my $newF = join ".", @spl, $last;
	    my $cmd = "mv $f $newF";
	    #say $cmd;
	    #
	    $cmd =~ s/\.o/.e/g;
	    $cmd =~ s/\.sim001//;
	    say $cmd;
	    system($cmd);
	}
	close $IN;
    }

    exit;
}

{
    my @a = 1..100;
    for (1..10) {
	say $a[rand(0+ @a)];
    }
    exit;
}

{
    sub h {
	my ($x, $n) = @_;
	return 0 if $x == 0;
	return -$x * log($x/$n);
    }

    sub HXY {
	my ($a, $b, $c, $d, $n) = @_;
	my $hA = h($a, $n);
	my $hB = h($b, $n);
	my $hC = h($c, $n);
	my $hD = h($d, $n);
	if ($hA + $hD >= $hB + $hC) {
	    return  $hA + $hB + $hC + $hD - 
		h($a+$d, $n) - h($b+$c, $n);
	} else {
	    return h($a+$b, $n) + h($c+$d, $n);
	}
    }

    say 5*HXY(80, 16, 0, 4, 100);
    say HXY(80, 0, 16, 4, 100);
    exit;
}

{
    # remove node1 and node2 from mcl cluster output

    my @files = readList('/home/glocke/DPiM/MCL_network-enrichment-PPI/ppis-mcl-dpim3.09-25-2015.nrBait.77.44.network/mcl.clusters/short.list');
    my $outDir = '/home/glocke/DPiM/MCL_network-enrichment-PPI/ppis-mcl-dpim3.09-25-2015.nrBait.77.44.network/mcl.clusters/clean';
    
    for my $f (@files) {
	my $new = (split/\//, $f)[-1];
	my $cmd = "grep -v node $f > $outDir/$new";
	say $cmd;
	system($cmd);
    }
    exit;
}

{
    # concatenate CORUM and GO clusters
    my $goFile = '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4';
    my $corumFile = '/home/glocke/DPiM/oldDpim/dpim3.1/corumClusters.tab';
    my $both = '/home/glocke/DPiM/oldDpim/dpim3.1/corum.plus.go.complexes.tab.raw';
    my (@corum, @go);
    readColsRef(\@corum, $corumFile, ['proteins'], 'line');
    readColsRef(\@go,  $goFile, ['proteins'], 'line');
    open my $OUT, ">", $both or die "can't write to $both. $!";
    say $OUT "# joined $goFile and $corumFile";
    say $OUT join "\t", qw(id      name    proteins);
    say $OUT $_->{line} for @go;
    say $OUT $_->{line} for @corum;
    close $OUT;
    exit;
}

{
    # verify that try2StepHarness's coCluster function works for the true set
    my $listFile = '/home/glocke/DPiM/oldDpim/dpim3.1/mcode1_10-12-2015/tmp.log';
    my $truthFile = '/home/glocke/DPiM/oldDpim/dpim3.1/corumClusters.tab';

    sub coClusterMap {
	my ($clusterFile) = @_;
	my @read = readCol($clusterFile, 'proteins', "\t");
	my %ret;
	for my $row (@read) {
	    my @spl = split /,/, $row;
	    for my $prot1 (@spl) {
		for my $prot2 (@spl) {
		    $ret{$prot1}{$prot2} = 1;
		}
	    }
	}

	return %ret;
    }
    say "calculate map";
    my %coCluster = coClusterMap($truthFile);

    my @pairs;
    for my $p (readList($listFile)) {
	chomp $p;
	push @pairs, [ split ',', $p ];
    }

    say "count positives";
    my $conditionPositive = 0;
    for my $p (@pairs) {
	$conditionPositive++ if exists $coCluster{$p->[0]}{$p->[1]};
    }
    say "condition positive = $conditionPositive";
    
    exit;
}

{
    my $ncbiFile = '/home/glocke/DPiM/corum/ncbi_corumUniProtID_to_entrez.tab';
    
    my %ncbi = readMultiMap($ncbiFile);
    for my $u (keys %ncbi) {
	say $u if 1 < keys %{ $ncbi{$u} };
    }
    exit;
}


{
    my $ncbiFile = '/home/glocke/DPiM/corum/ncbi_corumUniProtID_to_entrez.tab';
    my $uniFile = '/home/glocke/DPiM/corum/uniprot_corumUniProtID_to_entrez.tab';
    
    my %ncbi = readMultiMap($ncbiFile);
    my %uni = readMultiMap($uniFile);
    
    for my $u (keys %ncbi) {
	if (! exists $uni{$u}) {
	    say "ncbi has a match for $u but uni does not";
	} else {
	    for my $e (keys %{ $ncbi{$u} }) {
		say "ncbi maps $u to $e but uni does not" 
		    unless exists $uni{$u}{$e};
	    }
	}
    }

    for my $u (keys %uni) {
	if (! exists $ncbi{$u}) {
	    say "uni has a match for $u but ncbi does not";
	} else {
	    for my $e (keys %{ $uni{$u} }) {
		say "uni maps $u to $e but ncbi does not" 
		    unless exists $ncbi{$u}{$e};
	    }
	}
    }
    exit;
}

{
    my $in = '/home/glocke/DPiM/corum/allComplexes.csv';
    open my $IN, "<", $in or die "can't read from $in. $1";
    <$IN>; # skip header

    my %species;
    my $spCol = 3;
    while (<$IN>) {
	my @spl = split /;/;
	my $sp = $spl[$spCol];
	$species{$sp} = 1;
    }

    say join ", ", sort keys %species;
    exit;
}

{
    my $in = '/home/glocke/DPiM/corum/allComplexes.csv';
    open my $IN, "<", $in or die "can't read from $in. $1";
    <$IN>; # skip header

    my %uniprot;
    my $uniCol = 4;
    while (<$IN>) {
	my @spl = split /;/;
	my @ids = split /,/, $spl[$uniCol];
	$_ =~ s/\(//g for @ids;
	$_ =~ s/\)//g for @ids;
	$uniprot{$_}=1 for @ids;
    }

    say $_ for sort keys %uniprot;
    exit;
}

{
    my @txt = `ls $ENV{PWD}/*txt`;
    my @a;
    for my $i (0..$#txt) { 
	$txt[$i] =~ m/i([\d\.]+)\.txt/ or die "Cant parse $txt[$i]"; 
	push @a, [0+ $1, $i];
    }
    @a = sort { $a->[0] <=> $b->[0] } @a;
    print $txt[$_] for map { $_->[1] } @a;
    exit;
}


{
    my @a = (0,1);
    my ($x,$y,$z) = @a[0..2];
    say $_ for ($x,$y,$z);
    exit;
}

{
    # test if any lines in a GOList are subsets of one another
    my $in = '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4';

    my %list;
    open my $IN, "<", $in or die "Can't read from $in. $!";
    while (<$IN>) {
	my $line = $_;
	next if /^#/;
	next if /^term/;
	my ($term, $name, $prot) = split /\t/;
	$list{$term} = {map {$_ => 1} split /,/, $prot};
    }
    #die Dumper(\%list);

    sub testSubset {
	my ($h1, $h2) = @_;
	my @k1 = keys %$h1;
	my @k2 = keys %$h2;
	my ($test1, $test2) = (1,1);
	for my $k (@k1) {
	    if (! exists $h2->{$k}) {
		$test1 = undef;
		last;
	    }
	}
	for my $k (@k2) {
	    if (! exists $h1->{$k}) {
		$test2 = undef;
		last;
	    }
	}

	my $ret = 0;
	if ($test1 && $test2) {
	    $ret = -1;
	} elsif ($test1) {
	    $ret = 1;
	} elsif ($test2) {
	    $ret = 2;
	}

	return $ret
    }

    my @terms = keys %list;
    for my $i (0..($#terms-1)) {
	my $t1 = $terms[$i];
	for my $j (($i+1)..$#terms) {
	    my $t2 = $terms[$j];
	    my $test = testSubset($list{$t1}, $list{$t2});
	    if ($test == 1) {
		say "$t1 is a subset of $t2";
	    } elsif ($test == 2) {
		say "$t2 is a subset of $t1";
	    } elsif ($test == -1) {
		say "$t1 and $t2 list the same genes";
	    }
	}
    }

    exit;
}



{
    # select clusters that have at least four members
    my $in = '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes';
    my $out = '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4';

    my $header;
    readHeader($in, \$header);
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "# selected entries with at least four proteins from $in";
    say $OUT $header;

    open my $IN, "<", $in or die "Can't read from $in. $!";
    while (<$IN>) {
	my $line = $_;
	next if /^#/;
	next if /^term/;
	my @spl = split /\t/;
	my $prot = $spl[2];
	my $cnt=0;
	$cnt++ while $prot =~ /,/g;
	print $OUT $line if $cnt >= 3;
    }

    exit;
}


{
    # count the overlap between GOList sets
    my $in = '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network.GOComplexes.min4';

    my %list;
    open my $IN, "<", $in or die "Can't read from $in. $!";
    while (<$IN>) {
	my $line = $_;
	next if /^#/;
	next if /^term/;
	my ($term, $name, $prot) = split /\t/;
	$list{$term} = {map {$_ => 1} split /,/, $prot};
    }
    #die Dumper(\%list);

    sub countOverlap {
	my ($h1, $h2) = @_;
	my @k1 = keys %$h1;
	my $ret = 0;
	for my $k (@k1) {
	    if (exists $h2->{$k}) {
		$ret++; 
	    }
	}
	return $ret
    }

    my @terms = keys %list;
    for my $i (0..($#terms-1)) {
	my $t1 = $terms[$i];
	my $size1 = 0+ keys %{$list{$t1}};
	for my $j (($i+1)..$#terms) {
	    my $t2 = $terms[$j];
	    my $size2 = 0+ keys %{$list{$t2}};
	    my $lap = countOverlap($list{$t1}, $list{$t2});
	    say "$lap\t- $t1 ($size1) and $t2 ($size2) overlap by $lap genes" if $lap > 0;
	}
    }

    exit;
}

{
    my $in = '/home/glocke/DPiM/flybase/Drosophila_melanogaster.gene_info';
    my $out = '/home/glocke/DPiM/flybase/fbgn_locusTag.tab';
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    say $OUT "fbgn locusTag";
    
    my @cols = (5,3);
    while (<$IN>) {
	next if /^#/;
	my @spl = split /\t/;
	my ($fbgn, $locus) = @spl[@cols];
	$fbgn =~ s/FLYBASE://;
	say $OUT "$fbgn\t$locus";
    }
    exit;
}



{
    # HG tests on various GOHists

    my %files = (
	base   => '/home/glocke/DPiM/dpim4/tscCut/base.GOHist.list',
	high   => '/home/glocke/DPiM/dpim4/tscCut/toohigh.GOHist.list',
	low    => '/home/glocke/DPiM/dpim4/tscCut/toolow.GOHist.list',
	filter => '/home/glocke/DPiM/dpim4/tscCut/tscFilter.GOHist.list',
	);

    #my %lists;
    #for my $k (keys %files) {
    #my @list = readList($files{$k});
    #$lists{$k} = \@list;
    #}
    my $baseCmd = "~/utilScripts/applyScriptToList.pl -s /home/glocke/DPiM/scripts/secondary/compareGOHists.pl -ext vsBase -in in1 -l2arg in2";

    my @comp = qw(high low filter);
    for my $vs (@comp) {
	my $cmd = "$baseCmd -l $files{$vs} -l2 $files{base}";
	say $cmd;
	system($cmd);
    }

    exit;
}


{
    # count proteins in multiple files

    my $scr = "$ENV{DPSCR}/secondary/countProtein.pl";
    sub apply {
	my ($in) = @_;
	
	system("$scr -in $in -out tscCut/$in.baitApp");
	system("$scr -in $in -out tscCut/$in.preyApp -mode prey");
	system("$scr -in $in -out tscCut/$in.baitTsc -scalebytsc");
	system("$scr -in $in -out tscCut/$in.preyTsc -scalebytsc -mode prey");
    }

#all.combined.10-21-2015.rmDupes.sumIso
    my @files = qw(
all.combined.10-21-2015.rmDupes.sumIso.0_05.filtered.toolow
all.combined.10-21-2015.rmDupes.sumIso.0_05.filtered.toohigh
all.combined.10-21-2015.rmDupes.sumIso.0_05.tscFilter
);
    chdir("/home/glocke/DPiM/dpim4");

    for my $f (@files) {
	say $f;
	apply($f);
    }
    exit;
}

{
    my $hi = '/home/glocke/DPiM/dpim4/all.combined.10-15-2015.0_05.filtered.toohigh';
    my $lo = '/home/glocke/DPiM/dpim4/all.combined.10-15-2015.0_05.filtered.toolow';

    {
	my @read;
	readColsRef(\@read, $hi, [qw(sample_date search_id)]);
	my %parse = map { $_->{search_id} => $_->{sample_date}} @read;
	my %count;
	$count{$_}++ for values %parse;
	say "hi:\t$_\t$count{$_}" for sort keys %count;
    }

    {
	my @read;
	readColsRef(\@read, $lo, [qw(sample_date search_id)]);
	my %parse = map { $_->{search_id} => $_->{sample_date}} @read;
	my %count;
	$count{$_}++ for values %parse;
	say "lo:\t$_\t$count{$_}" for sort keys %count;
    }
    exit;
}

{
    my $hi = '/home/glocke/DPiM/dpim4/all.combined.10-15-2015.0_05.filtered.toohigh';
    my $lo = '/home/glocke/DPiM/dpim4/all.combined.10-15-2015.0_05.filtered.toolow';
    {
	my @read = readCol($hi, 'search_id');
	my %parse = map { $_ => 1 } @read;
	my $found = 0+ keys %parse;
	say "hi has $found search_id's"
    }
    
    {
	my @read = readCol($lo, 'search_id');
	my %parse = map { $_ => 1 } @read;
	my $found = 0+ keys %parse;
	say "lo has $found search_id's"
    }
    exit;
}

{
    my $low = -5;
    my $high = -2;
    my @a = ($low..$high);
    my @b = ($low..-1);
    say join ", ", @a;
    say join ", ", @b;
    exit;
}


{
    # DPiM3_Raw_data_2015_text.txt has lines with _TAG peptides
    # these lines lack a "gene" column
    # this mucks up my getLine subroutine
    # so add a dummy column

    my $in = '/home/glocke/DPiM/dpim4/DPiM3_Raw_data_2015_text.txt';
    my $out = '/home/glocke/DPiM/dpim4/DPiM3_Raw_data_2015_text.colFix';
    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";

    while (<$IN>) {
	if (/TAG/) {
	    chomp;
	    print $OUT $_, "\tTAG\n";
	} else {
	    print $OUT $_;
	}
    }
    exit;
}

{
    my $s = 'junkGood123';
    $s =~ s/\.*(?=Good)\d+//;
    say $s;
    exit;
}

{
    my @a = 1..10;
    ($a[1], $a[0]) = @a[0..1];
    say join ", ", @a; # did work
    exit;
}

{
    my @a = 10..1;
    say join ", ", @a; # didn't work
    exit;
}


{
    # find the number of experiments in a commContam file
    # this is (supposed to be) in the comment at the head of the file
    sub totalExpts {
	my ($file) = @_;

	open my $IN, "<", $file or die "Can't read from $file. $!";
	my $ret;
	while (<$IN>) {
	    if (/totalExpts = (\d+)/) {
		$ret = $1;
	    }
	}
	die "can't find totalExpts in $file." unless defined $ret;
	return $ret;
    }
    
    my @readCols = qw( fbgn Count );
    my $minCount = 10; # must have at least $x Count to register
    # ret{fbgn} = number of experiments in which ~ appeared
    sub appearanceHash {
	my ($file) = @_;
	
	my @read;
	readColsRef(\@read, $file, \@readCols);
	
	my %ret;
	for my $row (@read) {
	    next unless $row->{Count} >= $minCount;
	    $ret{$row->{fbgn}} = $row->{Count};
	}

	return %ret;
    }

    my $R = Statistics::R->new();
    $R->startR;

    # here's how the arguments are interpreted
    #         smiling  not-smiling
    # people   $g1_c1   $g1_c2     
    # rocks    $g2_c1   $g2_c2
    sub fisherTest {
	my ($group1_count1, $group1_count2, $group2_count1, 
	    $group2_count2) = @_;

	$R->send("conTab <- matrix(c($group1_count1, $group2_count1, 
                  $group1_count2, $group2_count2), nrow=2)");
	$R->send("fisher.test(conTab)\$p.value");
	my $read = $R->read;
	$read =~ s/\[1\] //;
	return $read;
    }

    sub fisherTests {
	my ($tests) = @_;
	
	my $nTests = 0+ @$tests;
	
	my (@b, @bX, @n, @nX);
	for my $t (@$tests) {
	    push @b,  $t->{badPresent};
	    push @bX, $t->{badAbsent};
	    push @n,  $t->{notPresent};
	    push @nX, $t->{notAbsent};
	}

	$R->send("bVec <- c(".(join ", ", @b).")");
	$R->send("bXVec <- c(".(join ", ", @bX).")");
	$R->send("nVec <- c(".(join ", ", @n).")");
	$R->send("nXVec <- c(".(join ", ", @nX).")");
	$R->send("matrify <- function(i, m11, m12, m21, m22) {
    return(matrix(c(m11[i], m12[i], m21[i], m22[i]), nrow=2))
}");
	$R->send("conTabs <- lapply(1:$nTests, matrify, m11=bVec, m12=nVec, m21=bXVec, m22=nXVec)");
	$R->send("pp <- vapply(conTabs, function(m) fisher.test(m)\$p.value, 1.1)");
	$R->send("pp");
	my @p = readR($R);
	$R->send("p.adjust(pp, method='fdr')");
	my @q = readR($R);

	my $it = each_array(@$tests, @p, @q);	
	while (my ($t, $p, $q) = $it->()) {
	    $t->{p} = $p;
	    $t->{q} = $q;
	}

	return;
    }

    sub readR {
	my $R = shift;
	my @read = split /\n/, $R->read;

	s/^\s*\[\d+\]\s*// for @read; # strip ' [123] ' at start of line
	
	my @ret;
	for (@read) {
	    push @ret, split;
	}

	return @ret;
    }

    my @writeCols = qw(fbgn badPresent badAbsent notPresent notAbsent ratio p q);
    # ratio = (badP/(badP+badA)) / (notP/(notP+notA))
    my $header = join "\t", @writeCols;
    my $format = join "\t", qw(%s %d %d %d %d %.4e %s %s);

    my @badFiles = qw( dpim2_all.120123.plusUniqMsInst.cp.badYear.commContam  dpim2_nrtap.120123.badYear.commContam );
    my @notFiles = qw( dpim2_all.120123.plusUniqMsInst.cp.notBadYear.commContam  dpim2_nrtap.120123.notBadYear.commContam );
    my @outFiles = qw( dpim2_all.120123.plusUniqMsInst.cp.badYear.commContam.fisherTests2 dpim2_nrtap.120123.notBadYear.commContam.fisherTests2 );
    my $it = each_array(@badFiles, @notFiles, @outFiles);

    while (my ($badF, $notF, $outF) = $it->()) {
	say $badF;
	my $badTot = totalExpts($badF);
	my %bad = appearanceHash($badF);
	my $notTot = totalExpts($notF);
	my %not = appearanceHash($notF);
	
	my @tests; # test[i] = map { $_ => ~ } writeCols
	my @absentFromNot; # all proteins in bad but not 'not-bad'
	my @absentFromBad; # contrariwise
	
	for my $fbgn (keys %bad) {
	    if (! exists $not{$fbgn}) {
		push @absentFromNot, $fbgn;
		next;
	    }
	    my $b = $bad{$fbgn};   # number of experiments including $fbgn
	    my $n = $not{$fbgn};
	    my $bX = $badTot - $b; # number of experiments not including $fbgn
	    my $nX = $notTot - $n; 
	    my $p = fisherTest($b, $bX, $n, $nX);

	    my $ratio = ($b/$badTot) / ($n/$notTot);

	    push @tests, { fbgn => $fbgn, badPresent => $b, badAbsent => $bX,
			   notPresent => $n, notAbsent => $nX, ratio => $ratio};
	}
	fisherTests(\@tests);

	@absentFromBad = grep { ! exists $bad{$_} } keys %not;

	@tests = sort {$a->{ratio} <=> $b->{ratio}} @tests;
	my @data;
	for my $col (@writeCols) {
	    push @data, [ map { $_->{$col} } @tests ];
	}

	my $preComments = "# the following were found in $badF but not $notF: ".
	    (join ", ", sort @absentFromNot)."\n";
	$preComments.= "# the following were found in $notF but not $badF: ".
	    (join ", ", sort @absentFromBad)."\n";
	$preComments.= "# q values derived using Benjamini-Hochberg correction";

	say "\twriting to $outF";
	writeCols($outF, \@data, $header, $format, $preComments);
    }

    exit;
}

{
    my @x = map { "cut$_" } 0..9;
    s/cut// for @x;
    say for @x;
    exit;
}

{
    # select all runs from the year of pain
    my $in =  '/home/glocke/DPiM/dpim3.1/dpim2_nrtap.120123';
    my $out = "$in.badYear";
    my $out2 = "$in.notBadYear";

    die "DateTime not working";
    my ($startDate, $endDate, $parseDate);
    #my $startDate = DateTime->new(
    #year => 2009,
#	month => 6,
#	day => 1);
 #   my $endDate = DateTime->new(
#	year => 2010,
#	month => 6,
#	day => 1);

 #   my $parseDate = DateTime::Format::Strptime->new(
  #      pattern   => '%Y-%m-%d');
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't write to $out. $!";
    open my $OUT2, ">", $out2 or die "Can't write to $out2. $!";
    my %row;
    while (getLineAPMS(\%row, $IN, 'line')) {
	my $date = $parseDate->parse_datetime($row{sample_date});
	#if (DateTime->compare($date, $startDate) >= 0 &&
	 #   DateTime->compare($date, $endDate) < 0) {
	  #  print $OUT $row{line};
	#} else {
	 #   print $OUT2 $row{line};	    
	#}
    }

    exit;
}

{
    # select all runs from the year of pain from raw data
    my $in = '/home/glocke/DPiM/dpim3_raw/dpim2_all.120123.plusUniqMsInst.cp';
    my $out = "$in.badYear";
    my $out2 = "$in.notBadYear";

    die "DateTime not working";
    my ($startDate, $endDate, $parseDate);
    #my $startDate = DateTime->new(
#	year => 2009,
#	month => 6,
#	day => 1);
 #   my $endDate = DateTime->new(
#	year => 2010,
#	month => 6,
#	day => 1);
#
 #   my $parseDate = DateTime::Format::Strptime->new(
  #      pattern   => '%Y-%m-%d');
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT,  ">", $out  or die "Can't write to $out. $!";
    open my $OUT2, ">", $out2 or die "Can't write to $out2. $!";
    my %row;
    while (getLineRawAPMS(\%row, $IN, 'line')) {
#	my $date = $parseDate->parse_datetime($row{sample_date});
#	if (DateTime->compare($date, $startDate) >= 0 &&
#	    DateTime->compare($date, $endDate) < 0) {
#	    print $OUT $row{line};
#	} else {
#	    print $OUT2 $row{line};	    
#	}
    }

    exit;
}

{
    # find uniq names in dpim2_raw
    #my $dp2_raw = '/home/glocke/DPiM/dpim3_raw/dpim2_all.120123.plusUniqMsInst.cp';
    my $dp3 = '/home/glocke/DPiM/dpim3.1/DPiM_Data_Summary_2014_edited.reformatted.rmDupes.tab';
    open my $IN, "<", $dp3 or die "can't read from $dp3. $!";
    my (%names, %row);
    while (getLineRawAPMS(\%row, $IN)) {
	$names{$row{user}} = 1;
    }

    die Dumper(\%names);
    
}

{
    # count uniq prey
    my $dp2_clean = '/home/glocke/DPiM/dpim3.1/dpim2_nrtap.120123';
    my $dp2_raw = '/home/glocke/DPiM/dpim3_raw/dpim2_all.120123.plusUniqMsInst.cp';
    my $dp3 = '/home/glocke/DPiM/dpim3.1/DPiM_Data_Summary_2014_edited.reformatted.rmDupes.tab';
    
    my %appear; # appear{$prey}{$experiment} = 1 if prey appeared in experiment
    sub findPrey {
	my ($inFile, $count, $exper) = @_;
	my $reader = \&getLineRawAPMS;
	$reader = \&getLineAPMS if $inFile =~ /nrtap/;
	
	open my $IN, "<", $inFile or die "Can't read from $inFile. $!";
	my %row;
	while ($reader->(\%row, $IN)) {
	    $count->{$row{prey_ref}}{$exper} = 1;
	}
    }

    findPrey($dp2_clean, \%appear, 'dp2_clean');
    findPrey($dp2_raw, \%appear, 'dp2_raw');
    findPrey($dp3, \%appear, 'dp3');
    
    for my $a (values %appear) {
	$a->{dp2_clean_dp3} = 1 if exists $a->{dp2_clean} || exists $a->{dp3};
	$a->{dp2_raw_dp3} = 1 if exists $a->{dp2_raw} || exists $a->{dp3};
    }

    my @experiments = qw(dp2_clean dp2_raw dp3 dp2_clean_dp3 dp2_raw_dp3);
    for my $e (@experiments) {
	my $count = 0+ grep { exists $_->{$e} } values %appear;
	say "$e\t$count";
    }
    exit;
}

{
    # count the number of times a given bait appears
    my $dp2_clean = '/home/glocke/DPiM/dpim3.1/dpim2_nrtap.120123.baitPullDown.tab';
    my $dp2_raw = '/home/glocke/DPiM/dpim3_raw/dpim2_all.120123.plusUniqMsInst.cp.baitPullDown.tab';
    my $dp3 = '/home/glocke/DPiM/dpim3.1/DPiM_Data_Summary_2014_edited.baitPullDown.noGFP.tab';

    my %count; # count{fbgn} = { dp2_clean => number of appearances }
    
    sub countBaits {
	my ($inFile, $count, $exper) = @_;

	my @baits = readCol($inFile, 'bait_ref');
	for my $b (@baits) {
	    $count->{$b}{$exper}++;
	}
    }

    countBaits($dp2_clean, \%count, 'dp2_clean');
    countBaits($dp2_raw, \%count, 'dp2_raw');
    countBaits($dp3, \%count, 'dp3');
    for my $c (values %count) {
	$c->{dp2_clean} //= 0;
	$c->{dp2_raw} //= 0;
	$c->{dp3} //= 0;

	$c->{dp2_clean_dp3} = $c->{dp2_clean} + $c->{dp3};
	$c->{dp2_raw_dp3} = $c->{dp2_raw} + $c->{dp3};
    }
    
    my @experiments = qw(dp2_clean dp2_raw dp3 dp2_clean_dp3 dp2_raw_dp3);
    my %uniqs;
    for my $e (@experiments) {
	$uniqs{$e} = 0+ grep { $_->{$e} > 0 } values %count;
    }
    say Dumper("unique baits", \%uniqs);

    my %experCount;
    for my $e (@experiments) {
	for my $i (1..2) {
	    $experCount{$e}{$i} = 0+ grep { $_->{$e} == $i } values %count;
	}
	$experCount{$e}{'3plus'} = 0+ grep { $_->{$e} > 2 } values %count;
    }

    for my $e (@experiments) {
	print $e;
	for my $size (qw(1 2 3plus)) {
	    print "\t", $experCount{$e}{$size};
	}
	print "\n";
    }


    exit;
}

{
    # find which experiments are in dpim2_clean but not raw
    my $dp2_clean = '/home/glocke/DPiM/dpim3.1/dpim2_nrtap.120123';
    my $dp2_raw = '/home/glocke/DPiM/dpim3_raw/dpim2_all.120123.plusUniqMsInst.cp';
    my (%clean, %raw); # x{sid} = {bait_ref, sample_date}
    sub readSid {
	my ($inFile, $count) = @_;
	my $reader = \&getLineRawAPMS;
	$reader = \&getLineAPMS if $inFile =~ /nrtap/;
	
	open my $IN, "<", $inFile or die "Can't read from $inFile. $!";
	my %row;
	while ($reader->(\%row, $IN)) {
	    $count->{$row{search_id}} //= 
		map { $_ => $row{$_} } qw(bait_ref sample_date);
	}
    }

    readSid($dp2_clean, \%clean);
    readSid($dp2_raw, \%raw);
    for my $sid (keys %clean) {
	say Dumper($clean{$sid}) if ! exists $raw{$sid};
    }
    exit;
}




{
    # find rows where the same prey appears more than once for the same bait
    my $in = '/home/glocke/DPiM/dpim3.0/DPiM_Data_Summary_2014_edited.reformatted.rmDupes.tab';    

    my %uniq; # uniq{sid}{prey_ref}

    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    while (getLineDPIM3(\%row, $IN, 'line')) {
	if (exists $uniq{$row{search_id}}{$row{prey_ref}}) {
	    print $row{line};
	} else {
	    $uniq{$row{search_id}}{$row{prey_ref}} = 1;
	}
    }

    exit;
}

{
    my $in = '/home/glocke/DPiM/dpim3.0/doubles.log';

    my %prey;
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    while (getLineAPMS(\%row, $IN, 'line')) {
	$prey{$row{prey_ref}}++;
    }

    for my $p (sort { $prey{$a} <=> $prey{$b} } keys %prey) {
	say "$p\t$prey{$p}";
    }
    exit;
}

{
    my %h;
    say "exist h{x}{y}? ", isTrue(exists $h{x}{y});
    $h{x}{y} = {a=>1};
    say "exist h{x}{y}? ", isTrue(exists $h{x}{y});
    exit;
}

{
    my $baitCol = 4;
    my $preyCol = 5;

    # print out all lines in AP-MS file where the bait matches $baits
    sub selectBaits {
	my ($in, $out, $baits) = @_;
	
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my %OUT;
	for my $b (keys %$baits) {
	    open $OUT{$b}, ">", "$out.$b" or die "can't write to $out.$b. $!";
	}
	
	while (my $line = <$IN>) {
	    $_ = $line;
	    my @spl = split;
	    next unless @spl > $preyCol;
	    next unless exists $OUT{$spl[$baitCol]};
	    next unless $spl[$preyCol] eq $spl[$baitCol];
	    print {$OUT{$spl[$baitCol]}} $line;
	}
	close $IN;
	close $_ for values %OUT;
    }
    
    # report all unique baits in AP-MS file
    sub enumerateBaits {
	my ($in) = @_;

	my %baits;

	my $lBaitCol = $baitCol+1;

	open my $IN, "<", $in or die "Can't read from $in. $!";
	
	while (my $line = <$IN>) {
	    $_ = $line;
	    my @spl = split;
	    next unless @spl > $lBaitCol;
	    $baits{$spl[$lBaitCol]} = 1;
	}

	return %baits;
    }

    my $logFile = '/home/glocke/DPiM/dpim3.0/tmp.log';
    my %baits = enumerateBaits($logFile);

    my $ccNrbait = '/home/glocke/DPiM/dpim3.0/dpim3.12232014.nrBait';
    my $gpNrbait = '/home/glocke/DPiM/dpim3.0/tmp.nrBait';
    my $redundant = '/home/glocke/DPiM/dpim3.0/dpim3.12232014';

    selectBaits($ccNrbait, "/home/glocke/DPiM/dpim3.0/nrBaitTest/ccRuns", 
		\%baits);
    selectBaits($gpNrbait, "/home/glocke/DPiM/dpim3.0/nrBaitTest/gpRuns", 
		\%baits);
    selectBaits($redundant, "/home/glocke/DPiM/dpim3.0/nrBaitTest/redundantRuns", 
		\%baits);

    exit;
}

{
    my $edgeFile = '/home/glocke/DPiM/analysis/nrDpim2/edgesMissingFromDpim3.edge.GOHist.binomTest.tab';
    my $nodeFile = '/home/glocke/DPiM/analysis/nrDpim2/edgesMissingFromDpim3.node.GOHist.binomTest.tab';

    sub readGO {
	my $file = shift;
	my @cols = qw(term p test);
	my @read = readCols($edgeFile, \@cols);
	return map { $_->{term} => { p => $_->{p}, test => $_->{test} } } @read;
    }
    
    my %edgeGO = readGO($edgeFile);
    my %nodeGO = readGO($nodeFile);
    for my $term (keys %edgeGO) {
	next unless exists $nodeGO{$term};
	next unless $edgeGO{$term}->{p} < 0.1 && $nodeGO{$term}->{p} < 0.1;
	say "$term" if $edgeGO{$term}->{test} ne $nodeGO{$term}->{test};
    }
    say "done";
    exit;
}

{
    my $edgeFile = '/home/glocke/DPiM/analysis/nrDpim2/edge.log';
    my $nodeFile = '/home/glocke/DPiM/analysis/nrDpim2/node.log';
    sub readLog {
	my $in = shift;
	open my $IN, "<", $in or die "Can't read from $in. $!";
	my $ret = {};
	while (my $line = <$IN>) {
	    if ($line =~ m/'(GO:\d+)' => (\d+)/) {
		$ret->{$1} = $2;
	    }
	}
	return $ret;
    }
   
    my $edge = readLog($nodeFile);
    my $node = readLog($edgeFile);

    for my $k (keys %$edge) {
	say "edge->{$k} / node->{$k} = $edge->{$k} / $node->{$k}"
	    if $edge->{$k} != $node->{$k};
    }
    for my $k (keys %$node) {
	say "edge->{$k} / node->{$k} = $edge->{$k} / $node->{$k}"
	    if $edge->{$k} != $node->{$k};
    }


    exit;
}

{
    my $in = '/home/glocke/DPiM/data/dpim2_nrtap.120123.badYear';
    my $out = '/home/glocke/DPiM/data/dpim2_nrtap.120123.badYear.short';
    

    my %uniq; # uniq{sid} = date

    die "DateTime not working";
    my ($startDate, $endDate, $parseDate);
#    my $parseDate = DateTime::Format::Strptime->new(
 #       pattern   => '%Y-%m-%d');
    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    while (getLineAPMS(\%row, $IN, 'line')) {
	$uniq{$row{search_id}} = $parseDate->parse_datetime($row{sample_date});
    }

    
    my @sID = sort { $uniq{$a} <=> $uniq{$b} || $a <=> $b } keys %uniq;
    my @date = map { $parseDate->format_datetime($uniq{$_}) } @sID;

    my $header = join "\t", qw(search_id sample_date);
    my $format = join "\t", qw(%d %s);
    writeCols($out, [\@sID, \@date], $header, $format);

    exit;
}

{
    my $f = '/home/glocke/DPiM/analysis/edgesMissingFromDpim3.goHist.tab';
    open my $IN, "<", $f or die "can't read from $f. $!";
    my $i = 0;
    while (my $line = <$IN>) {
	$line = " $line";
	chomp $line;
	my @spl = quotewords('\s+', 'keep', $line);
	@spl = map { "*$_*" } @spl;
	say join "\t", @spl;
	last if $i > 10;
	$i++;
    }
    exit;
}


{
    ### how many non-unique baits are there in dpim3 data alone?
    my $in = '/home/glocke/DPiM/data/apply_lc.out';
    open my $IN, "<", $in or die "can't read from $in. $!";
    my %row;
    my %sID; # sID{search_id}=bait ref
    while (getLineAPMS(\%row, $IN)) {
	$sID{$row{search_id}} = $row{bait_ref};
    }
    my %baitCnt;
    for my $b (values %sID) {
	$baitCnt{$b}++;
    }
    my @redun = grep { $baitCnt{$_} > 1 } keys %baitCnt;
    for my $bait (sort { $baitCnt{$a} <=> $baitCnt{$b} } @redun) {
	say "$bait\t$baitCnt{$bait}";
    }

    my @nonR = grep { $baitCnt{$_} == 1 } keys %baitCnt;
    say "".(0+ @redun)." redundant baits";
    say "".(0+ @nonR)." non-redundant baits";
    
    exit;
}

{
    ## pull out redundant baits
    my @top10 = qw(FBgn0028980 
FBgn0000084
GFP     
FBgn0011725
FBgn0086134
FBgn0025803
FBgn0031799
FBgn0043070
FBgn0000042
FBgn0011327
FBgn0023175
FBgn0001624
FBgn0262872 );
    # it so happens that the baits with 4 or more runs number 10, plus GFP

    my $in = '/home/glocke/DPiM/data/apply_lc.out';
    my $baseOut = '/home/glocke/DPiM/data/lcRedundant/dpim3.redundant';
    open my $IN, "<", $in or die "can't read from $in. $!";
    my %OUT;
    for my $bait (@top10) {
	my $out = "$baseOut.$bait.apms";
	open $OUT{$bait}, ">", $out or die "can't write to $out. $!";
    }
    my %row;

    while (getLineAPMS(\%row, $IN, 'line')) {
	my $bait = $row{bait_ref};
	if (exists $OUT{$bait}) {
	    print {$OUT{$bait}} $row{line};
	}
    }

    exit;
}

{
    my %h = (
	a => {aa => 1, ab => 2},
	b => {ba => 10, bb => 20},
	c => {ca => 100, cb => 200},
	);

    say max(values %{$h{a}});
    
    
    for my $n1 (sort {max(values %{$h{$b}}) <=> max(values %{$h{$a}})}
		keys %h) 
    {
	for my $n2 (sort {$h{$n1}{$b} <=> $h{$n1}{$a}}
		    keys %{$h{$n1}}) 
	{
	    say "h{$n1}{$n2} = $h{$n1}{$n2}";
	}
    }
    exit;
}

{
    my $in = '/home/glocke/DPiM/data/testingNewCode/sameBait.nrBait';
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my %row;
    my %preyCount;
    while (getLineAPMS(\%row, $IN, 'line') ) {
	$preyCount{$row{prey_ref}}++;
    }
    say "$_\t$preyCount{$_}" for sort {$preyCount{$a} <=> $preyCount{$b}} keys %preyCount;
    exit;
}

{
    # get first 10 prey from first ten baits
    # replace all baits with first bait
    my $in = '/home/glocke/DPiM/data/dpim3.12232014.nrBait';
    my $out = '/home/glocke/DPiM/data/testingNewCode/sameBait.nrBait';

    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't read from $out. $!";
    my %row;
    my %sampleCount;
    my $maxPerSample = 10;
    my $maxSamples = 10;
    my $firstBait;
    while (getLineAPMS(\%row, $IN, 'line') ) {
	$firstBait //= $row{bait_ref};
	my $newBait = $row{bait_ref};
	my $sID = $row{search_id};
	$sampleCount{$sID} //= 0;
	if ($newBait ne $firstBait) {
	    $row{line} =~ s/$newBait/$firstBait/g;
	}
	print $OUT $row{line} if $sampleCount{$sID} < $maxPerSample;
	$sampleCount{$sID}++;
	#say "maxSamples = $maxSamples; keys %sampleCount = "
	 #   , (0+ keys %sampleCount), "; sampleCount{$sID} = $sampleCount{$sID}"#
#	    , ";maxPerSample = $maxPerSample";
	last if ($maxSamples <= (keys %sampleCount)) && 
	    $sampleCount{$sID} >= $maxPerSample;
    }
    close $IN;
    close $OUT;
    exit;
    
}

{
    # get first 10 prey from first ten baits
    my $in = '/home/glocke/DPiM/data/dpim3.12232014.nrBait';
    my $out = '/home/glocke/DPiM/data/testingNewCode/super.short.nrBait';

    open my $IN, "<", $in or die "Can't read from $in. $!";
    open my $OUT, ">", $out or die "Can't read from $out. $!";
    my %row;
    my %sampleCount;
    my $maxPerSample = 10;
    my $maxSamples = 10;
    while (getLineAPMS(\%row, $IN, 'line') ) {
	my $sID = $row{search_id};
	$sampleCount{$sID} //= 0;
	print $OUT $row{line} if $sampleCount{$sID} < $maxPerSample;
	$sampleCount{$sID}++;
	#say "maxSamples = $maxSamples; keys %sampleCount = "
	 #   , (0+ keys %sampleCount), "; sampleCount{$sID} = $sampleCount{$sID}"#
#	    , ";maxPerSample = $maxPerSample";
	last if ($maxSamples <= (keys %sampleCount)) && 
	    $sampleCount{$sID} >= $maxPerSample;
    }
    close $IN;
    close $OUT;
    exit;
    
}

{
    exit;
    # grep out edges with a one-hit wonder
    my $oneFile = '/home/glocke/DPiM/analysis/allAPMSProteins.list';
    open my $IN, "<", $oneFile or die "can't read from $oneFile. $!";
    my @ones = <$IN>;
    close $IN;
    chomp @ones;
    @ones = sort @ones;

    say "found ", (0+ @ones), " OHW's";

    my $netFile = '/home/glocke/DPiM/data/dpim3.12232014.nrBait.network';

    my %found;
    {
	my $outFile = '/home/glocke/DPiM/analysis/allAPMSProteins.edges.tab';
	#open my $OUT, ">", $outFile or die "Can't write to $outFile. $!";
	#say $OUT "node1\tnode2\tscore";

	for my $o (@ones) {
	    my @grep = `grep $o $netFile`;
	    $found{$o} = 0+ @grep;
	    #print $OUT @grep;
	} 
	#close $OUT;
    }
    
    my $out2 = '/home/glocke/DPiM/analysis/allAPMSProteins.acceptance.tab';
    open my $OUT, ">", $out2 or die "can't write to $out2. $!";
    say $OUT "protein\tedges";
    for my $o (@ones) {
	say $OUT "$o\t$found{$o}";
    }
    close $OUT;
    exit;
}

{
        my %h = (
	a => {aa => 1, ab => 2},
	b => {ba => 10, bb => 20},
	c => {ca => 100, cb => 200},
	    );
	for my $k (keys %h) {
	    $_ = 0 for values %{$h{$k}};
	}

	say Dumper(\%h);
	exit;
}

{

    sub makeChecker {
	my @protList = @_;

	if (@protList == 0) {
	    return sub { return 1 };
	}

	my %protHash = map {$_ => 1} @protList;

	return sub {
	    my $p = shift;
	    return exists $protHash{$p};
	};
    }

    my $checker = makeChecker(qw(a b c));
    say "checking '$_': ", isTrue($checker->($_)) for qw(a b c d);
    $checker = makeChecker();
    say "checking '$_': ", isTrue($checker->($_)) for qw(a b c d);
    
    exit
}

{
    my @a = qw('string' 'strong' 'strange');
    say $_ for @a;
    map { $_ =~ s/'//g } @a;
    say $_ for @a;
    exit;
}

{
    my %hash = (1 => 2);
    say "Successful is no, unsuccessful is error: ".
	isTrue(exists $hash{owl}{bowling});
    exit;
}


{
    if (glob "*") {
	say "I found something";
    } else {
	say "I found nothing?";
    }
    if (glob "you didn't really find this did you?") {
	say "hopefully didn't find this weird thing";
    } else {
	say "I found what??";
    }

    exit;
}

{
    sub reader {
	my ($IN) = @_;

	return undef if eof($IN);
	
	my $line;
	do {
	    $line = <$IN>;
	} while (!eof($IN) && $line =~ /ak/);
	return undef if $line =~ /ak/;
	return $line;
    }

    my $in = "test.txt";
    open my $IN, "<", $in or die "can't read from $in. $!";
    
    say $_ while $_ = reader($IN);
    
    exit;
}


{
    my $log = '/home/glocke/DPiM/data/tmp.out.r';
    open my $IN, '<', $log or die "Can't read from $log. $!";
    my @in = <$IN>;
    my @hasMap = grep /searchId2IDx/, @in;
    chomp @hasMap;
    my %map;
    for my $line (@hasMap) {
	$line =~ m/\[(\d+)] = (\d+)/ or die "can't parse $line";
	$map{$1} = $2;
    }

    my @k = sort {$map{$a} <=> $map{$b}} keys %map;

    say "map{$_} = $map{$_}" for @k;
    
    my %revMap;
    for my $k (keys %map) {
	my $v = $map{$k};
	if (exists $revMap{$v}) {
	    die "multiple keys for $v";
	}
	$revMap{$v} = $k;
    }

    my @completeTest = min(values %map)..max(values %map);

    for my $v (@completeTest) {
	die "can't find revMap$v" unless exists $revMap{$v};
    }
    say "found keys for all values from ", $completeTest[0], " to "
	, $completeTest[-1];
    
    exit;
}

exit;
{
    my $f = '/home/glocke/DPiM/data/DPiM_Data_Summary_2014_edited.reformatted.out';
    my @read;
    readColsRef(\@read, $f, qw(unique_peptides total_peptides));
    for my $row (@read) {
	#if ($read{unique_peptides
    }
    exit;
}


    # for every protein in every edge, count GO terms
sub old_countByEdge {
    my ($goHist, $in, $goDB) = @_;
    my $MAXCNT = 20;    
    open my $IN, "<", $in or die "Can't read from $in. $!";
    my $protCount = 0;
    my %test;
    while (my $line = <$IN>) {
	while ($line =~ /(FBgn\d+)/g) {
	    $protCount++;
	    my $go = $goDB->{$1};
	    warn ">\tno GO terms for $1\n" unless defined $go;
	    $go //= ["UNKNOWN"];
	    for my $term (@$go) {
		if ($term eq 'GO:0043161') {
		    $test{$1} = 1;
		}
		die "NOT line = $line\n", Dumper($go) if $term =~ /NOT/;
		$goHist->{$term}++;
	    }
	}
	last if $protCount >= $MAXCNT;
    }
    close $IN;

    print Dumper($goHist);
    exit;
    
    # DEBUG open my $OUT, ">", 'totalEdgeLog';
    #say $OUT $_ for sort keys %test;
    
    
    return $protCount;
}

# for every unique protein in the network, count GO terms
sub countByNode {
    my ($goHist, $in, $goDB) = @_;

    my $MAXCNT=20;

    my %proteins;
    
    my $cnt = 0;
    open my $IN, "<", $in or die "Can't read from $in. $!";
    while (my $line = <$IN>) {
	while ($line =~ /(FBgn\d+)/g) {
	    $proteins{$1}++;
	    $cnt++;
	}
	last if $cnt >= $MAXCNT;
    }
    close $IN;

    for my $p (sort keys %proteins) {
	my $go = $goDB->{$p};
	warn ">\tno GO terms for $p\n" unless defined $go;
	$go //= ["UNKNOWN"];
	for my $term (@$go) {
	    if ($term eq 'GO:0043161') {
		## DEBUG say $p;
	    }
	    $goHist->{$term}++;
	}
    }
    print Dumper($goHist);
    exit;
    
    my $protCount = 0+ keys %proteins;
    return $protCount;
}


sub isTrue {
    my $arg = shift;
    return "yes" if $arg;
    return "no";
}



# one to many map
sub readMultiMap {
    my $in = shift;
    my @cols = qw(uniprot entrez );
    my @read;
    readColsRef(\@read, $in, \@cols);
    my %ret;
    for my $row (@read) {
	$ret{$row->{uniprot}}{$row->{entrez}} = 1;
    }
    return %ret;
}

#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(max min);
use HomeBrew::IO qw(checkExist readList readCol readColsRef readColsHashRef);
use DpimLib qw(getLineBiopAPMS networkHashFromEdgeList);

{
    ## finally translate 8467 symbol network to entrez network
    my $transFile = '/home/glocke/DPiM/human/CompPASS/symb2entrez.8467.manual.txt';
    my %trans;
    readColsHashRef(\%trans, $transFile, [qw(symbol entrez)]);

    my $in = "/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.txt";
    my $out = "/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.net";
    
    open my $IN, "<", $in or die "can't read from $in. $!";
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 converted gene symbol to entrez";
    say $OUT join "\t", "protein1 protein2 score";
    while (<$IN>) {
	next if /^Source/;
	chomp;
	my ($s1, $s2, $score) = split;
	$s1 = dateExcel($s1);
	$s2 = dateExcel($s2);
	my ($e1, $e2) = ($trans{$s1} // -1, $trans{$s2} // -1);
	if (! exists $trans{$s1}) {
	    warn "can't translate $s1\n";
	}
	if (! exists $trans{$s2}) {
	    warn "can't translate $s2\n";
	}
	say $OUT join "\t", $e1, $e2, $score;
    }

    close $IN;
    close $OUT;
    
    exit;
}


{
    # take uniprot's output and pull out everything with an entrez id
    # use id's present in the previous BioPlex
    # in cases of a symbol with multiple possible entrez id's,
    #   check to see if any of them are present in the raw data
    # also check for genes not mapped

    my %refSeq;
    my $transFile = '/home/glocke/DPiM/human/nsaf/symbol2entrez.tsv';
    readColsHashRef(\%refSeq, $transFile, [qw(symbol entrez)]);

    my %prev;
    my $prevFile = '/home/glocke/DPiM/human/CompPASS/symb2entrez.5884.txt';
    readColsHashRef(\%prev, $prevFile, [qw(symbol entrez)]);

    my $rawFile = '/home/glocke/DPiM/human/biopStyleData/gpsRes6_biop_05-13-2016.countBaitPrey';
    my %raw;
    readColsHashRef(\%raw, $rawFile, [qw(protein count)]);
    
    my %cols = (symbol => 4, species=> 5, entrez=>6);
    my $in = '/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.uniprot.tab';
    open my $IN, "<", $in or die "can't read from $in. $!";
    <$IN>; # skip header
    my %genes = %prev;
    while (<$IN>) {
	chomp;
	my @spl = split /\t/;
	next unless defined $spl[$cols{entrez}];
	next unless $spl[$cols{species}] eq 'Homo sapiens (Human)';
	
	my $enString = $spl[$cols{entrez}];
	my @en = split /;/, $enString;
	@en = grep {length($_)} @en;
	my @symb = split /\s+/, $spl[$cols{symbol}];
	for my $s (@symb) {
	    next if exists $prev{$s};
	    if (exists $refSeq{$s}) {
		$genes{$s} = $refSeq{$s};
		next;
	    }
	    for my $e (@en) {
		if (exists $genes{$s} && $genes{$s} != $e) {
		    warn "genes{$s} => $genes{$s} and $e\n";
		    if (exists $raw{$e}) {
			$genes{$s} = $e;
			warn "\t$e in 8400, $prev{$s} in 5900"
			    if exists $prev{$s} && $prev{$s} ne $e;
		    } else {
			$genes{$s} = min($genes{$s}, $e);
		    }
		} else {
		    $genes{$s} = $e;
		}
	    }
	}
    }

    my $out = '/home/glocke/DPiM/human/CompPASS/symb2entrez.8467.txt';
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 sent $in to uniprot.org, got these entrez id's back";
    say $OUT join "\t", qw(symbol entrez);
    say $OUT join "\t", $_, $genes{$_} for sort keys %genes;
    close $OUT;
    
    my $in2 = '/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.geneList';
    my %need = map {$_ => 1} readList($in2);
    for my $n (sort keys %need) {
	if (exists $genes{$n}) {
	    #say "yes $n";
	} else {
	    say $n;
	}
	
    }
    exit;
}

{
    # double check that I've gotten them all covered now
    my $in1 = '/home/glocke/DPiM/human/CompPASS/symb2entrez.8467.manual.txt';
    my %genes;
    readColsHashRef(\%genes, $in1, [qw(symbol entrez)]);

    my $missing;
    my $in2 = '/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.geneList';
    my %need = map {$_ => 1} readList($in2);
    for my $n (sort keys %need) {
	if (exists $genes{$n}) {
	    #say "yes $n";
	} else {
	    $missing = 1;
	    say $n;
	}
	
    }
    say "is anything missing? ", isTrue($missing);
    exit;
}

{
    sub fill  {
	my ($genes, $s, $e) = @_;

	if (exists $genes->{$s} && $genes->{$s} != $e) {
	    warn "genes->{$s} => $genes->{$s} and $e\n";
	    return;
	}
	$genes->{$s} = $e;
	return;
    }
    # find symbol to entrez matches in the 5884 file
    my $in = '/home/glocke/DPiM/human/CompPASS/BioPlex.v5884_HCIs.txt';

    my @cols = qw(source target GeneA   GeneB);
    my @read;
    readColsRef(\@read, $in, \@cols, undef, "\t");
    my %genes;
    for my $row (@read) {
	next if $row->{GeneB} eq 'NA';
	fill(\%genes, $row->{source}, $row->{GeneA});
	fill(\%genes, $row->{target}, $row->{GeneB});
    }
    
    my $out = '/home/glocke/DPiM/human/CompPASS/symb2entrez.5884.txt';
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 retrieved gene symbol -> entrez mapping from $in";
    say $OUT join "\t", qw(symbol entrez);
    say $OUT join "\t", $_, $genes{$_} for sort keys %genes;
    close $OUT;

    exit;
}


{
    # make a list of all gene symbols in the bioplex file

    my $in = "/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.txt";
    my $out = "/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.geneList";
    open my $IN, "<", $in or die "can't read from $in. $!";
    my %genes;
    while (<$IN>) {
	next if /^Source/;
	chomp;
	my ($s1, $s2, $score) = split;
	$s1 = dateExcel($s1);
	$s2 = dateExcel($s2);
	$genes{$s1} = 1;
	$genes{$s2} = 1;
    }
    close $IN;
   
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT $_ for sort keys %genes;
    exit;
}

{
    # special case: P01891 = HLA-A on some "as-is" annotation not in ref-seq
    # convert all such edges to entrez id 3015
    # this requires removing duplicate edges.  grumble grumble.
    my $HLAentrez = 3105;
    my $HLAuni = 'P01891';
    
    my $in = '/home/glocke/DPiM/human/CompPASS/BioPlex.v5884_HCIs.txt';
    my $out = '/home/glocke/DPiM/human/CompPASS/BioPlex.v5884_HCIs.net';

    my @cols = qw(GeneA   GeneB pInt UniprotB);
    my @read;
    readColsRef(\@read, $in, \@cols, undef, "\t");
    open my $OUT, ">", $out or die "can't write to $out. $!";
    say $OUT "# $0 pulled edges from $in";
    say $OUT join "\t", qw(protein1 protein2 score);
    my %HLAedges;
    for my $row (@read) {
	my ($p1, $p2, $s, $uB) = map { $row->{$_} } @cols;
	if ($uB eq $HLAuni) {
	    $p2 = $HLAentrez;
	}
	#next unless length($p2);
	if ($p1 == $HLAentrez) {
	    ($p1, $p2) = ($p2, $p1);
	}
	if ($p2 == $HLAentrez) {
	    $HLAedges{$p1} //= $s;
	    $HLAedges{$p1} = max($HLAedges{$p1}, $s);
	} else {
	    say $OUT join "\t", $p1, $p2, $s;
	}
    }
    for my $p1 (keys %HLAedges) {
	say $OUT join "\t", $p1, $HLAentrez, $HLAedges{$p1};
    }
    close $OUT;
    exit;
}

{
    # given
    my $P01891 = join '', 
    qw(MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRF
DSDAASQRMEPRAPWIEQEGPEYWDRNTRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ
MMYGCDVGSDGRFLRGYRQDAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQW
RAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT
WQRDGEDQTQDTELVETRPAGDGTFQKWVAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP
SSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSL
TACKV);
    my $P04439 = join '',
    qw(MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF
DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ
IMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQL
RAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLT
WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEL
SSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSL
TACKV);
    if ($P01891 eq $P04439) {
	say "equal"; 
    } else { say "not"; }
    exit;
}

{
    my $in = '/home/glocke/DPiM/human/CompPASS/NA.5884.log';
    open my $IN, "<", $in or die "can't read from $in. $!";
    my %uni;
    while (<$IN>) {
	my @spl = split;
	my ($ent, $uni) = ($spl[2], $spl[4]);
	die "found something not NA" if $ent ne 'NA';
	$uni{$uni} = 1;
    }
    say for sort keys %uni;
    exit;
}


sub isTrue {
    my $arg = shift;
    return "yes" if $arg;
    return "no";
}

sub dateExcel {
    # excel has helpfully changed SEPT11 to 11-Sep
    # and MARCH3 to 3-Mar
    my $symbol = shift;
    $symbol =~ s/(\d+)-Sep/SEPT$1/;
    $symbol =~ s/(\d+)-Mar/MARCH$1/;
    return $symbol;
}


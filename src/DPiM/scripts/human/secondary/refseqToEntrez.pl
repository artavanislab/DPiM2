#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw(tempfile);
use Statistics::R;
use HomeBrew::IO qw(checkExist readList readFastaRef);
#use DpimLib qw(getLineDP4APMS);

# make lists of NP_xxx and XP_xxx, feed them to biomart, make a table

my %opts = getCommandLineOptions();

{
    my $faList = $opts{in};
    my $out = $opts{out};

    my @faFiles = readList($faList);
    checkExist('f', $_) for @faFiles;

    my %YP = YPEntrez();
    
    my (@np, @xp, @yp);
    for my $f (@faFiles) {
	my @seqs;
	my @names;
	readFastaRef(\@seqs, $f, \@names, 'nocheck');
	for my $name (@names) {
	    my $rsPep = (split /\|/, $name)[3];
	    $rsPep =~ s/\.\d+$//;
	    if ($rsPep =~ /^XP/) {
		push @xp, $rsPep;
	    } elsif ($rsPep =~ /^NP/) {
		push @np, $rsPep;
	    } elsif ($rsPep =~ /^YP/) {
		die "my handmade hash doesn't include $rsPep"
		    unless exists $YP{$rsPep};
		push @yp, $rsPep;
	    } else {
		die "can't parse $f -> $name -> $rsPep";
	    }
	    last if @np >= 10;
	}
    }

    my (@entrez, @symbol, @refseq);
    my $R = Statistics::R->new();
    $R->startR;
    $R->send(qq{library(biomaRt)});
    $R->send(qq{hsap = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")});
    {
	my ($ent, $sym) = biomaRt($R, \@np, 'refseq_peptide');
	push @entrez, @$ent;
	push @symbol, @$sym;
	push @refseq, @np;
    }
    {
	my ($ent, $sym) = biomaRt($R, \@xp, 'refseq_peptide_predicted');
	push @entrez, @$ent;
	push @symbol, @$sym;
	push @refseq, @xp;
    }
    {
	my @ent = map { $YP{$_} } @yp;
	my ($dum, $sym) = biomaRt($R, \@ent, 'entrezgene');
	push @entrez, @ent;
	push @symbol, @$sym;
	push @refseq, @yp;
    }
    $R->stopR();
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

    my $usage = "usage: $0 -in input -out out.table\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "out=s");
    die $usage unless exists $opts{in} && exists $opts{out};

    for my $k (keys %defaults) {
	$opts{$k} //= $defaults{$k};
    }

    checkExist('f', $opts{in});

    return %opts;
}

# biomart not working for these, so I looked them up by hand
sub YPEntrez {

    return (
	YP_003024026 => 4535,
	YP_003024027 => 4536,
	YP_003024028 => 4512,
	YP_003024029 => 4513,
	YP_003024030 => 4509,
	YP_003024031 => 4508,
	YP_003024032 => 4514,
	YP_003024033 => 4537,
	YP_003024034 => 4539,
	YP_003024035 => 4538,
	YP_003024036 => 4540,
	YP_003024037 => 4541,
	YP_003024038 => 4519,
	
	);
}

# convert refseq_peptide(_predicted) to entrez and symbol
sub biomaRt {
    my ($R, $ids, $inputVocab) = @_;

    my $idString = join ",", map { "'$_'" } @$ids;
    say $idString;
    $R->send(qq{theseGenes <- c($idString)});
    my $attr = qq{"$inputVocab", "entrezgene", "external_gene_name"};
    $attr = qq'"entrezgene", "external_gene_name"'
	if $inputVocab eq 'entrezgene';
    $R->send(qq{getBM(attributes=c($attr), filters= "$inputVocab", values=theseGenes, mart=hsap)});
    my $read = $R->read;
    die Dumper($read);
}

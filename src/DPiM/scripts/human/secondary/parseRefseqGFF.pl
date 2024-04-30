#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw(tempfile);
use Statistics::R;
use Bio::SeqIO::genbank;
use HomeBrew::IO qw(checkExist readList writeCols);
#use DpimLib qw(getLineDP4APMS);

# read in gpff files for human proteins
# output two tables
# * entrez symbol refseq_peptide
# * entrez avg_exp_tryptic avg_length

my %opts = getCommandLineOptions();

{
    my $gpffList = $opts{in};
    my $symbolOut = $opts{sym};
    my $lengthOut = $opts{len};

    my @gpffFiles = readList($gpffList);
    checkExist('f', $_) for @gpffFiles;

    my %YP = YPEntrez();
    
    my (@entrez, @symbol, @refSeq);
    my (%sumLength, %sumN);
    for my $f (@gpffFiles) {
	say $f;
	my $gffIO = Bio::SeqIO->new(-file => $f, -format=> 'GenBank');
	while (my $pep = $gffIO->next_seq()) {
	    my $ref = $pep->display_id; # peptide ID
	    my $len = $pep->length();
	    my @features = $pep->get_SeqFeatures;
	    my ($sym, $trz);
	    for my $feat ( @features ) {
		next unless 'CDS' eq $feat->primary_tag();
		if ($feat->has_tag('gene')) {
		    my $sym0 = ($feat->get_tag_values('gene'))[0];
		    die "multiple symbols for display_id = $ref"
			if defined $sym && $sym ne $sym0;
		    $sym = $sym0;
		}
		if ($feat->has_tag('db_xref')) {
		    my @xref = $feat->get_tag_values('db_xref');
		    for my $xr (@xref) {
			$xr =~ /GeneID:(\d+)/ or next;
			my $trz0=$1;
			die "multiple entrez for display_id = $ref"
			    if defined $trz && $trz ne $trz0;
			$trz = $trz0;
		    }
		}
	    }
	    if (! defined $sym) {
		warn "can't find symbol for $ref";
		next;
	    } elsif (! defined $trz) {
		warn "can't find entrez for $ref";
		next;
	    }
	    #say "$ref $trz $sym $len";
	    push @entrez, $trz;
	    push @symbol, $sym;
	    push @refSeq, $ref;
	    $sumLength{$trz} += $len;
	    $sumN{$trz}++;
	}
    }
    my @sortI = sort {$entrez[$a] <=> $entrez[$b]} 0..$#entrez;
    @entrez = @entrez[@sortI];
    @symbol = @symbol[@sortI];
    @refSeq = @refSeq[@sortI];
    
    my $preComments = "# $0";

    {
	my $header = join "\t", qw(entrez symbol refseq_peptide);
	my $format = join "\t", qw(%d %s %s);
	writeCols($symbolOut, [\@entrez, \@symbol, \@refSeq], $header, $format,
		  $preComments);
    }
    {
	my @e = sort { $a <=> $b } keys %sumLength;
	my @avgLen = map { $sumLength{$_} / $sumN{$_} } @e;
	my @dummyTryptic = (-1) x (0+ @e);
	my $header = join "\t", qw(entrez avg_exp_tryptic avg_length);
	my $format = join "\t", qw(%s %d %.3f);
	writeCols($lengthOut, [\@e, \@dummyTryptic, \@avgLen], $header, 
		  $format, $preComments);
    }
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

    my $usage = "usage: $0 -in gpff.list -sym genesymbol.out -len ".
	"lengthTable.out\n";

    my %opts = ();
    GetOptions(\%opts, "in=s", "sym=s", "len=s");
    die $usage unless exists $opts{in} && exists $opts{sym} && 
	exists $opts{len};

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

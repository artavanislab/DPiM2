#!/usr/bin/env perl

use feature ':5.10'; 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max);
use HomeBrew::IO qw(checkExist readList readFastaRef readColsHashRef);
use HomeBrew::Bio qw(hamming);
use List::MoreUtils qw(each_array);
use Text::Levenshtein qw(distance);
use Math::Interpolate qw(linear_interpolate);
use POSIX;
use File::Slurp;
use JSON;
use DpimLib qw(networkHashFromEdgeList);


{
    my $coappearFile = '/home/glocke/DPiM/augRemap/apmsData/notchome/FBgn0022023-FBgn0034237.coappear.json';
    my @slurp = grep { $_ !~ /^#/ } read_file($coappearFile);
    ##
    ##unshift @slurp, "{\n";
    ##push @slurp, "}";
    ##my $json  = join "", @slurp;
    ##print $json;
    ##my $coap = decode_json($json);

    my $sumTogether = 0;
    my $nTogether = 0;
    for my $line (@slurp) {
	chomp $line;
	$line =~ s/^(\d+):/"$1":/;
	$line = "{$line}";
	my $coap = decode_json($line);
	my ($sid) = keys %$coap;
	my @v = values %{ $coap->{$sid} };
	if ($line =~ /238058/) {
	    #say $line;
	    #die Dumper($coap);
	}
	next unless @v > 1;
	$nTogether++;
	$sumTogether+= min(@v);
	say $sid;
    }

    exit;
    
}

{
    my $in = '/home/glocke/DPiM/augRemap/nrBait_yRej_nBrand_11-16-2016/qdir/test1-12-19-2016/dum.json';
    my $json = read_file($in);
    my $read = decode_json( $json );
    say join ", ", keys %$read;
    my @k = keys %{ $read->{ids} };
    my %fb2i = map { $read->{ids}{$_} => $_ } @k;

    my @pairs = ([qw(FBgn0001215 FBgn0259139)],
		 [qw(FBgn0266720 FBgn0267849)],
		 [qw(FBgn0022023 FBgn0034237)]);
    my @pairs2;
    for my $p (@pairs) {
	my ($p1, $p2) = @$p;
	my ($i1, $i2) = map {$fb2i{$_}} @$p;
	say "$p1 - $p2 -> ", $read->{score}{$i1}{$i2} // 
	    $read->{score}{$i2}{$i1} // "NO";
    }
    exit;
}


{
    ## change multi-line quotes into single-line quotes
    
    # pull complex members
    # do not try and 
    my $in = '/home/glocke/DPiM/prevDPIM/dpim1Net/publishedClusters/cell_5871_mmc4.csv';
    my $out = '/home/glocke/DPiM/prevDPIM/dpim1Net/publishedClusters/cell_5871_mmc4.sanitize.csv';
    my $slurp = read_file( $in ) ;
    my @spl = split '"', $slurp;

    open my $OUT, ">", $out or die "Can't write to $out: $!";
    my $inQuotes = undef;
    for my $elem (@spl) {
	if ($inQuotes) {
	    if ($elem =~ /;/) {
		$elem =~ s/\n//g;
	    } else {
		$elem =~ s/\n/;/g;
	    }
	    $elem = qq("$elem");
	    print $OUT $elem;
	    $inQuotes = undef;
	} else {
	    print $OUT $elem;
	    $inQuotes = 1;
	}
    }

    exit;
}

{
    ## select the data from dpim1 era
    my $dp1StatsFile = '/home/glocke/DPiM/augRemap/apmsData/dpim1/augDpim1.dp4.newBait.pepFDR.sumIso.applyLC.trans.statsBySID';
    my %dp1SID;
    readColsHashRef(\%dp1SID, $dp1StatsFile, [qw(search_id search_id)]);

    my $apmsFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.pepFDR.sumIso.applyLC.trans';
    open my $IN, "<", $apmsFile or die "Can't read from $apmsFile. $!";
    while (<$IN>) {
	my @spl = split;
	print if exists $dp1SID{$spl[0]};
    }
    
    exit;
}

{
    ## select all the data that was previously known, i.e. DPIM4
    exit;
    my $baitKeyFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable4.tap.correct1.tsv';
    my %prevBait;
    readColsHashRef(\%prevBait, $baitKeyFile, [qw(search_id prevbait)]);
    my %remove;
    for my $sid (keys %prevBait) {
	$remove{$sid} = 1 if $prevBait{$sid} eq 'FBgn0000000';
    }

    my $apmsFile = '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.pepFDR.sumIso.applyLC.trans';
    open my $IN, "<", $apmsFile or die "Can't read from $apmsFile. $!";
    while (<$IN>) {
	my @spl = split;
	print unless exists $remove{$spl[0]};
    }
    
    exit;
}

{
    my $dmel = 'MNISALLNSYYDLSGEQMNHINEFVARAVLHVVHHYILSVTPSLVLTLCCRSNHTCNFYNKMMSTLFREWGLAPLQIVNV
LRGVPWHPVPGRRHFNVIFTDSFAAFEEIRMEYYSREYNYNEHYFIFLQARDRLLQGEMRLIFDYCWRYRLIHCSIQVQK
SNGDILFYSYYPFGEHGCSDMEPQLINRYNGSMLVEPDLFPRKLRNFFGCPLRCALWDVPPFLTLDEDQEEVLRVNGGYE
GRLLLALAEKMNFTIAVRKVHVNMRDEALEMLRRDEVDLTLGGIRQTVARGMVATSSHNYHQTREVFGVLASSYELSSFD
ILFYPYRLQIWMGILGVVALSALIQLIVGRMLRERMGSRFWLNLELVFVGMPLLECPRSHTARLYCVMLMMYTLIIRTIY
QGLLYHLIRTHQLNRWPQTIESLVQKNFTVVLTPIVQEVLDEIPSVQHMRFRLLEANSELDPLYFLEANHQLRQHVTASA
LDIFIHFNRLSADKVHQRGEQGSGAHFEIVPEDIISMQLTMYLAKHSFLIDQLNEEIMWMRSVGLLSVWSRWELSESYLR
NEQSFQVLGTMELYAIFLMVLVGLIVGLLVFILELVSMRSIYLRKLFT';
    $dmel =~ s/\s//g;
    my $other = 'MNITELLRSYQNLAGEQASHINEYVARSLLHVVHTHIMSMTSSLVLILCCRNQNTCNFYNEMMIILFRSWGIAPLQIVNV
MQGMRGRHIPGRRHFNLIFTDSYAAFAEIEVEAYSRDYAYNEHYYIFLQARDHLMFAEMGRIFEHCWRHQLINCNVQVQR
ANGDILVYTYQPFGAHSCANMTPQLINRYNGSHMLHPQLYPRKLANFFGCPLRIALWHMPPYWFLSDSDPARITGGIEGR
LLQTISQRLNLSIQVRKPPPGRLPVRKRMLQMLHRDEADLTVGAIRPTVAYNQLATSSHNYYQTSIRLAVLRSHYELSSL
DILLYPYPRIIWLGISLVCALSLLLLLAMDRVLASDPGSLWLNVQLIFVGMPMEQIPRFRCKRLYCIMLMVYTLLIRTAY
QGMLYHLIRTHQLNRLPQSIEELVADNYTVISPIGMYNVLSGIPSMEHMRFHTLETETTLMQPLIYIAEHPRVKRQVVAS
VADVFHLFNHMRSYNADVNEERNCNQFEIISQDVINLQMSMYLRKHSFLVDEFNEQIMWMRSVGLLSIWTHWELDDRVSY
IAKPSSVRSLDLLELQIIFAVVLFGFMLSTIIFGLELLSLRFVRLRIWFHR';
    $other =~ s/\s//g;
    
    say isTrue($dmel eq $other);
    say hamming($dmel, $other);
    say distance($dmel, $other);
    exit;
}


{
    ## check how close dsec's src64B is to Dmel's src64B
    my $dsec = 'MGNKCCSKRQDQELALAYPTGGYKKSDYTFGQTHINSSGGGNMGGVLGQKHNNGGSLDSRYTPDPNHRGPLKIGGKGGVD
IIRPRTTPTGVPGVVLKRVVVALYDYKSRDESDLSFMKGDRMEVIDDTESDWWRVVNLTTRQEGLIPLNFVAEERSVNSE
DWFFENVLRKEADKLLLAEENPRGTFLVRPSEHNPNGYSLSVKDWEDGRGYHVKHYRIKPLDNGGYYIATNQTFPSLQAL
VMAYSKNALGLCHILSRPCPKPQPQMWDLGPELRDKYEIPRSEIQLLRKLGRGNFGEVFYGKWRNSIDVAVKTLREGTMS
TAAFLQEAAIMKKFRHNRLVALYAVCSQEEPIYIVQEYMSKGSLLDFLREGDGRYLHFEDLIYIATQVASGMEYLESKQL
IHRDLAARNVLIGENNVAKICDFGLARVIADDEYCPKQGSRFPVKWTAPEAIIYGKFSIKSDVWSYGILLMELFTYGQVP
YPGMHSREVIENIERGFRMPKPTNHYFPDNIYQLLLQCWDAVPEKRPTFEFLNHYFESFSVTSEVPYREVQD';
    $dsec =~ s/\s//g;

    my @dmel;
    my $dmFile = '/home/glocke/DPiM/nsaf/Src64B.AA.fa';
    readFastaRef(\@dmel, $dmFile, undef, 'nocheck');
    for my $dm (@dmel) {
	say "found a match" if $dm eq $dsec;
    }
    
    exit;
}

{
    my @x = 1..10;
    my @y = map { $_ * $_ } @x;
    my @x2 = map { $_/10 } 0..110;
    my @y2 = map { (linear_interpolate($_, \@x, \@y))[0] } @x2;
    say join ",", map { sprintf("%.3f", $_) } @y2;
    
    exit;
}

{
    my $netFile = '/home/glocke/DPiM/dpim4/witInstr/bestNet_02-01-2016.tsv';
    my %net;
    networkHashFromEdgeList(\%net, $netFile, 0, 1, 1);

    my %nodes;
    my $edges = 0;
    for my $k1 (keys %net) {
	$nodes{$k1} = 1;
	for my $k2 (keys %{$net{$k1}}) {
	    $edges++;
	    $nodes{$k2} = 1;
	}
    }

    say $edges;
    say 0+ keys %nodes;
    exit;
}

{
    #my $dp3File = '/home/glocke/DPiM/oldDpim/dpim3.1/dpim3.09-25-2015.nrBait.77.44.network';
    my $dp3File = '/home/glocke/DPiM/dpim4/witInstr/cons1_01-30-2016.tsv';
    my $dp4List = '/home/glocke/DPiM/dpim4/witInstr/consensus1_01-30-2016/qdir/tmp.list';
    # '/home/glocke/DPiM/dpim4/witInstr/percentile50/cons1_01-30-2016_50Percent.network';
    my %dp3;
    networkHashFromEdgeList(\%dp3, $dp3File, 0, 1, 1);

    my $minScore = 200;
    
    my ($dp3Tot, $dp4Tot) = (0, 0);
    for my $k1 (keys %dp3) {
	for my $k2 (keys %{$dp3{$k1}}) {
	    next unless $dp3{$k1}{$k2} > $minScore;
	    $dp3Tot++;
	}
    }

    
    my (@match, @tot);
    for my $dp4File (readList($dp4List)) {
	my $dp4Head = "$dp4File.head";
	system("head -12000 $dp4File > $dp4Head");
	my %dp4;
	networkHashFromEdgeList(\%dp4, $dp4Head, 0, 1, 1);
	
	my $dp4M = 0;
	my $dp4T = 0;
	for my $k1 (keys %dp4) {
	    for my $k2 (keys %{$dp4{$k1}}) {
		next unless $dp4{$k1}{$k2} > $minScore;
		$dp4T++;
		$dp4M++ if ($dp3{$k1}{$k2} //0) > $minScore;
	    }
	}

	push @match, $dp4M;
	push @tot, $dp4T;

	last if @match > 100;
	system("rm $dp4Head");
    }

    my $it = each_array(@match, @tot);
    while (my ($m, $t) = $it->()) {
	say join "\t", $m, $t, $dp3Tot;
    }
    exit;
}


{
    my @jsonFiles = readList('/home/glocke/DPiM/dpim4/witInstr/consensus1_01-30-2016/json.list');
    say "about to read...";
    my %fbgn;
    for my $f (@jsonFiles) {
	my $json = read_file($f);
	say "read $f";
	my $read = decode_json( $json );
    
	my $ids = $read->{ids};
	for my $fb (values %$ids) {
	    $fbgn{$fb} = 1;
	}
    }
    say "nKeys = ", 0+ keys %fbgn;
    say "$_\n" for sort keys %fbgn;
    exit;
}

{
    #my $jsonFile = '/home/glocke/DPiM/dpim4/consTest0_01-11-2016/HGConsensus000001.json';
    my $jsonFile = '/tmp/fileFXfdki';
    my $json = read_file($jsonFile);
    my $decoded = decode_json( $json );
    
    my $score = $decoded->{score};
    my $ids = $decoded->{ids};

    my @k = keys %$ids;
    my %indexByFBgn;
    for my $i (0..($#k-1)) {
	my $k1 = $k[$i];
	my $fb1 = $ids->{$k1};
	for my $j (($i+1)..$#k) {
	    my $k2 = $k[$j];
	    my $fb2 = $ids->{$k2};
	    $indexByFBgn{$fb1}{$fb2} = $indexByFBgn{$fb2}{$fb1} = $score->{$k1}{$k2} // 0;
	}
    }
    exit;
}


{
    
    my $a = [1..40];
    say "min(a) = ", min($a); #no
    say "max(a) = ", max($a); #no
    exit;
}

{
    my $jsonFile = '/home/glocke/DPiM/oldDpim/dpim3.1/test.01-06-2016.json';
    say "about to read...";
    my $json = read_file($jsonFile);
    say "read it...";
    my $read = decode_json( $json );
    say "decoded! hurray! ref = ", ref($read);

    my $ids = $read->{ids};
    my @k = keys %$ids;
    for my $k (@k[0..9]) {
	say "ids->{$k} = $ids->{$k}";
    }
    exit;
}

{
    my @a = qw(a b c);
    my @b = qw(1 2 3);
    say to_json(\@a, {pretty=>1});
    say to_json([\@a, \@b], {pretty=>1});
    
    my %x = ( (1+0) => 'b',
	      (2+0) => 'c',
	      (3+0) => 'd');
    say to_json(\%x);
    
    exit;
}


{
    #pdhyper (double x, double NR, double NB, double n, int log_p)
    sub pdhyper {
	my ($x, $NR, $NB, $n) = @_;
	
	my $sum = 0;
	my $term = 1;
	while ($x > 0 && $term >= DBL_EPSILON * $sum) {
	    $term *= $x * ($NB - $n + $x) / ($n + 1 - $x) / ($NR + 1 - $x);
	    $sum += $term;
	    $x--;
	}
 	
	return $sum;
    }
    
    my ($found, $red, $blue, $pulls) = (9999, 10000, 10000000, 10000);
    say "pdhyper($found, $red, $blue, $pulls) = ". pdhyper($found, $red, $blue, $pulls);
}




sub isTrue {
    my $arg = shift;
    return "yes" if $arg;
    return "no";
}

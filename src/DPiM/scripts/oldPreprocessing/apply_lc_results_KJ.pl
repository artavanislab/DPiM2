#!/usr/bin/perl -w 

sub process_cmd_args();
process_cmd_args();

my %experiments;
my $totalToRemove = 0;


open(INPUT, "<$ARGV[0]") or die "Cannot open $ARGV[0]\n";
while($buf = <INPUT>){
    @tokens = split("\t", $buf);
	chomp $buf;
	next if ($buf eq '');
    if ($buf =~ "^\t") {
        $prey = $tokens[2];
    } else {
        $instr_run_id = $tokens[0];
        $bait = $tokens[1];
        if ($bait ne $prey) {
			my %peptides;
			if (defined $experiments{$instr_run_id}) {
				%peptides = %{$experiments{$instr_run_id}};
			}
            $peptides{$prey} = '';
            $totalToRemove += 1;
            $experiments{$instr_run_id} = \%peptides;
        }
    }
}

open(DATA, "<$ARGV[1]") or die "Cannot open $ARGV[1]\n";
open(OUT, ">apply_lc.removed.out");
my $header = <DATA>; # defined but not used
#print $header;
print "tap_id\tsearch_id\tsample_date\ttotal_peptides\tbait_ref\tpref_ref\n";
while($buf = <DATA>){
    chomp $buf;
    @tokens = split("\t", $buf);
	$instr_run_id = $tokens[1];
    #my $experiment = $tokens[7];
    $peptide = $tokens[8];

    # Check to see if the experiment needs to have peptides removed.
    if(exists $experiments{$instr_run_id} ){
		my %peptides = % {$experiments{$instr_run_id}};
		if(exists $peptides{$peptide}){
			print OUT $buf."\n";
		    next;
		}
		else{
		    #print $buf . "\n";
			print $tokens[0]."\t".$tokens[3]."\t".$tokens[4]."\t".$tokens[5]."\t".$tokens[7]."\t".$tokens[8]."\n";
		}
    } else {
		#print $buf . "\n";
		print $tokens[0]."\t".$tokens[3]."\t".$tokens[4]."\t".$tokens[5]."\t".$tokens[7]."\t".$tokens[8]."\n";
	}
}
close(DATA);
close(OUT);

print STDERR "$totalToRemove to remove\n";

sub process_cmd_args(){
    if($#ARGV != 1){
	die "Usage apply_lc_results.pl [find_sequential_contam.out orignal_dpim_file]\n";
    }
}


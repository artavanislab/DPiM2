#!/usr/bin/env perl
use Bio::Seq;
$seq_obj = Bio::Seq->new(-seq => 'aaaatgggggggggggccccgtt',
			 -alphabet => 'dna' );
print "$seq_obj\n";

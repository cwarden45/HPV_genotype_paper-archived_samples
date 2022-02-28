use warnings;
use strict;
use Bio::SeqIO;

my $in_FA = "PAVE_download_220224.fasta";
my $out_FA = "PAVE_download_220224-reformat.fasta";

#code copied and modified from another lab's project

open(FA,">  $out_FA")||die("Cannot open $out_FA\n");

my $seqio_obj = Bio::SeqIO->new(-file => $in_FA, 
								-format => "fasta" );							
while (my $seq_obj = $seqio_obj->next_seq){
	my $ref_id = $seq_obj->id;
	my $ref_description = $seq_obj->desc;
	my $ref_seq = $seq_obj->seq;
	#print "$ref_id...\n";
	#print "$ref_description...\n";
	my $HPV_type = "PROBLEM";
	if ($ref_description =~ /\((HPV\d+)\)/){
		($HPV_type) = ($ref_description =~ /\((HPV\d+)\)/);
	}elsif ($ref_description =~ /HPV-m16031680A/){
		$HPV_type = "HPV.m16031680A";
	}elsif ($ref_description =~ /HPV-mICB2/){
		$HPV_type = "HPV.mICB2";
	}else{
		print "Write code to parse |$ref_description|\n";
		exit;
	}#end else
	print "|$HPV_type|\n";
	print FA ">$HPV_type\n$ref_seq\n";
}#end while (my $seq_obj = $seqio_obj->next_seq)

close(FA);
		
exit;
#!/usr/bin/perl -s
use warnings;
use strict;
use Bio::SeqIO;

my $seqio = Bio::SeqIO->new(-file => "sample.fa", '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) 
{
my $string = $seq->seq;
my @candidates = ();
my @prunedcandidates;
my $genename = "testgene";
open (my $log, '>', 'runlog.txt');
print $log "siRNA designer run log".localtime()."\n";
close $log;
#The fasta records are now read one at a time and can now be acted on

	#first we need to find siRNA candidate matches and load them into an array
	
	if ($string =~ m/(AA.{19,23}TT)/g) 
		{
		@candidates = ($string =~ m/(AA.{19,23}TT)/g) ;
#		print "$_\n" for @candidates;
		}
	else {print "no matches found"};
	
	#At this point we have identified any matching candidates, which now need to be pruned for stretches of 4 identical nucleotides or N bases.

	foreach (@candidates)
		{
		if (($_ =~ m/AAAA/g) or ($_ =~ m/TTTT/g) or ($_ =~ m/CCCC/g) or ($_ =~ m/GGGG/g) or ($_ =~ m/N/g)) { }
		else 
			{
#			print "candidate siRNA = $_\n";
			push @prunedcandidates, $_;
#			print "$_\n" for @prunedcandidates;
			}
		}
	#Any candidates which do not meet the selection criteria were removed and now we can focus on these to check for cross hybridisation with blast
	my $count = 0;
	foreach (@prunedcandidates)
		{
		#Create a temporary fa file to blast, and add it to the run log
		open(my $fa, '>', 'temp.fa');
		$count++;
		print $fa "\n>$genename\_candidate_siRNA_\#$count\n$_";
		system "blastn -db cho_mrna -word_size 16 -evalue 1000 -query temp.fa > temp.bln";
		
		#Parse the blast output and retain any useful information	

		system "cat runlog.txt temp.fa > runlogtemp.txt";
		system "mv runlogtemp.txt runlog.txt";
		
		}
}


#Should add the ability to take a remote fasta sequence as entry (or GI)
#Next thing to do is to make it write a new fasta containing the candidates, and to record the name of the input mRNA
#Add the ability to check thermostability, GC content, binding to loop regions of mRNA

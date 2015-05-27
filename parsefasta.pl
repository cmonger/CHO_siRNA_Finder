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
		#lazy use of system grep
		system("grep -A 1 \">\" temp.bln > tmp.txt");
		open(fi, "<tmp.txt");
		my @text = <fi>;
		close fi;
		my @titles = ();
		for (my $i=0; $i < scalar(@text); $i++) 
			{
			my $line = $text[$i];
			chomp($line);
			my $tmptitle;
			if ($line =~ />/) 
				{
				$tmptitle .= $line;	
				$line = $text[$i+1];
				chomp($line);
				$tmptitle .= " ".$line;
				push(@titles, $tmptitle);
				}						
       			}
		#array @titles now contains each blast hit for the siRNA candidate, can look for genes now		 
		my @blasthits;
        	foreach (@titles)
                #count number of ()
               		 {
                	if ($_ =~ m/(\(.+\))/g)
                        	{
                        	push (@blasthits, $1) ;
                        	}
                	}
        	my @uniqblasthits= uniq(@blasthits);
#       	print @uniqblasthits,"\n\n\n";

        	if (scalar @uniqblasthits == 1)
                	{
                	print "\n>$genename\_candidate_siRNA_\#$count has no off-target effects\n\n"
                	}

#               print @uniqblasthits,"\n\n\n";
		system "cat runlog.txt temp.fa > runlogtemp.txt";
		system "mv runlogtemp.txt runlog.txt";	

		}	
#		print scalar @titles, "\n";
		system("rm tmp.txt");
		#array @titles now contains each blast hit for the siRNA candidate, can look for genes now
}	
		


#Subroutines

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


#Should add the ability to take a remote fasta sequence as entry (or GI)
#Next thing to do is to make it write a new fasta containing the candidates, and to record the name of the input mRNA
#Add the ability to check thermostability, GC content, binding to loop regions of mRNA


#Ensure that only the coding region is searched (Look for an ORF)
#Score candidates using rational rnai guidelines http://www.protocol-online.org/prot/Protocols/Rules-of-siRNA-design-for-RNA-interference--RNAi--3210.html

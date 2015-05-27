#!/usr/bin/perl -s
use warnings;
#use strict;
use Bio::SeqIO;

#create log file 
open (my $log, '>', 'runlog.txt');
print $log "siRNA designer run log".localtime()."\n";
close $log;

#Begin to parse the multifasta file using bioperl and act on each record indiviually
my $seqio = Bio::SeqIO->new(-file => "sample.fa", '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) 
{
my $string = $seq->seq;
my @candidates = ();
my $candidateinfo;
my @prunedcandidates;
my $genename = $seq->desc; #correctly returns the gene information (not accession)
my $geneid = $seq->display_id; #correctly returns the GI
#The fasta records are now read one at a time and can now be acted on

	#first we need to find siRNA candidate matches and load them into an array, remembering the position of the match
	
	while ($string =~ m/(AA.{19,23}TT)/g) 
		{
		push (@candidates, $1) ;
		$candidateinfo->{$1}->{"start"} = $-[1];
		$candidateinfo->{$1}->{"end"} = $+[1];	
		}
	
	#At this point we have identified any matching candidates, which now need to be pruned for stretches of 4 identical nucleotides or N bases.

	foreach (@candidates)
		{
		if (($_ =~ m/AAAA/g) or ($_ =~ m/TTTT/g) or ($_ =~ m/CCCC/g) or ($_ =~ m/GGGG/g) or ($_ =~ m/N/g)) { }
		else 
			{
#			print "candidate siRNA = $_\n"; #debug line
			push @prunedcandidates, $_;
#			print "$_\n" for @prunedcandidates; #debug line
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
		my $candidatesequence = $_;

		#Parse the blast output and retain any useful information	
		#lazy use of system grep
		system("grep -A 1 \">\" temp.bln > tmp.txt");
		open("fi", "<tmp.txt");
		my @text = <fi>;
		close "fi";
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
#       	print @uniqblasthits,"\n\n\n"; #debug line

        	if (scalar @uniqblasthits == 1)
                	{
                	print $fa "\n$genename\_candidate_siRNA_\#$count has no off-target effects\n\n"
			}
		else {print $fa "\n$genename\_candidate_siRNA_\#$count has off-target effects in ".((scalar @uniqblasthits)-1)." genes\n\n"}
		
		#Now have a list of which genes have no off target effects. Now need to access GC content & thermodynamic stability
		my $nuccount=0;
		my $gccount=0;
		if (scalar @uniqblasthits == 1)
			{
			#compute GC
			while ($candidatesequence =~ m/(.)/g)
				{
				my $nuc = $1;
				chomp $nuc;
				$nuccount++;
				if (($nuc eq "G") or ($nuc eq "C")) {$gccount++;}
				}
			$candidatesequence->{"$_"}->{"gc"} = ((100 / $nuccount) * $gccount);
			}
		#The candidate sequence is now read into a data structure with its gc content and offtargets (already filtered so all should be 0)
		foreach my $keys (keys %$candidatesequence)
			{
#			print $candidatesequence->{$keys}->{"gc"}."\n"; #example syntax line
			$candidatesequence->{$keys}->{"offtargets"}= ((scalar @uniqblasthits)-1);
			}

		#At this point we finish any analysis on the current gene and log any data before starting on the next
#               print @uniqblasthits,"\n\n\n"; #debug line
		system "cat runlog.txt temp.fa > runlogtemp.txt";
		system "mv runlogtemp.txt runlog.txt";	

		}	
		
		#Remove any temporary files
#		print scalar @titles, "\n"; #debug line
		system("rm tmp.txt");
		system("rm temp.fa");
#		system("rm temp.bln");
}	
		


#Subroutines

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


#Should add the ability to take a remote fasta sequence as entry (or GI)
#Add the ability to check thermostability, GC content, binding to loop regions of mRNA


#Ensure that only the coding region is searched (Look for an ORF)
#Score candidates using rational rnai guidelines http://www.protocol-online.org/prot/Protocols/Rules-of-siRNA-design-for-RNA-interference--RNAi--3210.html

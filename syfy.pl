#!/usr/bin/perl -s
#use warnings;
#use strict;
use Bio::SeqIO;

my $filename = shift || 'sample.fa';
if (-e $filename) {}
else {print "\n File $filename not found, please check file exists\n"; exit;}

my $outfilename = shift || $filename.'_syfyresults.txt';
open ($syfyout, '>', $outfilename);


if (defined $hits) {$hits = $hits - 1}
else {$hits = 2;}

if (($help == 1) or ($h == 1)) 
	{
	print "\nCHO siRNA finder\n\nUsage: \"./syfy.pl <arguments> <inputfile>\" \n\nOptions: -hits=<number> (return X results (Default: 3)) -long (Use all motifs (Default: AA[N19]TT)) -advanced (enable refined scoring) -multilength  -help / -h (This helpful message)\n\n";
	print "Long mode enables searching for the motifs N2[CG]N8[AT]N8[AT]N2 and NA[N21] (more candidates but takes long)\n";
	print "Multilength enables searching for siRNA of length 23-27 (default 23 only)\n";
	print "Advanced enables additional scoring options\n\n";
	print "Input file must be a fasta file, can be a single or multi fasta.\n\n";
 	exit;
	}
print "\nCHO siRNA finder\n\nUsage: \"./syfy.pl <arguments> <inputfile>\" -h or -help for help/arguments\n\n";

#create log file 
open (my $log, '>', 'runlog.txt');
print $log "siRNA designer run log".localtime()."\n";
close $log;

#Begin to parse the multifasta file using bioperl and act on each record indiviually
my $seqio = Bio::SeqIO->new(-file => "$filename", '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) 
{
my $string = $seq->seq;
my @candidates = ();
my $candidateinfo;
my @prunedcandidates;
my $candidatesequence;
my $genename = $seq->desc; #correctly returns the gene information (not accession)
my $geneid = $seq->display_id; #correctly returns the GI
#The fasta records are now read one at a time and can now be acted on
	
	print "Finding siRNA candidates for $genename\n\n";

	#first we need to find siRNA candidate matches and load them into an array, remembering the position of the match
	if ($multilength == 1) {
        if ($string !~ m/(AA.{19,23}TT)/g) {print "No candidates found for $genename";}  
        while ($string =~ m/(AA.{19,23}TT)/g)
                {
                push (@candidates, $1) ;
                $candidateinfo->{$1}->{"start"} = $-[1];
                $candidateinfo->{$1}->{"end"} = $+[1];
                } }
	if ($string !~ m/(AA.{19}TT)/g) {print "No candidates found for $genename";}	 
	while ($string =~ m/(AA.{19}TT)/g) 
		{
		push (@candidates, $1) ;
		$candidateinfo->{$1}->{"start"} = $-[1];
		$candidateinfo->{$1}->{"end"} = $+[1];	
		}
	if ($long == 1) {
	while ($string =~ m/(..[CG].{8}[AT].{8}[AUT]..)/g)
                {
                push (@candidates, $1) ;
                $candidateinfo->{$1}->{"start"} = $-[1];
                $candidateinfo->{$1}->{"end"} = $+[1];
                }
        while ($string =~ m/(.A.{21})/g)
                { 
                push (@candidates, $1) ;
                $candidateinfo->{$1}->{"start"} = $-[1];
                $candidateinfo->{$1}->{"end"} = $+[1];
                }}

	
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
	#opening a file handle for later use in scoring
	open (my $sc, '>tmpscore.txt');

	#Any candidates which do not meet the selection criteria were removed and now we can focus on these to check for cross hybridisation with blast
	my $count = 0;
	foreach (@prunedcandidates)
		{
		#Create a temporary fa file to blast, and add it to the run log
		open(my $fa, '>', 'temp.fa');
		$count++;
		print $fa "\n>$genename\_candidate_siRNA_\#$count\n$_";

		system "blastn -db cho_mrna -word_size 16 -evalue 1000 -query temp.fa > temp.bln";
		$candidatesequence = $_;
		print $fa "\nCandidate match at position ".$candidateinfo->{$candidatesequence}->{"start"}." to ".$candidateinfo->{$candidatesequence}->{"end"};
		

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
                	print $fa "\n$genename\_candidate_siRNA_\#$count has no off-target effects\n"
			}
		else {print $fa "\n$genename\_candidate_siRNA_\#$count has off-target effects in ".((scalar @uniqblasthits)-1)." genes\n"}
		
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
#			print $candidateinfo->{$keys}->{"start"}."\n"; #as the 2 candidate data structures share key hashes, the start and end position can be returned from the original one
			$candidatesequence->{$keys}->{"offtargets"}= ((scalar @uniqblasthits)-1);
			
			#Calculate the thermodynamic values, we want a high negative value
			my $sense = $keys;
			my $antisense = reverse_complement($keys);
			$candidatesequence->{$keys}->{"stability"}= thermodynamic_stability($sense,$antisense);
			print $fa "The GC content of ",$keys," is ",$candidatesequence->{$keys}->{"gc"},"\n";
			print $fa "The thermodynamic stability of ",$keys," is ",$candidatesequence->{$keys}->{"stability"},"\n"; 
			}
		
		#Now to work out the best targets
		foreach my $keys (keys %$candidatesequence)
			{
			#Need to check the GC content, thermodynamic stability, offtarget effects, 100bp outside the start/end
			my $score= 0;
		
			#GC
			if (($candidatesequence->{$keys}->{"gc"} > 30) and ($candidatesequence->{$keys}->{"gc"} < 52))
				{
				$score = $score + 5;
				}
			else {$score = $score - 5};
			#Stability
			if ($candidatesequence->{$keys}->{"stability"} < -3)
				{
				$score = $score + 5;
				}
                        elsif ($candidatesequence->{$keys}->{"stability"} < -2.5)
                                {
                                $score = $score + 4;
                                }
                        elsif ($candidatesequence->{$keys}->{"stability"} < -2)
                                {
                                $score = $score + 3;
                                }
                        elsif ($candidatesequence->{$keys}->{"stability"} < -1.5)
                                {
                                $score = $score + 2;
                                }
                        elsif ($candidatesequence->{$keys}->{"stability"} < -1)
                                {
                                $score = $score + 1;
                                }
                        elsif ($candidatesequence->{$keys}->{"stability"} < -0.5)
                                {
                                $score = $score + 0.5;
                                }
                        elsif ($candidatesequence->{$keys}->{"stability"} > 1)
                                {
                                $score = $score - 5;
                                }
                        elsif ($candidatesequence->{$keys}->{"stability"} > 0.5)
                                {
                                $score = $score - 3;
                                }
			#UTR
			if ($candidateinfo->{$keys}->{"start"} > 150 and ($candidateinfo->{$keys}->{"end"} < (length($string) - 150)))
				{
				$score = $score + 5;
				}
                        elsif ($candidateinfo->{$keys}->{"start"} > 100 and ($candidateinfo->{$keys}->{"end"} < (length($string) - 100)))
                                {
                                $score = $score + 3;
                                }
                        elsif ($candidateinfo->{$keys}->{"start"} > 50 and ($candidateinfo->{$keys}->{"end"} < (length($string) - 50)))
                                {
                                $score = $score + 1;
                                }
			else {$score = $score - 5};
			
			#Advanced mode "Rational RNAi Design" http://www.protocol-online.org/prot/Protocols/Rules-of-siRNA-design-for-RNA-interference--RNAi--3210.html
			if ($advanced== 1)
				{
				if ($keys =~ /^.{14}[AT]/)
					{
					$score = $score +1;
					}
				if ($keys =~ /^.{15}[AT]/)
					{
					$score = $score	+1;
					}
                                if ($keys =~ /^.{16}[AT]/)
                                        {
                                        $score = $score +1;
                                        }
                                if ($keys =~ /^.{17}[AT]/)
                                        {
                                        $score = $score +1;
                                        }
				if ($keys =~ /^.{18}[AT]/)
					{
					$score = $score +1;
					}
				if ($keys =~ /^.{2}A/)
					{
					$score = $score +1;
					}
				if ($keys =~ /^.{9}U/)
					{
					$score = $score +1;
					}
				if ($keys =~ /^.{18}[AT]/)
					{
					$score = $score -1;
					}
				if ($keys =~ /^.{12}[ACT]/)
					{
					$score = $score -1;
					}
				}
	
			print $fa "$keys has a score of $score\n";
			
			#Now the canidates have been scored, they must be sorted based on this number and returned		
			print $sc "$keys\t$score\n";
			
			}
		
			

		#At this point we finish any analysis on the current gene and log any data before starting on the next
#               print @uniqblasthits,"\n\n\n"; #debug line
		system "cat runlog.txt temp.fa > runlogtemp.txt";
		system "mv runlogtemp.txt runlog.txt";	

		}	
		
		#generate score report for current gene
		system "cat tmpscore.txt | sort -rnk 2 > tmporderedscore.txt";
		#Maybe write to results/log file here (and repeat 2 system lines above)
		my @finalkeyname;
		my @finalkeyscore;
		
		open ("fks", '<tmporderedscore.txt');
		my @fks = <fks>;
		close "fks";
		
		#split the @fks into a name and score
		foreach (@fks)
			{
			(my $tseq, $tscore) = split(/\t/);
			push (@finalkeyname, $tseq);
			push (@finalkeyscore, $tscore) ;	
			}

		#returning the top 3 hits
		print $syfyout "\n$genename\n";	
		if (exists $finalkeyscore[0])
			{
			print $syfyout "sequence\tscore\tstart\tend\n";	
			foreach my $i (0..$hits)
				{
				chomp $finalkeyscore[$i];
				if (exists $finalkeyname[$i]) 
					{
					print $syfyout "$finalkeyname[$i]\t$finalkeyscore[$i]\t".$candidateinfo->{$finalkeyname[$i]}->{"start"}."\t".$candidateinfo->{$finalkeyname[$i]}->{"end"}."\n";
					}				
				}
			}
		else {print "No candidate siRNAs found for $genename, try long or multi mode.\n\n";}			 

		
		#Remove any temporary files
#		print scalar @titles, "\n"; #debug line
		system("rm tmp*");
		if (-e "temp.fa"){system("rm temp.fa")};
		if (-e "temp.bln"){system("rm temp.bln")};
}	
		


#Subroutines

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub thermodynamic_stability {
#Usage: thermodynamic_stability(sense,antisense);
#Returns the thermodynamic difference of the 5' end of the sense and antisense sirna transcript in -kcal/mol
#A more stable 5' sense transcript is desirable (A more negative return from this subfunct)
        my $s = shift;
        my $as = shift;
        my $energy;
        $energy->{"AA"} = 1.1;
        $energy->{"AC"} = 2.4;
        $energy->{"AG"} = 1.9;
        $energy->{"AT"} = 1.1;
        $energy->{"CA"} = 2.2;
        $energy->{"CC"} = 3.3;
        $energy->{"CG"} = 2.2;
        $energy->{"CT"} = 1.9;
        $energy->{"GA"} = 2.7;
        $energy->{"GC"} = 3.8;
        $energy->{"GG"} = 3.3;
        $energy->{"GT"} = 2.4;
        $energy->{"TA"} = 1.4;
        $energy->{"TC"} = 2.6;
        $energy->{"TG"} = 2.2;
        $energy->{"TT"} = 1.1;

        if ($s =~ m/(.{6})/) {$s = $1};
        if ($as =~ m/(.{6})/) {$as = $1};

        my $scount = 0;
        my $ascount = 0;

        for ($i = 0; $i < length($s); $i++) {
                if ($i != length($s)-1) {
                        my $nucpair= substr($s,$i,2);
                        $scount = $scount - $energy->{$nucpair};
                }
        }


        for ($i = 0; $i < length($as); $i++) {
                if ($i != length($as)-1) {
                        my $nucpair= substr($as,$i,2);
                        $ascount = $ascount - $energy->{$nucpair};
                 }
        }
        return $scount - $ascount;

}

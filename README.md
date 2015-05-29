# SYFY - A CHO siRNA design tool
SYFY, or si finder, is a tool for the optimal design of candidate siRNAs for the Chinese Hamster.<BR>This tool uses information on the optimal design of siRNAs from a variety of resources including<br>siRNA rational design guidelines (http://www.protocol-online.org/prot/Protocols/Rules-of-siRNA-design-for-RNA-interference--RNAi--3210.html) and thermodynamics calculations (Khvorova et al. 2003 - doi:10.1016/S0092-8674(03)00801-8)<br>in combination with a Blast search to identify potential off-target effects.<br>Default options use the motif AA[N19]TT to identify potential candidates


#Usage
./syfy.pl [options] [inputfasta] [outputfilename]<br>
Input file must be in the form of a fasta file (multifasta support included)<br>
Output filename generated if not specified

##Options
-h/-help	Display help information<br>
-hits=N		Return N number of siRNA candidates for each mRNA (Default = 3) 
-long		Perform a longer search to include sequence motifs N2[CG]N8[AUT]N8[AUT]N2 and NAN21<br>
-advanced	Incorporate additional scoring options<br>
-multilength	Find siRNAs of length 23-27 (Default = 23)<br>

##Disclaimer
siRNA predictions made with this tool are PREDICTIONS and not guaranteed to work. <br>
Author claims no responsibility for unsucessful gene knockdowns.

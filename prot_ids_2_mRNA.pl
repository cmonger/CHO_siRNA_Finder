#!/usr/bin/perl

$infile = shift;
$outfile = shift;

open (ids, $infile) or die;
@ids= <ids>;
close ids;

open ($out, '>', $outfile) or die "Could not open file '$outfile'";

foreach (@ids)
	{
	if ($_ =~ /(\d+)/)
		{
		@output = `lynx -dump http://www.ncbi.nlm.nih.gov/protein/$1 | grep 'http://www.ncbi.nlm.nih.gov/nuccore/'`;
		if ($output[-1] =~ /(\d+)$/)
			{print $out "$1\n";}		
		}
	}

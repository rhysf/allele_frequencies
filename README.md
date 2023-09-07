# allele_frequencies
Calculate allele frequencies from mpileup format. It will output several summary files that can be plotted or examined to help determine ploidy (usually along with any script or tool to calculate and compare depth of coverage).

perl Pileup-allele-frequencies.pl 

Usage: Pileup-allele-frequencies.pl -p <Pileup>

Optional: -m	Min depth to consider mutation [4]
          -r	Reference sequence (FASTA) for order, otherwise from pileup []
	  -i	Restrict to contig/supercontig []

Outputs:  -o	Output file name for percents [opt_p-test-for-polyploid-percents-MD-opt_m.tab]
          -a	Output file name for percents summary [opt_p-test-for-polyploid-percents-summary-MD-opt_m.tab]
	  -b	Output file name for bins [opt_p-test-for-polyploid-bins-MD-opt_m.tab]
	  -s	Output file name for summary [opt_p-test-for-polyploid-summary-MD-opt_m.tab]

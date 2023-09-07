#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use read_Pileup;
use read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: $0 -p <Pileup>\n
Optional: -m\tMin depth to consider mutation [4]
          -r\tReference sequence (FASTA) for order, otherwise from pileup []
	  -i\tRestrict to contig/supercontig []\n
Outputs:  -o\tOutput file name for percents [opt_p-test-for-polyploid-percents-MD-opt_m.tab]
          -a\tOutput file name for percents summary [opt_p-test-for-polyploid-percents-summary-MD-opt_m.tab]
	  -b\tOutput file name for bins [opt_p-test-for-polyploid-bins-MD-opt_m.tab]
	  -s\tOutput file name for summary [opt_p-test-for-polyploid-summary-MD-opt_m.tab]\n";
our($opt_a, $opt_p, $opt_m, $opt_r, $opt_o, $opt_b, $opt_s, $opt_i);
getopt('apmrobsi');
die $usage unless ($opt_p);
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_o) { $opt_o = ($opt_p . '-test-for-polyploid-percents-MD' . $opt_m . '.tab'); }
if(!defined $opt_a) { $opt_a = ($opt_p . '-test-for-polyploid-percents-summary-MD' . $opt_m . '.tab'); }
if(!defined $opt_b) { $opt_b = ($opt_p . '-test-for-polyploid-Bins-MD' . $opt_m . '.tab'); }
if(!defined $opt_s) { $opt_s = ($opt_p . '-test-for-polyploid-summary-MD' . $opt_m . '.tab'); }

# Save order of reference
my $fasta_struct;
if($opt_r) { $fasta_struct = fastafile::fasta_to_struct($opt_r); }

# Define percentile cut-offs
my %percentile_cut_offs;
foreach(qw(0-5 22-28 30-36 47-53 63-69 72-78 95-100)) {
	$percentile_cut_offs{$_}=0; 
}

# Store info
my (%pc_agree, %pc_agree_per_chr, %pc_agree_per_chr_bins, %pc_agree_per_chr_bins2);
my (%length_covered_by_minor_alleles, %length_covered_by_major_alleles);
my %length_covered_by_MD;

# Parse the pileup, saving counts in percentiles
warn "Parsing $opt_p...\n";
open my $fh, '<', $opt_p or die "Can't open $opt_p : $!";
PILEUP: while(my $line=<$fh>) {
	chomp $line;
	my @bits = split /\t/, $line;
	my ($contig, $pos, $ref_base, $depth, $aligned, $qual) = (@bits);

	# Ignore ambiguous bases or unwanted chr
	next PILEUP if($ref_base eq 'N');
	next PILEUP if($depth < $opt_m);
	next PILEUP if((defined $opt_i) && ($opt_i ne $contig));

	# get values
	my ($base_counts) = pileupfile::count_bases($ref_base, $aligned);

	# Summary
	$length_covered_by_minor_alleles{$contig}++; 
	$length_covered_by_MD{$contig}++;

	# Define order if no reference fasta provided
	if(!defined $opt_r) {
		# doesn't work anyway!
		#if(@{$$order} eq 0) { push @{$$order}, $contig; }
		#if(($$order[(scalar(@{$$order} -1))] ne $contig)) { push @{$$order}, $contig; }
		my $count_in_order = 0;
		if(defined $$fasta_struct{'order'}) { 
			$count_in_order = scalar(@{$$fasta_struct{'order'}}); 
			#warn "count in order = $count_in_order\n";
		}
		if($count_in_order eq 0) { 
			warn "$contig...\n";
			push @{$$fasta_struct{'order'}}, $contig; 
		}
		#die "my first thing should now be here: $$fasta_struct{'order'}[$count_in_order] which should match $contig\n";
		if($$fasta_struct{'order'}[($count_in_order - 1)] ne $contig) { 
			#die "$$fasta_struct{'order'}[($count_in_order - 1)] doesn't match $contig... what now?\n";	
			warn "$contig...\n";
			push @{$$fasta_struct{'order'}}, $contig; 
		}
	}
	
	# Percents and bins
	BASES: foreach my $base(keys %{$base_counts}) {
		my $tally = $$base_counts{$base};
		next BASES if($tally eq 0);

		my $percent_read_agree = (($tally / $depth) * 100);
		my $Binned_number = sprintf("%.0f", $percent_read_agree);
		$pc_agree_per_chr{$contig}{$Binned_number}++;
		$pc_agree{$Binned_number}++;

		# Put in bins for haploid, diploid, triploid.
		foreach my $bins(keys %percentile_cut_offs) {
			my @min_max = split /-/, $bins;
			if(($Binned_number >= $min_max[0]) && ($Binned_number <= $min_max[1])) { $pc_agree_per_chr_bins{$contig}{$bins}++; }
		}
	}
}
close $fh;

# Print headers
open my $ofh1, '>', $opt_o or die "Cannot open $opt_o: $!\n";
open my $ofh2, '>', $opt_b or die "Cannot open $opt_b: $!\n";
open my $ofh3, '>', $opt_s or die "Cannot open $opt_s: $!\n";
open my $ofh4, '>', $opt_a or die "Cannot open $opt_a: $!\n";
print $ofh1 "Chr\tPcReadAgree\tTally\n";
print $ofh2 "Chr\tBins\tTally\n";
print $ofh3 "Chr\tLength_Covered\tLength_Covered_Minor_Allele\tLength_Covered_Major_Allele\n";
print $ofh4 "PcReadAgree\tTally\n";

# Print files
#CONTIGS: foreach my $chr(@{$order}) {
CONTIGS: foreach my $chr(@{$$fasta_struct{'order'}}) {
	
	# Summary 
	print $ofh3 "$chr";
	my $length_covered = $length_covered_by_MD{$chr};
	if(!defined $length_covered) { $length_covered = 'no info'; }
	my $length_minor = $length_covered_by_minor_alleles{$chr};
	# default as major for whole contig
	my $length_major = ($length_covered);
	if(!defined $length_minor) {
		$length_minor = 'no info';
		if(!defined $length_covered) { $length_major = 'no info'; }
	} else { $length_major = ($length_covered - $length_minor); }
	print $ofh3 "\t$length_covered\t$length_minor\t$length_major\n";

	# Per contig/chromosome (only opt_p)
	foreach my $pc_agree(sort { $a <=> $b } keys %{$pc_agree_per_chr{$chr}}) {
		my $tally = $pc_agree_per_chr{$chr}{$pc_agree};
		print $ofh1 "$chr\t$pc_agree\t$tally\n";
	}

	# Per bins of haploid, diploid, triploid (only opt_p)
	foreach my $pc_agree_bins(sort keys %{$pc_agree_per_chr_bins{$chr}}) {
		my $tally = $pc_agree_per_chr_bins{$chr}{$pc_agree_bins};
		print $ofh2 "$chr\t$pc_agree_bins\t$tally\n";
	}
}

foreach my $read_agree_bin(sort { $a <=> $b } keys %pc_agree) {
	my $tally = $pc_agree{$read_agree_bin};
	print $ofh4 "$read_agree_bin\t$tally\n";
}

close $ofh1;
close $ofh2;
close $ofh3;
close $ofh4;

package pileupfile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use Data::Dumper;

### rfarrer@broadinstitute.org

sub count_bases {
	my ($ref_base, $aligned) = @_;
	my (%bases);

	# Count insertions (deletions are for subsequent positions)
	while($aligned =~ m/([\+|\-])(\d+)([ACGTNacgtn]+)/ig) {
		my $length_of_indel = $2;
		my $indel = substr $3, 0, $2;
		$indel =~ tr/a-z/A-Z/;
		#if($1 =~ m/\-/) { $bases{('-' . $indel)}++; }
		if($1 !~ m/\-/) { $bases{$indel}++; }
	}

	# Remove insertions and deletions (count '*' when they appear) to stop it being confused with mismatches
	while($aligned =~ m/[\-\+](\d+)([ACGTNacgtn]+)/ig) {
		my $length_of_insertion = $1;
		my $newaligned = substr $aligned, 0, length($`);
		my $offset = (length($`) + length($length_of_insertion) + $length_of_insertion + 1);
		$newaligned .= substr $aligned, $offset, (length($aligned) - $offset);
		$aligned = $newaligned;
	}

	# Start of read segments are given along with Quality scores as ASCII values -33 if -s option is not applied.
	while($aligned =~ m/\^./ig) { 
		my @quality_parts = split //, $&;
		#warn "looking at $aligned which has a $quality_parts[1]\ncomparing to:\n";
		#print Dumper($ascii);
		#my $quality = ($$ascii{$quality_parts[1]} -33);
		$aligned =~ s/\^.//g;
	}
	my $end_of_read_segment = ($aligned =~ tr/$//);

	# Number of matches and mismatches
	$bases{$ref_base} = ($aligned =~ tr/\.|\,//);
	$bases{'A'} += ($aligned =~ tr/A|a//);
	$bases{'C'} += ($aligned =~ tr/C|c//);
	$bases{'T'} += ($aligned =~ tr/T|t//);
	$bases{'G'} += ($aligned =~ tr/G|g//);
	$bases{'N'} += ($aligned =~ tr/N|n//);
	#$bases{'-' . $ref_base} += ($aligned =~ tr/\*//);
	#foreach my $insertion(keys %insertions) { $bases{$insertion} = $insertions{$insertion}; }
	foreach my $base(keys %bases) {
		delete $bases{$base} if($bases{$base} eq 0);
	}
	return \%bases;
}

1;

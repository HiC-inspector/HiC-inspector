#!/usr/bin/perl -w

=head1 NAME

	generate_contact_matrix.pl

=head1 SYNOPSIS
	
  	perl generate_contact_matrix.pl [-i infile] [-o outfile] [-cf chrfile] [-b bin] [-debug] [-h help]

=head1 DESCRIPTION

	Typical usage is as follows:

  	% perl generate_contact_matrix.pl -i interactions.txt -cf chromosome_sizes.txt -o contact_matrix.txt -b 1000000 -debug

=head2 Options

	The following options are accepted:
 
	-i, 	-infile=<file>	Specify an input file containing all identified interactions in the following format: 
				reads name, chr read1, start read1, end read1, score read1, strand read1, chr read2, start read2, end read2, score read2, strand read2 
				(example row is: 5_1_10000_10520 chr12 69090940 69090990 1000 + chr16 67878380 67878430 1000 -)
	-o, 	-outfile=<file>	Specify the output file where to write the contact matrix  
	-b, 	-bin=<num>		Specify the genomic window (alias "bin") to count the chromatin interactions
	-cf, 	-chrfile=<file>	Specify an input file containing chromosome sizes (tab separated text: 1-column is chr, 2-column is length) 
	-debug					Debug while running. Prints out comments during execution
	--help					This documentation

=head1 AUTHORS

	Guglielmo Roma <guglielmo.roma@crg.es>.
	Giancarlo Castellano <giancarlo.castellano@crg.es>.

=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;

my $USAGE = "perl generate_contact_matrix.pl [-i infile] [-o outfile] [-cf chrfile] [-b bin] [-debug] [-h help]\n";
my ($chrfile,$file,$matrixfile,$bin,$debug,$show_help);

&GetOptions(		'infile|i=s'		=> \$file,
			'outfile|o=s'		=> \$matrixfile,
			'chrfile|cf=s'		=> \$chrfile,
			'bin|b=s'		=> \$bin,
			'debug'			=> \$debug,
			'help|h'        	=> \$show_help
	   )
  or pod2usage(-verbose=>2);
pod2usage(-verbose=>2) if $show_help;

die "You must specify an input file with chromosome sizes\n Use -h for help"
  if !$chrfile;

die "You must specify an output file name for the contact matrix\n Use -h for help"
  if !$matrixfile;

die "You must specify an input file containing all identified interactions\n Use -h for help"
  if !$file;

# Set bin to 1 megabase, if not specified 
$bin = 1000000 unless $bin;

# Parse input file and count chromatin interactions
open(FILE, "$file");
my %interactions;
while (my $line = <FILE>) {
	my ($name,$chr1,$start1,$end1,$score1,$strand1,$chr2,$start2,$end2,$score2,$strand2) = split (/ /, $line);	
	# calculate chromosomal bins interacting each other
	my $bin1 = sprintf("%0d", $start1/$bin)+1;
	my $bin2 = sprintf("%0d", $start2/$bin)+1;	
	# add a new interaction for the specific chromosomal bins
	$interactions{$chr1}{$bin1}{$chr2}{$bin2}++;
}

# Parse input file with chromosome sizes and calculate max-bin
open(CHRFILE, $chrfile) || die("Cannot open file"); 
my %maxbin;
while (defined (my $line=<CHRFILE>)) {
	chomp($line);
	my ($chr,$size) = split /\t/, $line;
	next unless $chr;
	my $maxbin = sprintf("%.0f", $size/$bin) + 1;		
	print "chr $chr of size $size\t$maxbin\n";
	$maxbin{$chr}{'maxbin'}=$maxbin;
}
close (CHRFILE);

# Convert interactions into contact matrix
# Set variables to initial values
my $matrix_header = " ";
my $matrix;
my $chr_num = 0;
my $chr_first_col = 1;
my $chr_cols = "";

my %row_count;
my %row_labels;
my $rows = 1;

foreach my $chr (sort keys %maxbin) {
	my $maxbin = $maxbin{$chr}{'maxbin'};
	$debug && print STDOUT "Current chr is $chr has a max bin of $maxbin\n";
	
	# Calculate first and last columns occupied by current chr in the matrix
	my $chr_last_col = $chr_first_col + $maxbin - 1;
	$chr_cols .= "$chr\t$chr_first_col\t$chr_last_col\n";
	$chr_first_col = $chr_last_col + 1;

	# Retrieve all genome-wide interactions for current chr
	my $i;
	for ($i = 1; $i <= $maxbin; $i++) {
		$matrix .= $chr."_".$i;
		$row_labels{$rows} = $chr."_".$i;
		
		my $row_count;
		foreach my $chr2 (sort keys %maxbin) {
			my $maxbin2 = $maxbin{$chr2}{'maxbin'};
			my $k;
			
			for ($k = 1; $k <= $maxbin2; $k++) {
				# Extract interactions for the given bins
				# counting when chr at bin i interacts 
				# with chr2 at bin k, and viceversa
				# Divide the count by 2 on the diagonal
                                my $count = 0;
                                my $count2 = 0;
                                if ($interactions{$chr}{$i}{$chr2}{$k}) {
                                        $count = $interactions{$chr}{$i}{$chr2}{$k};
                                }
                                if ($interactions{$chr2}{$k}{$chr}{$i}) {
                                        $count2 = $interactions{$chr2}{$k}{$chr}{$i};
                                }
				my $final_count;
				if ($chr eq $chr2 && $i == $k) {
					$final_count = ($count + $count2) / 2;
				} else {
					$final_count = $count + $count2;
				}
				# In case of gaps, count, count2, and final_count
				# are not defined. Therefore final count is 0
				$final_count = 0 unless ($final_count);
				
				# Add current chr2 and bin to matrix header
				$matrix_header .= "\t".$chr2."_".$k if ($chr_num == 0);
				# Add final count to matrix
				$matrix .= "\t".$final_count;
				
				# Increment row count
				$row_count += $final_count;
				
				$debug && print STDOUT "$chr $i $chr2 $k $final_count $row_count\n";
			}
		}
		
		$row_count{$rows} = $row_count;
		$rows++;
		
		$matrix .= "\n";
		$chr_num++;
	}
}

# create cvg-corrected matrix
my ($matrix_cvg_corrected,$i,$j);
for ($i = 1; $i < $rows; $i++) {
	$matrix_cvg_corrected .= $row_labels{$i};
	for ($j = 1; $j < $rows; $j++) {
		$matrix_cvg_corrected .= "\t".$row_count{$i} * $row_count{$j};
	}
	$matrix_cvg_corrected .= "\n";
}

# Export matrix
open (MATRIXFILE, ">$matrixfile");
print MATRIXFILE $matrix_header."\n".$matrix;
close (MATRIXFILE);

# Export cvg-corrected matrix
open (CVGMATRIXFILE, ">$matrixfile.cvg");
print CVGMATRIXFILE $matrix_header."\n".$matrix_cvg_corrected;
close (CVGMATRIXFILE);

# Export info on the columns occupied by each chromosome in the matrix
open (CHRINFOFILE, ">$matrixfile.columns.txt");
print CHRINFOFILE $chr_cols;
close (CHRINFOFILE);

# export chromosome names in JSON format
open (CHRJSONFILE,">chrom.json");
print CHRJSONFILE "[\"".join("\", \"", sort keys %maxbin)."\"]";
close (CHRJSONFILE);

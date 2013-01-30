#!/usr/bin/perl -w

=head1 NAME

  hic-inspector.pl

=head1 SYNOPSIS
	
  perl hic-inspector.pl [-n missmatches] [-m multiplemappings] [-cf chrsizefile] [-df designfile] [-sf selectfile] [-dd datadir] [-rd restrictiondir] [-dfo dataformat] [-pd projectdir] [-g genome] [-fs fragmentsize] [-b bin] [-s step] [-t test] [-u utils] [-pr processors] [-h help]

=head1 DESCRIPTION

Typical usage is as follows:

  % perl hic-inspector.pl -df design.txt -dd /data/qseq -pd project -utils hg19

=head2 Options

The following options are accepted:

 -n, 	-missmatches=<int> 		Bowtie option -n: Max mismatches in seed (can be 0-3, default: -n 2)
 -m, 	-multiplemappings=<int> Bowtie option -m: Suppress all alignments if > <int> exist (def: no limit)
 -pd, 	-projectdir=<dir> 		Directory where to write all the results. This folder is created by the pipeline, if it does not exist
 -df, 	-designfile=<file>		Input file describing the experimental design (tab separated text: 1st-column is sample_name, 2nd-column is read1_file, 3rd-column is read2_file, 4th-column is restriction_enzyme_file)
 -dd, 	-datadir=<dir> 			Directory containing data to be analysed. These can be raw reads in qseq or fastq format(compressed or not), or mapped reads in BED format
 -rd, 	-restrictiondir=<dir> 	Directory containing restriction enzyme sites to be considered in the analysis. Should be provided in BED format
 -dfo, 	-dataformat=<dir> 		Format of sequencing reads. Valid options are: qseq (default), fastq, and bed
 -g, 	-genome=<file> 			Indexed genome file for the reads alignment
 -sf, 	-selectfile=<file>		Input file with user-defined genomic regions of interest (BED format)
 -fs, 	-fragment_size=<num>	Maximum expected fragment size
 -cf, 	-chrsizefile=<file>		Input file with chromosome sizes (tab separated text: 1st-column is chr, 2nd-column is size)
 -b, 	-bin=<list>    			Genomic windows, or "bins", to count for chromatin interactions. Several bins can be provided as a comma separated string (e.g. -b 100000,1000000)
 -s, 	-step=<list>   			Analysis steps to be performed. More steps can be provided either as comma-separated list [e.g. 1,2,3] or as dash-separated range [e.g. 1-3]. 
 								  Available steps:
									 STEP 1: copying files to local directory; 
									 STEP 2: converting qseq file to FASTQ format; 
									 STEP 3: mapping reads to genome; 
									 STEP 4: converting mapping output to BED format; 
									 STEP 5: filtering reads by proximity to restriction sites; 
									 STEP 6: filtering reads for regions of interest;
									 STEP 7: combining mate filtered outputs; 
									 STEP 8: calculating distances distribution between mate pairs; 
									 STEP 9: generating contact matrix; 
									 STEP 10: analyzing contact matrix.
 -t, 	-test         			Test mode. Prints out commands without executing them
 -u, 	-utils         			Specify the genome release to be used among those already provided by HiC-pipe (e.g. hg19)
 -pr,	-processors=<int>		Number of processors to be used -for allowing parallelization (default: 2)
 -help              			This documentation.


=head1 AUTHORS

Guglielmo Roma <guglielmo.roma@crg.es>
Giancarlo Castellano <giancarlo.castellano@crg.es>
Toni Hermoso Pulido <toni.hermoso@crg.cat>
=cut

BEGIN {
	# Calculate absolute path
	use Cwd 'abs_path';
	use File::Basename;
	our ($name,$path,$suffix) = fileparse(abs_path($0));
	require $path."./conf.pl";
}

my $path=$path;
my $scriptsdir = $path."/scripts";

use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use POSIX 'strftime';
use Benchmark;
use Parallel::ForkManager;

my $USAGE = "perl hic-inspector.pl [-n missmatches] [-m multiplemappings] [-cf chrsizefile] [-df designfile] [-sf selectfile] [-dd datadir] [-rd restrictiondir] [-dfo dataformat] [-pd projectdir] [-g genome] [-fs fragmentsize] [-b bin] [-s step] [-t test] [-u utils] [-pr processors] [-h help]\n";
my ($missmatches, $multiplemappings, $chrsizefile, $designfile, $selectfile, $restrictiondir, $datadir, $dataformat, $projectdir, $genome, $fragmentsize, $bin, $step, $test, $procs, $utils, $show_help);

&GetOptions(
			'missmatches|n=s'				=> \$missmatches,
			'multiplemappings|mm=s'			=> \$multiplemappings,
			'chrsizefile|cf=s'				=> \$chrsizefile,
			'designfile|df=s'				=> \$designfile,
			'selectfile|sf=s'				=> \$selectfile,
			'datadir|dd=s'					=> \$datadir,
			'restrictiondir|rd=s'			=> \$restrictiondir,
			'dataformat|dfo=s'				=> \$dataformat,
			'projectdir|pd=s'				=> \$projectdir,
			'genome|g=s'					=> \$genome,
			'fragmentsize|fs=s'				=> \$fragmentsize,
			'bin|b=s'						=> \$bin,
			'step|s=s'						=> \$step,
			'utils|u=s'						=> \$utils,
			'processors|pr=s'					=> \$procs,
			'test|t'						=> \$test,
			'help|h'        				=> \$show_help
	   )
  or pod2usage(-verbose=>2);
pod2usage(-verbose=>2) if $show_help;

###########################################
# Check analysis options and steps to run #
###########################################

my %conf =  %::conf;
my $debug = $conf{'debug'};
my $email = $conf{'email'};
my $rdir = $conf{'rdir'};
my $bowtiedir = $conf{'bowtiedir'};
my $bedtoolsdir = $conf{'bedtoolsdir'};

if ($utils) {
	$chrsizefile = $path."/utils/$utils/chromsizes.$utils";
	$genome = $path."/utils/$utils/$utils";
	if (!$restrictiondir) {
		#If restrictiondir is not defined, let's assume it's here
		$restrictiondir = $path."/utils/$utils/";
	}
}

die "You must specify a file describing the experimental design \n Use -h for help" if !$designfile;
die "You must specify the directory containing data to be analysed \n Use -h for help" if !$datadir;
die "You must specify a file with chromosome sizes \n Use -h for help" if !$chrsizefile;
die "You must specify the directory containing the restriction enzyme in BED format \n Use -h for help" if !$restrictiondir;
die "You must specify the path to local installation of R \n Use -h for help" if !$rdir;
die "You must specify the path to local installation of Bowtie \n Use -h for help" if !$bowtiedir;
die "You must specify the path to local installation of BedTools \n Use -h for help" if !$bedtoolsdir;

if ($dataformat) {
	die "Invalid -dataformat argument ($dataformat) specified\n Use -h for help" if ($dataformat ne 'qseq' && $dataformat ne 'fastq' && $dataformat ne 'bed');
} else {
	$dataformat = 'qseq';
}

# If unmapped reads are provided (e.g. qseq or fastq format, but not bed), mapping algorithm needs indexed genome file for the reads alignment
unless ($dataformat eq 'bed') {
	die "You must specify a genome file\n Use -h for help" if !$genome;
	
	# Setting mapping parameters
	$missmatches = "-n $missmatches" if ($missmatches);
	$multiplemappings = "-m $multiplemappings" if ($multiplemappings);
	$missmatches='' unless($missmatches);
	$multiplemappings='' unless($multiplemappings);	
}

# Some default parameters
$fragmentsize = 500 unless $fragmentsize;
$bin = 1000000 unless $bin;
$procs = 2 unless $procs; # 2 parallel processes


# Create an array with all bins
my @bins = split(/,/, $bin);

# Calculate analysis-steps to run
my (%validsteps, $allsteps_valid) = ('', 0);
if ($step) {
	my @step = split (/,/, $step);
    for my $i (0 .. @step-1) {
	    if ($step[$i] =~ m/^(\d+)-(\d+)$/) {
    		for my $nextstep ($1 .. $2) {
    	    	$validsteps{$nextstep}++;
    		}	
    	} elsif ($step[$i] =~ m/^(\d+)$/) {
			$validsteps{$1}++;
   		} else 	{
    		pod2usage ("Invalid -step argument ($step) specified. \n Use -h for help");
    	}
	}
} else {
	$allsteps_valid = 1;
}

# Get current process ID
my $pid = $$;

#######################################
# Let's start with the Hi-C analysis! #
#######################################

# Benchmark running time of current analysis 
my $t0 = Benchmark->new();

# Send email, if email address is provided
my $subject = "HIC-INSPECTOR running with process-id $pid [".strftime('%m-%d-%Y_%H:%M:%S', localtime)."]"; 
$subject .= "- TEST MODE " if $test;
print STDOUT "# $subject \n\n";
&sendEmail ($email, 'hic-inspector', $subject, '') if ($email);

# Preparing directories for storing results
print STDOUT "# Preparing directories for storing results \n\n";
$projectdir = "/$pid" if !$projectdir;
my $filteredbeddir = $projectdir."/filteredbed";
my $interactionsdir = $projectdir."/interactions";

#Configuration HTML files
my $confdir = $projectdir."/conf";
my $confchrom = $confdir."/chrom.json";
my $confsamples = $confdir."/samples.json";


&createdirs($projectdir, $filteredbeddir, $interactionsdir, $confdir);

# Read restdir content
my $resthashdir = &hashdir($restrictiondir);
my %resthashdir = %{$resthashdir};

# Read datadir content
my $datahashdir = &hashdir($datadir);
my %datahashdir = %{$datahashdir};

# First read design file
# Parse design file to store information into hashes
open (FH1, $designfile) or die "Cannot open $designfile: $!";
my (%samples, %mate);
while (<FH1>) {
	chomp;
	unless ($_=~/^\s*$/) {
	my @cas = split("\t", $_);
	my $sample = $cas[0];
	my $read1 = $cas[1];
	my $read2 = $cas[2];
	my $refile = $cas[3];
		
		die "Input files not found in specified datadir ($datadir) \n Use -h for help" if (!exists($datahashdir{$read1}) || !exists($datahashdir{$read2}) || !exists($resthashdir{$refile}));
		$samples{$sample}{$read1}{'enzyme'} = $refile;
		$samples{$sample}{$read2}{'enzyme'} = $refile;
		$samples{$sample}{$read1}{'sample'} = $sample;
		$samples{$sample}{$read2}{'sample'} = $sample;
		$mate{$sample}{$read1} = $read2;
	}
}
close (FH1);

#Iterate for every sample
foreach my $sample (keys %samples) {

	#Results dir with sample information
	my $resultsdir = $projectdir."/results.".$sample;
	my $matrixdir = $projectdir."/results.".$sample."/matrices";
	my $statisticsdir = $projectdir."/results.".$sample."/statistics";
	#Create dirs with sample info
	&createdirs($resultsdir, $matrixdir, $statisticsdir);
	
	my ($qseqdir, $fastqdir, $mapdir, $beddir);
	if ($dataformat eq 'bed') {
		$beddir = $datadir;
	} elsif ($dataformat eq 'qseq' || $dataformat eq 'fastq') {
		if ($dataformat eq 'qseq') {
			$qseqdir = $projectdir."/qseq";
			&createdirs($qseqdir);
		}
		$fastqdir = $projectdir."/fastq";
		$mapdir = $projectdir."/map";
		$beddir = $projectdir."/bed";
		&createdirs($fastqdir, $mapdir, $beddir);
	}
	
	# Max num processes
	my $pm = new Parallel::ForkManager($procs); 
	
	# Process each read-file separately
	foreach my $file (keys %{$samples{$sample}}) {
		
		$pm->start and next; # do the fork
		
		my $restrictionfile = $samples{$sample}{$file}{'enzyme'};
		print STDOUT "# Processing file: $file \n# Restriction file: $restrictionfile \n\n";
		
		# Check if file is compressed
		my $file_gz;
		if ($file =~ /(.+)\.(gz)$/) {
			$file_gz = $file;
			$file = $1;
		}
		
		# Defining filenames
		my ($qseqfile, $fastqfile, $mapfile, $bedfile);
		if ($dataformat eq 'bed') {
			$bedfile = $file;
		} else {
			if ($dataformat eq 'qseq') {
				$qseqfile = $file;
				$fastqfile = "$file.fastq";
			} elsif ($dataformat eq 'fastq') {
				$fastqfile = "$file";
			}
			$mapfile = "$file.map";
			$bedfile = "$file.bed";
		}
		my $filteredbedfile = "$file.filtered.bed";
		my $selectedbedfile = "$file.selected.bed";
		
		if (($validsteps{1} || $allsteps_valid) && $dataformat ne 'bed') {
			my $msg = "# STEP 1 => copying files to project directory:";
	
			my $destinationfile;
			if ($dataformat eq 'qseq') {
				$destinationfile = "$qseqdir/$qseqfile";
			} elsif ($dataformat eq 'fastq') {
				$destinationfile = "$fastqdir/$fastqfile";
			} elsif ($dataformat eq 'bed') {
				$destinationfile = "$beddir/$bedfile";
			}
	
			if ($file_gz) {
				my $cmd = "zcat $datadir/$file_gz > $destinationfile 2>&1";
				&executeCmd ($cmd, $msg, $debug, $test);
			} else {
				my $cmd = "cp $datadir/$file $destinationfile 2>&1";
				&executeCmd ($cmd, $msg, $debug, $test);
			}
		} else {
			print STDOUT "# Skipping STEP 1 => NOT copying files to project directory!\n\n";
		}
		
		if (($validsteps{2} || $allsteps_valid) && $dataformat eq 'qseq') {
			my $msg = "# STEP 2 => converting qseq file to FASTQ format:";
			my $cmd = "awk -f $scriptsdir/awk/qseq2fastq.awk $qseqdir/$qseqfile > $fastqdir/$fastqfile 2>&1";
			&executeCmd ($cmd, $msg, $debug, $test);
		} else {
			print STDOUT "# Skipping STEP 2 => NOT converting qseq file to FASTQ format!\n\n";
		}
		
		if (($validsteps{3} || $allsteps_valid) && $dataformat ne 'bed') {
			my $msg = "# STEP 3 => mapping reads to genome:";
			my $cmd = $bowtiedir."bowtie $missmatches $multiplemappings $genome $fastqdir/$fastqfile $mapdir/$mapfile 2>&1";
			
			# Execute command and save mapping statistics into file
			my $statistics = &executeCmd ($cmd, $msg, $debug, $test);
			if ($statistics) {
				open (FH, ">$statisticsdir/$mapfile");
				print FH $statistics;
				close (FH);
			}
		} else {
			print STDOUT "# Skipping STEP 3 => NOT mapping reads to genome.\n\n";
		}
		
		if (($validsteps{4} || $allsteps_valid) && $dataformat ne 'bed') {
			my $msg = "# STEP 4 => converting mapping output to BED format:";
			my $cmd = "awk -f $scriptsdir/awk/bowtie2bed.awk $mapdir/$mapfile > $beddir/$bedfile 2>&1";
			&executeCmd ($cmd, $msg, $debug, $test);
		} else {
			print STDOUT "# Skipping STEP 4 => NOT converting mapping output to BED format.\n\n";
		}
	
		if ($validsteps{5} || $allsteps_valid) {
			my $msg2 = "# STEP 5 => filtering reads by proximity to restriction sites:";
			my $cmd2 = $bedtoolsdir."windowBed -a $beddir/$bedfile -b $restrictiondir/$restrictionfile -l 0 -r $fragmentsize -sw -u | sort -k 4b,4 | sed -e 's/\\/[1-2]//g' > $filteredbeddir/$filteredbedfile 2>&1";
			&executeCmd ($cmd2, $msg2, $debug, $test);
		} else {
			print STDOUT "# Skipping STEP 5 => NOT filtering reads by proximity to restriction sites.\n";
		}
		
		print STDOUT "###\n\n";
		
		$pm->finish; # do the exit in the child process
	}
	
	$pm->wait_all_children;
	
	# Filter mate reads for regions of interest
	if (($validsteps{6} || $allsteps_valid) && $selectfile) {
		foreach my $read1 (keys %{$mate{$sample}}) {
			my $read2 = $mate{$sample}{$read1};
			my $sample_name = $sample.".selected";
			my $interactions_file = "$sample_name.interactions";
			
			$read1 =~ s/\.gz$//;
			$read2 =~ s/\.gz$//;
			
			my $filteredbedfile1 = "$read1.filtered.bed";
			my $selectedbedfile1 = "$read1.selected.bed";
			
			my $filteredbedfile2 = "$read2.filtered.bed";
			my $selectedbedfile2 = "$read2.selected.bed";
			
			my $msg1 = "# STEP 6 => filtering mate reads for user-specified regions of interest:\n";
			$msg1 .= "# 6.1: selecting read1 for regions of interest";
			my $cmd1 = $bedtoolsdir."intersectBed -a $filteredbeddir/$filteredbedfile1 -b $selectfile -u | sort -k 4 > $filteredbeddir/$selectedbedfile1 2>&1";
			&executeCmd ($cmd1, $msg1, $debug, $test);
			
			my $msg2 = "# 6.2: selecting read2 for regions of interest";
			my $cmd2 = $bedtoolsdir."intersectBed -a $filteredbeddir/$filteredbedfile2 -b $selectfile -u | sort -k 4 > $filteredbeddir/$selectedbedfile2 2>&1";
			&executeCmd ($cmd2, $msg2, $debug, $test);
			
			my $msg3 = "# 6.3: combining read1 by proximity to restriction sites and read2 selected for regions of interest";
			my $cmd3 = "join -j 4 -o 1.4 -o 1.1 -o 1.2 -o 1.3 -o 1.5 -o 1.6 -o 2.1 -o 2.2 -o 2.3 -o 2.5 -o 2.6 $filteredbeddir/$filteredbedfile1 $filteredbeddir/$selectedbedfile2 | sort -k 2b,3  > $filteredbeddir/$filteredbedfile1.$selectedbedfile2.txt 2>&1";
			&executeCmd ($cmd3, $msg3, $debug, $test);
			
			my $msg4 = "# 6.4: combining read2 by proximity to restriction sites and read1 selected for regions of interest";
			my $cmd4 = "join -j 4 -o 1.4 -o 1.1 -o 1.2 -o 1.3 -o 1.5 -o 1.6 -o 2.1 -o 2.2 -o 2.3 -o 2.5 -o 2.6 $filteredbeddir/$filteredbedfile2 $filteredbeddir/$selectedbedfile1 | sort -k 2b,3  > $filteredbeddir/$filteredbedfile2.$selectedbedfile1.txt 2>&1";
			&executeCmd ($cmd4, $msg4, $debug, $test);
			
			my $msg5 = "# 6.5: merge outputs ";
			my $cmd5 = "cat $filteredbeddir/$filteredbedfile1.$selectedbedfile2.txt $filteredbeddir/$filteredbedfile2.$selectedbedfile1.txt | awk -f $scriptsdir/awk/removeduplicates.awk > $interactionsdir/$interactions_file.txt 2>&1";
			&executeCmd ($cmd5, $msg5, $debug, $test);
			
			$mate{$sample}{"$read1.selected"} = "$read2.selected";
			$samples{$sample}{"$read1.selected"}{'sample'} = $sample_name;
			$samples{$sample}{"$read2.selected"}{'sample'} = $sample_name;
		}
	} else {
		print STDOUT "# Skipping STEP 6 => NOT filtering reads for regions of interest.\n";
	}
	 
	# Work on mate files
	foreach my $read1 (keys %{$mate{$sample}}) {
		my $read2 = $mate{$sample}{$read1};
		my $sample_name = $samples{$sample}{$read1}{'sample'};
		$read1 =~ s/\.gz$//;
		$read2 =~ s/\.gz$//;
		my $interactions_file = "$sample_name.interactions";
		my $matrix_file = "$sample_name.matrix";
		
		print STDOUT "# Generating contact matrices for mate pairs: $read1 and $read2 #\n\n";
		
		if (($validsteps{7} || $allsteps_valid) && $sample_name !~ /selected$/) {
			my $msg = "# STEP 7 => combining mate filtered outputs:";
			my $cmd = "join -j 4 -o 1.4 -o 1.1 -o 1.2 -o 1.3 -o 1.5 -o 1.6 -o 2.1 -o 2.2 -o 2.3 -o 2.5 -o 2.6 $filteredbeddir/$read1.filtered.bed $filteredbeddir/$read2.filtered.bed | awk -f $scriptsdir/awk/removeduplicates.awk > $interactionsdir/$interactions_file.txt 2>&1";
			&executeCmd ($cmd, $msg, $debug, $test);
		} else {
			print STDOUT "# Skipping STEP 7 => combining mate filtered outputs.\n\n";
		}
		
		if ($validsteps{8} || $allsteps_valid) {
			my $msg = "# STEP 8 => calculating distances distribution between mate pairs:";
			my $cmd = "awk -f $scriptsdir/awk/computedistances.awk $interactionsdir/$interactions_file.txt > $matrixdir/$interactions_file.distances.txt 2>&1";
			&executeCmd ($cmd, $msg, $debug, $test);
		} else {
			print STDOUT "# Skipping STEP 8 => NOT calculating distances distribution between mate pairs.\n\n";	
		}
		
		if ($validsteps{9} || $validsteps{10} || $allsteps_valid) {

			my $pmb = new Parallel::ForkManager($procs);

			foreach my $bin (@bins) {

				$pmb->start and next; # do the fork

				if ($validsteps{9} || $allsteps_valid) {
					my $msg1 = "# STEP 9 => generating contact matrix:";
					my $cmd1 = "perl $scriptsdir/perl/generate_contact_matrix.pl -i $interactionsdir/$interactions_file.txt -cf $chrsizefile -o $matrixdir/$matrix_file.$bin.txt -b $bin 2>&1";
					&executeCmd ($cmd1, $msg1, $debug, $test);
				} else {
					print STDOUT "# Skipping STEP 9 => NOT generating contact matrix\n\n";
				}
				
				if ($validsteps{10} || $allsteps_valid) {
					my $msg2 = "# STEP 10 => analyzing contact matrix:";
					my $cmd2 = $rdir."Rscript $scriptsdir/R/analyze_contact_matrix.R $matrix_file.$bin.txt $bin $interactions_file.distances.txt $matrix_file.$bin.txt.columns.txt $matrix_file.$bin.txt.cvg $resultsdir 2>&1";
					&executeCmd ($cmd2, $msg2, $debug, $test);
				} else {
					print STDOUT "# Skipping STEP 10 => NOT analyzing contact matrix.\n\n";
				}

				$pmb->finish; # do the exit in the child process
			}

			$pmb->wait_all_children;

		}
		
		print STDOUT "\n###\n\n";
	}

}

#########
# TBD: run final script to generate HTML file for browsing all the pictures and matrices
############
#Create conf files
#Chromosomes
&list_chromosomes($chrsizefile, $confchrom);
#Samples
&list_samples(\%samples, \@bins, $confsamples);

#Copy index.html, CSS & JS
system("cp -rf $scriptsdir/html/* $projectdir");


# Compute analysis time
my $t1 = Benchmark->new;
my $td = timediff ($t1, $t0);

## Send mail, if email address is provided
my $subject2 = "HIC-INSPECTOR with process-id $pid completed after ".timestr($td)." [".strftime('%m-%d-%Y_%H:%M:%S', localtime)."]";
$subject2 .= "- TEST MODE " if $test;
print STDOUT "# $subject2 \n\n";
&sendEmail ($email, 'hic-inspector', $subject2, '') if ($email);

###################
### subroutines ###
###################

sub sendEmail {
 	my ($to, $from, $subject, $message) = @_;
 	my $sendmail = '/usr/lib/sendmail';
 	open (MAIL, "|$sendmail -oi -t");
 	print MAIL "From: $from\n";
 	print MAIL "To: $to\n";
 	print MAIL "Subject: $subject\n\n";
 	print MAIL "$message\n";
 	close (MAIL);
}

sub executeCmd {
 	my ($cmd, $msg, $debug, $test) = @_;
 	
	print STDOUT "$msg\n";
	$debug && print STDOUT "$cmd\n";
	my $out = `$cmd` unless $test;
	print STDOUT "# done!\n\n";
	
	return $out;
}

sub hashdir {
    my ($dir) = @_;
    opendir my $dh, $dir or die $!;
    my %hash;
    while (my $file = readdir($dh)) {
        next if $file =~ m[^\.{1,2}$];
        $hash{$file} = $dir."/".$file;
    }
    
    return \%hash;
}

sub createdirs {
	my (@dirs) = @_;
	foreach my $dir (@dirs) {
		mkdir ($dir, 0755) or die "cannot create $dir" unless (-d $dir);
		print STDOUT "# mkdir $dir\n";
	}
}

sub list_chromosomes {
	
	my $filein = shift;
	my $fileout = shift;
	
	my @chromlist;
	open(FILEIN, $filein) || die "cannot open $filein";
	while (<FILEIN>) {
		
		my ($chrom) = $_=~ /^(\S+)\s/;
		push(@chromlist, "\"".$chrom."\"");
		
	}
	close(FILEIN);
	
	open(FILEOUT, ">$fileout") || die "cannot write";
	
	my $chromstr = "[".join(", ", @chromlist)."]";
	print FILEOUT $chromstr;
	close(FILEOUT);
	
}

sub list_samples {
	
	my $samples = shift;
	my $bins = shift;
	my $fileout = shift;
	
	my $samplestr = "[\n";
	my @binarr = @{$bins};
	my @binstr;
	
	foreach my $bin (@binarr) {
		my $binn = "\"bin.".$bin."\"";
		push(@binstr, $binn);
	}
	
	my $binstr = join(", ", @binstr);
	
	my @numarr = keys %{$samples};
	my $iter =0;
	
	foreach my $sample (keys %{$samples}) {
		
		$iter++;
		$samplestr.= "\t{\n";
		$samplestr.= "\t\t\"id\": \"".$sample."\"\,\n";
		$samplestr.= "\t\t\"value\": [".$binstr."]\n";
		if ($iter == scalar(@numarr)) {$samplestr.= "\t}\n";}
		else { $samplestr.= "\t},\n";}

	}

	$samplestr.="]";
	
	open(FILEOUT, ">$fileout") || die "cannot write";

	print FILEOUT $samplestr;
	close(FILEOUT);
}

# File: qseq2fastq.awk
#?? QSEQ2FASTQ.AWK
#
#?v    Written 2011-03-19
#?v    By      Giancarlo Castellano and Guglielmo Roma
#
# This is an AWK program to convert Illumina QSEQ files to FASTQ format 
#
#? Usage:
#?       AWK -f qseq2fastq.awk file

BEGIN { 
	FS="\t"; 
	OFS=""; 
}

{
	print "@",$3,"_",$4,"_",$5,"_",$6,"\n",$9,"\n+",$3,"_",$4,"_",$5,"_",$6,"\n",$10;
}
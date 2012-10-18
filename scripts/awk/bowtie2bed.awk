# File: bowtie2bed.awk
#?? BOWTIE2BED.AWK
#
#?v    Written 2011-03-19
#?v    By      Giancarlo Castellano and Guglielmo Roma
#
# This is an AWK program to convert bowtie output files into bed format file. 
#
#? Usage:
#?       AWK -f bowtie2bed.awk mapfile

BEGIN { 
	FS="\t"; 
	OFS="\t"; 
}

{
    chr = $3;
    strand = $2;
    start = $4;
    end = start + length($5);
    name = $1;
    score = 1000;
    print chr, start, end, name, score, strand;
}

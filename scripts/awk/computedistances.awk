# File: computedistances.awk
#?? COMPUTEDISTANCES.AWK
#
#?v    Written 2011-03-19
#?v    By      Giancarlo Castellano and Guglielmo Roma
#
# This is an AWK program to convert bowtie output files into bed format file. 
#
#? Usage:
#?       AWK -f computedistances.awk combinedfile

BEGIN { 
	OFS="\t"; 
}

{
	if ($2=$7) print $8-$3
}
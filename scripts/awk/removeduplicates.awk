# File: removeduplicates.awk
#?? REMOVEDUPLICATES.AWK
#
#?v    Written 2011-03-19
#?v    By      Giancarlo Castellano and Guglielmo Roma
#
# This is an AWK program to remove read mates that have the same positions. 
#
#? Usage:
#?       AWK -f removeduplicates.awk

BEGIN { 
	OFS="\t"; 
}

{
	if ($3!=$8 && $4!=$9) print $0
}
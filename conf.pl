BEGIN {
    package main;

    %conf = (
    	# Email address
		'email' 		=> 'email@example.com',
    	# Path to local tools
		'rdir' 			=> "/soft/general/R-2.15.0/bin/",
		'bowtiedir' 	=> "/soft/molbio/bowtie-0.12.7/",
		'bedtoolsdir' 	=> "/soft/molbio/bedtools-2.15.0/bin/",
		# Debug while running 
		'debug' 		=>1 
	);
}

1;

BEGIN {
    package main;

    %conf = (
    	# Email address
		'email' 		=> 'email@example.com',
    	# Path to local tools
		'rdir' 			=> "/usr/bin/",
		'bowtiedir' 	=> "/data/hic-inspector/soft/bowtie-1.0.0/",
		'bedtoolsdir' 	=> "/data/hic-inspector/soft/bedtools-2.17.0/bin/",
		# Debug while running 
		'debug' 		=>1 
	);
}

1;

YACPC: Yet another CLIP-seq peak caller
=========

This peak caller was specifically written to handle low quality datasets.

Works on .bam files as input.

Requires one library file - a .gtf of the genome.

Will do a local Poisson on the CLIP data for the assigned gene.
That is, the signal across the assigned gene is modeled by Poisson.

If given a negative control, will model the negative control as a negative binomial.
This is also modelled for the signal across the assigned gene.

For yeast ribosomal loci, it shifts to a local mode (1 kbp around peak).

Example: :: 

	$ python call_peaks.py -c clip.bam -t genome.gtf

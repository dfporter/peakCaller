YACPC: Yet another CLIP-seq peak caller
=========

This peak caller was specifically written to handle low quality datasets.

Works on .bam files as input.

Requires one library file - a .gtf of the genome.

This code is now deprecated.
A new version has been written from scratch, which fixes some of the problems in this version (poor code organization, little useful documentation, poorly written sections of code, inability to call multiple peaks in long sequences of high coverage, few options for different analysis, ect.), and will form its own repository.

Will do a local Poisson on the CLIP data for the assigned gene.
That is, the signal across the assigned gene is modeled by Poisson.

If given a negative control, will model the negative control as a negative binomial.
This is also modelled for the signal across the assigned gene.

For yeast ribosomal loci, it shifts to a local mode (1 kbp around peak).

Example: :: 

	$ python call_peaks.py -c clip.bam -t genome.gtf

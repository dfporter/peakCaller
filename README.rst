Yet another CLIP-seq peak caller.
=========

Works on .bam files as input.

Will do a local Poisson on the CLIP data.

If given a negative control, will model the negative control as a ZTNB.

Example: :: 

	$ python call_peaks.py -c clip.bam

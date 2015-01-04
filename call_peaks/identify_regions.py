import sys
import os
import glob

def identify_regions(globstr, cutoff=10):
	"""Identifies regions above cutoff in bam files.
	Input argument specifies bam files to work on.
	"""
	print globstr
	for f in glob.glob(globstr):
		print ">%s" % f
	for filename in glob.glob(globstr):
		bash_cmd = """f=%s
""" % filename
		bash_cmd += """
	BAM=$f
	b=$(basename "$f" .bam)
	# Output positive strand.
	samtools view -F 0x10 -b $BAM > ${BAM}.pos
	
	# Extract all regions above cutoff on the positive strand.
	samtools mpileup ${BAM}.pos > pos.mpileup
	rm ${BAM}.pos
	perl extractAllRegionsAboveCutoff.pl pos.mpileup + > ${BAM}.peaks

	# Output negative strand.
	samtools view -f 0x10 -b $BAM  > ${BAM}.neg

	# Extract all regions above cutoff on the negative strand.
	samtools mpileup ${BAM}.neg > neg.mpileup
	rm ${BAM}.neg
	perl extractAllRegionsAboveCutoff.pl neg.mpileup - > neg.peaks

	# Combine files.
	cat neg.peaks >> ${BAM}.peaks

	# Clean up.
	rm neg.peaks
	rm neg.mpileup
	rm pos.mpileup

	# This leaves one output file, $BAM.peaks, which has both strands.

	# Sort by peak height.
	sort -k5nr ${BAM}.peaks > ${BAM}.regions
	rm ${BAM}.peaks
"""
		os.system(bash_cmd)


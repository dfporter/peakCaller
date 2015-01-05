import sys
import os
import glob

def identify_regions(globstr, src_dir, cutoff=1):
	"""Identifies regions above cutoff in bam files.
	Input argument specifies the bam files to work on.
	Outputs .regions files. Returns nothing.
	src_dir is the directory with the required perl script.
	"""
	print "Identifying regions above %i coverage in %s..." % (
                cutoff, globstr)
	perl_script_path = src_dir + '/extractAllRegionsAboveCutoff.pl'
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
	perl %s pos.mpileup + > ${BAM}.peaks

	# Output negative strand.
	samtools view -f 0x10 -b $BAM  > ${BAM}.neg

	# Extract all regions above cutoff on the negative strand.
	samtools mpileup ${BAM}.neg > neg.mpileup
	rm ${BAM}.neg
	perl %s neg.mpileup - > neg.peaks

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
""" % (perl_script_path, perl_script_path)
		os.system(bash_cmd)


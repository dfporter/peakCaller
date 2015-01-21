import sys
import os
import time
import re
import subprocess
import argparse
from peak import peak
import clip_tools
import identify_regions


def parse_args():
    src_dir = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(description=
    """Call peaks in CLIP-seq data.""")
    parser.add_argument('-p', '--peaks', default="create",
                        help="""File of regions of interest.""")
    parser.add_argument('-c', '--clip_reads',
                        help=""".bam file of clip-seq.""")
    parser.add_argument('-b', '--background_reads', default=False,
                        help=""".bam file of background.""")
    parser.add_argument('-a', '--annotation_bed',
                        default="%s/lib/cerevisiae_genes.bed" % (src_dir),
                        help=""".bed file of gene locations.""")
    parser.add_argument('-g', '--gain', type=float, default=1.0, 
                        help="""Coefficient to inflate the background.
Normalizing a CLIP-seq dataset to a much larger (~20-fold) RNA-seq dataset
tends to inflate the CLIP-seq peak out of proportion, so that the ZTNB
statistics fail to eliminate any peaks. Adding a -g 4 option will multiply
the background by 4 vs CLIP-seq when producing ZTNB p values.
Of course, this replaces the statistical intrepetation of the ZTNB
value with the ZTNB p value being merely a measure of relative enrichment.""")
    parser.add_argument('-o', dest='output_folder', default='.',
                        help="""Output directory.""")
    args=parser.parse_args()
    return args

if __name__ == '__main__':
    src_dir = os.path.dirname(os.path.realpath(__file__))
    args = parse_args(src_dir)
    # Variable initilization.
    clip_bam_filename = args.clip_reads
    regions_above_cutoff_filename = args.peaks
    if regions_above_cutoff_filename == "create":
        identify_regions.identify_regions(clip_bam_filename, src_dir)
        regions_above_cutoff_filename = "%s.regions" % clip_bam_filename
    control_bam_filename = args.background_reads
    normalization_coefficient = 1
    peakHeights = dict()
    # Set the normalization coefficient by dividing total read number
    # of CLIP reads by total read number of RNA-seq. Uses only mapped reads.
    if(args.background_reads):
        normalization_coefficient = clip_tools.normalize(
            clip_bam_filename, control_bam_filename)
        normalization_coefficient = args.gain * normalization_coefficient
        sys.stderr.write("Normalization coefficient is %.4f (Gain: %.3f)\n" % (
            normalization_coefficient, float(args.gain)))
    # Count the number of peaks.
    num_peaks = 0
    with open(regions_above_cutoff_filename, 'r') as f:
        for line in f:
            num_peaks += 1
    print "\n\nRegions above cutoff in .regions file: %i" % num_peaks
    # File io initialization.
    peaksBasename = re.match(
        r'(.+)\.regions',
        os.path.basename(regions_above_cutoff_filename)).group(1)
    peaksBasename = re.sub(r'.bam', '', peaksBasename)
    results_folder = args.output_folder + '/' + peaksBasename + '/'
    if(not os.path.exists(results_folder)):
        print "Creating output folder %s..." % results_folder
        os.system("mkdir %s" % results_folder)    
    peaks_ranges_filename = results_folder + 'init_ranges' 
    # Start the clock.
    startTime = time.time()
    #skip = r'''
    # This is the first main loop of the program:
    # For each peak, get the reads mapping in the general region from
    # peaks .bam and RNA-seq .bam and output to init_ranges.
    clip_tools.define_ranges(results_folder,
                             regions_above_cutoff_filename,
                             clip_bam_filename,
                             #control_bam_filename,
                             normalization_coefficient,
                             peaks_ranges_filename,
                             peakHeights, startTime, num_peaks)
    
    clip_tools.remove_duplicate_ranges(
        results_folder + 'init_ranges',
        results_folder + '/ranges_no_dups')
    #'''
    peak_stats = dict()
    print "Processing ranges in %s..." % peaks_ranges_filename
    clip_tools.process_ranges(results_folder,
                              "%s/ranges_no_dups" % results_folder,
                              clip_bam_filename,
                              control_bam_filename,
                              normalization_coefficient,
                              peak_stats)

    # Call R and reformat the results, which are put in a file r.out.
    # callR() returns a dict with key = peak number, value = pvalue.
    clip_tools.call_R(results_folder, peak_stats, src_dir)
    clip_tools.add_p_value("%s/ranges_no_dups" % results_folder,
                           results_folder + 'ranges_with_stats',
                           peak_stats)
    os.system("mv %s/ranges_with_stats %s/ranges" % (results_folder,
                                                       results_folder))
    print "Converting ranges with stats to a .peaks file..."
    clip_tools.ranges_with_stats_to_peaks('ranges', results_folder,
                                          annotation_file=args.annotation_bed)

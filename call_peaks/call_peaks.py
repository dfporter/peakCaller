"""
Calls peaks in .bam files of CLIP-seq data.

This script needs a .gtf for the whole genome (given with -t).
It will create a .bed file of annotations from the .gtf. The
.bed file can also be given with -a, if it already exists.

A control dataset, such as RNA-seq, can also be provided.

"""
import sys
import os
import time
import re
import subprocess
import argparse
import peak
import clip_tools
import identify_regions
import pysam
import annotation_bed_from_gtf

def parse_args():
    src_dir = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(description=
    """Call peaks in CLIP-seq data.""")
    parser.add_argument('-p', '--peaks', default="create",
                        help="""File of regions of interest.""")
    parser.add_argument('-c', '--clip_reads',
                        help="""(Required) .bam file of clip-seq.""")
    parser.add_argument('-b', '--background_reads', default=False,
                        help=""".bam file of background.""")
    parser.add_argument('-a', '--annotation_bed',
                        help=""".bed file of gene locations.""")
    parser.add_argument('-t', '--gtf',
        default="%s/lib/Saccharomyces_cerevisiae.EF4.70.gtf" % src_dir,
        help="(Required) .gtf file for the given genome build.")
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
    args = parse_args()
    
    # Variable initilization.
    if not args.annotation_bed or (not os.path.exists(args.annotation_bed)):
        annotation_bed_from_gtf.create_bed_from_gtf(args.gtf, src_dir + "/lib/tmp_genes.bed")
        args.annotation_bed = src_dir + "/lib/tmp_genes.bed"
        os.system(
            f"bedtools sort -i {args.annotation_bed} > {args.annotation_bed}.tmp " + \
            f"&& mv {args.annotation_bed}.tmp {args.annotation_bed}")
    clip_bam_filename = args.clip_reads
    regions_above_cutoff_filename = args.peaks
    
    if regions_above_cutoff_filename == "create":
        identify_regions.identify_regions(clip_bam_filename, src_dir)
        regions_above_cutoff_filename = "%s.regions" % clip_bam_filename
        
    control_bam_filename = args.background_reads
    normalization_coefficient = 1.0
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
    with open(regions_above_cutoff_filename, 'r') as f:
        for num_peaks, line in enumerate(f):
            pass
    num_peaks += 1
    
    print("\n\nRegions above cutoff in .regions file: %i" % num_peaks)
    
    # File io initialization.
    peaksBasename = re.match(r'(.+)\.regions',
        os.path.basename(regions_above_cutoff_filename)).group(1)
    peaksBasename = re.sub(r'.bam', '', peaksBasename)
    results_folder = args.output_folder + '/' + peaksBasename + '/'
    os.makedirs(results_folder, exist_ok=True)
    peaks_ranges_filename = results_folder + 'init_ranges' 
    
    # Start the clock.
    startTime = time.time()
    
    # This is the first main loop of the program:
    # For each peak, get the reads mapping in the general region from
    # peaks .bam and RNA-seq .bam and output to init_ranges.
#     clip_tools.define_ranges(regions_above_cutoff_filename,
#                              clip_bam_filename,
#                              normalization_coefficient,
#                              peaks_ranges_filename,
#                              peakHeights, startTime, num_peaks)
    
    print(f"Removing duplicate ranges.")
    clip_tools.remove_duplicate_ranges(
        results_folder + 'init_ranges',
        results_folder + '/ranges_no_dups')
    
    print(f"Assigning to genes.")
    clip_tools.assign_to_gene(
        results_folder + '/ranges_no_dups',
        results_folder + '/ranges_with_ann',
        annotation_file=args.annotation_bed)
    
    print("Processing ranges in %s..." % peaks_ranges_filename)
    # binned_genes is returned to correct significances for gene length.
    (peaks, binned_genes) = clip_tools.process_ranges(
        "%s/ranges_with_ann" % results_folder,
        clip_bam_filename,
        control_bam_filename,
        normalization_coefficient,
        gtf_filename=args.gtf,
        annotation_file=args.annotation_bed)
    print("Finished processing ranges. Calling R...")
    
    # Call R and reformat the results, which are put in a file r.out.
    # callR() returns a dict with key = peak number, value = pvalue.
    clip_tools.call_R(results_folder, peaks, src_dir, which='background')
    clip_tools.call_R(results_folder, peaks, src_dir, which='clip')
    #clip_tools.add_p_value("%s/ranges_with_ann.cor_height" % results_folder,
    #                       results_folder + 'ranges_with_stats',
    #                       peak_stats)
    #os.system("mv %s/ranges_with_stats %s/ranges" % (results_folder,
    #                                                   results_folder))
    
    print("Converting ranges with stats to a .peaks file...")
    clip_tools.ranges_with_stats_to_peaks(
        results_folder, binned_genes, peaks,
        annotation_file=args.annotation_bed)
    
    print("Finished with file %s." % clip_bam_filename)

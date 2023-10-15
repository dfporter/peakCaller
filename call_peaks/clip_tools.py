  # int/int = float
import sys
import os
import pybedtools
import HTSeq
from scipy.stats import poisson
import scipy
import numpy as np
import time
import re
import subprocess
import math
import random
import argparse
from .peak import peak
from .gene import gene
from .gtf_data import gtf_data


def normalize(clipReadsFname, background_bam_filename):
    # Get total read number from CLIP-seq .bam
    cmdl = ["samtools", "flagstat", str(clipReadsFname)]
    stats = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
    stp = re.compile(r'\\n(\d+) \+ (\d+) mapped(.+)')
    parseout = stp.search(str(stats.communicate()))
    if parseout:
        clipTotal = int(parseout.group(1)) + int(parseout.group(2))
    # Get total read number from RNA-seq .bam
    cmdl = ["samtools", "flagstat", str(background_bam_filename)]
    stats = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
    stp = re.compile(r'\\n(\d+) \+ (\d+) mapped(.+)')
    parseout = stp.search(str(stats.communicate()))
    if parseout:
        backgroundTotal = int(parseout.group(1)) + int(parseout.group(2))
        print("clip reads:%i rna seq reads:%i" % (clipTotal, backgroundTotal))
    return float(clipTotal)/float(backgroundTotal)
    #clip peak/clip total cf background peak/background total
    #is the same as
    #clip peak cf normalCoef * background peak


def dataset_size(clipReadsFname):
    cmdl = ["samtools", "flagstat", str(clipReadsFname)]
    stats = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
    stp = re.compile(r'\\n(\d+) \+ (\d+) mapped(.+)')
    parseout = stp.search(str(stats.communicate()))
    if parseout:
        clipTotal = int(parseout.group(1)) + int(parseout.group(2))
        return clipTotal
    else:
        print("Error: expected output missing.")


def guess_bam_file(peaksFname, silent=False):
    t = os.path.dirname(os.path.abspath(peaksFname))
    foldername = os.path.dirname(peaksFname).rstrip('/')
    #foldername = re.match(r'.+/([^/]+)$', t).group(1)
    if(not silent):
        print("\tGuessing bam file for foldername %s..." % t)
    bamfile = r'/home/dp/Desktop/bams/' + foldername + '.bam'
    if not silent:
        print("\t...Guessing bam file is %s." % bamfile)
    return bamfile


def remove_duplicate_ranges(bedfilename, outfilename):
    ranges = dict()
    with open(bedfilename, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            a_range = (s[1], s[2], s[3], s[6])
            if(a_range in ranges):
                if(ranges[a_range]['height'] > int(s[5])):
                    continue
            ranges[a_range] = {'height': int(s[5]), 'name': int(s[4]), 'line': li}
            #print int(s[4])
    with open(outfilename, 'w') as f:
        for a_range in ranges:
            f.write(ranges[a_range]['line'])


def merge_peaks(bedfile, results_folder):
    cmdl = "sort -k1,1 -k2,2n %s/%s > %s/%s.formerge" % (
        results_folder, bedfile, results_folder, bedfile)
    print(cmdl)
    os.system(cmdl)
    # The output of bedtools merge will contain the highest peak height in column 5
    # of the merged peaks, the the peak numbers in column 4 separated by ;.
    # Use -d 5 option to merge peaks within 5 nucleotides
    cmdl = "bedtools merge -s -scores max -nms -i %s/%s.formerge > %s/%s.merged.bed" % (
                results_folder, bedfile, results_folder, bedfile)
    print(cmdl)
    os.system(cmdl)
    cmdl = "rm %s/%s.formerge" % (results_folder, bedfile)
    print(cmdl)
    #os.system(cmdl)

    
def call_R(results_folder, peaks, src_path, which='background'):
    # R cannot handle huge files, so we'll chunk the input file.
    bin_filenames = {}
    bin_files = {}
    peak_stats_by_name = {}
    for index, p in enumerate(peaks, start=0):
        if index % 500:  # Not opening/closing a file.
            if which == 'background':
                bin_files[file_num].write(p.background_line_for_R)
            if which == 'clip':
                bin_files[file_num].write(p.clip_line_for_R)
        else:
            file_num = int(index/500)
            bin_filenames[file_num] = results_folder + "peaksForR.%s.%i.txt" % (which, file_num)
            bin_files[file_num] = open(bin_filenames[file_num], 'w')
            if which == 'background':
                bin_files[file_num].write(p.background_line_for_R)
            if which == 'clip':
                bin_files[file_num].write(p.clip_line_for_R)
            if index:
                # Close the last file, if there was one.
                bin_files[file_num - 1].close()
    try:
        bin_files[file_num].close()
    except:
        pass
    for index in bin_filenames:
        filename = bin_filenames[index]
        cmd = r"""Rscript %s/calculate_ztnb_p_value.r %s""" % (
            src_path, filename)
        os.system(cmd)
        with open('r.out', 'r') as f:
            for li in f:
                s = li.rstrip('\n').split('\t')
                peak_num = str(re.match('\s*(\d+)$', s[0]).group(1))
                peak_stats_by_name[peak_num] = float(s[3])
    for p in peaks:
        try:
            if which=='background':
                p.background_nb_p_value = peak_stats_by_name[p.name]
            if which=='clip':
                p.clip_nb_p_value = peak_stats_by_name[p.name]
        except:
            print("keys are %s" % str(list(peak_stats_by_name.keys())))
      
    
def ranges_with_stats_to_peaks(results_folder, binned_genes, peaks, annotation_file=''):
    filename = results_folder + '/ranges.peaks'
    not_by_gene_filename = filename.rstrip('peaks') + 'pvalue_by_bin.peaks'
    header = "Chr\tstart\tend\tname\theight\tstrand\tclip_poisson\tbackground_nb\tclip_nb\tgene\n"
    by_bin_f = open(not_by_gene_filename, 'w')
    by_bin_f.write(header)
    peaks = sorted(peaks, key=lambda x: float(x.height), reverse=True)
    for p in peaks:
        if not hasattr(p, 'background_nb_p_value'):
            print("Peak %s has no back p value." % p.gene_name)
        if not hasattr(p, 'clip_nb_p_value'):
            print("Peak %s has no clip p value." % p.gene_name)    
    with open(filename, 'w') as f:
        f.write(header)
        for p in peaks:
            if p.gene_name in binned_genes:
                f.write(p.write_as_peaks_format(
                    multiply_p_values=float(len(binned_genes[p.gene_name].bins))
                    ))
            by_bin_f.write(p.write_as_peaks_format())
    # The PMA1 5'UTR peak is sometimes assigned to the upstream gene, so we
    # just fix it by hand with sed.
    #print "Re-assigning YGL007C-A peaks to YGL008C (autofixing a common error)."
    #os.system("sed -i'.bak' 's/YGL007C-A/YGL008C/g' %s" % (annotated_results_file))
    return True

    
def assign_to_gene(in_file, out_file, annotation_file=''):
    bed_format_filename = os.path.dirname(os.path.realpath(in_file)) + '/ranges_no_dups.bed'
    print(bed_format_filename)
    bed_format = open(bed_format_filename, 'w')
    lines_from_ranges_no_dups = {}
    with open(in_file, "r") as f:
        for li in f:
            s = li.rstrip("\n").split('\t')
            lines_from_ranges_no_dups[s[4]] = li.rstrip('\n')
            bed_li = "{chrm}\t{start}\t{end}\t{name}\t{value}\t{strand}\n".format(
                chrm=s[1], start=s[2], end=s[3], name=s[4], value=s[5], strand=s[6]) 
            bed_format.write(bed_li + "\n")            
    bed_format.close()
    cmdl = "closestBed -s -a %s -b %s > %s_raw" % (
            bed_format_filename, annotation_file, out_file)
    
    os.system(cmdl)
    os.system("sed -i'' 's/YGL007C-A/YGL008C/g' %s" % (
            out_file + '_raw'))
    pat1 = r'ID=(\w+):([\w\-\(\)\d]+):'
    pat2 = r'ID=([\w\-\(\)\d]+);'
    out_lines = []
    known_ranges = set()
    with open(out_file + '_raw', 'r') as f:
        for li in f:
            s = li.rstrip("\n").split("\t")
            if(len(s) > 9):
                iv = (s[0], s[1], s[2], s[3], s[4], s[5])
                if iv in known_ranges:
                    continue  # Don't count the same range twice.
                known_ranges.add(iv)
                m = re.search(pat1, li)
                if(m is not None):
                    out_lines.append("%s\t%s\n" % (
                        lines_from_ranges_no_dups[s[3]], m.group(2)))
                    continue
                m = re.search(pat2, li)
                if(m is not None):
                    out_lines.append("%s\t%s\n" % (
                        lines_from_ranges_no_dups[s[3]], m.group(1)))
                    continue
    with open(out_file, 'w') as f:
        for line in out_lines:
            f.write(line)
            

def define_ranges(regions_above_cutoff_filename, clipReadsFname,
                  normalCoef, peaks_ranges_filename,
                  peak_heights, start_time, num_peaks):
    results_folder = os.path.dirname(regions_above_cutoff_filename) + '/'
    if(not os.path.exists(results_folder)):
        cmdl = r'mkdir ' + results_folder
        os.system(cmdl)
    peaksNewRanges = open(peaks_ranges_filename, 'w')
    # File io for diagnostic .wig file (aka bedgraph) for iv and
    # .bed file for stats regions.
    ivBedgraph = open(results_folder + 'iv.wig', 'w')
    statsIvBed = open(results_folder + 'statsIv.bed', 'w')
    lineO = "type=bedGraph\n"  # Required header line for bedgraph.
    ivBedgraph.write(lineO)
    f = open(regions_above_cutoff_filename, 'r')
    for num_processed, line in enumerate(f):
        s = line.rstrip("\n").split('\t')
        iv = [s[0], int(s[1]), int(s[2]), s[5]]  # [chr,start,stop,strand]
        p = peak(str(num_processed), iv)
        p.add_reads_and_adjust_range(clipReadsFname)
        peak_heights[p.name] = s[4]
        peaksNewRanges.write("%i\t%s\n" % (num_processed, p.write_range()))
        p.write_bedgraphs(ivBedgraph, statsIvBed)
        # Write current timing.
        if(num_processed > 0 and not (num_processed % 100)):
            sys.stderr.write("Defining range %i of %i:\t" % (
                num_processed, num_peaks))
            elapsedTime = (time.time() - start_time)
            timeToFinish = (num_peaks-num_processed)/(num_processed/elapsedTime)
            outline = "Time elapsed: %.1f m," % (elapsedTime/60)
            outline += " Time remaining: %.1f m\n" % (timeToFinish/60)
            sys.stderr.write(outline)
    f.close()
    peaksNewRanges.close()
    ivBedgraph.close()
    statsIvBed.close()
    

def process_ranges(peaks_ranges_filename, clip_bam_filename, background_bam_filename,
                   normalCoef,
                   gtf_filename='lib/Saccharomyces_cerevisiae.EF4.70.gtf',
                   annotation_file=''):
    """Find the true heights and statistics for an init_ranges file.
    The input file (init_ranges) contains peak numbers, locations,
    and, in the last two columns: height and position of max height.
    However, those last two columns are only approximations.
    """
    results_folder = os.path.dirname(peaks_ranges_filename)
    peak_heights = dict()
    #peak_stats = {}
    peak_ranges_file = open(peaks_ranges_filename, 'r')
    f_corrected_height = open(peaks_ranges_filename + '.cor_height', 'w')
    gtf_info = gtf_data(gtf_filename)
    binned_genes = {}
    use_local_for_these_genes = ['RDN18-1', 'RDN37-1', 'RDN58-1',
                                 'RDN18-2', 'RDN37-2', 'RDN58-2',
                                 'RDN5-1', 'RDN5-1',]
                                 #'YLR154W-A', 'YLR154W-B', 'YLR154W-C']
    read_ends_clip_peak = HTSeq.GenomicArray(chroms="auto", stranded=True)
    read_ends_clip_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    bins_clip_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    bins_clip_peak = HTSeq.GenomicArray(chroms="auto", stranded=True)
    bins_background_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    read_ends_background_peak = HTSeq.GenomicArray(chroms="auto", stranded=True)
    read_ends_background_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    (peaks, peak_heights) = peak_objects_from_ranges_file(peak_ranges_file)
    (peaks, peak_heights) = take_only_highest_peak_per_gene(peaks, peak_heights)
    for index, p in enumerate(peaks, start=1):
        if not (index % 100):
            print("Processing range %i." % index)
        p.add_reads_and_adjust_range(clip_bam_filename)
        if p.gene_name in gtf_info.txpt_ranges:
            if p.gene_name in binned_genes and (p.gene_name not in use_local_for_these_genes):
                p.put_clip_reads_in_bins(use_bins_from_gene=binned_genes[p.gene_name])
                p.calculate_poisson(a_binned_gene=binned_genes[p.gene_name])
            else:
                gene_iv = [p.iv[0], int(gtf_info.txpt_ranges[p.gene_name][0]),
                           int(gtf_info.txpt_ranges[p.gene_name][1]), p.iv[3]]
                if p.gene_name not in use_local_for_these_genes:
                    binned_genes[p.gene_name] = gene(
                        name=p.gene_name, gene_iv=gene_iv, gtf_info=gtf_info)
                    binned_genes[p.gene_name].add_clip_reads_in_gene_and_bin(clip_bam_filename)
                if background_bam_filename and (p.gene_name not in use_local_for_these_genes):
                    binned_genes[p.gene_name].add_background_reads_in_gene_and_bin(
                        background_bam_filename)
                if p.gene_name not in use_local_for_these_genes:
                    p.put_clip_reads_in_bins(use_bins_from_gene=binned_genes[p.gene_name])
                    p.calculate_poisson(a_binned_gene=binned_genes[p.gene_name])
                else:
                    p.put_clip_reads_in_bins(use_bins_from_gene=False)
                    p.calculate_poisson(a_binned_gene=False)
        else:
            print("Unknown gene %s..." % p.gene_name)
            p.put_clip_reads_in_bins(use_bins_from_gene=False)
            print(p.calculate_poisson())
        #peakStats[p.name] = {'Poisson': p.pvalue}
        if(background_bam_filename):
            if p.gene_name in gtf_info.txpt_ranges and (p.gene_name not in use_local_for_these_genes):
                p.set_line_for_R(
                    normalCoef=normalCoef, a_binned_gene=binned_genes[p.gene_name])
            else:
                # Add reads near the peak.
                p.addBackgroundReads(background_bam_filename)
                p.set_line_for_R(normalCoef=normalCoef) 
        else:
            p.set_line_for_R(normalCoef=normalCoef, output_zeros=True)
        p.add_reads_to_bedgraph(read_ends_clip_peak)
        p.add_background_reads_to_bedgraph(read_ends_background_peak)
        p.add_bins_to_bedgraph(bins_clip_peak)
        f_corrected_height.write("%s\t%s\n" % (p.name, p.write_range()))
        # Remove all the info from the .bam files to conserve space.
        p.remove_read_info()
    for a_gene in binned_genes:
        binned_genes[a_gene].add_reads_to_bedgraph(read_ends_clip_gene)
        binned_genes[a_gene].add_background_reads_to_bedgraph(read_ends_background_gene)
        binned_genes[a_gene].add_bins_to_bedgraph(bins_clip_gene, which_set="clip")
        binned_genes[a_gene].add_bins_to_bedgraph(
            bins_background_gene, which_set="background")
    outdir = os.path.dirname(os.path.realpath(peaks_ranges_filename))
    bins_clip_gene.write_bedgraph_file(outdir + '/bins_clip_gene_plus.wig', strand="+")
    bins_clip_gene.write_bedgraph_file(outdir + '/bins_clip_gene_minus.wig', strand="-")
    bins_clip_peak.write_bedgraph_file(outdir + '/bins_clip_peak_plus.wig', strand="+")
    bins_clip_peak.write_bedgraph_file(outdir + '/bins_clip_peak_minus.wig', strand="-")
    bins_background_gene.write_bedgraph_file(outdir+ '/bins_background_gene_plus.wig', strand="+")
    bins_background_gene.write_bedgraph_file(outdir+ '/bins_background_gene_minus.wig', strand="-") 
    read_ends_clip_peak.write_bedgraph_file(outdir+ '/read_ends_clip_peak_plus.wig', strand="+")
    read_ends_clip_gene.write_bedgraph_file(outdir+ '/read_ends_clip_gene_plus.wig', strand="+")
    read_ends_background_peak.write_bedgraph_file(outdir+ '/read_ends_background_peak_plus.wig', strand="+")
    read_ends_background_gene.write_bedgraph_file(outdir+ '/read_ends_background_gene_plus.wig', strand="+")
    read_ends_clip_peak.write_bedgraph_file(outdir+ '/read_ends_clip_peak_minus.wig', strand="-")
    read_ends_clip_gene.write_bedgraph_file(outdir+ '/read_ends_clip_gene_minus.wig', strand="-")
    read_ends_background_peak.write_bedgraph_file(outdir+ '/read_ends_background_peak_minus.wig', strand="-")
    read_ends_background_gene.write_bedgraph_file(outdir+ '/read_ends_background_gene_minus.wig', strand="-")
    f_corrected_height.close()
    return (peaks, binned_genes)

def take_only_highest_peak_per_gene(peaks, peak_heights):
    new_peaks_list = []
    new_peak_heights = {}
    peaks_by_gene = {}
    for p in peaks:
        if(p.gene_name not in peaks_by_gene):
            peaks_by_gene[p.gene_name] = p
        else:
            if peaks_by_gene[p.gene_name].height < p.height:
                peaks_by_gene[p.gene_name] = p
    for p in peaks:
        if p.name == peaks_by_gene[p.gene_name].name:
            new_peaks_list.append(p)
            new_peak_heights[p.name] = p.height
    print("Total peaks: %i" % len(peaks))
    print("Taking only the highest per gene: %i" % len(new_peaks_list))
    return (new_peaks_list, new_peak_heights)


def peak_objects_from_ranges_file(peak_ranges_file):
    peak_heights = {}
    peaks = []
    for line in peak_ranges_file:  # For each peak.
        s = line.rstrip("\n").split('\t')
        iv = [s[1], int(s[2]), int(s[3]), s[6]]
        pos_of_peak = int(s[7])
        a_gene = s[-1]
        p = peak(str(s[0]), iv, height=int(s[5]), pos_of_peak=pos_of_peak,
                 gene_name=s[-1],gtf_info=False)
        peak_heights[p.name] = s[4]
        peaks.append(p)
    return (peaks, peak_heights)

from __future__ import division  # int/int = float
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
from peak import peak
from gene import gene
from gtf_data import gtf_data


def normalize(clipReadsFname, backgroundReadsFname):
    # Get total read number from CLIP-seq .bam
    cmdl = ["samtools", "flagstat", str(clipReadsFname)]
    stats = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
    stp = re.compile(r'\\n(\d+) \+ (\d+) mapped(.+)')
    parseout = stp.search(str(stats.communicate()))
    if parseout:
        clipTotal = int(parseout.group(1)) + int(parseout.group(2))
    # Get total read number from RNA-seq .bam
    cmdl = ["samtools", "flagstat", str(backgroundReadsFname)]
    stats = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
    stp = re.compile(r'\\n(\d+) \+ (\d+) mapped(.+)')
    parseout = stp.search(str(stats.communicate()))
    if parseout:
        backgroundTotal = int(parseout.group(1)) + int(parseout.group(2))
    print "clip reads:%i rna seq reads:%i" % (clipTotal, backgroundTotal)
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
        print "Error: expected output missing."


def guess_bam_file(peaksFname, silent=False):
    t = os.path.dirname(os.path.abspath(peaksFname))
    foldername = os.path.dirname(peaksFname).rstrip('/')
    #foldername = re.match(r'.+/([^/]+)$', t).group(1)
    if(not silent):
        print "\tGuessing bam file for foldername %s..." % t
    bamfile = r'/home/dp/Desktop/bams/' + foldername + '.bam'
    if not silent:
        print "\t...Guessing bam file is %s." % bamfile
    return bamfile


def call_R(results_folder, peakStats, src_path):
    # R cannot handle huge files, so we'll chunk the input file.
    initialFile = open(results_folder + 'peaksForR.txt', 'r')
    numLinesInitialFile = 0
    numLines = 0
    binFiles = dict()
    binFilenames = dict()
    i = 0 #bin
    binFilenames[i] = results_folder + "peaksForR.%i.txt" % i
    binFiles[i] = open(binFilenames[i], 'w')
    for li in initialFile:
        numLinesInitialFile += 1
        numLines += 1
        binFiles[i].write(li)
        if(not(numLines % 500)):
            binFiles[i].close()
            i += 1
            binFilenames[i] = results_folder + "peaksForR.%i.txt" % i
            binFiles[i] = open(binFilenames[i], 'w')
    binFiles[i].close()
    initialFile.close()
    numLinesInSeparateFiles = 0
    # Call R and also check the sizes are correct.
    for i in binFilenames:
        binF = open( binFilenames[i], 'r')
        for li in binF:
            numLinesInSeparateFiles += 1
        binF.close()
        cmdl = r"""Rscript %s/calculate_ztnb_p_value.r %s""" % (
            src_path, binFilenames[i])
        os.system(cmdl)
        rf = open('r.out', 'r')
        for li in rf:
            s = li.rstrip("\n").split("\t")
            peakNum = re.match('\s*(\d+)$', s[0]).group(1)
            a_peakHeight = s[1]
            a_peakFdr = float(s[3])
            peakStats[str(peakNum)]['ZTNB'] = a_peakFdr
        rf.close()
    print "initial file had %i lines; separate files had %i lines" % (
        numLinesInitialFile, numLinesInSeparateFiles)
    return peakStats


def retrieve_p_values(premerge_fname, range_to_ann_line):
    peak_name_to_merged_peak = dict()
    for a_range in range_to_ann_line:
        s = range_to_ann_line[a_range]['line'].rstrip('\n').split('\t')
        peak_names = s[3].split(';')
        for a_name in peak_names:
            peak_name_to_merged_peak[a_name] = a_range
    with open(premerge_fname, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            # s[3] is the peak number.
            a_range = (s[0], s[1], s[2], s[5])
            peak_name = s[3]
            if(peak_name in peak_name_to_merged_peak):
                merged_peak_range = peak_name_to_merged_peak[peak_name]
            else:
                continue
            if(merged_peak_range in range_to_ann_line):
                known_s = range_to_ann_line[merged_peak_range]['line'].rstrip('\n').split('\t')
                known_s[6] = str(min(float(known_s[6]), float(s[6])))  # Set poisson
                known_s[7] = str(min(float(known_s[7]), float(s[7])))  # Set nb
                range_to_ann_line[merged_peak_range]['line'] = "\t".join(known_s) + "\n"
    return range_to_ann_line


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
    print cmdl
    os.system(cmdl)
    # The output of bedtools merge will contain the highest peak height in column 5
    # of the merged peaks, the the peak numbers in column 4 separated by ;.
    # Use -d 5 option to merge peaks within 5 nucleotides
    cmdl = "bedtools merge -s -scores max -nms -i %s/%s.formerge > %s/%s.merged.bed" % (
                results_folder, bedfile, results_folder, bedfile)
    print cmdl
    os.system(cmdl)
    cmdl = "rm %s/%s.formerge" % (results_folder, bedfile)
    print cmdl
    #os.system(cmdl)

    
def ranges_with_stats_to_peaks(bedfile, annotation_file=''):
    results_folder = os.path.dirname(bedfile)
    annotated_results_file = "%s/%s.ann" %(results_folder, os.path.basename(bedfile))
    cmdl = "closestBed -s -a %s -b %s > %s" % (bedfile, annotation_file,
                                               annotated_results_file)
    # The PMA1 5'UTR peak is sometimes assigned to the upstream gene, so we
    # just fix it by hand with sed.
    print "Re-assigning YGL007C-A peaks to YGL008C (autofixing a common error)."
    os.system("sed -i'.bak' 's/YGL007C-A/YGL008C/g' %s" % (annotated_results_file))
    print cmdl
    os.system(cmdl)
    range_to_ann_line = read_annotations_take_highest_peak_per_gene(annotated_results_file)
    # Sets the poisson and ztnb in the range_to_ann_line to be the minimum of the merged peaks.
    retrieve_p_values(bedfile, range_to_ann_line)
    write_ranges_with_gene_id_and_stats("%s/%s.peaks" % (
        results_folder, os.path.basename(bedfile)), range_to_ann_line)
    return True

    
def assign_to_gene(in_file, out_file, annotation_file=''):
    bed_format_filename = os.path.dirname(os.path.realpath(in_file)) + '/ranges_no_dups.bed'
    print bed_format_filename
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
            
    
def write_ranges_with_gene_id_and_stats(
    filename, range_to_ann_line):
    # Sort by height
    sorted_keys = sorted(range_to_ann_line.keys(),
                         key=lambda x: float(
                             range_to_ann_line[x]['line'].split("\t")[4]),
                         reverse=True)
    with open(filename, 'w') as f:
        header = "Chr\tstart\tend\tname\theight\tstrand\tpoisson\tztnb\tgene\n"
        f.write(header)
        for a_range in sorted_keys:
            if a_range in range_to_ann_line:
                li = range_to_ann_line[a_range]['line'].rstrip('\n')
                s = li.split('\t')
                s[3] = '.'
                li = "\t".join(s) + "\t"
                li += range_to_ann_line[a_range]['gene']
                li += "\n"
                f.write(li)


def read_annotations_take_highest_peak_per_gene(annotation_filename):
    pfile = open(annotation_filename, 'r')
    peaks_by_id = dict()
    pat1 = r'ID=(\w+):([\w\-\(\)\d]+):'
    pat2 = r'ID=([\w\-\(\)\d]+);'
    for li in pfile:
        s = li.rstrip("\n").split("\t")
        if(len(s) > 9):
            m = re.search(pat1, li)
            if(m is not None):
                if(m.group(2) in peaks_by_id):
                    if(float(peaks_by_id[m.group(2)][4]) >= float(s[4])):
                        continue
                peaks_by_id[m.group(2)] = s[0:6]+ [m.group(2)]
                continue
            m = re.search(pat2, li)
            if(m is not None):
                if(m.group(1) in peaks_by_id):
                    if(float(peaks_by_id[m.group(1)][4]) >= float(s[4])):
                        continue
                peaks_by_id[m.group(1)] = s[0:6] + [m.group(1)]
                continue
            gene_name = s[11]
            if(gene_name in peaks_by_id):
                if(float(peaks_by_id[gene_name][4]) >= float(s[4])):
                        continue
            peaks_by_id[gene_name] = s[0:6] + [gene_name]
    pfile.close()
    peaks_by_range = dict()
    for a_gene in peaks_by_id:
        s = peaks_by_id[a_gene]
        a_range = (s[0], s[1], s[2], s[5])
        line = "\t".join(str(x) for x in s[0:6]) + "\t"
        line += "\t".join(str(x) for x in [1.0, 1.0]) + "\n"
        peaks_by_range[a_range] = {'gene': a_gene, 'line': line}
    return peaks_by_range


def define_ranges(regions_above_cutoff_filename, clipReadsFname,
                  normalCoef, peaks_ranges_filename,
                  peak_heights, start_time, num_peaks):
    results_folder = os.path.dirname(regions_above_cutoff_filename)
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
        p.writeBedgraphs(ivBedgraph, statsIvBed)
        # Write current timing.
        if(num_processed > 0 and not (num_processed % 100)):
            sys.stderr.write("processing peak %i of %i:\t" % (
                num_processed, num_peaks))
            elapsedTime = (time.time() - start_time)
            timeToFinish = (num_peaks-num_processed)/(num_processed/elapsedTime)
            outline = "Time elapsed: %f.1 m," % (elapsedTime/60)
            outline += " Time remaining: %f.1 m\n" % (timeToFinish/60)
            sys.stderr.write(outline)
    f.close()
    peaksNewRanges.close()
    ivBedgraph.close()
    statsIvBed.close()
    

def process_ranges(peaks_ranges_filename, clipReadsFname, backgroundReadsFname,
                   normalCoef, peakStats,
                   gtf_filename='lib/Saccharomyces_cerevisiae.EF4.70.gtf',
                   annotation_file=''):
    """Find the true heights and statistics for an init_ranges file.
    The input file (init_ranges) contains peak numbers, locations,
    and, in the last two columns: height and position of max height.
    However, those last two columns are only approximations.
    """
    results_folder = os.path.dirname(peaks_ranges_filename)
    peak_heights = dict()
    f = open(peaks_ranges_filename, 'r')
    f_corrected_height = open(peaks_ranges_filename + '.cor_height', 'w')
    inputToR = open(results_folder + 'peaksForR.txt', 'w')
    gtf_info = gtf_data(gtf_filename)
    binned_genes = {}
    use_local_for_these_genes = ['RDN18-1', 'RDN37-1', 'RDN58-1',
                                 'RDN18-2', 'RDN37-2', 'RDN58-2',
                                 'RDN5-1', 'RDN5-1',]
                                 #'YLR154W-A', 'YLR154W-B', 'YLR154W-C']
    read_ends_clip_peak = HTSeq.GenomicArray(chroms="auto", stranded=True)
    read_ends_clip_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    bins_clip_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    bins_background_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    read_ends_background_peak = HTSeq.GenomicArray(chroms="auto", stranded=True)
    read_ends_background_gene = HTSeq.GenomicArray(chroms="auto", stranded=True)
    for line in f:
        s = line.rstrip("\n").split('\t')
        iv = [s[1], int(s[2]), int(s[3]), s[6]]
        pos_of_peak = int(s[7])
        a_gene = s[-1]
        p = peak(str(s[0]), iv, height=int(s[5]), pos_of_peak=pos_of_peak,
                 gene_name=s[-1],gtf_info=False)
        p.add_reads_and_adjust_range(clipReadsFname)
        peak_heights[p.name] = s[4]
        p.put_clip_reads_in_bins()
        if a_gene in gtf_info.txpt_ranges:
            if a_gene in binned_genes:
                p.calculate_poisson(a_binned_gene=binned_genes[a_gene])
            else:
                gene_iv = [s[1], int(gtf_info.txpt_ranges[a_gene][0]),
                           int(gtf_info.txpt_ranges[a_gene][1]), s[6]]
                binned_genes[a_gene] = gene(name=s[-1], gene_iv=gene_iv, gtf_info=gtf_info)
                if a_gene not in use_local_for_these_genes:
                    binned_genes[a_gene].add_clip_reads_in_gene_and_bin(clipReadsFname)
                if backgroundReadsFname and a_gene not in use_local_for_these_genes:
                    binned_genes[a_gene].add_background_reads_in_gene_and_bin(
                        backgroundReadsFname)
                if a_gene not in use_local_for_these_genes:
                    p.calculate_poisson(a_binned_gene=binned_genes[a_gene])
                else:
                    p.calculate_poisson(a_binned_gene=False)
        else:
            print "Unknown gene %s..." % a_gene
            print p.calculate_poisson()
        peakStats[p.name] = {'Poisson': p.pvalue}
        if(backgroundReadsFname):
            if a_gene in gtf_info.txpt_ranges and a_gene not in use_local_for_these_genes:
                lineO = p.write_background_bins(normalCoef=normalCoef,
                                                    a_binned_gene=binned_genes[a_gene])
                if(lineO):
                    inputToR.write(lineO)
            else:
                # Add reads near the peak in RNA-seq data
                p.addBackgroundReads(backgroundReadsFname)
                # Write bins of background data so that R can be called
                # to do a statistical test
                lineO = p.write_background_bins(normalCoef=normalCoef) 
                if(lineO):
                    inputToR.write(lineO)
        else:
            lineO = p.write_background_bins(normalCoef=normalCoef, output_zeros=True)
            if(lineO):
                inputToR.write(lineO)
        p.add_reads_to_bedgraph(read_ends_clip_peak)
        p.add_background_reads_to_bedgraph(read_ends_background_peak)
        f_corrected_height.write(
            "%s\t%s\n" % (p.name, p.write_range()))
    for a_gene in binned_genes:
        binned_genes[a_gene].add_reads_to_bedgraph(read_ends_clip_gene)
        binned_genes[a_gene].add_background_reads_to_bedgraph(read_ends_background_gene)
        binned_genes[a_gene].add_bins_to_bedgraph(bins_clip_gene, which_set="clip")
        binned_genes[a_gene].add_bins_to_bedgraph(
            bins_background_gene, which_set="background")
    outdir = os.path.dirname(os.path.realpath(peaks_ranges_filename))
    bins_clip_gene.write_bedgraph_file(outdir + '/bins_clip_gene_plus.wig', strand="+")
    bins_clip_gene.write_bedgraph_file(outdir + '/bins_clip_gene_minus.wig', strand="-")
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
    inputToR.close()

not_used = '''
def read_gff(filename):
    gene_ranges = {}
    pat = r'ID=([\w\-\(\)\d]+);'
    with open(filename, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            m = re.search(pat, li)
            if m is not None:
                gene_ranges[m.group(1)] = (int(s[1]), int(s[2]))
    return gene_ranges
'''

def add_p_value(ranges_without_stats_fname,
                ranges_with_stats_fname,
                peakStats):
    without_stats_f = open(ranges_without_stats_fname, 'r')
    with_stats_f  = open(ranges_with_stats_fname, 'w')
    numProcessed = -1
    for line in without_stats_f:
        numProcessed += 1
        s = line.rstrip("\n").split('\t')
        name = str(s[0])
        if name in peakStats:
            # We leave out the last column of peaksNewRanges.txt because
            # it is just the center of the peak, which is now unnecessary
            # information
            lineOut = "\t".join(s[1:7])
            lineOut += "\t%e\t%e\n" % (peakStats[name]['Poisson'], peakStats[name]['ZTNB'])
            with_stats_f.write(lineOut)
        else:
            print "Error: Did not find peak %s in peakStats." % str(numProcessed)
    without_stats_f.close()
    with_stats_f.close()

not_used = '''
def separate_by_p_value(results_folder, peakStats):
    newRangesF = open(results_folder + 'peaksNewRanges.txt', 'r')
    past_pvalF = open(results_folder + 'peaksPastFdr', 'w')
    allPeaksF  = open(results_folder + 'allpeaks', 'w')
    numProcessed = -1
    for line in newRangesF:
        numProcessed += 1
    #    if(numProcessed > 100):
    #        break
        s = line.rstrip("\n").split('\t')
        name = str(s[0])
        if name in peakStats:
            # We leave out the last column of peaksNewRanges.txt because
            # it is just the center of the peak, which is now unnecessary
            # information
            lineOut = "\t".join(s[1:7])
            lineOut += "\t%e\t%e\n" % (peakStats[name]['Poisson'], peakStats[name]['ZTNB'])
            allPeaksF.write(lineOut)
            if(peakStats[name]['ZTNB'] < 0.05 and peakStats[name]['Poisson'] < 0.001):
                past_pvalF.write(lineOut)
        else:
            print "error! did not find peak %s in peakStats" % str(numProcessed)
    newRangesF.close()
    past_pvalF.close()
    allPeaksF.close()


def normalize_peak_height_by_dataset_size(peaks_fname,
                                          out_peaks_fname, dataset_size):
    out_f = open(out_peaks_fname, 'w')
    peaks_lines = dict()
    coefficient = 1000000  # Per million
    header_line = True
    with open(peaks_fname, 'r') as f:
        for li in f:
            if(header_line):
                out_f.write(li)
                header_line = False
                continue
            s = li.rstrip('\n').split('\t')
            height = s[4]
            new_height = coefficient * float(height)/float(dataset_size)
            new_li = "\t".join(s[0:4])
            new_li += "\t%f\t" % new_height
            new_li += "\t".join(s[5:len(s)])
            new_li += "\n"
            out_f.write(new_li)
    out_f.close()


def normalize_peak_height_by_gene_abundance(peaks_fname,
                                            out_peaks_fname,
                                            geneAbundance):
    out_f = open(out_peaks_fname, 'w')
    peaks_lines = dict()
    coefficient = 1
    abundance_list = list()
    for gene in geneAbundance:
        abundance_list.append(geneAbundance[gene])
    median_abundance = np.median(abundance_list)
    header_line = True
    with open(peaks_fname, 'r') as f:
        for li in f:
            if(header_line):
                out_f.write(li)
                header_line = False
                continue
            s = li.rstrip('\n').split('\t')
            gene = s[8]
            if(gene in geneAbundance and geneAbundance[gene] > 0):
                genes_rel_abundance = float(geneAbundance[gene])/float(median_abundance)
                genes_rel_abundance = max(0.2, genes_rel_abundance)
                genes_rel_abundance = min(5, genes_rel_abundance)
                height_per_abundance = coefficient * float(s[4])/genes_rel_abundance
                if(re.search(r'YGL008C', li)):
                    print "norm by gene abundance():"
                    print "YBR115C mRNA has %f abundance." % genes_rel_abundance
                    print "CLIP-seq height is %f" % float(s[4])
                    print "new height is %f" % height_per_abundance
            else:
                continue
            new_li = "\t".join(s[0:4])
            new_li += "\t%f\t" % height_per_abundance
            new_li += "\t".join(s[5:len(s)])
            new_li += "\n"
            out_f.write(new_li)
    out_f.close()
'''

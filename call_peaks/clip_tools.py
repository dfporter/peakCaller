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


def moments(a): #takes averaged peak as argument
    squaredDiff = 0
    cubedDiff = 0
    sumOfObservations = 0
    summedXvalues = 0
    fourthPower = 0
    mean = 0
    #calculate array size
    arrsize = 0
    for step in list(a.steps()):
        if((step[0].start>11000) or (step[0].end>11000)):
            pass
        else:
            for x in range(step[0].start, step[0].end):
                arrsize+=1
    
    #alternative method: create a massive list
    observations = scipy.zeros((arrsize,2))
    obsList = list()
    for step in list(a.steps()):
        if((step[0].start>11000) or (step[0].end>11000)):
            pass
        else:
            #print step[0].start
            #print step[0].end
            for x in range(step[0].start, step[0].end):
                summedXvalues += step[1] * (x)
                sumOfObservations += float(step[1])
                observations[ x ] = (x, step[1])
                obsList.extend([x] * int(step[1]) )
    mean = float(summedXvalues)/sumOfObservations
    print "stats:"
    print observations[:,0]
    print "describe:", scipy.stats.describe(obsList)
    print scipy.stats.tmean(obsList, limits=(1,1000000))
    print " test", scipy.stats.skewtest(obsList)
    #print scipy.stats.tmin(observations,axis=0)
    for step in list(a.steps()):
        if((step[0].start>11000) or (step[0].start<9000) or (step[0].end>11000)):
            pass
        else:        
                squaredDiff += step[1] * ((x-mean)**2)
                cubedDiff += step[1] * ((x-mean)**3)
                fourthPower += step[1] * ((x-mean)**4)
                
    variance = float(squaredDiff)/float(sumOfObservations)
    standardDev = variance ** (0.5)
    skew = (float(cubedDiff)/float(sumOfObservations)) * (variance**(-3/2))
    kurtosis = float(fourthPower)/float(sumOfObservations) / (variance **2)
    print "variance: %f standard_deviation: %f skew: %f mean: %f kurtosis: %f " % (variance, standardDev, skew, mean, kurtosis)


def call_R(resultsFolder, peakStats):
    # R cannot handle huge files, so we'll chunk the input file.
    initialFile = open(resultsFolder + 'peaksForR.txt', 'r')
    numLinesInitialFile = 0
    numLines = 0
    binFiles = dict()
    binFilenames = dict()
    i = 0 #bin
    binFilenames[i] = resultsFolder + "peaksForR.%i.txt" % i
    binFiles[i] = open(binFilenames[i], 'w')
    for li in initialFile:
        numLinesInitialFile += 1
        numLines += 1
        binFiles[i].write(li)
        if(not(numLines % 500)):
            binFiles[i].close()
            i += 1
            binFilenames[i] = resultsFolder + "peaksForR.%i.txt" % i
            binFiles[i] = open(binFilenames[i], 'w')
    binFiles[i].close()
    initialFile.close()
    numLinesInSeparateFiles = 0
    # Call R and also check the sizes are correct
    for i in binFilenames:
        binF = open( binFilenames[i], 'r')
        for li in binF:
            numLinesInSeparateFiles += 1
        binF.close()
        cmdl = r"""R CMD BATCH '--args """ + binFilenames[i]
        cmdl += """' calculate_ztnb_p_value.r"""
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


def create_fasta(bedfname, resultsFolder):
    """Ignores strand. Not usable."""
    cmdl = r'fastaFromBed -fi ~/Desktop/ensemblEF4/ef4.fa '
    cmdl = cmdl + "-bed %s/%s -fo %s/%s.fa" % (resultsFolder, bedfname, resultsFolder, bedfname)
    print cmdl
    os.system(cmdl)
    analyzeFasta("%s/%s.fa" % (resultsFolder, bedfname), resultsFolder)


def count_motifs_in_fasta(fasta_filename):
    num_seqs = 0.0
    num_with_dual_site = 0.0
    num_with_any_site = 0.0
    seq_sizes = list()
    obs_chrs = list()
    sequences = dict((p.name, p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for region in sequences:
        num_seqs += 1.0
        seq_sizes.append(len(sequences[region]))
        obs_chrs.append(region.split(':')[0])
        if(len(re.findall(r'TAAT', sequences[region])) > 0):
            num_with_any_site += 1.0
        if(len(re.findall(r'TAAT\w{0,30}TAAT', sequences[region])) > 0):
            num_with_dual_site += 1.0
    zero_sites = num_seqs - num_with_any_site
    if(num_seqs == 0):
        return False
    fraction_any_site = num_with_any_site/num_seqs
    fraction_dual_site = num_with_dual_site/num_seqs
    print "**********"
    print "Sequences\t%.0f" % num_seqs 
    print "With a UAAU: %f (%f)" % (num_with_any_site, fraction_any_site)
    print "With a dual UAAU: %f (%f)" % (num_with_dual_site, fraction_dual_site)
    return (num_seqs, num_with_any_site, num_with_dual_site, seq_sizes, obs_chrs)

    
def analyze_fasta(fasta_filename, resultsFolder, do_random=True,
                  genome_fasta='~/Desktop/ensemblEF4/ef4.fa',
                  chr_sizes='/home/dp/Desktop/ensemblEF4/chr_lengths'):
    # Process information from the given fasta file
    print "CLIP:"
    (n_seqs_clip, n_any_site_clip, n_dual_site_clip, seq_sizes, obs_chrs
     ) = count_motifs_in_fasta(fasta_filename)
    if(not do_random):
        return
    # Generate a negative control dataset
    backup = r''''
    chrSizeTable = {
    'IX':439888,
    'VII':1090940,
    'XII':1078177,
    'II':813184,
    'VI':270161,
    'X':745751,
    'IV':1531933,
    'XIII':924431,
    'V':576874,
    'Mito':85779,
    'I':230218,
    'III':316620,
    'XI':666816,
    'XV':1091291,
    'XIV':784333,
    'XVI':948066,
    'VIII':562643}'''
    chrSizeTable = dict()
    with open(chr_sizes, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            chrSizeTable[s[0]] = int(s[1])
    randomFa = open(resultsFolder +'/random.bed', 'w')
    print "Writing to %s" % resultsFolder +'/random.bed'
    for r in range(1, 500):
        rlen = random.choice(seq_sizes)
        rchr = random.choice(obs_chrs)
        rstart = random.choice(range(1, chrSizeTable[rchr]))
        rend = rstart + rlen
        rstrand = random.choice( [ '+', '-' ] )
        lineO = "%s\t%i\t%i\t.\t.\t%s\n" % (rchr, rstart, rend, rstrand)
        randomFa.write(lineO)
    randomFa.close()
    cmdl = "fastaFromBed -fi %s -bed " % genome_fasta
    cmdl += resultsFolder + r'random.bed -fo ' + resultsFolder + r'random.fa'
    os.system(cmdl)
    # Process information from the negative control fasta file
    print "Random seqs:"
    (n_seqs_rand, n_any_site_rand, n_dual_site_rand, _, _) = count_motifs_in_fasta(
        resultsFolder + r'random.fa')
    oddsRatioDualSites, pvalueDualSite = scipy.stats.fisher_exact(
        [[n_seqs_clip - n_dual_site_clip, n_dual_site_clip],
        [n_seqs_rand - n_dual_site_rand, n_dual_site_rand]])
    oddsRatioAnySite, pvalueAnySite = scipy.stats.fisher_exact(
        [[n_seqs_clip - n_any_site_clip, n_any_site_clip],
        [n_seqs_rand - n_any_site_rand, n_any_site_rand]])
    print "p value Dual site\t%e\tp value any site\t%e\n" % (pvalueDualSite, pvalueAnySite)


def retrieve_p_values(resultsFolder, bedfile, range_to_ann_line):
    premerge_fname = "%s/%s" % (resultsFolder, bedfile)
    peak_name_to_merged_peak = dict()
    for a_range in range_to_ann_line:
        s = range_to_ann_line[a_range]['line'].rstrip('\n').split('\t')
        peak_names = s[3].split(';')
        for a_name in peak_names:
            peak_name_to_merged_peak[a_name] = a_range
    with open(premerge_fname, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            # s[3] is the peak number
            a_range = (s[0], s[1], s[2], s[5])
            peak_name = s[3]
            if(peak_name in peak_name_to_merged_peak):
                merged_peak_range = peak_name_to_merged_peak[peak_name]
            else:
                continue
            if(merged_peak_range in range_to_ann_line):
                known_s = range_to_ann_line[merged_peak_range]['line'].rstrip('\n').split('\t')
                known_s[6] = str(min(float(known_s[6]), float(s[6])))  # Set poisson
                known_s[7] = str(min(float(known_s[7]), float(s[7])))  # Set ztnb
                range_to_ann_line[merged_peak_range]['line'] = "\t".join(known_s) + "\n"
            #peakStats[s[3]]= {'Poisson': float(s[6]),
            #                  'ZTNB': float(s[7])}
            #range_to_ranges_line[a_range] = {'Poisson': float(s[6]),
            #                          'ZTNB': float(s[7]),
            #                           'line': li.rstrip('\n')}
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


def merge_peaks(bedfile, resultsFolder):
    cmdl = "sort -k1,1 -k2,2n %s/%s > %s/%s.formerge" % (
        resultsFolder, bedfile, resultsFolder, bedfile)
    print cmdl
    os.system(cmdl)
    # The output of bedtools merge will contain the highest peak height in column 5
    # of the merged peaks, the the peak numbers in column 4 separated by ;.
    # Use -d 5 option to merge peaks within 5 nucleotides
    cmdl = "bedtools merge -s -scores max -nms -i %s/%s.formerge > %s/%s.merged.bed" % (
                resultsFolder, bedfile, resultsFolder, bedfile)
    print cmdl
    os.system(cmdl)
    cmdl = "rm %s/%s.formerge" % (resultsFolder, bedfile)
    print cmdl
    #os.system(cmdl)

    
def ranges_with_stats_to_peaks(bedfile, resultsFolder,
                               annotation_file='gff2clean.bed', use_merged=False):
    if(use_merged):
        cmdl = "closestBed -s -a %s/%s.merged.bed -b %s > %s/%s.ann" % (
            resultsFolder, bedfile, annotation_file, resultsFolder, bedfile)
        # The PMA1 5'UTR peak is sometimes assigned to the upstream gene, so we
        # just fix it by hand with sed
        os.system("sed -i'.bak' 's/YGL007C-A/YGL008C/g' %s/%s.ann" % (
            resultsFolder, bedfile))
        print cmdl
        os.system(cmdl)
        range_to_ann_line = read_annotations_take_highest_peak_per_gene(
            "%s/%s.ann" % (resultsFolder, bedfile))
        # Sets the poisson and ztnb in the range_to_ann_line to be the minimum of the merged peaks
        retrieve_p_values(resultsFolder, bedfile, range_to_ann_line)
        write_ranges_with_gene_id_and_stats("%s/%s.peaks" % (
            resultsFolder, bedfile), range_to_ann_line)
        return True
    cmdl = "closestBed -s -a %s/%s -b %s > %s/%s.ann" % (
            resultsFolder, bedfile, annotation_file, resultsFolder, bedfile)
    os.system("sed -i'' 's/YGL007C-A/YGL008C/g' %s/%s.ann" % (
            resultsFolder, bedfile))
    os.system(cmdl)
    range_to_ann_line = read_annotations_take_highest_peak_per_gene("%s/%s.ann" % (
             resultsFolder, bedfile))
    retrieve_p_values(resultsFolder, bedfile, range_to_ann_line)
    write_ranges_with_gene_id_and_stats("%s/%s.peaks" % (
            resultsFolder, bedfile), range_to_ann_line)

    
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
            #print "%s did not match regex..." % li
            gene_name = s[11]
            if(gene_name in peaks_by_id):
                if(float(peaks_by_id[gene_name][4]) >= float(s[4])):
                        continue
            peaks_by_id[gene_name] = s[0:6] + [gene_name]
    print peaks_by_id
    pfile.close()
    peaks_by_range = dict()
    for gene in peaks_by_id:
        s = peaks_by_id[gene]
        a_range = (s[0], s[1], s[2], s[5])
        line = "\t".join(str(x) for x in s[0:6]) + "\t"
        line += "\t".join(str(x) for x in [1.0, 1.0]) + "\n"
        peaks_by_range[a_range] = {'gene': gene, 'line': line}
        #print "adding gene %s at the end" % gene
    print peaks_by_range
    return peaks_by_range

old = r'''
def take_only_highest_peak(resultsFolder, bedfile):
    peaks_fname = "%s/%s.pvalue.sort" % (resultsFolder, bedfile)
    outfname = "%s/%s.pvalue.highest" % (resultsFolder, bedfile)
    pfile = open(peaks_fname, 'r')
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
                peaks_by_id[m.group(2)] = s[0:6] + [s[12]] + [m.group(2)]
                continue
            m = re.search(pat2, li)
            if(m is not None):
                if(m.group(1) in peaks_by_id):
                    if(float(peaks_by_id[m.group(1)][4]) >= float(s[4])):
                        continue
                peaks_by_id[m.group(1)] = s[0:6] + [s[12]] + [m.group(1)]
                continue
            print "%s did not match regex..." % li
    pfile.close()
    outf = open(outfname, 'w')
    for geneid in peaks_by_id:
        li = "\t".join(peaks_by_id[geneid]) + "\n"
        outf.write(li)
    outf.close()
'''

def define_ranges(resultsFolder,
                  regions_above_cutoff_filename,
                  clipReadsFname,
                  backgroundReadsFname,
                  normalCoef,
                  peaks_ranges_filename,
                  peakHeights, startTime, numPeaks):
    if(not os.path.exists(resultsFolder)):
        cmdl = r'mkdir ' + resultsFolder
        print cmdl
        os.system(cmdl)
    peaksNewRanges = open(peaks_ranges_filename, 'w')
    # file io for diagnostic wig file (aka bedgraph) for iv and
    # bed file for stats regions
    ivBedgraph = open(resultsFolder + 'iv.wig', 'w')
    statsIvBed = open(resultsFolder + 'statsIv.bed', 'w')
    lineO = "type=bedGraph\n" #required header line for bedgraph
    ivBedgraph.write(lineO)
    # For each region of interest in the input .bed file,
    # write the peak number, CLIP height and bins from the given region to
    # a peaksForR.txt file for statistical testing in R
    f = open(regions_above_cutoff_filename, 'r')
    numProcessed = -1
    for line in f:
        numProcessed += 1
        s = line.rstrip("\n").split('\t')
        iv = [s[0], int(s[1]), int(s[2]), s[5]]
        # .bed files are CHR start   stop    name    score   strand,
        # so iv is [chr,start,stop,strand]
        p = peak(str(numProcessed), iv)
        peakHeights[p.name] = s[4]
        # Write current timing
        if(numProcessed > 0 and not (numProcessed % 100)):
            sys.stderr.write("processing peak %i of %i:\t" % (
                numProcessed, numPeaks))
            elapsedTime = (time.time() - startTime)
            timeToFinish = (numPeaks-numProcessed)/(numProcessed/elapsedTime)
            outline = "Time elapsed: %f.1 m," % (elapsedTime/60)
            outline += " Time remaining: %f.1 m\n" % (timeToFinish/60)
            sys.stderr.write(outline)
        #if(numProcessed > 100):
        #   break
        p.add_reads_in_peak(clipReadsFname) # Add reads under the peak from CLIP-seq data
        p.set_center_add_reads_in_region(clipReadsFname) # Add reads in the general region
        p.adjust_range() # Changes self.iv[1]/.iv[2] start and stop positions
        peaksNewRanges.write(
            "%i\t%s\n" % (numProcessed, p.write_range()))
        p.writeBedgraphs(ivBedgraph, statsIvBed)
    f.close()
    peaksNewRanges.close()
    ivBedgraph.close()
    statsIvBed.close()
    

def process_ranges(resultsFolder,
                   peaks_ranges_filename,
                   clipReadsFname,
                   backgroundReadsFname,
                   normalCoef,
                   peakStats):
    peakHeights = dict()
    f = open(peaks_ranges_filename, 'r')
    inputToR = open(resultsFolder + 'peaksForR.txt', 'w')
    for line in f:
        s = line.rstrip("\n").split('\t')
        iv = [s[1], int(s[2]), int(s[3]), s[6]]
        center = int(s[7])
        p = peak(str(s[0]), iv, height=s[5], center=center)
        peakHeights[p.name] = s[4]
        p.add_reads_in_peak(clipReadsFname)
        if(not p.set_center_add_reads_in_region(clipReadsFname)):
            print "Hit a fail condition on line %s..." % line
            continue
        p.set_clip_bins()
        p.calculate_poisson()
        peakStats[p.name] = {'Poisson': p.pvalue}
        if(backgroundReadsFname):
            # Add reads near the peak in RNA-seq data
            p.addBackgroundReads(backgroundReadsFname)
            # Write bins of background data so that R can be called
            # to do a statistical test
            lineO = p.write_background_bins(normalCoef) 
            if(lineO):
                inputToR.write(lineO)
        else:
            lineO = p.write_background_bins(1.0, output_zeros=True)
            if(lineO):
                inputToR.write(lineO)
    inputToR.close()
    print "Finished processing ranges."


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
            print "error! did not find peak %s in peakStats" % str(numProcessed)
    without_stats_f.close()
    with_stats_f.close()


def separate_by_p_value(resultsFolder, peakStats):
    newRangesF = open(resultsFolder + 'peaksNewRanges.txt', 'r')
    past_pvalF = open(resultsFolder + 'peaksPastFdr', 'w')
    allPeaksF  = open(resultsFolder + 'allpeaks', 'w')
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

import sys
import os
import pybedtools
import HTSeq
from scipy.stats import poisson
import pysam
import re
import subprocess
from gene import gene

class peak:
    """A CLIP-seq peak.
    Initially set with only a genomic range and name.
    Contains methods to add reads from a bam file of CLIP reads
    and from a bam file of control reads.
    """
    def __init__(self, name, iv, height=0, pos_of_peak=0, gene_name='NA'):
        self.name = name
        self.iv = iv
        self.height = height
        self.bin_size = 50
        self.halfWidth = 1000  # For viewing bams
        self.range_for_bins = 500  # For the poisson and ZTNB test
        self.pos_of_peak = pos_of_peak
        self.gene_name = gene_name
        if(self.iv[3] == "+"):
            self.strandSelect = "-F 0x10"
        if(self.iv[3] == "-"):
            self.strandSelect = "-f 0x10"
        #define the 1kb range around the given peak
        print "Created peak %s with range %s and height %i" % (str(name), str(iv), height)
        self.center_of_region = int(float(self.iv[1] + self.iv[2])/2.0)
        self.set_boundries_for_viewing_bams(use_pos_of_peak=False)
        #print "Center set at %i" % self.center_of_region

    def add_reads_and_adjust_range(self, clipReadsFname):
        """Find the position of the peak, its extent and the read distribution.
        Called immediately after initializing a peak to define its shape.
import pysam
import peak
p = peak.peak("0", ['XVI', 96581, 96582, '-'])
clipReadsFname = '/home/dp/Desktop/bams/PUF_domain.bam'
p.add_reads_in_peak('/home/dp/Desktop/bams/PUF_domain.bam')

reload(peak)
        """
        self.add_reads_in_peak(clipReadsFname) # Add reads under the peak from CLIP-seq data
        self.set_center_add_reads_in_region(clipReadsFname)
        self.adjust_range_around_peak()  # Changes self.iv[1]/.iv[2] start and stop positions
        
    def set_boundries_for_viewing_bams(self, use_pos_of_peak=False):
        if use_pos_of_peak:
            diff = self.pos_of_peak - self.halfWidth
            self.lower_bound = diff if (diff > 0) else 1
            self.upper_bound = self.pos_of_peak + self.halfWidth
        else:
            diff = self.center_of_region - self.halfWidth
            self.lower_bound = diff if (diff > 0) else 1
            self.upper_bound = self.center_of_region + self.halfWidth
            
    def add_reads_in_peak(self, clipReadsFname):
        """Adds all reads in the genomic interval
        of the peak. Passed a bam file name. Only reads in the specific interval,
        not in the larger self.upper_bound/lowerBound or bins range.
        """
        bamfile = pysam.AlignmentFile(clipReadsFname, "rb")
        s = bamfile.fetch(self.iv[0], self.iv[1], self.iv[2])
        self.reads = list()
        self.ga_true_coverage = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in s:
            if((self.iv[3]=="+" and not r.is_reverse) or (self.iv[3]=="-" and r.is_reverse)):
                riv = HTSeq.GenomicInterval(
                     self.iv[0], r.reference_start, r.reference_end, self.iv[3])
                self.ga_true_coverage[riv] += 1
        bamfile.close()

    def find_position_of_max_coverage(self):
        """Finds the position of the highest coverage.
        Sets self.pos_of_peak and self.ga_true_coverage.
        Also sets self.height, the height of the peak.
        These are not the values used for statistics, which use only the
        5' end of the reads (self.ga_read_starts), not read depth by position.
        The only peak function that calls find_positon_of_max_coverage is
        set_center_add_reads_in_region(),
        which is not called by any other peak function (only called externally).
        self.ga_true_coverage is first set by add_reads() viewing the bam file with
        the initial self.iv range, and then later reset by viewing the bam file
        with the self.upper/lower_bound range.
        """
        # Find the maximum of the peak.
        highestValue = 0
        for s in self.ga_true_coverage.steps():
            if(int(s[1]) > highestValue):
                highestValue = int(s[1])
                highestValueIV = s[0]
                #print s
        #print highestValueIV.start, highestValueIV.end, highestValue
        if highestValue == 0:
            print "Error: why is this range empty? %s" % str(self.ga_true_coverage.steps())
            self.pos_of_peak = float(self.iv[1] + self.iv[2])/2.0
            self.height = 0
            return False
        highestValueCoord = float(highestValueIV.end + highestValueIV.start)/float(2)
        self.pos_of_peak = int(highestValueCoord)
        print "Setting peak to (pos=%i, height=%i), previously height=%i." % (
            self.pos_of_peak, highestValue, self.height)
        self.height = highestValue
        self.set_boundries_for_viewing_bams(use_pos_of_peak=True)
        #print "find_center set center as %i" % self.pos_of_peak
        #print "and height as %f" % float(self.height)
        return True

    def set_center_add_reads_in_region(self, clipReadsFname):
        """Sets and fills in the genomic array of the peak's area, reading from
        the bam file and using the self.lower_bound/upperBound range.
        Assumes self.ga exists. Finds the peak position of self.ga.
        Uses self.lower_bound, which was initialized as 1000 kb before
        the center of the peak.
        """
        #print "In determineGA function, making a call to findCenter..."
        if(not self.find_position_of_max_coverage()):
            return False
        self.add_reads_in_general_region(clipReadsFname)
        return True
    
    def add_reads_in_general_region(self, clipReadsFname):
        """Add reads around the position of highest coverage.
        When a peak is first created, reads are added to the given range
        (everywhere continuously above 10 reads), and then
        set_center_add_reads_in_region() is called, which finds the position
        of max coverage and sets self.pos_of_peak, then resets self.lower_bound
        and self.upper_bound around that position. Then this function is called
        to add reads in that region.
        """
        self.ga_true_coverage = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        self.ga_read_starts = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        bamfile = pysam.AlignmentFile(clipReadsFname, "rb")
        s = bamfile.fetch(self.iv[0], self.lower_bound, self.upper_bound)
        #self.reads = list()
        for r in s:
            if((self.iv[3]=="+" and not r.is_reverse) or (self.iv[3]=="-" and r.is_reverse)):
                riv = HTSeq.GenomicInterval(
                     self.iv[0], r.reference_start, r.reference_end, self.iv[3])
                self.ga_true_coverage[riv] += 1
                r_pos = HTSeq.GenomicPosition(
                    self.iv[0], r.reference_start, self.iv[3])
                self.ga_read_starts[r_pos] += 1
        bamfile.close()

    def addBackgroundReads(self, backgroundReadsFname):
        self.background_ga = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        bamfile = pysam.AlignmentFile(backgroundReadsFname, "rb")
        s = bamfile.fetch(self.iv[0], self.lower_bound, self.upper_bound)
        for r in s:
            if((self.iv[3]=="+" and not r.is_reverse) or (self.iv[3]=="-" and r.is_reverse)):
                r_pos = HTSeq.GenomicPosition(
                    self.iv[0], r.reference_start, self.iv[3])
         #           str(s[2]), int(s[3]), self.iv[3])
                self.background_ga[r_pos] += 1
        bamfile.close()
        nm = """
        region = "%s:%s-%s" % (self.iv[0], self.lower_bound, self.upper_bound)
        cmdl = ["samtools", "view", self.strandSelect, backgroundReadsFname, region]
        readsRegion = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
        # Set self.backgroundReads
        self.background_reads = list()
        background_samtools_result = readsRegion.communicate()[0].split('\n')
        for r in background_samtools_result:
            #read_name \t number \t chr \t locus \t 255 \t 30M
            if(r is not ""):
                self.background_reads.append(r)
        # Set self.background_ga
        self.background_ga = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in self.background_reads:
            s = r.split('\t')
            matchLen = 0
            if(len(s) > 5):
                for segment in re.findall(r'(\d+)[MD]', s[5]):
                    matchLen += int(segment)
                read_end = int(s[3])
                if(self.iv[3] == "-"):
                    read_end = int(s[3]) + matchLen - 1
                r_pos = HTSeq.GenomicPosition(
                    str(s[2]), int(s[3]), self.iv[3])
                self.background_ga[r_pos] += 1
            else:
                print "Read with less than 5 columns in samtools output: %s" % str(r)            
"""
    def writeBedgraphs(self, ivBedgraph, statsIvBed):
        #format=chromA  chromStartA  chromEndA  dataValueA
        line_newRegions = "%s\t%i\t%i\t%i\n" % (self.iv[0], self.iv[1], self.iv[2], self.height)
        ivBedgraph.write(line_newRegions)
        line_statsIv = "%s\t%i\t%i\t%s\t%i\t%s\n" % (
            self.iv[0], self.lower_bound, self.upper_bound,
            self.name, self.height, self.iv[3])
        statsIvBed.write(line_statsIv)

    def adjust_range_around_peak(self):
        """Sets self.center_of_region, and self.iv[0] and self.iv[1],
        based on the point at which the peak drops below 20% of its height.
        Expects self.pos_of_peak to be known and does not reset it.
        """
        cutoff = int(0.2 * self.height)
        right_edge = self.pos_of_peak + 100000
        left_edge = self.pos_of_peak - 100000
        found_edge = False
        start_of_search_region = self.pos_of_peak
        while(not found_edge):
            end_of_search_region = start_of_search_region + 100
            k = list(self.ga_true_coverage[HTSeq.GenomicInterval(self.iv[0],
                                                     start_of_search_region,
                                                     end_of_search_region,
                                                     self.iv[3])].steps())
            # Find the range below cutoff with the lowest start position,
            # in this bin.
            for j in k:
                if ((j[1] <= cutoff) and (j[0].start < right_edge)):
                        right_edge = j[0].start
                        found_edge = True
            start_of_search_region += 100
        found_edge = False
        start_of_search_region = self.pos_of_peak - 100
        end_of_search_region = self.pos_of_peak
        if(start_of_search_region < 1):
            left_edge = 1
            found_edge = True
        while(not found_edge):
            k = list(self.ga_true_coverage[HTSeq.GenomicInterval(
                    self.iv[0],
                     start_of_search_region,
                     end_of_search_region,
                     self.iv[3])].steps())
            # Find the range below cutoff with the highest start position.
            for j in k:
                if ((j[1] <= cutoff) and (j[0].end > left_edge)):
                        left_edge = j[0].end
                        found_edge = True
            start_of_search_region -= 100
            if(start_of_search_region < 1):
                left_edge = 1
                found_edge = True
        print "Left edge set as %i" % left_edge
        if(left_edge >= 0):
            self.iv[1] = left_edge
        else:
            self.iv[1] = 0
        self.iv[2] = right_edge
        if(abs(self.iv[2] - self.iv[1]) < 40):
                self.iv[1] = self.pos_of_peak - 20
                if(self.iv[1] < 1):
                    self.iv[1] = 1
                self.iv[2] = self.pos_of_peak + 20
        self.center_of_region = int(float(self.iv[1] + self.iv[2])/2.0)
        self.set_boundries_for_viewing_bams()
        #print "adjust_range(): center set as %i based on range extremes of 0.2 height" % self.center
        return True

    def write_range(self):
        # .bed is chr start stop name score strand
        lineO = "%s\t%i\t%i\t%s\t%i\t%s\t%i" % (
            self.iv[0], self.iv[1], self.iv[2], str(self.name),
            int(self.height), self.iv[3], self.pos_of_peak)
        print "write_range(): %s" % lineO
        return lineO

    def set_bin_locations(self):
        self.num_bins = int(float(self.range_for_bins)/float(self.bin_size))
        # Define bin_lower_bound so that the center bin countains the center of the peak.
        self.bin_lower_bound = self.pos_of_peak - self.range_for_bins - int(float(self.bin_size)/2.0) 
        if(self.bin_lower_bound < 1):
            self.bin_lower_bound = 1

    def max_value_in_range_read_ends(self, chrm, start, end, strand):
        # Calculate the maximimum value in the bin.
        highV = 0
        k = list(self.ga_read_starts[HTSeq.GenomicInterval(
            chrm, start, end, strand)].steps() )
        for j in k:
            if (int(j[1]) > highV):
                highV = int(j[1])
        return highV

    def set_clip_bins(self):
        self.set_bin_locations()
        self.clip_bin = list()
        for i in range(self.bin_lower_bound,
                       self.pos_of_peak + self.range_for_bins - self.bin_size,
                       self.bin_size):
            highV = self.max_value_in_range_read_ends(self.iv[0], i, i+self.bin_size,
                                                      self.iv[3])
            if(highV >= 0):
                self.clip_bin.append(float(highV))
        if(len(self.clip_bin) < 1):
            self.clip_bin = [0]

    def set_background_bins(self, normalCoef):
        background_bin = list()
        verb = False
        for i in range(self.bin_lower_bound,
                       self.pos_of_peak + self.range_for_bins - self.bin_size,
                       self.bin_size):
            #bin is i to (i+l)
            #calculate the maximimum value in the bin
            #if(verb):
            #    print "i=%i" % i
            highV= 0
            k = list(self.background_ga[HTSeq.GenomicInterval(
                self.iv[0], i, i+self.bin_size, self.iv[3])].steps())
            for j in k:
                #if(verb):
                #    print j
                if (j[1] > highV):
                    highV = int(j[1])
            if((highV * normalCoef) > 0):
                background_bin.append(float(highV * normalCoef))
                #if(verb):
                #    print "adding %i normalized by %f to %f" % (
                #        highV, normalCoef, float(highV * normalCoef))
        background_bin.sort(reverse=True)
        if(len(background_bin) < 1):
            background_bin = [0]
        #if(verb):
        #    self.ga.write_bedgraph_file("peak91_peak.wig", "+")
        #    self.background_ga.write_bedgraph_file("peak91_background.wig", "+")
        return background_bin
            
    def write_background_bins(self, normalCoef, a_binned_gene=False, output_zeros=False):
        #suppose file is of the format:
        #peakNumber     height  1,4,2,5,2
        #the R code takes a file in which one column is a comma separated list of integers and sets y equal to the set of integers
        #t = read.table(file="", sep="\t")
        #y =as.numeric(strsplit(as.character(t[1,3]), ',', fixed=TRUE)[[1]])
        #fit a negative binomial to the bins
        #ff = fitdistr(y, "Negative Binomial")
        #est_size = ff$estimate[1]
        #est_mu = ff$estimate[2]
        #this code then gives the p value for getting the observed height from a negative binomial distribution built from the given bins
        #1 - pnbinom(q=t[1,2], size=est_size, mu=est_mu)
        
        # If there is no background, make sure there is something
        # for R to model so R doesn't crash
        if a_binned_gene:
            background_bin = a_binned_gene.background_bins
            lineO = "%s\t%s\t" % (self.name, self.reads_in_peak_bin)
            for j in background_bin:
                lineO += "%f," % float(j)
            lineO += "\n"
            print "*(************"
            print lineO
            return lineO
        if(output_zeros):
            background_bin = [0.0, 0.0, 0.0, 0.0, 1.0]
        else:
            background_bin = self.set_background_bins(normalCoef)
        if(len(background_bin) < 3):
            background_bin = [0.0, 0.0, 0.0, 0.0, 1.0]
        lineO = "%s\t%s\t" % (self.name, self.reads_in_peak_bin)
        for j in background_bin:
            lineO += "%f," % float(j)
        lineO += "\n"
        return lineO

    def find_reads_in_peak_bin(self):
        """Reads in bins are counted from the 5' end, while the center of the peak is counted
        as the true max height of the peak. So counting the "height" of the peak
        as the number of reads in the center bin may give a lower number than to the side.
        We will therefore allow some flexibility in where the "center bin" is.
        """
        #self.reads_in_peak_bin = self.clip_bin[int(float(self.num_bins)/2.0)]
        left_bins = []
        right_bins = []
        # What bins overlap the peak range?
        print "peak range %s, bin range %s" % (str(self.iv),
                                               str((self.bin_lower_bound,
                       self.pos_of_peak + self.range_for_bins - self.bin_size,
                       self.bin_size)))
        for i in range(self.bin_lower_bound-1000,
                       self.pos_of_peak + self.range_for_bins - self.bin_size+1000,
                       self.bin_size):
            if(i <= self.iv[1] <= i+self.bin_size):
                # Left peak boundary overlaps this bin
                overlap_left = max([(i+self.bin_size-self.iv[1]), (self.iv[1]-i)])
                print "Comparing left boundary from %s with bin %s. Overlap %i." % (
                    str(self.iv), str((i, i+self.bin_size)), overlap_left)
                if float(overlap_left) > (0.1 * float(self.bin_size)):
                    print "bin ."
                    left_bins.append(i)
                else:
                    # We are guaranteed that all further bins are past the left edge.
                    # So we will include the next bin.
                    # In the case of an edge, two bins are included - we take the lower one.
                    left_bins.append(i+self.bin_size)
            if(i <= self.iv[2] <= i+self.bin_size):
                # Right peak boundary overlaps this bin
                overlap_right = max([(i+self.bin_size-self.iv[2]), (self.iv[2]-i)])
                print "Comparing right boundary from %s with bin %s. Overlap %i." % (
                    str(self.iv), str((i, i+self.bin_size)), overlap_right)
                if float(overlap_right) > (0.1 * float(self.bin_size)):
                    print "Including bin."
                    right_bins.append(i)
                else:
                    right_bins.append(i+self.bin_size)
        print "left bins %s" % str(left_bins)
        print "right bins %s" % str(right_bins)
        self.leftmost_bin = min(left_bins)
        self.rightmost_bin = max(right_bins)
        if self.leftmost_bin == self.rightmost_bin:
            self.rightmost_bin = self.leftmost_bin + self.bin_size
        values_in_bins = [0]
        for i in range(self.leftmost_bin, self.rightmost_bin, self.bin_size):
            values_in_bins.append(self.max_value_in_range_read_ends(
                self.iv[0], i, i+self.bin_size, self.iv[3]))
            print "bin %i has value %f" % (i, float(values_in_bins[-1]))
        self.reads_in_peak_bin = max(values_in_bins)
        
    def calculate_poisson(self, a_binned_gene=False, return_boolean=False):
        """
        """
        if a_binned_gene:
            print "Using binned clip for poisson..."
            gene_bins = a_binned_gene.bins
            print "Bins from clip on this gene (%s): %s" % (self.gene_name, str(gene_bins))
            mu = float(sum(gene_bins))/float(len(gene_bins))
            pois_dist = poisson(mu)
            self.find_reads_in_peak_bin()
            self.pvalue = 1-pois_dist.cdf(self.reads_in_peak_bin)
            return "%s\tRange=%s,mu=%e\treads_in_peak_bin=%e\tp_value=%e\tbins_around_peak=%s\n" % (
                self.name, str(self.iv), mu, self.reads_in_peak_bin,
                self.pvalue, str(self.clip_bin))
        #var = scipy.stats.tvar(self.clip_bin)
        #mean = scipy.stats.tmean(self.clip_bin)
        pois_dist = poisson(float(sum(self.clip_bin))/float(len(self.clip_bin)))
        self.find_reads_in_peak_bin()
        self.pvalue = 1-pois_dist.cdf(self.reads_in_peak_bin)
        # Poisson with smaller range.
        middle_bin_start = int((len(self.clip_bin) * 0.5))
        middle_bin_end = int(len(self.clip_bin) * 1.5)
        pois_dist = poisson(float(sum(self.clip_bin[middle_bin_start:middle_bin_end+1]))/float(
            len(self.clip_bin[middle_bin_start:middle_bin_end+1])))
        middle_region_pvalue = 1-pois_dist.cdf(self.reads_in_peak_bin)
        #if middle_region_pvalue > self.pvalue:
        #    self.pvalue = middle_region_pvalue
        if(return_boolean):
            if(self.pvalue < 0.001):
                return True
            else:
                return False
        mu = float(sum(self.clip_bin))/float(len(self.clip_bin))
        return "%s\tRange=%s,mu=%e\treads_in_peak_bin=%e\tp_value=%e\tbins=%s\n" % (
            self.name, str(self.iv), mu, self.reads_in_peak_bin,
            self.pvalue, str(self.clip_bin))
    
    def rescalePeak(self):
        self.relga = HTSeq.GenomicArray( ['chrI'], stranded=True)
        for s in self.ga.steps():
            relStart = int(10000 + s[0].start - int(self.center+self.halfWidth))
            if(relStart < -1):
                 next
            relEnd = int(10000 + s[0].end - int(self.center+self.halfWidth))
            if(relEnd < 0):
                 next        
            #print "start/stop= %i, %i" % (relStart, relEnd)
            if(relEnd < relStart):
                print "uh oh! backwards %i %i"
                exit
            if(abs(relEnd - relStart) > 10000):
                next
            else:
                if(relEnd == relStart):
                    r = HTSeq.GenomicPosition('chrI', relStart, s[0].strand)
                else:
                    r = HTSeq.GenomicInterval('chrI', relStart, relEnd, s[0].strand) 
                self.relga[ r ] += 1000* float(s[1])/float(self.height)

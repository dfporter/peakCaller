import sys
import os
import pybedtools
import HTSeq
from scipy.stats import poisson
import re
import subprocess

class peak:
    """A CLIP-seq peak.
    Initially set with only a genomic range and name.
    Contains methods to add reads from a bam file of CLIP reads
    and from a bam file of control reads.
    """
    def __init__(self, name, iv, height=0, center=0):
        self.name = name
        self.iv = iv
        self.height = height
        self.halfWidth = 1000  # For viewing bams
        self.range_for_bins = 500 # For the poisson and ZTNB test
        if(self.iv[3] == "+"):
            self.strandSelect = "-F 0x10"
        if(self.iv[3] == "-"):
            self.strandSelect = "-f 0x10"
        #define the 1kb range around the given peak
        #print "Created peak %s with range %s" % (str(name), str(iv))
        if(not center):
            #print "center argument not given..."
            self.center = int(float(self.iv[1] + self.iv[2])/2.0)
        else:
            #print "center given as argument..."
            self.center = center
        self.set_boundries_for_viewing_bams()
        #print "Center set at %i" % self.center
        
    def set_boundries_for_viewing_bams(self):
        if(self.center - self.halfWidth >= 0):
            self.lowerBound = self.center - self.halfWidth
        else:
            self.lowerBound = 0
        self.upperBound = self.center + self.halfWidth
            
    def add_reads_in_peak(self, clipReadsFname): #reads is a bedtool object
        """Adds all reads in the genomic interval
        of the peak. Passed a bam file name. Only reads in the specific interval,
        not in the larger self.upperBound/lowerBound or bins range.
        """
        region = "%s:%s-%s" % (self.iv[0], self.iv[1], self.iv[2])
        cmdl = ["samtools", "view", self.strandSelect,
                clipReadsFname, region]
        reads = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
        self.reads = list()
        s = reads.communicate()[0].split('\n')
        for r in s:
            # Format is:
            # read_name \t number \t chr \t locus \t 255 \t 30M
            if(r is not ""):
                self.reads.append(r)
        #print "add_reads(): Added %i reads" % len(self.reads)
                
    def findCenter(self):
        """Finds the center of the peak.
        Called when determineGA is called.
        Sets self.center and self.height based on the reads in self.ga,
        which is first set by add_reads() viewing the bam file with
        the initial self.iv range, and then later re-set by viewing the bam file
        with the self.upper/lowerBound range.
        """
        self.ga = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in self.reads:
            s = r.split('\t')
            matchLen = 0
            if(len(s)<6):
                print "Error in read ||%s|| in peak %s iv=%s" % (
                    str(r), str(self.name), str(self.iv))
                return False
                continue
            for segment in re.findall(r'(\d+)M', s[5]):
                matchLen += int(segment)
            riv = HTSeq.GenomicInterval(
                str(s[2]), int(s[3]), int(s[3]) + matchLen, self.iv[3])
                #chrom, start, stop, strand
            self.ga[riv] += 1 
        # Find the maximum of the peak
        highestValue = 0
        highestValueIV = 0
        for s in self.ga.steps():
            if(s[1] > highestValue):
                highestValue = s[1]
                highestValueIV = s[0]
        self.height = highestValue
        #print highestValueIV.start, highestValueIV.end, highestValue
        highestValueCoord = float(highestValueIV.end + highestValueIV.start)/float(2)
        self.center = int(highestValueCoord)
        self.height = highestValue
        self.set_boundries_for_viewing_bams()
        #print "findCenter set center as %i" % self.center
        #print "and height as %f" % float(self.height)
        return True

    def set_center_add_reads_in_region(self, clipReadsFname):
        """Sets and fills in the genomic array of the peak's area, reading from
        the bam file and using the self.lowerBound/upperBound range.
        Assumes self.ga exists. Finds the center of self.ga.
        Uses self.lowerBound, which was initialized as 1000 kb before
        the center of the peak.
        """
        #print "In determineGA function, making a call to findCenter..."
        if(not self.findCenter()):
            return False
        self.add_reads_in_general_region(clipReadsFname)
        return True
    
    def add_reads_in_general_region(self, clipReadsFname):
        region = "%s:%s-%s" % (self.iv[0], self.lowerBound, self.upperBound)
        cmdl = ["samtools", "view", self.strandSelect, clipReadsFname, region]
        readsRegion = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
        # Reset self.reads
        self.reads = list()
        clip_samtools_result = readsRegion.communicate()[0].split('\n')
        for r in clip_samtools_result:
            #read_name \t number \t chr \t locus \t 255 \t 30M
            if(r is not ""):
                self.reads.append(r)
        #print "determineGA: added %i reads" % len(self.reads)
        # Reset self.ga
        self.ga = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in self.reads:
            s = r.split('\t')
            matchLen = 0
            if(len(s) > 5):
                for segment in re.findall(r'(\d+)M', s[5]):
                    matchLen += int(segment)
                riv = HTSeq.GenomicInterval(
                    str(s[2]), int(s[3]), int(s[3]) + matchLen, self.iv[3])
                self.ga[riv] += 1
            else:
                print "Read with less than 5 columns in samtools output: %s" % str(r)

    def addBackgroundReads(self, backgroundReadsFname):
        region = "%s:%s-%s" % (self.iv[0], self.lowerBound, self.upperBound)
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
                for segment in re.findall(r'(\d+)M', s[5]):
                    matchLen += int(segment)
                riv = HTSeq.GenomicInterval(
                    str(s[2]), int(s[3]), int(s[3]) + matchLen, self.iv[3])
                self.background_ga[riv] += 1
            else:
                print "Read with less than 5 columns in samtools output: %s" % str(r)            

    def writeBedgraphs(self, ivBedgraph, statsIvBed):
        #format=chromA  chromStartA  chromEndA  dataValueA
        line_newRegions = "%s\t%i\t%i\t%i\n" % (self.iv[0], self.iv[1], self.iv[2], self.height)
        ivBedgraph.write(line_newRegions)
        line_statsIv = "%s\t%i\t%i\t%s\t%i\t%s\n" % (
            self.iv[0], self.lowerBound, self.upperBound,
            self.name, self.height, self.iv[3])
        statsIvBed.write(line_statsIv)

    def adjust_range(self):
        """Sets self.center and self.iv[0] and self.iv[1],
        based on the point at which the peak drops below 20% of its height.
        """
        cutoff = int(0.2 * self.height)
        rightEdge = self.center + 5
        leftEdge = self.center - 100000
        foundEdge = False
        startOfSearchRegion = self.center
        edgeOfSearchRegion = 100
        while(not foundEdge):
            k = list(self.ga[HTSeq.GenomicInterval(self.iv[0],
                                                     startOfSearchRegion,
                                                     startOfSearchRegion + edgeOfSearchRegion,
                                                     self.iv[3])].steps())
            #self.ga is the CLIP-seq genomic interval
            for j in k:
                if (j[1] <= cutoff):
                        rightEdge = j[0].start
                        foundEdge = True
                        break
            startOfSearchRegion += 100
        foundEdge = False
        startOfSearchRegion = self.center
        edgeOfSearchRegion = 100
        while(not foundEdge):
            try:
                k = list(self.ga[HTSeq.GenomicInterval(
                    self.iv[0],
                     startOfSearchRegion - edgeOfSearchRegion,
                     startOfSearchRegion,
                     self.iv[3])].steps())
            except:
                print str([self.iv[0],
                     startOfSearchRegion - edgeOfSearchRegion,
                     startOfSearchRegion,
                     self.iv[3]])
            #self.ga is the CLIP-seq genomic interval
            #can't seem to sort this list, so we'll have to cheese this
            #k = k.sort(key= lambda x: x[0].start, reverse=True)
            for j in k:
                #print "looking at %s" % str(j)
                if ((j[1] <= cutoff) and (j[0].end > leftEdge)):
                        #print "Found a new leftEdge: %s" % str(j)
                        leftEdge = j[0].end
                        foundEdge = True #stop adjusting range k.
                        #no break, look through this entire k list for the highest value of j[0].start
            startOfSearchRegion -= 100
            if(startOfSearchRegion < 1):
                leftEdge = 1
                foundEdge = True
            #except:
             #   print "Error in boundry: %s" % str(self.iv)
             #   foundEdge = True
        #print "Left edge set as %i" % leftEdge
        if(leftEdge >= 0):
            self.iv[1] = leftEdge
        else:
            self.iv[1] = 0
        self.iv[2] = rightEdge
        if(abs(self.iv[2] - self.iv[1]) < 40):
                self.iv[1] = self.center - 20
                if(self.iv[1] < 1):
                    self.iv[1] = 1
                self.iv[2] = self.center + 20
        self.center = int(float(self.iv[1] + self.iv[2])/2.0)
        self.set_boundries_for_viewing_bams()
        #print "adjust_range(): center set as %i based on range extremes of 0.2 height" % self.center
        return True

    def write_range(self):
        #.bed is chr start stop name score strand
        lineO = "%s\t%i\t%i\t%s\t%i\t%s\t%i" % (
            self.iv[0], self.iv[1], self.iv[2], str(self.name),
            int(self.height), self.iv[3], self.center)
        #print "write_range(): %s" % lineO
        return lineO

    def set_clip_bins(self):
        binSize = 20
        self.clip_bin = list()
        lowerBound = self.center - self.range_for_bins
        if(lowerBound < 1):
            lowerBound = 1
        for i in range(lowerBound, self.center + self.range_for_bins - binSize, binSize):
            #bin is i to (i+l)
            #calculate the maximimum value in the bin
            highV= 0
            k = list(self.ga[HTSeq.GenomicInterval(
                self.iv[0], i, i+binSize, self.iv[3])].steps() )
            for j in k:
                if (j[1] > highV):
                    highV = int(j[1])
            if(highV >= 0):
                self.clip_bin.append(float(highV))
        if(len(self.clip_bin) < 1):
            self.clip_bin = [0]

    def set_background_bins(self, normalCoef):
        binSize = 20
        background_bin = list()
        lowerBound = self.center - self.range_for_bins
        if(lowerBound < 1):
            lowerBound = 1
        verb = False
        if(int(self.name) == 91):
            #verb = True
            print "Setting background bins for peak 91 with iv %s" % (
                str(self.iv))
            print "Center=%i, lowerBound=%i, upperBound=%i" % (
                self.center, lowerBound, self.center + self.range_for_bins)
        for i in range(lowerBound, self.center + self.range_for_bins - binSize, binSize):
            #bin is i to (i+l)
            #calculate the maximimum value in the bin
            #if(verb):
            #    print "i=%i" % i
            highV= 0
            k = list(self.background_ga[HTSeq.GenomicInterval(
                self.iv[0], i, i+binSize, self.iv[3])].steps())
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
            
    def write_background_bins(self, normalCoef, output_zeros=False):
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
        if(output_zeros):
            background_bin = [0.0, 0.0, 0.0, 0.0, 1.0]
        else:
            background_bin = self.set_background_bins(normalCoef)
        if(len(background_bin) < 3):
            background_bin = [0.0, 0.0, 0.0, 0.0, 1.0] 
        lineO = "%s\t%s\t" % (self.name, self.height)
        for j in background_bin:
            lineO += "%f," % float(j)
        lineO += "\n"
        return lineO


    def calculate_poisson(self, return_boolean=False):
        """
        """
        #var = scipy.stats.tvar(self.clip_bin)
        #mean = scipy.stats.tmean(self.clip_bin)
        pois_dist = poisson(float(sum(self.clip_bin))/float(len(self.clip_bin)))
        self.pvalue = 1-pois_dist.cdf(self.height)
        if(return_boolean):
            if(self.pvalue < 0.001):
                return True
            else:
                return False
        mu = float(sum(self.clip_bin))/float(len(self.clip_bin))
        return "%s\tmu=%e\theight=%e\tp_value=%e\tbins=%s\n" % (self.name, mu, self.height,
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

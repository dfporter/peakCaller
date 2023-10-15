import sys
import os
import pybedtools
import HTSeq
from scipy.stats import poisson
import pysam
import re
import subprocess
from . import gtf_data

class gene:
    
    def __init__(self, name='Unkown_gene', gene_iv=['I', 1, 20, "+"],
                 gtf_info=False):
        self.name = name
        self.iv = gene_iv
        self.bin_size = 50
        self.bin_lower_bound = self.iv[1] - self.bin_size
        self.bin_upper_bound = self.iv[2] + self.bin_size
        if(self.bin_lower_bound < 1):
            self.bin_lower_bound = 1
        self.gtf_info = gtf_info
        self.find_bins_in_introns()

    def add_clip_reads_in_gene_and_bin(self, clipReadsFname):
        self.add_clip_reads_in_gene(clipReadsFname)
        self.put_clip_reads_in_bins()
        
    def add_background_reads_in_gene_and_bin(self, backgroundReadsFname):
        self.add_background_reads_in_gene(backgroundReadsFname)
        self.put_background_reads_in_bins()

    def add_reads_to_bedgraph(self, bedgraph_obj):
        if not hasattr(self, 'ga_read_starts'):
            return False
        for iv, value in self.ga_read_starts.steps():
            bedgraph_obj[iv] += value

    def add_background_reads_to_bedgraph(self, bedgraph_obj):
        if not hasattr(self, 'ga_background_read_starts'):
            return False
        for iv, value in self.ga_background_read_starts.steps():
            bedgraph_obj[iv] += value

    def add_bins_to_bedgraph(self, bedgraph_obj, which_set="clip"):
        if which_set == "clip" and (not hasattr(self, 'bins')):
            return False
        if which_set == "background" and (not hasattr(self, 'background_bins')):
            return False
        adjust_index_for_introns = 0
        for index, pos in enumerate(range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size)):
            if pos in self.bins_in_introns:
                adjust_index_for_introns += 1
                continue
            iv = HTSeq.GenomicInterval(
                self.iv[0], pos, pos+self.bin_size, self.iv[3])
            if which_set == "clip":
                try:
                    bedgraph_obj[iv] += self.bins[index-adjust_index_for_introns]
                except:
                    print("Error in gene.add_bins_to_bedgraph:")
                    print("index %i not in self.bins %s" % (index-adjust_index_for_introns, str(self.bins)))
                    print("bins in introns %s" % str(self.bins_in_introns))
            if which_set == "background":
                try:
                    bedgraph_obj[iv] += self.background_bins[index-adjust_index_for_introns]
                except:
                    print("Error in gene.add_bins_to_bedgraph:")
                    print("index %i not in self.bins %s" % (index-adjust_index_for_introns, str(self.background_bins)))
                    print("bins in introns %s" % str(self.bins_in_introns)) 

    def add_clip_reads_in_gene(self, clipReadsFname):
        """Adds all reads in the genomic interval and bins.
        """
        bamfile = pysam.AlignmentFile(clipReadsFname, "rb")
        s = bamfile.fetch(self.iv[0], int(self.iv[1]) - 200, int(self.iv[2]) + 200)
        self.reads = list()
        self.ga_read_starts = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in s:
            if(self.iv[3]=="+" and not r.is_reverse):
                r_pos = HTSeq.GenomicPosition(
                     self.iv[0], r.reference_start, self.iv[3])
                self.ga_read_starts[r_pos] += 1
            if(self.iv[3]=="-" and r.is_reverse):
               r_pos = HTSeq.GenomicPosition(
                     self.iv[0], r.reference_end-1, self.iv[3])
               self.ga_read_starts[r_pos] += 1
        bamfile.close()
        
    def find_bins_in_introns(self):
        if not self.gtf_info:
            self.bins_in_introns = []
            return True
        if self.name not in self.gtf_info.has_introns:
            self.bins_in_introns = []
            return True
        self.bins_in_introns = []
        for i in range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size):
            for intron in self.gtf_info.intron_list[self.name]:
                if (intron[0] <= i <= intron[1]) and (
                    intron[0] <= i+self.bin_size <= intron[1]):
                    self.bins_in_introns.append(i)
                    
    def put_clip_reads_in_bins(self):
        self.bins = list()
        for i in range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size):
            if i in self.bins_in_introns:
                continue
            self.bins.append(
                self.total_reads_in_bin(
                    self.ga_read_starts,
                    self.iv[0], i, i+self.bin_size, self.iv[3])
            )
        if(len(self.bins) < 1):
            self.bins = [0]

    def add_background_reads_in_gene(self, readsFname):
        """Adds all reads in the genomic interval and bins.
        """
        bamfile = pysam.AlignmentFile(readsFname, "rb")
        s = bamfile.fetch(self.iv[0], int(self.iv[1]) - 200, int(self.iv[2]) + 200)
        self.reads = list()
        self.ga_background_read_starts = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in s:
            if(self.iv[3]=="+" and not r.is_reverse):
                r_pos = HTSeq.GenomicPosition(
                     self.iv[0], r.reference_start, self.iv[3])
                self.ga_background_read_starts[r_pos] += 1
            if (self.iv[3]=="-" and r.is_reverse):
               r_pos = HTSeq.GenomicPosition(
                     self.iv[0], r.reference_end-1, self.iv[3])
               self.ga_background_read_starts[r_pos] += 1
        bamfile.close()
        self.put_background_reads_in_bins()
        
    def put_background_reads_in_bins(self):
        self.background_bins = list()
        for i in range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size):
            if i in self.bins_in_introns:
                continue
            self.background_bins.append(
                self.total_reads_in_bin(
                    self.ga_background_read_starts,
                    self.iv[0], i, i+self.bin_size, self.iv[3])
                )
        if(len(self.background_bins) < 1):
            self.background_bins = [0]
            
    def max_value_in_range_read_ends(self, ga_obj, chrm, start, end, strand):
        # Calculate the maximimum value in the bin.
        highV = 0
        k = list(ga_obj[HTSeq.GenomicInterval(
            chrm, start, end, strand)].steps() )
        for j in k:
            if (int(j[1]) > highV):
                highV = int(j[1])
        return highV
    
    def total_reads_in_bin(self, ga_obj, chrm, start, end, strand):
        # Calculate the maximimum value in the bin.
        total_reads = 0
        for iv, value in list(ga_obj[HTSeq.GenomicInterval(
            chrm, start, end, strand)].steps()):
            total_reads += (iv.end - iv.start) * (value)
        return total_reads

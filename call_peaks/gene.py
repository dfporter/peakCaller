import sys
import os
import pybedtools
import HTSeq
from scipy.stats import poisson
import pysam
import re
import subprocess

class gene:
    
    def __init__(self, name='Unkown', gene_iv=['I', 1, 2, "+"]):
        self.name = name
        self.iv = gene_iv
        self.bin_size = 50
        self.bin_lower_bound = self.iv[1] - 50
        self.bin_upper_bound = self.iv[2] + 50
        if(self.bin_lower_bound < 1):
            self.bin_lower_bound = 1

    def add_clip_reads_in_gene_and_bin(self, clipReadsFname):
        self.add_clip_reads_in_gene(clipReadsFname)
        self.put_clip_reads_in_bins()
        
    def add_background_reads_in_gene_and_bin(self, backgroundReadsFname):
        self.add_background_reads_in_gene(backgroundReadsFname)
        self.put_background_reads_in_bins()
        
    def add_clip_reads_in_gene(self, clipReadsFname):
        """Adds all reads in the genomic interval and binds.
        """
        bamfile = pysam.AlignmentFile(clipReadsFname, "rb")
        s = bamfile.fetch(self.iv[0], int(self.iv[1]) - 200, int(self.iv[2]) + 200)
        self.reads = list()
        self.ga_read_starts = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in s:
            if((self.iv[3]=="+" and not r.is_reverse) or (
                self.iv[3]=="-" and r.is_reverse)):
                r_pos = HTSeq.GenomicPosition(
                     self.iv[0], r.reference_start, self.iv[3])
                self.ga_read_starts[r_pos] += 1
        bamfile.close()
        
    def put_clip_reads_in_bins(self):
        self.bins = list()
        for i in range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size):
            highV = self.max_value_in_range_read_ends(
                self.ga_read_starts,
                self.iv[0], i, i+self.bin_size, self.iv[3])
            if(highV >= 0):
                self.bins.append(float(highV))
        if(len(self.bins) < 1):
            self.bins = [0]

    def add_background_reads_in_gene(self, readsFname):
        """Adds all reads in the genomic interval and bins.
        """
        bamfile = pysam.AlignmentFile(readsFname, "rb")
        s = bamfile.fetch(self.iv[0], int(self.iv[1]) - 200,
                          int(self.iv[2]) + 200)
        self.reads = list()
        self.ga_background_read_starts = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        for r in s:
            if((self.iv[3]=="+" and not r.is_reverse) or (
                self.iv[3]=="-" and r.is_reverse)):
                r_pos = HTSeq.GenomicPosition(
                     self.iv[0], r.reference_start, self.iv[3])
                self.ga_background_read_starts[r_pos] += 1
        bamfile.close()
        
    def put_background_reads_in_bins(self):
        self.background_bins = list()
        for i in range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size):
            highV = self.max_value_in_range_read_ends(
                self.ga_background_read_starts,
                self.iv[0], i, i+self.bin_size, self.iv[3])
            if(highV >= 0):
                self.background_bins.append(float(highV))
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

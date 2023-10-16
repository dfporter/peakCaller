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
    def __init__(self, name, iv, height=0, pos_of_peak=0, gene_name='NA',
                 gtf_info=False):
        self.name = name
        self.iv = iv
        self.height = height
        self.bin_size = 50
        self.halfWidth = 1000  # For viewing bams.
        self.range_for_bins = 500  # For the poisson and NB test.
        self.pos_of_peak = pos_of_peak
        self.gene_name = gene_name
        self.center_of_region = int(float(self.iv[1] + self.iv[2])/2.0)
        self.set_boundries_for_viewing_bams(use_pos_of_peak=False)

    def add_reads_and_adjust_range(self, clipReadsFname):
        """Find the position of the peak, its extent and the read distribution.
        Called immediately after initializing a peak to define its shape.
        """
        self.add_reads_in_peak(clipReadsFname) # Add reads under the peak from CLIP-seq data.
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
        highest_value = 0
        for s in self.ga_true_coverage.steps():
            if(int(s[1]) > highest_value):
                highest_value = int(s[1])
                highest_value_iv = s[0]
        if highest_value == 0:
            self.pos_of_peak = float(self.iv[1] + self.iv[2])/2.0
            self.height = 0
            return False
        self.pos_of_peak = int(
            float(highest_value_iv.end + highest_value_iv.start)/float(2))
        self.height = highest_value
        self.set_boundries_for_viewing_bams(use_pos_of_peak=True)
        return True

    def set_center_add_reads_in_region(self, clipReadsFname):
        """Sets and fills in the genomic array of the peak's area, reading from
        the bam file and using the self.lower_bound/upperBound range.
        Assumes self.ga exists. Finds the peak position of self.ga.
        Uses self.lower_bound, which was initialized as 1000 kb before
        the center of the peak.
        """
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
        for r in s:
            if(self.iv[3]=="+" and not r.is_reverse):
                riv = HTSeq.GenomicInterval(
                     self.iv[0], r.reference_start, r.reference_end, self.iv[3])
                self.ga_true_coverage[riv] += 1
                r_pos = HTSeq.GenomicPosition(
                    self.iv[0], r.reference_start, self.iv[3])
                self.ga_read_starts[r_pos] += 1
            if(self.iv[3]=="-" and r.is_reverse):
                r_pos = HTSeq.GenomicPosition(
                    self.iv[0], r.reference_end-1, self.iv[3])
                self.ga_read_starts[r_pos] += 1                
        bamfile.close()

    def add_reads_to_bedgraph(self, bedgraph_obj):
        for iv, value in self.ga_read_starts.steps():
            bedgraph_obj[iv] += value
            
    def add_background_reads_to_bedgraph(self, bedgraph_obj):
        if not hasattr(self, 'ga_background_read_starts'):
            return False
        for iv, value in self.ga_background_read_starts.steps():
            bedgraph_obj[iv] += value

    def add_bins_to_bedgraph(self, bedgraph_obj):
        if not hasattr(self, 'leftmost_bin_index'):
            print("Error: Asked to report bins, but no bins are set for %s." % self.name)
        for i in range(self.leftmost_bin_index, self.rightmost_bin_index + 1, 1):
            _bin = self.clip_bins[i]
            iv = HTSeq.GenomicInterval(*_bin['iv'])
            bedgraph_obj[iv] += _bin['signal']
            
    def addBackgroundReads(self, backgroundReadsFname):
        self.ga_background_read_starts = HTSeq.GenomicArray([self.iv[0]], stranded=True)
        bamfile = pysam.AlignmentFile(backgroundReadsFname, "rb")
        s = bamfile.fetch(self.iv[0], self.lower_bound, self.upper_bound)
        for r in s:
            if(self.iv[3]=="+" and not r.is_reverse):
                r_pos = HTSeq.GenomicPosition(
                    self.iv[0], r.reference_start, self.iv[3])
                self.ga_background_read_starts[r_pos] += 1
            if(self.iv[3]=="-" and r.is_reverse):
                r_pos = HTSeq.GenomicPosition(
                    self.iv[0], r.reference_end-1, self.iv[3])
                self.ga_background_read_starts[r_pos] += 1
        bamfile.close()
        
    def write_bedgraphs(self, peak_iv_bedgraph, region_for_stats_bed):
        """Write the peak boundaries as .bedgraph and region used for
        stats as .bed. Called only during define_ranges().
        """
        #format=chromA  chromStartA  chromEndA  dataValueA
        line_peak_boundaries = "%s\t%i\t%i\t%i\n" % (self.iv[0], self.iv[1], self.iv[2], self.height)
        peak_iv_bedgraph.write(line_peak_boundaries)
        line_stats_iv = "%s\t%i\t%i\t%s\t%i\t%s\n" % (
            self.iv[0], self.lower_bound, self.upper_bound,
            self.name, self.height, self.iv[3])
        region_for_stats_bed.write(line_stats_iv)

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
        return True

    def write_range(self):
        """Return a .bed formatted line of the range of the peak.
        .bed format is chr start stop name score strand
        """
        lineO = "%s\t%i\t%i\t%s\t%i\t%s\t%i" % (
            self.iv[0], self.iv[1], self.iv[2], str(self.name),
            int(self.height), self.iv[3], self.pos_of_peak)
        return lineO

    def write_as_peaks_format(self, multiply_p_values=False):
        line = "{chrm}\t{start}\t{end}\t{name}\t{height}\t{strand}".format(
            chrm=self.iv[0], start=self.iv[1], end=self.iv[2],
            name=self.name, height=self.height, strand=self.iv[3])
        if multiply_p_values:
            line += "\t{clip_poisson}\t{nb_background}\t{nb_clip}".format(
                clip_poisson=min(1.0, self.poisson_pvalue_by_gene),
                nb_background=min(
                    1.0, float(self.background_nb_p_value * multiply_p_values)),
                nb_clip=min(1.0, float(self.clip_nb_p_value * multiply_p_values)))
        else:
            line += "\t{clip_poisson}\t{nb_background}\t{nb_clip}".format(
                clip_poisson=min(1.0, self.poisson_pvalue_by_gene),
                nb_background=min(1.0, self.background_nb_p_value),
                nb_clip=min(1.0, self.clip_nb_p_value))
            
        line += "\t%s\n" % self.gene_name
        return line
    
    def set_bin_locations(self, use_bins_from_gene=False):
        """Set bin locations so they can be used for statistics.
        Set the genomic positions of the left edges of the leftmost
        and rightmost bin.
        """
        #use_bins_from_gene = False
        if use_bins_from_gene:
            self.offset = self.pos_of_peak - use_bins_from_gene.bin_lower_bound
            self.offset = self.offset % self.bin_size  # Genes and peaks must have the same bin size.
        else:
            self.offset = int(float(self.bin_size)/2.0) 
        self.num_bins = int(float(2 * self.range_for_bins)/float(self.bin_size))
        self.bin_lower_bound = self.pos_of_peak - self.range_for_bins - self.offset
        if(self.bin_lower_bound < 1):
            self.bin_lower_bound = 1
        self.bin_upper_bound = self.bin_lower_bound + (2 * self.range_for_bins)

    def max_value_in_range_read_ends(self, chrm, start, end, strand):
        """Calculate the maximimum value in the bin.
        Not used, I think.
        """
        highV = 0
        k = list(self.ga_read_starts[HTSeq.GenomicInterval(
            chrm, start, end, strand)].steps() )
        for j in k:
            if (int(j[1]) > highV):
                highV = int(j[1])
        return highV

    def put_clip_reads_in_bins(self, use_bins_from_gene=False):
        self.set_bin_locations(use_bins_from_gene=use_bins_from_gene)
        self.clip_bins = list()
        for i in range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size):
            iv = (self.iv[0], i, i + self.bin_size, self.iv[3])
            value_in_bin = self.total_reads_in_bin(self.ga_read_starts, *iv)
            self.clip_bins.append({'iv': iv, 'signal': value_in_bin})
        if(len(self.clip_bins) < 1):
            self.clip_bins = [{'iv': self.iv, 'signal': 0.0}]
            print("No clip bins for peak %s?" % self.name)
            
    def total_reads_in_bin(self, ga_obj, chrm, start, end, strand):
        """Sum all reads in a given bin.
        """
        total_reads = 0
        for iv, value in list(ga_obj[HTSeq.GenomicInterval(
            chrm, start, end, strand)].steps()):
            total_reads += (iv.end - iv.start) * (value)
        return total_reads
    
    def put_background_reads_in_bins(self):
        self.background_bins = list()
        verb = False
        for i in range(self.bin_lower_bound,
                       self.bin_upper_bound,
                       self.bin_size):
            iv = (self.iv[0], i, i + self.bin_size, self.iv[3])
            value_in_bin = self.total_reads_in_bin(self.ga_background_read_starts, *iv)
            self.background_bins.append({'iv': iv, 'signal': value_in_bin})
        if(len(self.background_bins) < 1):
            self.background_bins = [{'iv': self.iv, 'signal': 0}]
            print("No background bins for peak %s?" % self.name)
        return self.background_bins
    
    def process_bins(self, bin_list):
        nonzero_bins = list()
        for _bin in bin_list:
            if _bin > 0:
                nonzero_bins.append(_bin)
        if len(nonzero_bins) > 0:
            bin_list = nonzero_bins
        return bin_list
    
    def compose_line_for_R(self, bin_list, normal_coef):
        line = "%s\t%s\t" % (self.name, self.reads_in_peak_bin)
        for j in bin_list:
            line += "%f," % float(j * normal_coef)
        line += "\n"
        return line
    
    def set_line_for_R(self, normalCoef=1.0, a_binned_gene=False,
                              output_zeros=False):
        """Write the background for R to model.
        If there is no background, make sure there is something
        (a near-zero background) for R to model so R doesn't crash.
        TO-DO: a better solution.
        """
        if a_binned_gene:
            # NB vs background.
            background_bin = a_binned_gene.background_bins
            # Remove bins with zero signal (probably not real bins).
            background_bin = self.process_bins(background_bin)
            self.background_line_for_R = self.compose_line_for_R(
                background_bin, normalCoef)
            # NB vs CLIP signal in gene.
            clip_in_gene_bins = a_binned_gene.bins
            # We allow zero-signal bins for CLIP.
            self.clip_line_for_R = self.compose_line_for_R(
                clip_in_gene_bins, 1.0)
            return self.background_line_for_R
        # Not using information from a gene, just local region.
        if(output_zeros):
            background_bin = [0.0, 0.0, 0.0, 0.0, 1.0]
        else:
            background_bin = self.put_background_reads_in_bins()
            background_bin = [normalCoef * x['signal'] for x in background_bin]
        if(len(background_bin) < 3):
            background_bin = [0.0, 0.0, 0.0, 0.0, 1.0]
        self.background_line_for_R = self.compose_line_for_R(
            background_bin, 1.0)
        self.clip_line_for_R = self.compose_line_for_R(
            [_bin['signal'] for _bin in self.clip_bins], 1.0)
        return self.background_line_for_R

    def find_reads_in_peak_bin(self):
        """Reads in bins are counted from the 5' end, while the center of the peak is counted
        as the true max height of the peak. So counting the "height" of the peak
        as the number of reads in the center bin may give a lower number than to the side.
        We will therefore allow some flexibility in where the "center bin" is.
        """
        left_bins = []
        right_bins = []
        # Set some defaults in case the peak is so long we can't find overlapping borders.
        self.leftmost_bin = sorted(self.clip_bins, key=lambda x: x['iv'][1])[0]
        self.rightmost_bin = sorted(self.clip_bins, key=lambda x: x['iv'][1])[-1]
        # What bins overlap the peak range?
        for _bin in self.clip_bins:
            if(_bin['iv'][1] <= self.iv[1] <= _bin['iv'][2]):
                # Left peak boundary overlaps this bin.
                left_bins.append(_bin)
            if(_bin['iv'][1] <= self.iv[2] <= _bin['iv'][2]):
                # Right peak boundary overlaps this bin.
                right_bins.append(_bin)
        if len(left_bins) > 0:
             self.leftmost_bin = sorted(left_bins, key=lambda x: x['iv'][1])[0]
        if len(right_bins) > 0:
            self.rightmost_bin = sorted(right_bins, key=lambda x: x['iv'][1])[-1]
        self.leftmost_bin_index = self.clip_bins.index(self.leftmost_bin)
        self.rightmost_bin_index = self.clip_bins.index(self.rightmost_bin)
        if self.leftmost_bin_index == self.rightmost_bin_index:
            if (self.leftmost_bin_index + 1) > (len(self.clip_bins) - 1):
                self.rightmost_bin_index = self.leftmost_bin_index
                print("Error: only one bin? Peak %i. Bins %s." % (self.name, str(self.clip_bins)))
            else:
                self.rightmost_bin_index = self.leftmost_bin_index + 1
        values_in_bins = [_bin['signal'] for _bin in self.clip_bins]
        for i in range(self.leftmost_bin_index, self.rightmost_bin_index + 1, 1):
            values_in_bins.append(self.clip_bins[i]['signal'])
        self.reads_in_peak_bin = max(values_in_bins)
        
    def calculate_poisson(self, a_binned_gene=False):
        """Calculate a p value by Poisson.
        """
        if a_binned_gene:
            gene_bins = a_binned_gene.bins
            #nonzero_bins = []
            #for _bin in gene_bins:
            #    if _bin > 0:
            #        nonzero_bins.append(_bin)
            #if len(nonzero_bins) > 0:
            #    gene_bins = nonzero_bins
            mu = float(sum(gene_bins))/float(len(gene_bins))
            pois_dist = poisson(mu)
            self.find_reads_in_peak_bin()
            self.poisson_pvalue = 1 - pois_dist.cdf(self.reads_in_peak_bin)
            # Correct for multiple hypothesis testing - that is,
            # get the p value for obtaining a significant bin somewhere in
            # the given gene. We assume this is exon + 200 bases, as an approximation.
            # An exact value is not used because an exact transcript length
            # is not generally known.
            self.poisson_pvalue_by_gene = self.poisson_pvalue * float(len(gene_bins) + 2)
            line = "%s\tRange=%s,mu=%e\treads_in_peak_bin=%e\tp_value=%e\tbins_around_peak=%s\n" % (
                self.name, str(self.iv), mu, self.reads_in_peak_bin,
                self.poisson_pvalue, str(self.clip_bins))
            return line
        total_across_bins = float(sum([_bin['signal'] for _bin in self.clip_bins]))
        mu = total_across_bins/float(len(self.clip_bins)) 
        pois_dist = poisson(mu)
        self.find_reads_in_peak_bin()
        self.poisson_pvalue = 1-pois_dist.cdf(self.reads_in_peak_bin)
        # No gene length information; use the length of the binned region.
        self.poisson_pvalue_by_gene = self.poisson_pvalue * float(self.num_bins)
        return "%s\tRange=%s,mu=%e\treads_in_peak_bin=%e\tp_value=%e\tbins=%s\n" % (
            self.name, str(self.iv), mu, self.reads_in_peak_bin,
            self.poisson_pvalue, str(self.clip_bins))

    def remove_read_info(self):
        del self.ga_true_coverage
        del self.ga_read_starts
        if hasattr(self, 'ga_background_read_starts'):
            del self.ga_background_read_starts

"""
Creates the annotation .bed file from a given .gtf needed for call_peaks.

.gtf files are 1-based.
.bed files are 0-based.
"""
import sys
import os
import csv
import argparse
import re


def create_bed_from_gtf(gtf_filename, bed_filename):
    exons_by_txpt = {}
    genes = {}
    with open(gtf_filename, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            m = re.search("transcript_id ([^;]+);", li)
            if m is None:
                print "Error on line %s" % li
            txpt_id = re.sub('''"''', '', m.group(1))
            iv = (s[0], int(s[3]), int(s[4]), s[6])
            if s[2] == 'CDS':
                genes.setdefault(txpt_id, {})
                genes[txpt_id].setdefault('CDS', [])
                genes[txpt_id]['CDS'].append(iv)
            if s[2] == "exon":
                genes.setdefault(txpt_id, {})
                genes[txpt_id].setdefault('exon', [])
                genes[txpt_id]['exon'].append(iv)
    for txpt in genes:
        genes[txpt]['longest'] = 0
        if 'CDS' in genes[txpt]:
            starts = [x[1] for x in genes[txpt]['CDS']]
            ends = [x[2] for x in genes[txpt]['CDS']]
            try:
                max_range = (genes[txpt]['CDS'][0][0], min(starts), max(ends),
                             genes[txpt]['CDS'][0][3])
                genes[txpt]['max_cds'] = max_range
            except:
                print "error with %s" % str(genes[txpt])
        if 'exon' in genes[txpt]:
            starts = [x[1] for x in genes[txpt]['exon']]
            ends = [x[2] for x in genes[txpt]['exon']]
            max_range = (genes[txpt]['exon'][0][0], min(starts), max(ends),
                         genes[txpt]['exon'][0][3])
            genes[txpt]['max_exon'] = max_range
        if 'max_cds' in genes[txpt] and ('max_exon' in genes[txpt]):
            if genes[txpt]['max_cds'] != genes[txpt]['max_exon'] and (
                (genes[txpt]['max_cds'][2] - genes[txpt]['max_cds'][1] + 3) != (
                    genes[txpt]['max_exon'][2] - genes[txpt]['max_exon'][1])):
                print "Error: CDS and exon lengths don't match for %s" % txpt
                print genes[txpt]
    with open(bed_filename, 'w') as f:
        for txpt in genes:
            if 'max_exon' not in genes[txpt]:
                print "Missing exon lengths for %s" % txpt
                continue
            iv = genes[txpt]['max_exon']
            f.write("{chrm}\t{start}\t{end}\t{name}\t1\t{strand}\n".format(
                chrm=iv[0], start=str(iv[1] - 1), end=str(iv[2] - 1),
                name="ID=%s;" % txpt, strand=iv[3]))


if __name__ == "__main__":
    gtf_filename = sys.argv[1]
    bed_filename = sys.argv[2]
    create_bed_from_gtf(gtf_filename, bed_filename)
    

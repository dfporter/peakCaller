import sys
import os
import csv
import argparse
import re


if __name__ == "__main__":
    gtf_filename = "/home/dp/Desktop/ensemblEF4/Saccharomyces_cerevisiae.EF4.70.gtf"
    exons_by_txpt = {}
    with open(gtf_filename, 'r') as f:
        for li in gtf_filename:
            s = li.rstrip('\n').split('\t')
            m = re.search("transcript_id ([^;]+);", li)
            if m is None:
                print "Error on line %s" % li
            txpt_id = m.group(1)
            if s[2] == "exon":
                exons_by_txpt.setdefault(txpt_id, set())
                exons_by_txpt[txpt_id].add({
                    'start': int(s[3]),
                    'end': int(s[4]),
                    'strand': s[6]})

    with open("exons_by_transcript", "w") as f:
        for name in exons_by_txpt:
            line = name
            for exon in exons_by_txpt[name]:
                line += ""

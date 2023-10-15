import re


class gtf_data:

    def __init__(self, gtf_filename='lib/Saccharomyces_cerevisiae.EF4.70.gtf'):
        self.has_introns = set()
        self.intron_list = {}
        self.txpt_ranges = {}
        self.read_file(gtf_filename)
        
    def read_file(self, gtf_filename):
        # GTF files are 1-based. HTSeq genomic arrays are 0-based.
        # We move all coordinates back by one so they will
        # be treated as 0-based.
        with open(gtf_filename, 'r') as f:
            self.exons_by_txpt = {}
            for li in f:
                s = li.rstrip('\n').split('\t')
                m = re.search("transcript_id ([^;]+);", li)
                if m is None:
                    print("Error on line %s" % li)
                    continue
                txpt_id = re.sub(r'''"''', '', m.group(1))
                if s[2] == "exon":
                    this_exon = {
                        'start': int(s[3]) - 1,
                        'end': int(s[4]) - 1,
                        'strand': s[6]}
                    self.exons_by_txpt.setdefault(txpt_id, [])
                    self.exons_by_txpt[txpt_id].append(this_exon)
        for txpt in self.exons_by_txpt:
            if len(self.exons_by_txpt[txpt]) > 1:
                self.has_introns.add(txpt)
                self.intron_list[txpt] = []
                last_end = 0 
                for exon in sorted(self.exons_by_txpt[txpt],
                                   key=lambda x: x['start']):
                    if last_end:
                        an_intron = [last_end, exon['start']]
                        self.intron_list[txpt].append(an_intron)
                    last_end = exon['end']
            first_exon_start = min(
                [x['start'] for x in self.exons_by_txpt[txpt]])
            last_exon_end = max(
                [x['end'] for x in self.exons_by_txpt[txpt]])
            self.txpt_ranges[txpt] = [first_exon_start, last_exon_end]

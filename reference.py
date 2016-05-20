import pysam

# Class representing a reference genome build
class Reference(object):

    # Constructor
    def __init__(self, fastafile):
        self.ref_tabix = pysam.Fastafile(fastafile)

    # Read reference sequence of a given genomic region
    def read(self, *arg):
        if len(arg) == 3: chrom, start, end = arg[0], arg[1], arg[2]
        elif len(arg) == 1:
            chrom, pos = arg[0].split(':')
            if '-' in pos: [start, end] = pos.split('-')
            else: start, end = pos, pos
        else: return None
        chrom, start, end = str(chrom), int(start), int(end)

        goodchrom = chrom
        if not goodchrom in self.ref_tabix.references:
            goodchrom = 'chr' + chrom
            if not goodchrom in self.ref_tabix.references:
                if chrom == 'MT':
                    goodchrom = 'chrM'
                    if goodchrom not in self.ref_tabix.references: return None
                else: return None

        if end < start: return ''
        if start < 0: start = 1

        if pysam.__version__ in ['0.7.7', '0.7.8', '0.8.0']: last = self.ref_tabix.getReferenceLength(goodchrom)
        else: last = self.ref_tabix.get_reference_length(goodchrom)

        if end > last: end = last
        seq = self.ref_tabix.fetch(goodchrom, start - 1, end)
        return str(seq.upper().decode("utf-8"))





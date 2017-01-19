import gzip
import reference
import sys
import os
from optparse import OptionParser
import datetime


#######################################################################################################################

# Class representing an exon sequence
class ExonSequence(object):

    # Constructor
    def __init__(self, index, seq):
        self.index = index
        self.seq = seq


# Class representing a single Ensembl transcript
class EnsemblTranscript(object):

    # Constructor
    def __init__(self, line):
        self.exons = []
        cols = line.split('\t')
        self.ID = cols[0]
        self.gene = cols[1]
        self.geneID = cols[3]
        self.chrom = cols[5]
        self.strand = int(cols[6])
        self.transcriptStart = int(cols[7])
        self.transcriptEnd = int(cols[8])
        self.codingStart = int(cols[9])
        self.codingStartGenomic = int(cols[10])
        self.codingEndGenomic = int(cols[11])
        # Initializing and adding exons
        for i in range(1, len(cols) - 12, 2):
            self.exons.append(EnsemblExon(int((i + 1) / 2), int(cols[11 + i]), int(cols[12 + i])))

    # Return CDS, UTR5 and UTR3 dictionaries (keys: exon index, value: sequence)
    def getSeqs(self, ref):
        cds_seqs = []
        utr5_seqs = []
        utr3_seqs = []
        # int_seqs = dict()

        for i in range(len(self.exons)):
            exon = self.exons[i]

            if self.strand == 1:
                if self.codingStartGenomic > exon.end:
                    utr5_seqs.append(ExonSequence(i+1, ref.read(self.chrom, exon.start+1, exon.end)))
                    continue
                if self.codingEndGenomic < exon.start+1:
                    utr3_seqs.append(ExonSequence(i+1, ref.read(self.chrom, exon.start+1, exon.end)))
                    continue
                if exon.start+1 <= self.codingStartGenomic <= exon.end:
                    utr5_seqs.append(ExonSequence(i+1,ref.read(self.chrom, exon.start+1, self.codingStartGenomic-1)))
                    cdsstart = self.codingStartGenomic-1
                else:
                    cdsstart = exon.start
                if exon.start+1 <= self.codingEndGenomic <= exon.end:
                    utr3_seqs.append(ExonSequence(i+1, ref.read(self.chrom, self.codingEndGenomic+1, exon.end)))
                    cdsend = self.codingEndGenomic
                else:
                    cdsend = exon.end
                cds_seqs.append(ExonSequence(i+1, ref.read(self.chrom, cdsstart+1, cdsend)))

            else:
                if self.codingStartGenomic < exon.start+1:
                    utr5_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, exon.start+1, exon.end))))
                    continue
                if self.codingEndGenomic > exon.end:
                    utr3_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, exon.start+1, exon.end))))
                    continue
                if exon.start+1 <= self.codingEndGenomic <= exon.end:
                    utr3_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, exon.start+1, self.codingEndGenomic-1))))
                    cdsstart = self.codingEndGenomic-1
                else:
                    cdsstart = exon.start
                if exon.start+1 <= self.codingStartGenomic <= exon.end:
                    utr5_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, self.codingStartGenomic+1, exon.end))))
                    cdsend = self.codingStartGenomic
                else:
                    cdsend = exon.end
                cds_seqs.append(ExonSequence(i+1, reverseComplement(ref.read(self.chrom, cdsstart+1, cdsend))))

        return cds_seqs, utr5_seqs, utr3_seqs


# Class representing a single exon of an Ensembl transcript
class EnsemblExon(object):

    # Constructor
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start


#######################################################################################################################


# Class representing a single Ensembl transcript
class RefSeqTranscript(object):

    # Constructor
    def __init__(self, line):
        cols = line.split('\t')
        self.ID = cols[0]
        self.version = cols[1]
        self.gene = cols[2]
        self.exons = []
        for exonstr in cols[3].split(','): self.exons.append(exonstr.split('-'))
        self.cds = cols[4].split('-')
        self.seq = cols[5]

    # Return CDS, UTR5 and UTR3 dictionaries (keys: exon index, value: sequence)
    def getSeqs(self):
        cds_seqs = []
        utr5_seqs = []
        utr3_seqs = []

        cds_start = int(self.cds[0])
        cds_end = int(self.cds[1])

        for i in range(1,len(self.exons)+1):
            exon = self.exons[i-1]
            if not int(exon[1]) < cds_start and not int(exon[0]) > cds_end:
                if int(exon[0]) < cds_start: start = cds_start
                else: start = int(exon[0])
                if int(exon[1]) > cds_end: end = cds_end
                else: end = int(exon[1])
                cds_seqs.append(ExonSequence(i, self.seq[start-1:end]))

            if not int(exon[0]) > cds_start:
                if cds_start < int(exon[1]): end = cds_start-1
                else: end = int(exon[1])
                utr5_seqs.append(ExonSequence(i, self.seq[int(exon[0])-1:end]))

            if not int(exon[1]) < cds_end:
                if cds_end > int(exon[0]): start = cds_end+1
                else: start = int(exon[0])
                utr3_seqs.append(ExonSequence(i, self.seq[start-1:int(exon[1])]))

        return cds_seqs, utr5_seqs, utr3_seqs


#######################################################################################################################


# Return reverse complement sequence
def reverseComplement(seq):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c", "n": "n", '-': '-'}
    ret = ''
    for base in seq[::-1]: ret += complement[base]
    return ret

# Read Ensembl transcript database
def readEnsemblData(filename, idkey):
    ret = dict()
    n_transcripts = 0
    for line in gzip.open(filename,'r'):
        line = line.strip()
        if line=='' or line.startswith('#'): continue
        t = EnsemblTranscript(line)
        n_transcripts += 1
        if idkey:
            ret[t.ID] = t
        else:
            if t.gene not in ret.keys(): ret[t.gene] = []
            ret[t.gene].append(t)
    return ret, n_transcripts

# Read RefSeq transcript database
def readRefSeqData(filename, idkey):
    ret = dict()
    n_transcripts = 0
    for line in gzip.open(filename,'r'):
        line = line.strip()
        if line=='' or line.startswith('#'): continue
        t = RefSeqTranscript(line)
        n_transcripts += 1
        if idkey:
            ret[t.ID] = t
        else:
            if t.gene not in ret.keys(): ret[t.gene] = []
            ret[t.gene].append(t)
    return ret, n_transcripts

# Read configuration file
def readConfigFile(dir):
    ret37 = None
    ret38 = None
    for line in open(dir+'/config.txt'):
        line = line.strip()
        if line == '': continue
        [key, val] = line.split('=')
        key = key.strip()
        val = val.strip()
        if key == 'GRCh37': ret37 = val
        if key == 'GRCh38': ret38 = val
    return ret37, ret38

# Read Ensembl info txt file
def readEnsemblInfo(path):
    for line in open(path[:-3]+'.txt'):
        line = line.strip()
        if line == '': continue
        if line[0] == '#':
            return line[line.find('GRCh3'):line.find('GRCh3')+6], line[line.find('Ensembl release'):]

# Read RefSeq info txt file
def readRefSeqInfo(filename):
    for line in gzip.open(filename,'r'):
        line = line.strip()
        return 'RefSeq data (downloaded: '+line[line.find('downloaded:')+12:]+')'

# Check type of data file (i.e. Ensembl or RefSeq)
def checkTypeOfDataFile(filename):
    for line in gzip.open(filename,'r'):
        line = line.strip()
        if line=='' or line.startswith('#'): continue
        if line.startswith('ENST'): return 'Ensembl'
        if line.startswith('NM_'): return 'RefSeq'
        return None

# Check type of input file (i.e. Ensembl or RefSeq)
def checkTypeOfInputFile(filename):
    for line in open(filename):
        line = line.strip()
        if line=='' or line.startswith('#'): continue
        if line.startswith('ENST'): return 'Ensembl'
        if line.startswith('NM_'): return 'RefSeq'
        return None

# Compare query transcript with a set of transcripts
def compare(transcript, comparewith, dataXtype, dataYtype, dataX_build, dataY_build, ref_GRCh37, ref_GRCh38, options, out):
    identical = False
    cds_identical = False
    for x in comparewith:
        flags = []
        if dataXtype == 'Ensembl' and dataYtype == 'Ensembl': flags = compareEnsemblVsEnsembl(transcript, x, dataX_build, dataY_build, ref_GRCh37, ref_GRCh38, options)
        elif dataXtype == 'Ensembl' and dataYtype == 'RefSeq': flags = compareEnsemblVsRefSeq(transcript, x, dataX_build, dataY_build, dataXtype, options)
        elif dataXtype == 'RefSeq' and dataYtype == 'Ensembl': flags = compareEnsemblVsRefSeq(x, transcript, dataX_build, dataY_build, dataXtype, options)

        if all([('CDS' not in flag) for flag in flags]): cds_identical = True

        if len(flags) == 0:
            identical = True
            if not options.discr: flags.append('.')
            else: continue

        out.write('\t'.join([transcript.ID, x.ID, ';'.join(flags)])+'\n')

    return identical, cds_identical

# Compare two Ensembl transcripts
def compareEnsemblVsEnsembl(enst_transcript1, enst_transcript2, dataX_build, dataY_build, ref_GRCh37, ref_GRCh38, options):
    flags = []

    # Extract CDS, UTR5 and UTR3 sequences for both releases using the appropriate genome build
    if dataX_build == 'GRCh37': cds_1, utr5_1, utr3_1 = enst_transcript1.getSeqs(ref_GRCh37)
    else: cds_1, utr5_1, utr3_1 = enst_transcript1.getSeqs(ref_GRCh38)
    if dataY_build == 'GRCh37': cds_2, utr5_2, utr3_2 = enst_transcript2.getSeqs(ref_GRCh37)
    else: cds_2, utr5_2, utr3_2 = enst_transcript2.getSeqs(ref_GRCh38)

    # Adding CDS flags
    if not len(cds_1) == len(cds_2):
        if options.simple: flags.append('CDS')
        else: flags.append('CDS:DN')
    else:
        discrv = []
        for i in range(len(cds_1)):
            if not cds_1[i].seq == cds_2[i].seq:
                discrv.append(str(cds_1[i].index))
        if len(discrv) > 0:
            if options.simple: flags.append('CDS')
            else: flags.append('CDS:'+','.join(discrv))

    # Adding UTR5 flags
    if not len(utr5_1) == len(utr5_2):
        if options.simple: flags.append('UTR5')
        else: flags.append('UTR5:DN')
    else:
        discrv = []
        for i in range(len(utr5_1)):
            if not utr5_1[i].seq == utr5_2[i].seq:
                discrv.append(str(utr5_1[i].index))
        if len(discrv) > 0:
            if options.simple: flags.append('UTR5')
            else: flags.append('UTR5:'+','.join(discrv))

    # Adding UTR3 flags
    if not len(utr3_1) == len(utr3_2):
        if options.simple: flags.append('UTR3')
        else: flags.append('UTR3:DN')
    else:
        discrv = []
        for i in range(len(utr3_1)):
            if not utr3_1[i].seq == utr3_2[i].seq:
                discrv.append(str(utr3_1[i].index))
        if len(discrv) > 0:
            if options.simple: flags.append('UTR3')
            else: flags.append('UTR3:'+','.join(discrv))

    return flags

# Compare an Ensembl and a RefSeq transcript
def compareEnsemblVsRefSeq(enst_transcript, refseq_transcript, dataX_build, dataY_build, dataXtype, options):
    flags = []

    # Extract CDS, UTR5 and UTR3 sequences for the Ensembl dataset
    if dataX_build == 'GRCh37' or dataY_build == 'GRCh37': cds_enst, utr5_enst, utr3_enst = enst_transcript.getSeqs(ref_GRCh37)
    else: cds_enst, utr5_enst, utr3_enst = enst_transcript.getSeqs(ref_GRCh38)

    # Extract CDS, UTR5 and UTR3 sequences for the RefSeq dataset
    cds_refseq, utr5_refseq, utr3_refseq = refseq_transcript.getSeqs()

    # Adding CDS flags
    if not len(cds_enst) == len(cds_refseq):
        if options.simple: flags.append('CDS')
        else: flags.append('CDS:DN')
    else:
        discrv = []
        for i in range(len(cds_enst)):
            if not cds_enst[i].seq == cds_refseq[i].seq:
                if dataXtype == 'Ensembl': idx = cds_enst[i].index
                else: idx = cds_refseq[i].index
                discrv.append(str(idx))
        if len(discrv) > 0:
            if options.simple: flags.append('CDS')
            else: flags.append('CDS:'+','.join(discrv))

    # Adding UTR5 flags
    if not len(utr5_enst) == len(utr5_refseq):
        if options.simple: flags.append('UTR5')
        else: flags.append('UTR5:DN')
    else:
        discrv = []
        for i in range(len(utr5_enst)):
            if not utr5_enst[i].seq == utr5_refseq[i].seq:
                if dataXtype == 'Ensembl': idx = utr5_enst[i].index
                else: idx = utr5_refseq[i].index
                discrv.append(str(idx))
        if len(discrv) > 0:
            if options.simple: flags.append('UTR5')
            else: flags.append('UTR5:'+','.join(discrv))

    # Adding UTR3 flags
    if not len(utr3_enst) == len(utr3_refseq):
        if options.simple: flags.append('UTR3')
        else: flags.append('UTR3:DN')
    else:
        discrv = []
        for i in range(len(utr3_enst)):
            if not utr3_enst[i].seq == utr3_refseq[i].seq:
                if dataXtype == 'Ensembl': idx = utr3_enst[i].index
                else: idx = utr3_refseq[i].index
                discrv.append(str(idx))
        if len(discrv) > 0:
            if options.simple: flags.append('UTR3')
            else: flags.append('UTR3:'+','.join(discrv))

    return flags


#######################################################################################################################


# Version
ver = '1.1'

# Directory of transcript_map
dir = os.path.dirname(os.path.realpath(__file__))

# Check if any command line arguments
if not len(sys.argv) > 1:
    print '\nUsage: python path/to/transcript_map/transcript_map.py <options>\n'
    quit()

# Command line argument parsing
descr = 'transcript_map v'+ver
parser = OptionParser(usage='python path/to/transcript_map/transcript_map.py <options>', description=descr)
parser.add_option('-x', default=None, dest='dataX', action='store', help="Ensembl database file, release 1")
parser.add_option('-y', default=None, dest='dataY', action='store', help="Ensembl database file, release 2")
parser.add_option('-i', default=None, dest='input', action='store', help="Input file containing ENST IDs of interest")
parser.add_option('-o', default='output.txt', dest='output', action='store', help="Output file name")
parser.add_option('-d',  default=False, dest='discr', action='store_true', help="Show only discrepant transcripts")
parser.add_option('-s',  default=False, dest='simple', action='store_true', help="Simple output")
(options, args) = parser.parse_args()

# Welcome message
print '\n'+'-'*100
print 'transcript_map v'+ver+' started: '+str(datetime.datetime.now())+'\n'

# Check if input and data files exist
if not os.path.isfile(options.input): print 'Error: Input file ('+options.input+') cannot be found.\n'; quit()
if not os.path.isfile(options.dataX): print 'Error: Dataset X file ('+options.dataX+') cannot be found.\n'; quit()
if not os.path.isfile(options.dataY): print 'Error: Dataset Y file ('+options.dataY+') cannot be found.\n'; quit()

# Check database file types (i.e. Ensembl or RefSeq)
dataXtype = checkTypeOfDataFile(options.dataX)
dataYtype = checkTypeOfDataFile(options.dataY)

# Checks if Ensembl TXT files exist if required
if dataXtype=='Ensembl' and not os.path.isfile(options.dataX[:-3]+'.txt'): print 'Error: Dataset X txt file ('+options.dataX[:-3]+'.txt) cannot be found.\n'; quit()
if dataYtype=='Ensembl' and not os.path.isfile(options.dataY[:-3]+'.txt'): print 'Error: Dataset Y txt file ('+options.dataY[:-3]+'.txt) cannot be found.\n'; quit()

# Check input file type (i.e. Ensembl or RefSeq)
inputtype = checkTypeOfInputFile(options.input)

# Read input file
inputlist = []
for line in open(options.input): inputlist.append(line.strip())
print 'Input list ['+options.input+'] contains '+str(len(inputlist))+' transcripts\n'

# Check if data file and input file types allowed
if dataXtype == 'RefSeq' and dataYtype == 'RefSeq':
    print 'Error: you cannot compare two RefSeq transcript databases'
    quit()
if not dataXtype == inputtype:
    print 'Error: transcript types in input file and dataset X do not agree'
    quit()

# Read transcript database X
sys.stdout.write('Dataset X ['+options.dataX+']: READING...')
sys.stdout.flush()
dataX_build = dataY_build = ''
if dataXtype == 'Ensembl':
    dataX, N = readEnsemblData(options.dataX, True)
    dataX_build, infoX = readEnsemblInfo(options.dataX)
else:
    dataX, N = readRefSeqData(options.dataX, True)
    infoX = readRefSeqInfo(options.dataX)
sys.stdout.write('\rDataset X ['+options.dataX+']: '+infoX+', '+str(N)+' transcripts')
print ''

# Read transcript database Y
sys.stdout.write('Dataset Y ['+options.dataY+']: READING...')
sys.stdout.flush()
if dataYtype == 'Ensembl':
    dataY, N = readEnsemblData(options.dataY, (dataYtype == dataXtype))
    dataY_build, infoY = readEnsemblInfo(options.dataY)
else:
    dataY, N = readRefSeqData(options.dataY, (dataYtype == dataXtype))
    infoY = readRefSeqInfo(options.dataY)
sys.stdout.write('\rDataset Y ['+options.dataY+']: '+infoY+', '+str(N)+' transcripts')
print '\n'

# Read configuration file
refpath37, refpath38 = readConfigFile(dir)

# Initialize GRCh37 reference genome
if refpath37 is not None:
    if not os.path.isfile(refpath37):
        print '\nError: GRCh37 reference genome file ('+refpath37+') cannot be found.\n'
        quit()
    ref_GRCh37 = reference.Reference(refpath37)
else: ref_GRCh37 = None

# Initialize GRCh38 reference genome
if refpath38 is not None:
    if not os.path.isfile(refpath38):
        print '\nError: GRCh38 reference genome file ('+refpath38+') cannot be found.\n'
        quit()
    ref_GRCh38 = reference.Reference(refpath38)
else: ref_GRCh38 = None

# Check if required reference genome files are specified
if ref_GRCh37 is None and 'GRCh37' in [dataX_build, dataY_build]:
    print '\nError: GRCh37 reference genome file needs to be specified in configuration file.\n'
    quit()
if ref_GRCh38 is None and 'GRCh38' in [dataX_build, dataY_build]:
    print '\nError: GRCh38 reference genome file needs to be specified in configuration file.\n'
    quit()

# Initialize output file
out = open(options.output, 'w')
out.write('# Data set X: '+infoX+'\n')
out.write('# Data set Y: '+infoY+'\n')

# Iterate through the input list of transcripts
n_inboth = 0
n_identical = 0
n_cds_identical = 0
i = 0
for transcriptID in inputlist:
    i+=1
    flags = []
    comparewith = []

    sys.stdout.write('\rAnalysing transcripts '+str(i)+'/'+str(len(inputlist)))
    sys.stdout.flush()

    if transcriptID not in dataX.keys(): flags.append('NFX')
    else: transcript = dataX[transcriptID]

    if dataYtype == dataXtype:
        if transcriptID not in dataY.keys(): flags.append('NFY')
        else: comparewith = [dataY[transcriptID]]
    else:
        if transcript.gene not in dataY.keys(): flags.append('NFY')
        else: comparewith = dataY[transcript.gene]

    if len(flags) > 0:
        out.write('\t'.join([transcriptID, '.', ';'.join(flags)]) + '\n')
        continue

    n_inboth += 1

    identical, cds_identical = compare(transcript, comparewith, dataXtype, dataYtype, dataX_build, dataY_build, ref_GRCh37, ref_GRCh38, options, out)
    if identical: n_identical += 1
    if cds_identical: n_cds_identical += 1

print ' - Done.'

# Close output file
out.close()

# Goodbye message
print '\nSummary:'
print '- '+str(n_inboth) + ' transcripts present in both datasets'
print '- '+str(n_identical)+' transcripts have identical match'
print '- '+str(n_cds_identical)+' transcripts have CDS-identical match'
print '\nOutput written to file: '+options.output
print '\nFinished: '+str(datetime.datetime.now())
print '-'*100+'\n'
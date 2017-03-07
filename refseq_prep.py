import urllib
import urllib2
import gzip
import os
import sys
import datetime
from optparse import OptionParser


###########################################################################################################################

# Version and date
ver = '1.2'
today = datetime.date.today()
datestr = today.strftime('%d %b, %Y')


# Process RefSeq data file
def processFile(outfile, errorfile, txtfile):
    record = []
    counter = 0
    for line in gzip.open('refseqdata.gz','r'):
        line = line.strip()
        if line.startswith('LOCUS'):
            if processRecord(outfile, errorfile, txtfile, record): counter += 1
            record = []
        record.append(line)

    if processRecord(outfile, errorfile, txtfile, record): counter += 1
    urllib.urlcleanup()
    return counter

# Process record
def processRecord(outfile, errorfile, txtfile, record):
    if len(record) == 0: return False
    try:
        exons = []
        dna = ''
        cds = ''
        gene = ''

        sections = splitToSections(record)

        locusline = sections['LOCUS'][0]
        ID = locusline.split()[1]
        if not ID.startswith('NM_'): return False

        versionline = sections['VERSION'][0]
        version = versionline.split()[1]
        version = version.split('.')[1]

        for line in sections['FEATURES']:
            try:
                if line.startswith('exon'):
                    coords = line.split()[1]
                    [start, end] = coords.split('..')
                    exons.append(start+'-'+end)
                if line.startswith('CDS'):
                    coords = line.split()[1]
                    [cds_start, cds_end] = coords.split('..')
                    cds = cds_start+'-'+cds_end
                if line.startswith('/db_xref="HGNC:') and gene == '':
                    gene = line[15:-1]
            except: pass

        for line in sections['ORIGIN']:
            if line.startswith('ORIGIN'): continue
            dna += line
        sequence = ''
        for x in dna:
            if x in ['a','c','g','t']:
                sequence += x.upper()

        if cds == '' or sequence == '' or version == '' or gene == '' or len(exons)==0:
            for line in record: errorfile.write(line+'\n')
            errorfile.write('\n')
            return False

        data = {'ID': ID, 'exons': exons, 'cds': cds, 'seq': sequence, 'ver': version, 'gene': gene}
        output(outfile, txtfile, data)

        return True

    except:
        for line in record: errorfile.write(line+'\n')
        errorfile.write('\n')
        return False

# Retrieve list of RefSeq data files
def getListOfFiles():
    ret = []
    i = 0
    while True:
        i += 1
        url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.'+str(i)+'.rna.gbff.gz'
        try: f = urllib2.urlopen(urllib2.Request(url))
        except: break
        ret.append(url)
    urllib.urlcleanup()
    return ret

# Write record to output file
def output(outfile, txtfile, data):
    if len(data['exons']) > 0: exonsstr = ','.join(data['exons'])
    else: exonsstr = '.'
    record = '\t'.join([data['ID'],data['ver'],data['gene'],exonsstr,data['cds'],data['seq']])
    outfile.write(record+'\n')
    txtfile.write(data['gene']+'\t'+data['ID']+'\n')

# Split record to LOCUS, VERSION, FEATURES and ORIGIN sections
def splitToSections(record):
    ret = {'LOCUS': [], 'VERSION': [], 'FEATURES': [], 'ORIGIN': []}
    keywords = ['LOCUS', 'DEFINITION', 'ACCESSION', 'VERSION', 'KEYWORDS', 'SOURCE', 'REFERENCE', 'FEATURES', 'ORIGIN']
    key = ''
    for line in record:
        for k in keywords:
            if line.startswith(k):
                key = k
                break
        if key in ret.keys(): ret[key].append(line)
    return ret


###########################################################################################################################

# Check if any command line arguments
if not len(sys.argv) > 1:
    print '\nUsage: python path/to/refseq_prep/refseq_prep.py <options>\n'
    quit()

# Command line argument parsing
descr = 'refseq_prep v'+ver
parser = OptionParser(usage='python path/to/refseq_prep/refseq_prep.py <options>', version=ver, description=descr)
parser.add_option('-o', "--out", default='output', dest='output', action='store', help="Output file name prefix [default value: %default]")
(options, args) = parser.parse_args()

# Welcome message
print '\n'+'-'*100
print 'refseq_prep v'+ver+' for transcript_map started: '+str(datetime.datetime.now())+'\n'

# Retrieve list of RefSeq data files
sys.stdout.write('Checking RefSeq FTP site ...')
sys.stdout.flush()
files = getListOfFiles()
sys.stdout.write('\r'+str(len(files))+' RefSeq data files will be downloaded and processed\n')
sys.stdout.flush()
print ''

# Initialize output files
outfile = gzip.open(options.output+'.gz','wb')
errorfile = open(options.output+'.error','w')
txtfile = open(options.output+'.txt','w')
outfile.write('#Created by refseq_prep v'+ver+'; RefSeq data downloaded: '+datestr+'\n')
outfile.write('#'+'\t'.join(['ID','VERSION','GENE','EXONS','CDS','SEQUENCE'])+'\n')
txtfile.write('# Created by refseq_prep '+ver+'\n')
txtfile.write('GENE\tID\n')

# Iterate through
total = 0
for i in range(len(files)):
    url = files[i]

    # Progress information
    sys.stdout.write('\rProcessing data file '+str(i+1)+'/'+str(len(files))+' ... ')
    sys.stdout.flush()

    # Download RefSeq data file
    urllib.urlretrieve(url, 'refseqdata.gz')

    # Process RefSeq data file
    total += processFile(outfile, errorfile, txtfile)

    # Remove RefSeq data file
    os.remove('refseqdata.gz')

sys.stdout.flush()
print ' - Done.'

# Close output files
outfile.close()
errorfile.close()
txtfile.close()

# Goodbye message
print '\n'+str(total)+' NM records retrieved: '+options.output+'.gz'
print '(Erroneous records: '+options.output+'.error)'
print ''
print 'Finished: '+str(datetime.datetime.now())
print '-'*100+'\n'

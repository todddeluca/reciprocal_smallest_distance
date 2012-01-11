#!/usr/bin/env python

'''
What is a FASTA format file/string?
This module follows the NCBI conventions: http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
'''

import cStringIO
import re
import math


TEST_FASTA = {'GOOD_TEST': '''>ns|id|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
''',
              'LONG_TEST': '''>ns|id|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
>ns|id2|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
''',
              'BLANK_IN_TEST': '''>ns|id|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
  
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
''',
              'BLANK_START_TEST': '''  
>ns|id|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
''',
              'DATA_START_TEST': '''CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
>ns|id|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
''',
              'TWO_NAMELINES_TEST': '''>ns|id|a long description
>ns|id2|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
''',
              'TRAILING_NAMELINE_TEST': '''>ns|id|a long description
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA
>ns|id2|a long description
'''}

def test():
    import traceback, StringIO
    for name, testFasta in TEST_FASTA.items():
        try:
            print
            print name
            print 'non-strict'
            for lines in readFastaLines(StringIO.StringIO(testFasta), strict=False):
                print lines
            print 'non-strict, all'
            for lines in readFastaLines(StringIO.StringIO(testFasta), strict=False, goodOnly=False):
                print lines
            print 'strict'
            for lines in readFastaLines(StringIO.StringIO(testFasta)):
                print lines
            print 'readFasta'
            for lines in readFasta(StringIO.StringIO(testFasta)):
                print lines
        except Exception:
            traceback.print_exc()


def idFromName(line):
    '''
    line: a fasta nameline
    returns: an id parsed from the fasta nameline.  The id is the first whitespace separated token after an optional namespace, etc.  See the examples below.
    This covers a lot of cases that a sane person would put on a nameline.  So needless to say it covers very few cases.
    Examples in the form nameline => return value:
    id => id
    id desc => id
    >id => id
    >id desc => id
    >ns|id => id
    >ns|id desc => id
    >ns|id| => id
    >ns|id|desc => id
    ns|id => id
    ns|id desc => id
    ns|id| => id
    ns|id|desc => id
    ns|id blah|desc => id
    Example namelines not covered:
    JGI-PSF GENOMES ftp://ftp.jgi-psf.org/pub/JGI_data/Nematostella_vectensis/v1.0/annotation/proteins.Nemve1FilteredModels1.fasta.gz
    >jgi|Nemve1|18|gw.48.1.1
    >jgi|Nemve1|248885|estExt_fgenesh1_pg.C_76820001
    '''
    # This could probably be done with one regex, but I am too stupid and this way I can read it.

    # remove the leading '>' if there is one.
    if line.startswith('>'):
        line = line[1:]

    # keep only everything after the first pipe.  will keep everything if there is no first pipe.
    pipe = line.find('|')
    if pipe > -1:
        line = line[line.find('|')+1:]

    # keep everything before the second pipe.  will keep everything if there is no second pipe.
    pipe = line.find('|')
    if pipe > -1:
        line = line[:pipe]

    # return the first token as the id.
    return line.split()[0]


def prettySeq(seq, n=60):
    '''
    seq: one long bare (no nameline) sequence. e.g. MASNTVSAQGGSNRPVRDFSNIQDVAQFLLFDPIWNEQPGSIVPWKMNREQALAERYPELQTSEPSEDYSGPVESLELLPLEIKLDIMQYLSWEQISWCKHPWLWTRWYKDNVVRVSAITFED
    n: maximum length of sequence lines
    returns: seq split over multiple lines, all terminated by newlines.
    '''
    if len(seq) == 0:
        raise Exception('zero-length sequence', seq)
    seq = ''.join(seq.strip().split())
    chunks = int(math.ceil(len(seq)/float(n)))
    pretty = ''
    for i in range(chunks):
        pretty += seq[i*n:(i+1)*n] + '\n'
    return pretty


def numSeqsInFastaDb(path):
    num = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                num += 1
    return num


def readIds(fastaFile, strict=True):
    '''
    fastaFile: a file-like object or a path to a fasta file
    yields: id in each nameline.
    '''
    for nameline in readNamelines(fastaFile, strict):
        yield idFromName(nameline)


def readNamelines(fastaFile, strict=True):
    '''
    fastaFile: a file-like object or a path to a fasta file
    yields: each nameline
    '''
    for nameline, seq in readFasta(fastaFile, strict):
        yield nameline
        

def readFasta(fastaFile, strict=True):
    '''
    fastaFile: a file-like object or a path to a fasta file
    yields: a tuple of (nameline, sequence) for each sequence in the fasta file.
    '''
    for lines in readFastaLines(fastaFile, strict):
        nameline = lines[0].strip()
        seq = ''.join((l.strip() for l in lines[1:]))
        yield nameline, seq


def readFastaLines(fastaFile, strict=True, goodOnly=True, filterBlankLines=False):
    '''
    fastaFile: a file-like object or a path to a fasta file
    yields: a seq of fasta sequence lines for each sequence in the fasta file.
    the first line is the nameline.  the other lines are the sequence data lines.  lines include newlines.
    '''
    if isinstance(fastaFile, basestring):
        with open(fastaFile) as fh:
            for lines in _fastaSeqIter(fh, strict, goodOnly, filterBlankLines):
                yield lines
    else:
        for lines in _fastaSeqIter(fastaFile, strict, goodOnly, filterBlankLines):
            yield lines

    
def splitSeq(seq):
    '''
    seq: a well-formed fasta sequence string containing a single nameline, including '>' and sequence data lines.
    returns: tuple of nameline, including '>', without a newline, and concatenated sequence lines, without newlines
    e.g. ['>blahname', 'AFADFDSAFAFAFAFFAFAF']
    '''
    lines = seq.splitlines()
    name = lines[0].strip()
    chars = ''.join([l.strip() for l in lines[1:]])
    return [name, chars]
    

def _fastaSeqIter(filehandle, strict=True, goodOnly=True, filterBlankLines=False):
    '''
    filehandle: file object containing fasta-formatted sequences.
    strict: if True, raise an exception when a malformed fasta sequence is encountered.
      A malformed sequence is a sequence with a blank line, a sequence line not preceded by a nameline, or a nameline not followed by a sequence line.
      E.g. a nameline directly after a nameline, like '>foo\n>bar\nCTAGCTAGGGCA\n'
    goodOnly: if True, only yield well-formed fasta sequences, ones with a nameline and one or more sequence datalines and no blank lines.
      if False and strict is False, all sequences, malformed or otherwise, will be yielded.  there will always be at least one (possibly blank) line.
    filterBlankLines: if True, no blank lines (lines only containing whitespace) will be yielded and blank lines will not raise an exception.
    Parses the filehandle, yielding one fasta sequence at a time.
    yields: a seq of fasta sequence lines.  the first line is the nameline.  the other lines are the sequence data lines.  
    '''
    for lines in _splitOnNamelines(filehandle, filterBlankLines):
        if not lines[0] or lines[0][0] != '>':
            if strict:
                raise Exception('FASTA error: sequence must start with a nameline.', lines)
            elif not goodOnly:
                yield lines
        elif len(lines) < 2:
            if strict:
                raise Exception('FASTA error: sequence must contain at least one sequence data line.', lines)
            elif not goodOnly:
                yield lines
        elif '' in (line.strip() for line in lines): # contains a blank line
            if strict:
                raise Exception('FASTA error: blank lines not allowed')
            elif not goodOnly:
                yield lines
        else: # a good sequence
            yield lines


def _splitOnNamelines(filehandle, filterBlankLines=False):
    '''
    split the lines in filehandle on namelines.
    yields: seq of lines, where the first line is a nameline (except if filehandle starts with a non-nameline) the other lines are lines until the next nameline
    or the end of the file.  lines include newlines.  length of yielded seq always contains at least one line.
    Said another way, the seq of lines will always contain at least one line.  only the first line will ever be a nameline.
    filterBlankLines: if True, no blank lines (lines only containing whitespace) will be yielded.
    '''
    lines = []
    for line in filehandle:
        if line and line[0] == '>': # a nameline
            if lines:
                yield lines # yield last sequence
            lines = [line] # start new sequence
        elif not filterBlankLines:
            lines.append(line)
        elif line.strip():
            lines.append(line)
    if lines:
        yield lines


def isNameLine(line):
    return line.startswith('>')


def _tern(comp, trueVal, falseVal):
    ''' a tawdry ternary operator '''
    if comp:
        return trueVal
    return falseVal


def head(query, n):
    '''returns the first n sequences in query.'''
    count = 0
    headstr = ''
    for line in query.splitlines(keepends=1):
        if line.startswith('>'):
            count += 1
            if count > n: break
        headstr += line
    return headstr


def dbSize(query):
    '''returns the number of sequence characters'''
    size = 0
    for line in query.splitlines():
        if isNameLine(line):
            continue
        size += len(line.strip())
    return size


def numChars(query):
    '''
    synonym for dbSize().  returns the number of character (e.g. bases or residues for nucleotide or protein sequences).
    '''
    return dbSize(query)
    

def numSeqs(query):
    '''
    synonym for size(), whose name is a little more specific as to what is being measured: the number of sequences.
    '''
    return size(query)


def size(query):
    '''
    query: string containing fasta formatted seqeunces
    returns: the number of sequences
    '''
    fh = cStringIO.StringIO(query)
    size = numSeqsInFile(fh)
    fh.close()
    return size


def numSeqsInFile(file):
    '''
    file: file like object containing fasta formatted sequences
    '''
    return sum([1 for line in file if isNameLine(line.strip())])


def numSeqsInPath(path):
    '''
    path: path to fasta formatted db
    returns: number of sequences in fasta db
    '''
    fh = open(path)
    size = numSeqsInFile(fh)
    fh.close()
    return size


def main():
    pass

        
if __name__ == '__main__':
    main()
                        

##################
# DECPRECATED CODE
##################


def fastaSeqIterOld(filehandle, strict=True):
    '''
    filehandle: file object containing fasta-formatted sequences.
    ignoreParseError: if True, parsing ignores namelines that do not have sequence character lines.  For example, '>foo\n>bar\nABCD\n'
    would yield the 'bar' sequence, ignoring the 'foo' sequence that has no sequence characters associated with it.
    In all cases blank lines are ignored, no matter where they occur.
    Generator function yielding a string representing a single fasta sequence (name line including '>' and sequence lines)
    for each fasta sequence in filehandle.
    returns: a generator.
    notes:
    This function was modified from fastaSeqIterStrict to handle bogus fasta input like this:
    >ORFP:20136 YOL048C, Contig c378 4079-4399 reverse complement
    MLFKVSNFTSLTLLSLIPIVGPILANQLMAPKRTFTYLQRYFLLKGFSKKQAKDFQYEHYASFICFGMSAGLLELIPFFTIVTISSNTVGAAKWCSSLLKGERKKD*
    >ORFP:18671 , Contig c238 100299-100300 reverse complement

    >ORFP:20137 , Contig c378 4878-5189 reverse complement
    MKVGIELISHSQTSHGTHVNSTVLAEKTPQPLEKPSKEHSISKESNINRWLKI

    LRRQFDIWFPETIPTMKVRYELLKKNFIKEIFNSRAFIYPFLVSILYYLY*
    The old function, even with error handling turned on, would not concatenate all the sequence characters of the 3rd sequence
    since they are separated by a blank line.
    '''
    # states: seeking_nameline, seeking_seqline, in_seq.  
    state = 'seeking_nameline'
    fasta = ''
    for line in filehandle:
        # ignore all blank lines
        if not line.strip():
            continue
        elif state == 'seeking_nameline' and line.startswith('>'):
            state = 'seeking_seqline'
            fasta = line
        elif state == 'seeking_nameline' and not ignoreParseError:
            raise Exception('FASTA parse error.  Looking for name line and found line which is neither blank nor nameline.  line=%s'%line)
        elif state == 'seeking_seqline' and line.startswith('>'):
            if ignoreParseError:
                # skip nameline without sequence and restart with this nameline.
                state = 'seeking_seqline'
                fasta = line
            else:
                raise Exception('FASTA parse error.  Looking for sequence line and found name line.  line=%s'%line)
        elif state == 'seeking_seqline':
            state = 'in_seq'
            fasta += line
        elif state == 'in_seq' and line.startswith('>'):
            yield fasta
            state = 'seeking_seqline'
            fasta = line
        elif state == 'in_seq':
            fasta += line
        else:
            raise Exception('FASTA parse error.  Unrecognized state.  state=%s, line=%s'%(state, line))

    if state == 'in_seq':
        yield fasta
    elif state == 'seeking_seqline' and not ignoreParseError:
        raise Exception('FASTA parse error.  Looking for sequence line and found end of file.')
    elif state == 'seeking_nameline':
        pass
    else:
        raise Exception('FASTA parse error.  Unrecognized state found at end of file.  state=%s'%state)








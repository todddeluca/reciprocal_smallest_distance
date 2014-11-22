#!/usr/bin/env python

'''
What is a FASTA format file/string?
This module follows the NCBI conventions: http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
'''

import cStringIO
import math


def idFromName(line):
    '''
    line: a fasta nameline
    returns: an id parsed from the fasta nameline.  The id is the first
    whitespace separated token after an optional namespace, etc.  See the
    examples below.  This covers a lot of cases that a sane person would put on
    a nameline.  So needless to say it covers very few cases.  Examples in the
    form nameline => return value:

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
    seq: one long bare (no nameline) sequence. e.g.
    MASNTVSAQGGSNRPVRDFSNIQDVAQFLLFDPIWNEQPGSIVPWKMNREQALAERYPELQTSEPSEDYSGPVESLELLPLEIKLDIMQYLSWEQISWCKHPWLWTRWYKDNVVRVSAITFED
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


def readIds(fastaFile):
    '''
    fastaFile: a file-like object or a path to a fasta file
    Yields each id in each nameline in each sequence in the fasta file.
    '''
    for nameline in readNamelines(fastaFile):
        yield idFromName(nameline)


def readNamelines(fastaFile):
    '''
    fastaFile: a file-like object or a path to a fasta file
    Yields each nameline in each sequence in the fasta file.
    '''
    for nameline, seq in readFasta(fastaFile):
        yield nameline


def readFasta(fastaFile):
    '''
    fastaFile: a file-like object or a path to a fasta file
    Yields a tuple of (nameline, sequence) for each sequence in the fasta file.
    Newlines are stripped from the nameline and sequence lines, and the sequence
    lines are concatenated into one long sequence string.
    Here is an examle (edited for length):
    ('>sp|P31946|1433B_HUMAN',
     'MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKT')
    '''
    for lines in readFastaLines(fastaFile):
        nameline = lines[0].strip()
        seq = ''.join((l.strip() for l in lines[1:]))
        yield nameline, seq


def readFastaLines(fastaFile):
    '''
    fastaFile: a file-like object or a path to a fasta file
    yields: a seq of fasta sequence lines for each sequence in the fasta file.
    the first line is the nameline.  the other lines are the sequence data lines.  lines include newlines.
    '''
    if isinstance(fastaFile, basestring):
        with open(fastaFile) as fh:
            for lines in relaxedFastaSeqIter(fh):
                yield lines
    else:
        for lines in relaxedFastaSeqIter(fastaFile):
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


def relaxedFastaSeqIter(filehandle):
    '''
    Parse the lines in filehandle, first removing any blank lines, and then
    yielding all well-formed fasta sequences and ignoring badly-formed seqs.  A
    well-formed sequence has exactly one nameline followed by one or more
    sequence lines.

    Well-formed example:

        >sp|P27348|1433T_HUMAN
        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR
        EGAEN

    Badly-formed example (no nameline):

        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR
        EGAEN

    Badly-formed example (no sequence lines):

        >sp|P27348|1433T_HUMAN
    '''
    # lines guaranteed to have no blank lines by filterBlanks()
    # lines guaranteed to have have exactly one nameline as the first
    # element (except possibly the first lines yielded, which might not
    # have a nameline if the filehandle starts with a sequence line).
    for lines in splitFastaOnNamelines(filterBlanks(filehandle)):
        if lines[0][0] == '>' and len(lines) >= 2:
            yield lines


def filterBlanks(lines):
    '''
    Yield each line in lines that contains non-whitespace characters.
    Used to remove blank lines from FASTA files.
    '''
    for line in lines:
        if line.strip():
            yield line


def splitFastaOnNamelines(filehandle):
    '''
    Split the lines in filehandle on namelines.  Yields a seq of lines, where
    the first line in the seq is a nameline (except if filehandle starts with a
    non-nameline) and the other lines are lines until the next nameline or the
    end of the file.  Lines include newlines.  The seq of lines will always
    contain at least one line.  Only the first line will ever be a nameline.

    Example input (note that this is not well-formed fasta, since it starts
    with a sequence line, has a nameline with no sequence lines, and has blank
    lines within a sequence):

        VLSSIEQKSNEEGSEEKGPEVREYREKVETELQGVCDTVLGLLDSHLIKEAGDAESRVFY
        >sp|P31947|1433S_HUMAN
        >sp|P27348|1433T_HUMAN
        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR

        EGAEN

        >sp|P63104|1433Z_HUMAN

        MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR
        MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE

    Example yielded output (note how every sequence except the first starts
    with a nameline):
        yield ['VLSSIEQKSNEEGSEEKGPEVREYREKVETELQGVCDTVLGLLDSHLIKEAGDAESRVFY\n']
        yield ['>sp|P31947|1433S_HUMAN\n']
        yield ['>sp|P27348|1433T_HUMAN\n', 
        'MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR\n',
        '\n',
        'EGAEN\n',
        '\n']
        yield ['>sp|P63104|1433Z_HUMAN\n',
        '\n',
        'MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR\n',
        'MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE\n']

    Well-formed example input (note how the first line of the input is a
    nameline):

        >sp|P27348|1433T_HUMAN
        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR
        EGAEN
        >sp|P63104|1433Z_HUMAN
        MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR
        MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE

    Well-formed example output (note how the first element of the first yielded
    list is a nameline):

        yield ['>sp|P27348|1433T_HUMAN\n', 
        'MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR\n',
        'EGAEN\n']
        yield ['>sp|P63104|1433Z_HUMAN\n',
        'MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR\n',
        'MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE\n']

    '''
    lines = []
    for line in filehandle:
        if line and line[0] == '>': # a nameline
            if lines:
                yield lines # yield current sequence
            lines = [line] # start new sequence
        else:
            lines.append(line) # add to current sequence

    if lines:
        yield lines # yield current sequence


def isNameLine(line):
    return line.startswith('>')


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


def readFastaLinesOld(fastaFile, strict=True, goodOnly=True, filterBlankLines=False):
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
    Split the lines in filehandle on namelines.  Example nameline: '>lcl|12345'
    Yields a seq of lines, where the first line is a nameline (except if filehandle starts with a non-nameline) and the other lines are lines until the next nameline
    or the end of the file.  Lines include newlines.  The seq of lines will always contain at least one line.  Only the first line will ever be a nameline.

    filterBlankLines: if True, no blank lines (lines only containing whitespace) will be yielded.


        >sp|P27348|1433T_HUMAN
        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR
        EGAEN

        >sp|P63104|1433Z_HUMAN
        MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR
        MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE
    '''
    lines = []
    for line in filehandle:
        if line and line[0] == '>': # a nameline
            if lines:
                yield lines # yield current sequence
            lines = [line] # start new sequence
        elif not filterBlankLines: # a blank or non-blank line
            lines.append(line) # add to current sequence
        elif line.strip(): # a non-blank line
            lines.append(line) # add to current sequence
    if lines:
        yield lines # yield current sequence


def _fastaSeqIter2(filehandle):
    '''
    Untested code.

    This is a state-machine parser for FASTA formatted files.  The page
    http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml describes this format.
    This parser is tolerant of blank lines that come before a nameline and/or
    after the sequence lines of a sequence, but not within a sequence.

    This is OK because it contains a blank line after the sequence lines of the
    first sequence and before the sequence lines of the second sequence:

        >sp|P27348|1433T_HUMAN
        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR
        EGAEN

        >sp|P63104|1433Z_HUMAN
        MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR
        MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE

    This is NOT ok because one sequence contains a blank line after then
    nameline and before the sequence lines and the other sequence contains a
    blank line between two sequence lines:

        >sp|P27348|1433T_HUMAN

        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR
        EGAEN
        >sp|P63104|1433Z_HUMAN
        MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR

        MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE

    Sequences must start with a nameline.  This is NOT ok:

        MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR
        EGAEN

    Sequences must contain at least one sequence line.  This is NOT ok:

        >sp|P63104|1433Z_HUMAN

    For every sequence in the fasta file, this parser yields the lines of each
    sequence in the fasta file.  Here is an example of what it yields:

        ['>sp|P27348|1433T_HUMAN\n',
        'MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWR\n',
        'EGAEN\n']

    This parser will raise an Exception
    '''
    # parsing states
    LFN = 'Looking for nameline'
    LFS = 'Looking for sequence line'
    LFA = 'Looking for any line' # a sequence or a nameline
    # line states: blank, nameline (starts with '>'), seq line (starts with
    # anything else), or EOF.
    BL = 'blank line'
    NL = 'nameline'
    SL = 'sequence line'

    state = LFN
    lines = []
    for line in filehandle:
        linestate = BL if not line.strip() else NL if line[0] == '>' else SL
        if state == LFN:
            if linestate == BL:
                continue
            elif linestate == NL:
                lines = [line] # Start a new sequence
                state = LFS # Found nameline.  Now look for seq lines.
            elif linestate == SL:
                raise Exception(PARSING_ERROR, 'Expecting a {} or {}.  Found a {}'.format(BL, NL, SL), line, state, linestate)
            else:
                raise Exception(PARSING_ERROR, 'Unrecognized line state', line, state, linestate)
        elif state == LFS:
            if linestate == BL:
                raise Exception(PARSING_ERROR, 'Expecting a {}.  Found a {}'.format(SL, BL), line, state, linestate)
            elif linestate == NL:
                raise Exception(PARSING_ERROR, 'Expecting a {}.  Found a {}'.format(SL, NL), line, state, linestate)
            elif linestate == SL:
                lines.append(line) # Add to the current sequence
                state = LFA # Found a seq line, Now look for more seq lines or a new sequence.
            else:
                raise Exception(PARSING_ERROR, 'Unrecognized line state', line, state, linestate)
        elif state == LFA:
            if linestate == BL:
                yield lines # Emit the current sequence
                lines = []
                state = LFN # Look for a new sequence.
            elif linestate == NL:
                yield lines # Emit the current sequence
                lines = [line] # Start a new sequence
                state = LFS # Found nameline.  Now look for seq lines.
            elif linestate == SL:
                lines.append(line) # Add to the current sequence.
            else:
                raise Exception(PARSING_ERROR, 'Unrecognized line state', line, state, linestate)
        else:
            raise Exception(PARSING_ERROR, 'Unrecognized parsing state', line, state, linestate)

    # EOF
    if state == LFN:
        pass # Done
    elif state == LFS:
        raise Exception(PARSING_ERROR, 'Expecting a {}.  Found an EOF.'.format(SL), state)
    elif state == LFA:
        yield lines # Emit the current sequence
        pass # Done
    else:
        raise Exception(PARSING_ERROR, 'Unrecognized parsing state', state)








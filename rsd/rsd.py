#!/usr/bin/env python2.7
# RSD: The reciprocal smallest distance algorithm.
#   Wall, D.P., Fraser, H.B. and Hirsh, A.E. (2003) Detecting putative orthologs, Bioinformatics, 19, 1710-1711.
# Original author: Dennis P. Wall, Department of Biological Sciences, Stanford University.
# Contributors: I-Hsien Wu, Computational Biology Initiative, Harvard Medical School
# Maintainer: Todd F. DeLuca, Center for Biomedical Informatics, Harvard Medical School
#

# This program is written to run on linux.  It has not been tested on Windows.
# To run this program you need to have installed on your system:
# Python 2.7
# NCBI BLAST 2.2.24  
# paml 4.4
# Kalign 2.04 (recommended) or clustalw 2.0.9 (deprecated)
# see README FOR FULL DETAILS

import argparse
import cStringIO
import glob
import logging
import os
import re
import shutil
import subprocess
import time

import fasta
import nested
import util


PAML_ERROR_MSG = 'paml_error'
FORWARD_DIRECTION = 0
REVERSE_DIRECTION = 1
DASHLEN_RE = re.compile('^(-*)(.*?)(-*)$')

JIKE_DEBUG = util.getBoolFromEnv('ROUNDUP_RSD_DEBUG', False)

MAX_HITS = 3
MATRIX_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'jones.dat')
CODEML_CONTROL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'codeml.ctl')

# Constants used when aligning seqs with clustalw.  Kalign does not need these.
USE_CLUSTALW = util.getBoolFromEnv('RSD_USE_CLUSTALW', False)
CLUSTAL_INPUT_FILENAME = 'clustal_fasta.faa'
CLUSTAL_ALIGNMENT_FILENAME = 'clustal_fasta.aln'


#################
# BLAST FUNCTIONS
#################
#
# Used to compute blast hits between two genomes, parse the results, and save the best hits to a file
#

def formatForBlast(fastaPath):
    # os.chdir(os.path.dirname(fastaPath))
    # cmd = 'formatdb -p -o -i'+os.path.basename(fastaPath)
    # cmd = 'formatdb -p -o -i'+fastaPath
    # redirect stdout to /dev/null to make the command quiter.
    cmd = 'makeblastdb -in {} -dbtype prot -parse_seqids >/dev/null'.format(fastaPath)
    subprocess.check_call(cmd, shell=True)


def getHitId(hit):
    return hit[0]


def getHitEvalue(hit):
    return hit[1]


def loadBlastHits(path):
    '''
    path: location of stored blast hits computed by computeBlastHits()
    returns: mapping object from query id to hits.  used to be a bsddb, now is a dict.
    '''
    return util.loadObject(path)


def getBlastHits(queryFastaPath, subjectIndexPath, evalue, limitHits=MAX_HITS, workingDir='.', copyToWorking=False):
    '''
    queryFastaPath: location of fasta file of query sequences
    subjectIndexPath: location and name of blast-formatted indexes.
    evalue: a string or float representing the maximum evalue threshold of hits to get.
    workingDir: creates, uses, and removes a directory under workingDir.
    copyToWorking: if True, copy query fasta path and subject index files to within the working directory and use the copies to blast.
      can improve performance if the working directory is on local disk and the files are on a slow network.
    blasts every sequence in query agaist subject, adding hits that are better than evalue to a list stored in a dict keyed on the query id.
    '''
    # work in a nested tmp dir to avoid junking up the working dir.
    with nested.NestedTempDir(dir=workingDir, nesting=0) as tmpDir:
        if copyToWorking:
            localFastaPath = os.path.join(tmpDir, 'query.fa')
            shutil.copyfile(queryFastaPath, localFastaPath)
            localIndexDir = os.path.join(tmpDir, 'local_blast')
            os.makedirs(localIndexDir, 0770)
            localIndexPath = os.path.join(localIndexDir, os.path.basename(subjectIndexPath))
            for path in glob.glob(subjectIndexPath+'*'):
                if os.path.isfile:
                    shutil.copy(path, localIndexDir)
            queryFastaPath = localFastaPath
            subjectIndexPath = localIndexPath
        blastResultsPath = os.path.join(tmpDir, 'blast_results')
        # blast query vs subject, using /opt/blast-2.2.22/bin/blastp
        cmd = 'blastp -outfmt 6 -evalue %s -query %s -db %s -out %s'%(evalue, queryFastaPath, subjectIndexPath, blastResultsPath)
        subprocess.check_call(cmd, shell=True)
        # parse results
        hitsMap = parseResults(blastResultsPath, limitHits)
    return hitsMap


def computeBlastHits(queryFastaPath, subjectIndexPath, outPath, evalue, limitHits=MAX_HITS, workingDir='.', copyToWorking=False):
    '''
    queryFastaPath: location of fasta file of query sequences
    subjectIndexPath: location and name of blast-formatted indexes.
    evalue: a string or float representing the maximum evalue threshold of hits to get.
    outPath: location of file where blast hits are saved.
    workingDir: creates, uses, and removes a directory under workingDir.  
    copyToWorking: if True, copy query fasta path and subject index files to within the working directory and use the copies to blast.
      can improve performance if the working directory is on local disk and the files are on a slow network.
    Runs getBlastHits() and persists the hits to outPath.
    '''
    hitsMap = getBlastHits(queryFastaPath, subjectIndexPath, evalue, limitHits, workingDir, copyToWorking)
    util.dumpObject(hitsMap, outPath)


def parseResults(blastResultsPath, limitHits=MAX_HITS):
    '''
    returns: a map from query seq id to a list of tuples of (subject seq id, evalue) for the top hits of the query sequence in the subject genome
    '''
    # parse tabular results into hits.  thank you, ncbi, for creating results this easy to parse.
    hitsMap = {}
    hitsCountMap = {}
    prevSeqId = None
    prevHitId = None
    fh = open(blastResultsPath)
    for line in fh:
        splits = line.split()
        try:
            seqId = fasta.idFromName(splits[0]) # remove namespace prefix, e.g. 'gi|'
            hitId = fasta.idFromName(splits[1])
            hitEvalue = float(splits[10])
        except Exception as e:
            logging.exception('parseResults(): prevSeqId: {}, prevHitId: {}, line: {}'.format(prevSeqId, prevHitId, line))
        # results table reports multiple "alignments" per "hit" in ascending order by evalue
        # we only store the top hits.
        if prevSeqId != seqId or prevHitId != hitId:
            prevSeqId = seqId
            prevHitId = hitId
            if seqId not in hitsCountMap:
                hitsCountMap[seqId] = 0
                hitsMap[seqId] = []
            if not limitHits or hitsCountMap[seqId] < limitHits:
                hitsCountMap[seqId] += 1                
                hitsMap[seqId].append((hitId, hitEvalue))
    fh.close()
    return hitsMap
    
    
###############
# RSD FUNCTIONS
###############


def pamlGetDistance(path):
    filename = '%s/2AA.t'%path
    
    # adding a pause on the off-chance that the filesystem might be lagging a bit, causing the open() to fail below.
    # I think it is more likely that codeml in runPaml_all() is failing before writing the file.
    if not os.path.isfile(filename):
        import time
        time.sleep(0.5)

    with open(filename) as rst:
        get_rst = rst.readlines()
    os.unlink(filename)
        
    if not get_rst:
        raise Exception(PAML_ERROR_MSG, path)
        		
    str = ''
    for line in get_rst[1:]:
        cd1 = line.split()
        if not len(cd1) > 1:
            str += "%s "%(line.split('\n')[0])
            continue
        if len(cd1) > 1:
            str+="%s %s"%(cd1[0], cd1[1])

    dist = float(str.split()[2])
    return dist


def alignFastaKalign(input):
    '''
    input: string containing fasta formatted sequences to be aligned.
    runs alignment program kalign
    Returns: fasta-formatted aligned sequences
    '''
    alignedFasta = util.run(['kalign', '-f', 'fasta'], input) # output clustalw format
    return alignedFasta.replace('\n\n', '\n') # replace fixes a bug in Kalign version 2.04, where if a seq is exactly 60 chars long, an extra newline is output.
    

def alignFastaClustalw(input, path):
    '''
    input: string containing fasta formatted sequences to be aligned.
    path: working directory where fasta will be written and clustal will write output files.
    runs alignment program clustalw
    Returns: fasta-formatted aligned sequences
    '''
    clustalFastaPath = os.path.join(path, CLUSTAL_INPUT_FILENAME)
    clustalAlignmentPath = os.path.join(path, CLUSTAL_ALIGNMENT_FILENAME)
    util.writeToFile(input, clustalFastaPath)
    try:
        subprocess.check_call('clustalw -output=fasta -infile=%s -outfile=%s 2>&1 >/dev/null'%(clustalFastaPath, clustalAlignmentPath), shell=True)
    except Exception:
        logging.exception('runClustal Error:  clustalFastaPath data = %s'%open(clustalFastaPath).read())
        raise
    alignedFasta = util.readFromFile(clustalAlignmentPath)
    return alignedFasta
    

def dashlen_check(seq):
    '''
    Objective: calculate the density of gaps in a sequence at 5' and 3' ends --  caused by poor alignment or by diff length seqs
    Arguments: sequence
    Result: the number of bases to be cut from the subjects 5' and 3' ends, and the divergence of the trimmed seq.
    '''
    seq = seq.strip()
    # trim the dashes from the front and end
    (frontDashes, trimmedSeq, endDashes) = DASHLEN_RE.search(seq).groups()
    # logging.debug('dashlen_check: seq=%s'%seq)
    # all dashes -- do not trim anything
    if not trimmedSeq:
        return (0, 0)

    # ignore trims < 10.
    frontTrim = len(frontDashes)
    if frontTrim < 10:
        frontTrim = 0
    endTrim = len(endDashes)
    if endTrim < 10:
        endTrim = 0

    trimmedSeqDivergence = (trimmedSeq.count('-') / float(len(trimmedSeq)))
    return (frontTrim, endTrim, trimmedSeqDivergence)


def makeGetSeqForId(genomeFastaPath):
    '''
    genomeFastaPath: location of fasta file.  also location/name of blast formatted indexes of the fasta file.
    '''
    # suck fasta file into memory, converting it into a map from id to sequence
    # in memory dict performs much better than on-disk retrieval with xdget or fastacmd.
    # and genome fasta files do not take much space (on a modern computer).
    fastaMap = {}
    for (seqNameline, seq) in fasta.readFasta(genomeFastaPath):
        seqId = fasta.idFromName(seqNameline)
        fastaMap[seqId] = seq
    def getSeqForIdInMemory(seqId):
        return fastaMap[seqId]
    return getSeqForIdInMemory
    

def makeGetHitsOnTheFly(genomeIndexPath, evalue, workingDir='.'):
    '''
    genomeIndexPath: location of blast formatted indexes.  usually same directory/name as genome fasta path
    returns: a function that returns that takes as input a sequence id and sequence and returns the blast hits
    workingDir: a directory in which to create, use, and delete temporary files and dirs.
    '''
    def getHitsOnTheFly(seqid, seq):
        with nested.NestedTempDir(dir=workingDir, nesting=0) as tmpDir:
            queryFastaPath = os.path.join(tmpDir, 'query.faa')
            # add 'lcl|' to make ncbi blast happy.
            util.writeToFile('{0}\n{1}\n'.format('>lcl|'+seqid, seq), queryFastaPath)
            hitsDb = getBlastHits(queryFastaPath, genomeIndexPath, evalue, workingDir=workingDir)
        return hitsDb.get(seqid)
    return getHitsOnTheFly


def makeGetSavedHits(filename):
    '''
    returns a function which can be used to get the hits
    from a file containing pre-computed blast results
    '''
    # in memory retrieval is faster than on-disk retrieval with bsddb, but this has a minor impact on overall roundup performance.
    hitsDb = loadBlastHits(filename)
    def getHitsInMemory(seqid, seq):
        return hitsDb.get(seqid)
    return getHitsInMemory


def getGoodEvalueHits(seqId, seq, getHitsFunc, getSeqFunc, evalue):
    '''
    returns: a list of pairs of (seqid, sequence, evalue) that have an evalue below evalue
    '''
    goodhits = []

    hits = getHitsFunc(seqId, seq)
    
    # check for 3 or fewer blast hits below evalue threshold
    if hits:
        hitCount = 0
        for hit in hits:
            if hitCount >= MAX_HITS:
                break
            hitSeqId = getHitId(hit)
            hitEvalue = getHitEvalue(hit)
            if hitEvalue < evalue:
                hitCount += 1
                hitSeq = getSeqFunc(hitSeqId)
                goodhits.append((hitSeqId, hitSeq, hitEvalue))

    if JIKE_DEBUG:
        for data in goodhits:
            print 'hit\t{0}\t{2}'.format(*data)
    return goodhits


def getDistanceForAlignedSeqPair(seqId, alignedSeq, hitSeqId, alignedHitSeq, workPath):

    # paranoid check: aligned and trimmed seqs need to be the same length.
    # if len(alignedSeq) != len(alignedHitSeq):
    #     raise Exception('getDistanceForAlignedSeqPairs: different lengths for seqs: '+str(((seqId, alignedSeq), (hitSeqId, alignedHitSeq))))

    dataFileName = 'datafile.seq'
    treeFileName = 'treefile.seq'
    outFileName = 'outfile.seq'
    dataFilePath = os.path.join(workPath, dataFileName)
    treeFilePath = os.path.join(workPath, treeFileName)
    outFilePath = os.path.join(workPath, outFileName)
    
    # heading is number of seqs and length of each seq (which all need to be the same len).
    heading = '2 %s\n'%len(alignedSeq)
    pamlData = heading + '%s\n%s\n'%(seqId, alignedSeq) + '%s\n%s\n'%(hitSeqId, alignedHitSeq)
    # logging.debug('pamlData=%s'%pamlData)
    util.writeToFile(pamlData, dataFilePath)
    
    # workPath is simply your folder that will contain codeml (Yang 2000), codeml.ctl (the codeml control file), and the jones.dat (Jones et. al, 1998)
    # write the codeml control file that will run codeml
    # run the codeml 
    
    try:
        subprocess.check_call('codeml >/dev/null', cwd=workPath, shell=True) # /dev/null to silence extraneous codeml output
        distance = pamlGetDistance(workPath)
        if JIKE_DEBUG:
            print 'dist\t{0}\t{1}'.format(hitSeqId, distance)
        return distance
    finally:
        for filePath in [dataFilePath, treeFilePath, outFilePath]:
            if os.path.exists(filePath):
                os.remove(filePath)
            

def getGoodDivergenceAlignedTrimmedSeqPair(seqId, seq, hitSeqId, hitSeq, workPath):
    '''
    aligns seq to hit.  trims aligned seq and hit seq.
    returns: pairs of pairs of id and aligned trimmed sequences for sequences in hits,
    and a predicate function that, given a divergence threshold, says if the divergence of the sequences exceeds the threshold.
    e.g. ((seqId, alignedTrimmedSeq), (hitSeqId, alignedTrimmedHitSeq), divergencePredicateFunc)
    '''
    # ALIGN SEQ and HIT
    # need to align the sequences so we'z can study the rate of evolution per site
    inputFasta = '>%s\n%s\n>%s\n%s\n'%(seqId, seq, hitSeqId, hitSeq)
    if USE_CLUSTALW:
        alignedFasta = alignFastaClustalw(inputFasta, workPath)
    else:
        alignedFasta = alignFastaKalign(inputFasta)
    try:
        (alignedIdAndSeq, alignedHitIdAndSeq) = ((fasta.idFromName(seqNameline), seq) for seqNameline, seq in fasta.readFasta(cStringIO.StringIO(alignedFasta)))
    except Exception as e:
        e.args += (inputFasta, alignedFasta)
        raise
    
    # CHECK FOR EXCESSIVE DIVERGENCE AND TRIMMING
    # find most diverged sequence
    # sort sequences by dash count.  why?
    divIdSeqs = []
    for id, seq in (alignedIdAndSeq, alignedHitIdAndSeq):
        dashCount = seq.count('-')
        div = dashCount / float(len(seq))
        g = (dashCount, div, id, seq)
        divIdSeqs.append(g)
    divIdSeqs.sort()

    if JIKE_DEBUG:
        for data in divIdSeqs:
            if data[2] != seqId:
                print 'div\t{}\t{}'.format(data[2], data[1])
        
    # check for excessive divergence
    leastDivergedDashCount, leastDivergedDiv, leastDivergedId, leastDivergedSeq = divIdSeqs[0]
    # check for excessive divergence and generate dashtrim.
    mostDivergedDashCount, mostDivergedDiv, mostDivergedId, mostDivergedSeq = divIdSeqs[1]
    # dashtrim = dashlen_check(mostDivergedSeq, divergence)
    startTrim, endTrim, trimDivergence = dashlen_check(mostDivergedSeq)
    # logging.debug('dashtrim='+str(dashtrim))
    # trim and add seqs to output
    def divergencePredicate(divergenceThreshold):
        '''Why this logic?  Ask Dennis.  Function closed over local variables that returns whether or not the alignment of the sequences is too diverged.'''
        if leastDivergedSeq and leastDivergedDiv > divergenceThreshold:
            return True
        if (startTrim or endTrim) and trimDivergence >= divergenceThreshold:
            return True
        return False
            
    alignedTrimmedIdAndSeq, alignedTrimmedHitIdAndSeq = [(id, seq[startTrim:(len(seq)-endTrim)]) for id, seq in (alignedIdAndSeq, alignedHitIdAndSeq)]
    return alignedTrimmedIdAndSeq, alignedTrimmedHitIdAndSeq, divergencePredicate


def minimumDicts(dicts, key):
    '''
    dicts: list of dictionaries.
    key: a key present in every dict in dicts.
    returns: list of d in dicts, s.t. d[key] <= e[key] for every d, e in dicts.
    e.g.: [{'a':4, 'b':1}, {'a':5, 'b':0}, {'b': 0, 'a': 3}], 'b' -> [{'a':5, 'b':0} and {'b': 0, 'a': 3}] (not necessarily in that order)
    '''
    if not dicts:
        return []
    sortedDicts = sorted(dicts, key=lambda x: x[key])
    minValue = sortedDicts[0][key]
    return [d for d in sortedDicts if d[key] == minValue]


def computeOrthologs(queryFastaPath, subjectFastaPath, divEvalues, getForwardHits, getReverseHits, querySeqIds=None, workingDir='.'):
    '''
    queryFastaPath: fasta file path for query genome.
    subjectFastaPath: fasta file path for subject genome.
    divEvalues: list of (div, evalue) tuples.  orthologs for the given div and evalue combinations are computed.
    getForwardHits: a function mapping a query seq id to a list of subject genome blast hits.  see makeGetSavedHits() and makeGetHitsOnTheFly().
    getReverseHits: a function mapping a subject seq id to a list of query genome blast hits.  see makeGetSavedHits() and makeGetHitsOnTheFly().
    querySeqIds: a list of sequence ids for the query genome.  orthologs are only computed for those sequences.
      If False, orthologs are computed for every sequence in the query genome.
    workingDir: under workingDir, a temp directory is created, worked in (files and dirs created and deleted), and removed.
    returns: a mapping from (div, evalue) tuples to lists of orthologs.
    '''
    # optimization: internally swap query and subject if subject has fewer sequences than query and no querySeqIds were given.
    #   compute orthologs and unswap results.
    #   roundup time complexity is roughly linear in the number of sequences in the query genome.
    genomeSwapOptimization = True
    if not querySeqIds and genomeSwapOptimization and fasta.numSeqsInFastaDb(subjectFastaPath) < fasta.numSeqsInFastaDb(queryFastaPath):
        # print 'roundup(): subject genome has fewer sequences than query genome.  internally swapping query and subject to improve speed.'
        isSwapped = True
        # swap query and subject, forward and reverse
        queryFastaPath, subjectFastaPath = subjectFastaPath, queryFastaPath
        getForwardHits, getReverseHits = getReverseHits, getForwardHits
    else:
        isSwapped = False

    # make functions to look up a sequence from a sequence id.
    getQuerySeqFunc = makeGetSeqForId(queryFastaPath)
    getSubjectSeqFunc = makeGetSeqForId(subjectFastaPath)

    # if no querySeqIds were specified, get orthologs for every query sequence
    if not querySeqIds:
        querySeqIds = list(fasta.readIds(queryFastaPath))
        
    # get orthologs for every (div, evalue) combination
    with nested.NestedTempDir(dir=workingDir, nesting=0) as tmpDir:
        divEvalueToOrthologs = _computeOrthologsSub(querySeqIds, getQuerySeqFunc, getSubjectSeqFunc, divEvalues, getForwardHits, getReverseHits, workingDir)

    # if swapped query and subject genome, need to swap back the ids in orthologs before returning them.
    if isSwapped:
        swappedDivEvalueToOrthologs = divEvalueToOrthologs
        divEvalueToOrthologs = {}
        for divEvalue, swappedOrthologs in swappedDivEvalueToOrthologs.items():
            orthologs = [(query, subject, distance) for subject, query, distance in swappedOrthologs]
            divEvalueToOrthologs[divEvalue] = orthologs

    return divEvalueToOrthologs

    
def _computeOrthologsSub(querySeqIds, getQuerySeqFunc, getSubjectSeqFunc, divEvalues, getForwardHits, getReverseHits, workingDir):
    '''
    querySeqIds: a list of sequence ids from query genome.  Only orthologs for these ids are searched for.
    getQuerySeqFunc: a function that takes a seq id and returns the matching sequence from the query genome.
    getSubjectSeqFunc: a function that takes a seq id and returns the matching sequence from the subject genome.
    divEvalues: a list of (div, evalue) pairs which are thresholds for finding orthologs.  All pairs are searched simultaneously.
    getForwardHits: a function that takes a query seq id and a query seq and returns the blast hits in the subject genome.
    getReverseHits: a function that takes a subject seq id and a subject seq and returns the blast hits in the query genome.
    find orthologs for every sequence in querySeqIds and every (div, evalue) combination.
    return: a mapping from (div, evalue) pairs to lists of orthologs.
    '''
    # Note: the divs and evalues in divEvalues are strings which need to be converted to floats at the appropriate times below.
    
    # copy config files to working dir
    shutil.copy(MATRIX_PATH, workingDir)
    shutil.copy(CODEML_CONTROL_PATH, workingDir)

    divEvalueToOrthologs = dict(((div, evalue), list()) for div, evalue in divEvalues)
    maxEvalue = max(float(evalue) for div, evalue in divEvalues)
    maxDiv = max(float(div) for div, evalue in divEvalues)

    # get ortholog(s) for each query sequence
    for queryId in querySeqIds:
        if JIKE_DEBUG:
            print
            print 'forward\t{0}'.format(queryId)
        querySeq = getQuerySeqFunc(queryId)
        # get forward hits, evalues, alignments, divergences, and distances that meet the loosest standards of all the divs and evalues.
        # get forward hits and evalues, filtered by max evalue
        idSeqEvalueOfForwardHits = getGoodEvalueHits(queryId, querySeq, getForwardHits, getSubjectSeqFunc, maxEvalue)
        hitDataList = [{'hitId': hitId, 'hitSeq': hitSeq, 'hitEvalue': hitEvalue} for hitId, hitSeq, hitEvalue in idSeqEvalueOfForwardHits]
        # get alignments and divergences
        for hitData in hitDataList:
            (queryId, alignedQuerySeq), (hitId, alignedHitSeq), tooDivergedPred = getGoodDivergenceAlignedTrimmedSeqPair(queryId, querySeq, hitData['hitId'], hitData['hitSeq'], workingDir)
            hitData['alignedQuerySeq'] = alignedQuerySeq
            hitData['alignedHitSeq'] = alignedHitSeq
            hitData['tooDivergedPred'] = tooDivergedPred
        # filter by max divergence.
        hitDataList = [hitData for hitData in hitDataList if not hitData['tooDivergedPred'](maxDiv)]
        # get distances of remaining hits, discarding hits for which paml generates no rst data.
        distancesHitDataList = []
        for hitData in hitDataList:
            try:
                hitData['distance'] = getDistanceForAlignedSeqPair(queryId, hitData['alignedQuerySeq'], hitData['hitId'], hitData['alignedHitSeq'], workingDir)
                distancesHitDataList.append(hitData)
            except Exception as e:
                if e.args and e.args[0] == PAML_ERROR_MSG:
                    continue
                else:
                    raise
                
        # filter hits by specific div and evalue combinations.
        divEvalueToMinimumDistanceHitDatas = {}
        minimumHitIdToDivEvalues = {}
        minimumHitIdToHitData = {}
        for divEvalue in divEvalues:
            div, evalue = divEvalue
            # collect hit datas that pass thresholds.
            goodHitDatas = []
            for hitData in distancesHitDataList:
                if hitData['hitEvalue'] < float(evalue) and not hitData['tooDivergedPred'](float(div)):
                    goodHitDatas.append(hitData)
            # get the minimum hit or hits.
            minimumHitDatas = minimumDicts(goodHitDatas, 'distance')
            divEvalueToMinimumDistanceHitDatas[divEvalue] = minimumHitDatas
            for hitData in minimumHitDatas:
                minimumHitIdToDivEvalues.setdefault(hitData['hitId'], []).append(divEvalue)
                minimumHitIdToHitData[hitData['hitId']] = hitData # possibly redundant, since if two divEvalues have same minimum hit, it gets inserted into dict twice.  
        
        # get reverese hits that meet the loosest standards of the divs and evalues associated with that minimum distance hit.
        # performance note: wasteful or necessary to realign and compute distance between minimum hit and query seq?
        for hitId in minimumHitIdToHitData:
            if JIKE_DEBUG:
                print 'reverse\t{0}'.format(hitId)
            hitData = minimumHitIdToHitData[hitId]
            hitSeq = hitData['hitSeq']
            # since minimum hit might not be associated with all divs and evalues, need to find the loosest div and evalue associated with this minimum hit.
            maxHitEvalue = max(float(evalue) for div, evalue in minimumHitIdToDivEvalues[hitId])
            maxHitDiv = max(float(div) for div, evalue in minimumHitIdToDivEvalues[hitId])
            # get reverse hits and evalues, filtered by max evalue
            idSeqEvalueOfReverseHits = getGoodEvalueHits(hitId, hitSeq, getReverseHits, getQuerySeqFunc, maxHitEvalue)
            revHitDataList = [{'revHitId': revHitId, 'revHitSeq': revHitSeq, 'revHitEvalue': revHitEvalue} for revHitId, revHitSeq, revHitEvalue in idSeqEvalueOfReverseHits]
            # if the query is not in the reverese hits, there is no way we can find an ortholog
            if queryId not in [revHitData['revHitId'] for revHitData in revHitDataList]:
                continue
            for revHitData in revHitDataList:
                values = getGoodDivergenceAlignedTrimmedSeqPair(hitId, hitSeq, revHitData['revHitId'], revHitData['revHitSeq'], workingDir)
                (hitId, alignedHitSeq), (revHitId, alignedRevHitSeq), tooDivergedPred = values
                revHitData['alignedHitSeq'] = alignedHitSeq
                revHitData['alignedRevHitSeq'] = alignedRevHitSeq
                revHitData['tooDivergedPred'] = tooDivergedPred
            # filter by max divergence.
            revHitDataList = [revHitData for revHitData in revHitDataList if not revHitData['tooDivergedPred'](maxHitDiv)]
            # if the query is not in the reverese hits, there is no way we can find an ortholog
            if queryId not in [revHitData['revHitId'] for revHitData in revHitDataList]:
                continue
            # get distances of remaining reverse hits, discarding reverse hits for which paml generates no rst data.
            distancesRevHitDataList = []
            for revHitData in revHitDataList:
                try:
                    revHitData['distance'] = getDistanceForAlignedSeqPair(hitId, revHitData['alignedHitSeq'], revHitData['revHitId'], revHitData['alignedRevHitSeq'], workingDir)
                    distancesRevHitDataList.append(revHitData)
                except Exception as e:
                    if e.args and e.args[0] == PAML_ERROR_MSG:
                        continue
                    else:
                        raise


            # if passes div and evalue thresholds of the minimum hit and minimum reverse hit == query, write ortholog.
            # filter hits by specific div and evalue combinations.
            for divEvalue in minimumHitIdToDivEvalues[hitId]:
                div, evalue = divEvalue
                # collect hit datas that pass thresholds.
                goodRevHitDatas = []
                for revHitData in distancesRevHitDataList:
                    if revHitData['revHitEvalue'] < float(evalue) and not revHitData['tooDivergedPred'](float(div)):
                        goodRevHitDatas.append(revHitData)
                # get the minimum hit or hits.
                minimumRevHitDatas = minimumDicts(goodRevHitDatas, 'distance')
                if queryId in [revHitData['revHitId'] for revHitData in minimumRevHitDatas]:
                    divEvalueToOrthologs[divEvalue].append((queryId, hitId, hitData['distance']))
                    if JIKE_DEBUG:
                        print 'ortholog\t{0}\t{1}\t{2}'.format(queryId, hitId, hitData['distance'])

    return divEvalueToOrthologs


def computeOrthologsUsingOnTheFlyHits(queryFastaPath, subjectFastaPath, divEvalues, querySeqIds=None, workingDir='.'):
    '''
    Convenience function around computeOrthologs()
    querySeqIds: a list of sequence ids from query genome to find orthologs for.  If empty/falsy, will compute orthologs for every sequence in query genome.
    queryFastaPath: location and name of of fasta file and blast indexes of the query genome. e.g. /groups/rodeo/roundup/genomes/current/Homo_sapiens.aa/Homo_sapiens.aa
    subjectFastaPath: location and name of of fasta file and blast indexes of the subject genome.
    workingDir: a directory in which to create, use, and delete temporary files and dirs.
    This computes blast hits on-the-fly, so it slower than rounduPrecompute() for computing orthologs for full genomes.
    '''
    if JIKE_DEBUG:
        print
        print 'query_genome\t'+os.path.basename(queryFastaPath)
        print 'subject_genome\t'+os.path.basename(subjectFastaPath)

    # get blast hits using the least stringent evalue from among all the evalues in divEvalues.
    maxEvalue = str(max(float(evalue) for div, evalue in divEvalues))
    getForwardHits = makeGetHitsOnTheFly(subjectFastaPath, maxEvalue, workingDir)
    getReverseHits = makeGetHitsOnTheFly(queryFastaPath, maxEvalue, workingDir)
    divEvalueToOrthologs = computeOrthologs(queryFastaPath, subjectFastaPath, divEvalues, getForwardHits, getReverseHits, querySeqIds, workingDir)
    return divEvalueToOrthologs


def computeOrthologsUsingSavedHits(queryFastaPath, subjectFastaPath, divEvalues, forwardHitsPath, reverseHitsPath, querySeqIds=None, workingDir='.'):
    '''
    Convenience function around computeOrthologs()
    returns: a mapping from (div, evalue) pairs to lists of orthologs.
    '''    
    getForwardHits = makeGetSavedHits(forwardHitsPath)
    getReverseHits = makeGetSavedHits(reverseHitsPath)
    divEvalueToOrthologs = computeOrthologs(queryFastaPath, subjectFastaPath, divEvalues, getForwardHits, getReverseHits, querySeqIds, workingDir)
    return divEvalueToOrthologs
    

def writeToOutfile(orthologs, outfile):
    '''
    orthologs: a list of tuples of (queryid, subjectid, distance).
    outfile: where to write the orthologs
    write the orthologs to the outfile in the canonical format: one ortholog per line.  each line is tab-separated query id subject id and distance.
    '''
    data = ''.join(['%s\t%s\t%s\n'%(query, subject, distance) for query, subject, distance in orthologs])
    with open(outfile, 'w') as fh:
        fh.write(data)


if __name__ == '__main__':
    # rsd from the command line
    # format fasta files?
    # precompute blast hits or run for a few sequences?
    # --precompute
    # --seqidfile
    # --format
    # --workdir
    # -o <outfile>
    # -d, --divergence
    # -e, --evalue
    parser = argparse.ArgumentParser(description='Compute orthologs using the reciprocal smallest distance (RSD) algorithm between the query genome and the subject genome.  See http://bioinformatics.oxfordjournals.org/content/19/13/1710 for a description of the algorithm.')
    parser.add_argument('-q', '--query-genome', required=True, help='FASTA format sequence file, with unique ids on each nameline either in the form ">id" or ">ns|id|...".')
    parser.add_argument('-s', '--subject-genome', required=True, help='FASTA format sequence file, with unique ids on each nameline either in the form ">id" or ">ns|id|...".')
    parser.add_argument('-o', '--outfile', required=True, help='File in which to write orthologs.  Orthologs are written: query_id subject_id distance.')
    parser.add_argument('-d', '--divergence', type=float, default='0.8', help='Theshold for the maximum divergence allowed between a query and subject sequence.  A number > 0 and < 1. e.g. 0.2, or 0.5.  Default is 0.8.')
    parser.add_argument('-e', '--evalue', type=float, default='1e-5', help='Theshold for the maximum BLAST e-value allowed between a query and subject sequence.  e.g. 1e-20, or 0.005.  Default is 1e-5.')
    parser.add_argument('--ids', help='Path to file containing seq ids (one per line) in query_genome for which to compute orthologs.  If you only have one or a few sequences of interest it can be much faster to limit computation to those sequences.  The default is to compute othologs for all sequences in query_genome.  The sequence ids in the file must correspond to ids on the fasta namelines of query_genome.')
    parser.add_argument('--no-blast-cache', default=False, action='store_true', help='If this option is given, blast hits will not be precomputed for every sequence in each genome.  Using this option Can be faster if computing orthologs for only a few sequences.  Consider using in conjunction with --ids.')
    parser.add_argument('--no-format', default=False, action='store_true', help='If this option is given, genome fasta files will not be formatted for blast.  This is useful if blast formatted indices already exist and are located in the same directory as the fasta files.')
    parser.add_argument('--workdir', default='.', help='Directory under which to work.  will create a subdirectory under this dir in which to write temporary files, etc.  This subdirectory will be removed when rsd finishes.  Default is "."')
    parser.add_argument('-v', '--verbose', default=False, action='store_true')

    args = parser.parse_args()

    assert args.divergence > 0.0 and args.divergence < 1.0
    assert args.evalue >= 0.0

    if args.ids:
        with open(args.ids) as fh:
            # one id per line.  ignore blank lines and comment lines
            ids = [i for i in (line.strip() for line in fh) if i and not i.startswith('#')] 
    else:
        ids = None
    if args.verbose:
        print 'query sequence ids:', ids

    div = str(args.divergence)
    evalue = str(args.evalue)
    outfile = os.path.abspath(args.outfile)
    if args.verbose:
        print 'outfile:', outfile
        print 'divergence:', div
        print 'evalue:', evalue
    
    with nested.NestedTempDir(dir=args.workdir, nesting=0) as tmpDir:
        # format fasta files if needed.  do it in tmpDir to be clean.
        if args.no_format:
            # assume blast formatted index files coexist with the fasta files
            queryFastaPath = args.query_genome
            subjectFastaPath = args.subject_genome
            if args.verbose:
                print 'fasta files:', queryFastaPath, subjectFastaPath
        else:
            if args.verbose:
                print 'copying fasta files'
            queryFastaPath = os.path.join(tmpDir, os.path.basename(args.query_genome))
            subjectFastaPath = os.path.join(tmpDir, os.path.basename(args.subject_genome))
            if args.verbose:
                print 'fasta files:', queryFastaPath, subjectFastaPath
            shutil.copy(args.query_genome, queryFastaPath)
            shutil.copy(args.subject_genome, subjectFastaPath)
            if args.verbose:
                print 'formatting fasta files'
            formatForBlast(queryFastaPath)
            formatForBlast(subjectFastaPath)

        # print queryFastaPath, subjectFastaPath, outfile, div, evalue
        divEvalues = [(div, evalue)]

        # precompute blast hits if needed.  do it in tmpDir to be clean.
        if args.no_blast_cache: 
            getForwardHits = makeGetHitsOnTheFly(subjectFastaPath, evalue, tmpDir)
            getReverseHits = makeGetHitsOnTheFly(queryFastaPath, evalue, tmpDir)
            if args.verbose:
                print 'computing orthologs'
            divEvalueToOrthologs = computeOrthologsUsingOnTheFlyHits(queryFastaPath, subjectFastaPath, divEvalues, ids, tmpDir)
        else:
            if args.verbose:
                print 'precomuting blast hits'
            forwardHitsPath = os.path.join(tmpDir, 'query_subject.blast.hits.pickle')
            reverseHitsPath = os.path.join(tmpDir, 'subject_query.blast.hits.pickle')
            computeBlastHits(queryFastaPath, subjectFastaPath, forwardHitsPath, evalue, workingDir=tmpDir)
            computeBlastHits(subjectFastaPath, queryFastaPath, reverseHitsPath, evalue, workingDir=tmpDir)
            if args.verbose:
                print 'computing orthologs'
            divEvalueToOrthologs = computeOrthologsUsingSavedHits(queryFastaPath, subjectFastaPath, divEvalues, forwardHitsPath, reverseHitsPath, ids, tmpDir)

        if args.verbose:
            print 'writing {0} orthologs to outfile'.format(len(divEvalueToOrthologs[divEvalues[0]]))
        writeToOutfile(divEvalueToOrthologs[divEvalues[0]], outfile)

    
#################
# DEPRECATED CODE
#################

      
# do not cross this line...or else.

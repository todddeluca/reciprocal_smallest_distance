#!/usr/bin/env python

'''
general, widely applicable utilities
to fill in some gaps in the python language

ONLY DEPENDENCIES ON STANDARD LIBRARY MODULES ALLOWED.
'''

import datetime
import math
import hashlib # sha
import itertools
import tempfile
import shutil
import time
import os
import sys
import subprocess


def humanBytes(num):
    '''
    http://blogmag.net/blog/read/38/Print_human_readable_file_size
    num: a number of bytes.
    returns a human-readable version of the number of bytes
    Byte (B), Kilobyte (KB), Megabyte (MB), Gigabyte (GB), Terabyte (TB), Petabyte (PB), Exabyte (EB), Zettabyte (ZB), Yottabyte (YB)
    '''
    for x in ['B', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB']:
        if num < 1024.0:
            return "%3.1f%s" % (num, x)
        num /= 1024.0


def run(args, stdin=None, shell=False):
    '''
    for python 2.7 and above, consider using subprocess.check_output().
    
    args: commandline string treated as a one element list, or list containing command and arguments.
    stdin: string to be sent to stdin of command.
    shell: defaults to False to avoid shell injection attacks
    
    Basically, if you want to run a command line, pass it as a string via args, and set shell=True.
    e.g. 'ls databases/fasta'
    If you do not want shell interpretation, break up the commandline and args into a list and set shell=False.
    e.g. ['ls', 'databases/fasta']
    Runs command, sending stdin to command (if any is given).  If shell=True, executes command through shell,
    interpreting shell characters in command and arguments.  If args is a string, runs args like a command line run
    on the shell.  If shell=False, executes command (the first item in args) with the other items in args as arguments.
    If args is a string, it is executed as a command.  If the string includes arguments, strange behavior will ensue.
    This is a convenience function wrapped around the subprocess module.

    returns: stdout of cmd (as string), if returncode is zero.
    if returncode is non-zero, throws Exception with the 'returncode' and 'stderr' of the cmd as attributes.
    '''
    p = subprocess.Popen(args, shell=shell, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate(stdin)
    if p.returncode != 0:
        e = Exception('Error running command.  args='+str(args)+' returncode='+str(p.returncode)+'\nstdin='+str(stdin)+'\nstderr='+str(error))
        e.returncode = p.returncode
        e.stderr = error
        raise e
    else:
        return output


def dispatch(name, args=None, keywords=None):
    '''
    name: name of callable/function including modules, etc., e.g. 'foo_package.gee_package.bar_module.wiz_func'
    that can be imported from the current sys.path
    args: a list of arguments for the callable/function.
    keywords: a dict of keyword parameters for the callable/function.
    using the fully qualifed name, loads module, finds and calls function with the given args and keywords
    returns: the return value of the called callable.
    '''
    if args is None:
        args = []
    if keywords is None:
        keywords = {}
    modname, attrname = name.rsplit(".", 1)
    __import__(modname)
    mod = sys.modules[modname]
    func = getattr(mod, attrname)
    return func(*args, **keywords)


# remove at will.
def testing(msg='Hello, World.'):
    print msg
    return msg


def strToBool(value):
    '''
    An arbitrary set of human-readable strings is mapped to False.  Everything else is true.
    What is false? Ingoring case, 'F', 'FALSE', '0', '0.0', 'None'
    '''
    return str(value).upper() not in ('F', 'FALSE', '0', '0.0', 'NO', 'N', 'NONE')


def getBoolFromEnv(key, default=True):
    '''
    looks in os.environ for key.  If key is set to 'F', 'False', '0', 'N', 'NO', or some other falsy value (case insensitive), returns false.
    Otherwise, if key is set, returns True.  Otherwise returns the default, which defaults to True.
    '''
    if os.environ.has_key(key):
        return strToBool(os.environ[key])
    else:
        return default


###########################
# CONTEXT MANAGER UTITLITES
###########################
# Generic context manager utilities.  Useful for turning objects or object factories into context managers.

class ClosingFactoryCM(object):
    '''
    context manager for creating a new obj from a factory function when entering a context an closing the obj when exiting a context.
    useful, for example, for creating and closing a db connection each time.
    Calls obj.close() when the context manager exits.
    '''
    def __init__(self, factory):
        self.factory = factory
        self.obj = None
        
    def __enter__(self):
        self.obj = None
        self.obj = self.factory()
        return self.obj

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.obj is not None:
            self.obj.close()


class FactoryCM(object):
    '''
    context manager for creating a new object from a factory function.  Might be useful for getting an object from a pool (e.g. a db connection pool).
    '''
    def __init__(self, factory):
        self.factory = factory

    def __enter__(self):
        return self.factory()

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


class NoopCM(object):
    '''
    context manager for getting an object.  Useful when a context manager is required instead of a simple object.  e.g. to reuse an existing db connection.
    '''
    def __init__(self, obj):
        self.obj = obj

    def __enter__(self):
        return self.obj

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


#######################################
# 
#######################################

def truePred(*args, **keywords):
    return True


def retryErrorExecute(operation, args=[], keywords={}, pred=truePred, numTries=1, delay=0, backoff=1):
    '''
    pred: function takes Exception as arg, returns True to retry operation, False otherwise.  Default is to
      always return True.
    numTries: number of times to try, including the first time. numTries < 0 means try an infinite # of times
      if numTries == 0, does not execute operation.  Simply returns.
    delay: pause (in seconds) between tries.  default = 0 (no delay)
    backoff: delay is multiplied by this factor after every retry, so the length of successive delays are
      delay, delay*backoff, delay*backoff*backoff, etc. default = 1 (no backoff)
    execute operation.  if an exception occurs, pass it to the predicate.  if the
    predicate returns true, retry the operation if there are any tries left.  Otherwise, raise the exception.  
    '''
    # could make backoff a function, so delay = backoff(delay), for more flexibility than just an exponential relationship.
    error = None
    for i in xrange(numTries):
        try:
            return operation(*args, **keywords)
        except Exception, e:
            # re-raise exception if pred fails
            if not pred(e): raise
            # re-raise exception if that was the last try
            if i == (numTries-1): raise
            # else retry
            time.sleep(delay)
            delay *= backoff

# example:
# myFuncReturnValue = retryErrorExecute(myfunc, [param1, param2, param3], pred=customPred, numTries=10, delay=10, backoff=1.4)


##########################
# TEMPORARY DIRS AND FILES
##########################

def withTempDir(**keywords):
    '''
    keywords: accepts 'suffix', 'prefix', and 'dir' keywords with string values that are passed to tempfile.mkdtemp
    creates temporary directory, yields it, and deletes it after the yield returns or if an exception occurs.
    yields: path to a temporary directory.
    '''
    tempDir = None
    try:
        tempDir = tempfile.mkdtemp(**keywords)
        yield tempDir
        shutil.rmtree(tempDir)
    except:
        if tempDir:
            shutil.rmtree(tempDir)
        raise

def withOpenFile(name, mode='r', bufsize=-1):
    '''
    name, mode, and bufsize are the same as in the builtin function open().
    opens the file, yields it, and closes it after the yield returns, or after an exception occurs.
    yields: a filehandle for the open file.
    This implementation will not actually work correctly until code is ported to python2.5, since exception passing (control flow)
    is broken for yields prior to 2.5.
    '''
    fh = None
    try:
        fh = open(name, mode, bufsize)
        yield fh
        fh.close()
    except:
        if fh:
            fh.close()
        raise

################
# DATES AND TIME
################

def lastMonth(thisMonth=None):
    '''
    thisMonth: datetime or date obj.  defaults to today.
    returns: a date obj from the month before the month in thisMonth.  e.g. if this month is 2006/01/31, then 2005/12/01 is returned.
    '''
    if thisMonth == None:
        thisMonth = datetime.date.today()
    try:
        last = datetime.date(thisMonth.year, thisMonth.month-1, 1)
    except:
        last = datetime.date(thisMonth.year-1, 12, 1)
    return last


########
# RANDOM
########


class SimpleNamespace(object):
    '''
    use this if you want to instantiate an object to serve as a namespace.
    e.g. foo = SimpleNamespace(); foo.bar = 1; print foo.bar; # prints '1'
    '''
    pass


def mergeListOfLists(lists):
    '''
    lists: a list of lists
    returns: a list containing all the elements each list within lists
    '''
    merge = []
    for l in lists:
        merge.extend(l)
    return merge


def groupsOfN(iterable, n):
    '''
    iterable: some iterable collection
    n: length of lists to collect
    Iterates over iterable, returning lists of the next n element of iterable, until iterable runs out of elements.
    Last list may have less than n elements.
    returns: lists of n elements from iterable (except last list might have [0,n] elements.)
    '''
    seq = []
    count = 0
    it = iter(iterable)
    try:
        while 1:
            seq.append(it.next())
            count += 1
            if count == n:
                count = 0
                yield seq
                seq = []
    except StopIteration:
        if seq:
            yield seq


def splitIntoN(input, n, exact=False):
    '''
    input: a sequence
    n: the number of evenly-sized groups to split input into.  must be an integer > 0.
    exact: if True, returns exactly n sequences, even if some of them must be empty.
      if False, will not return empty sequences, so will return < n sequences when n > len(input).
    split input into n evenly sized sequences.  some sequences might have one less element than other sequences.
    first sequences returned have more elements than later sequences.
    
    e.g. if input had 10 elements, [1,2,3,4,5,6,7,8,9,10] and n=3, input would be split into these sequences: [1,2,3,4],[5,6,7],[8,9,10]
    e.g. if input had 2 elements, [1,2] and n=3, input would be split into these sequences: [1], [2], []
    e.g. if input had 0 elements, [] and n=3, input would be split into these sequences: [], [], []
    e.g. if exact=False and input had 2 elements, [1,2] and n=3, input would be split into these sequences: [1], [2]
    e.g. if exact=False and input had 0 elements, [] and n=3, input would be split into no sequences.  I.e. no sequences would be yielded.
    yields: n evenly sized sequences, or if exact=False, up to n evenly sized, non-empty sequences.
    '''
    size = len(input) // n
    numExtra = len(input) % n
    start = 0
    end = size
    for i in range(n):
        if i < numExtra:
            end += 1
        if not exact and start == end: # only empty sequences left, so exit early.
            break
        yield input[start:end]
        start = end
        end = end + size


def isInteger(num):
    '''
    num: possibly an integer
    What is an integer?  In this case it is a thingy which can be converted to an integer and a float
    and when done so the two values are equal.
    Returns true if successful, false otherwise
    '''
    try:
        int(num)
        return int(num) == float(num)
    except (TypeError, ValueError):
        return False
    
def isNumber(num):
    '''
    num: possibly a number
    Attempts to convert num to a number.
    Returns true if successful, false otherwise
    '''
    try:
        float(num)
        return True
    except (TypeError, ValueError):
        return False
    

def makeCounter(n=0, i=1):
    '''
    n: number to start counting from
    i: increment
    This could be implemented the normal way using a closure, but python does not support writing to variables in a closure,
    forcing one to implement the function by wrapping the counter variable in a list.  I guess a generator is more pythonic.
    usage: counter = makeCounter(1); counter.next(); counter.next(); # etc.
    returns: a generator object which counts up from n by i, starting with n, when next() is called on the generator.
    '''
    while 1:
        yield n
        n = n + i

        
########################
# DESCRIPTIVE STATISTICS
########################

def mean(nums):
    return float(sum(nums))/len(nums)


def variance(nums):
    m = mean(nums)
    return sum([(n - m)**2 for n in nums]) / float(len(nums))


def stddev(nums):
    return math.sqrt(variance(nums))


#######################################
# COMBINATORICS FUNCTIONS
#######################################

def permute(items, n):
    '''
    deprecated: use itertools.permutations()
    returns a list of lists every permutation of n elements from items
    '''
    return list(itertools.permutations(items, n))


def choose(items, n):
    '''
    deprecated: use itertools.combinations()
    items: a list
    returns: a list of lists of every combination of n elements from items
    '''
    return list(itertools.combinations(items, n))


def every(pred, seq):
    """ returns True iff pred is True for every element in seq """
    for x in seq:
        if not pred(x): return False
    return True


def any(pred, seq):
    """ returns False iff pred is False for every element in seq """
    for x in seq:
        if pred(x): return True
    return False


################################
# SERIALIZATION HELPER FUNCTIONS
################################

def loadObject(pickleFilename, protocol=-1):
    import cPickle
    fh = open(pickleFilename)
    obj = cPickle.load(fh)
    fh.close()
    return obj


def dumpObject(obj, pickleFilename, protocol=-1):
    import cPickle
    fh = open(pickleFilename, 'w')
    cPickle.dump(obj, fh, protocol=protocol)
    fh.close()
    return obj


################################
# FILE COMPARISION
################################


def writeToFile(data, filename, mode='w'):
    '''
    opens file, writes data, and closes file
    flushing used to improve consistency of writes in a concurrent environment.
    '''
    fh = open(filename, mode)
    fh.write(data)
    fh.flush()
    fh.close()


def readFromFile(filename, mode='r'):
    '''
    opens file, reads data, and closes file
    returns: contents of file
    '''
    fh = open(filename, mode)
    data = fh.read()
    fh.close()
    return data


def differentFiles(filename1, filename2):
    '''
    compares the contents of the two files using the SHA digest algorithm.
    returns: True if the contents of the files are different.  False otherwise.
    throws: an exception if either file does not exist.
    '''
    file1 = open(filename1)
    file2 = open(filename2)

    s1 = hashlib.sha1() # sha.new()
    s2 = hashlib.sha1() # sha.new()
    for l in file1:
        s1.update(l)
    for l in file2:
        s2.update(l)
    isDifferent = (s1.hexdigest() != s2.hexdigest())

    file1.close()
    file2.close()
    return isDifferent


if __name__ == '__main__':
    pass


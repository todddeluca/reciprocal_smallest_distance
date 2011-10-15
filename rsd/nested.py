'''
Code for creating nested directory structures to avoid having a single directory with millions of files.
Code for creating nested temp files and dirs.
'''

import os
import hashlib # sha
import uuid
import shutil

# SET THE DEFAULT TMP DIR ROOT.
if os.environ.has_key('NESTED_TMP_DIR'):
    DEFAULT_TMP_DIR = os.environ['NESTED_TMP_DIR']
elif os.environ.has_key('TMPDIR'):
    DEFAULT_TMP_DIR = os.environ['TMPDIR']
elif os.environ.has_key('TEMP'):
    DEFAULT_TMP_DIR = os.environ['TEMP']
elif os.environ.has_key('TMP'):
    DEFAULT_TMP_DIR = os.environ['TMP']
else:
    DEFAULT_TMP_DIR = os.getcwd()
DEFAULT_NESTED_LEVELS = int(os.environ.get('NESTED_LEVELS', 0)) # should be >= 0
DEFAULT_TMP_PREFIX = 'tmp'
DEFAULT_DIRS_MODE = 0777


########################################################
# FUNCTIONS FOR CREATING (NESTED) TEMPFILES AND TEMPDIRS
########################################################

class NestedTempDir(object):
    '''
    context manager for creating, using and deleting a temp dir using the 'with' statement.
    the absolute path to the created temp dir is 'returned' by the context manager.
    '''
    def __init__(self, dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS, suffix='', prefix=DEFAULT_TMP_PREFIX, dirsMode=DEFAULT_DIRS_MODE, dirsUseUmask=False):
        self.dir, self.nesting, self.suffix, self.prefix, self.dirsMode, self.dirsUseUmask = dir, nesting, suffix, prefix, dirsMode, dirsUseUmask
        self.path = None
        
    def __enter__(self):
        self.path = makeTempDir(dir=self.dir, nesting=self.nesting, suffix=self.suffix, prefix=self.prefix, dirsMode=self.dirsMode, dirsUseUmask=self.dirsUseUmask)
        return self.path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.path:
            shutil.rmtree(self.path)


def makeTempFile(dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS, suffix='', prefix=DEFAULT_TMP_PREFIX, dirsMode=DEFAULT_DIRS_MODE, dirsUseUmask=False, mode='w+b'):
    '''
    Creates a unique file under dir and any intermediate directories that do not yet exist.
    The caller is responsible for deleting the file.
    warning: permissions are set as if the file was opened with the builtin function open().
    warning: the caller is responsible for closing the returned file object.
    returns: fo, a file object open for writing.  the path is fo.name.
    '''
    path = makeTempPath(dir=dir, nesting=nesting, suffix=suffix, prefix=prefix, dirsMode=dirsMode, dirsUseUmask=dirsUseUmask)
    return open(path, mode=mode)


def makeTempDir(dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS, suffix='', prefix=DEFAULT_TMP_PREFIX, dirsMode=DEFAULT_DIRS_MODE, dirsUseUmask=False):
    '''
    dirsMode: mode/permissions of created nested directories.  e.g. 0775.
    creates a unique directory nested under dir and any intermediate dirs that do not yet exist.
    warning: permissions on the directory are as if the dir was created with os.mkdir().
    returns: absolute path to the created dir.
    '''
    path = makeTempPath(dir=dir, nesting=nesting, suffix=suffix, prefix=prefix, dirsMode=dirsMode, dirsUseUmask=dirsUseUmask)
    os.mkdir(path)
    return path


def makeTempPath(dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS, suffix='', prefix=DEFAULT_TMP_PREFIX, dirsMode=DEFAULT_DIRS_MODE, dirsUseUmask=False):
    '''
    Returns: absolute path to a unique name, nested under the tmp dir.  This name has not yet been created as a file or dir.
    '''
    name = prefix + uuid.uuid4().hex + suffix
    path = makeNestedPath(name, dir=dir, nesting=nesting, dirsMode=dirsMode, dirsUseUmask=dirsUseUmask)
    return path


def makeNestedPath(name, dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS, makeDirs=True, dirsMode=DEFAULT_DIRS_MODE, dirsUseUmask=False):
    '''
    name: basename of file or dir to nest.  should not end with a slash.  E.g. not 'foo/'.
    nesting: between 0 and 20.  Useful for keeping the number of files in a temp dir small.
    makeDirs: if True, creates any missing sub directories between dir and name.
    dirsMode: mode/permissions of created nested directories.  e.g. 0775.  see also dirsUseUmask.
    dirsUseUmask: By default False.  If True, nested dirs are created using umask of current user, instead of dirsMode.
    Given a name and a directory, hashes name to a path nested under dir, creates any nested dirs in path that do not yet exist,
    and sets the perms of those subdirs created to 02777 by default.
    returns: absolute path to name, including name.
    '''
    basename = os.path.basename(name)
    if not basename:
        raise Exception('Error generating basename for %s.  Name arg should not end in a slash.')

    path = makeNestedSeedDir(basename, dir=dir, nesting=nesting, makeDirs=makeDirs, dirsMode=dirsMode, dirsUseUmask=dirsUseUmask)
    path = os.path.join(path, basename)
    return path


def makeNestedSeedDir(seed, dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS, makeDirs=True, dirsMode=DEFAULT_DIRS_MODE, dirsUseUmask=False):
    '''
    seed: basename of file or dir to nest.  should not end with a slash.  E.g. not 'foo/'.
    nesting: between 0 and 20.  Useful for keeping the number of files in a temp dir small.
    makeDirs: if True, creates any missing sub directories between dir and seed.
    dirsMode: mode/permissions of created nested directories.  e.g. 0775.  see also dirsUseUmask.
    dirsUseUmask: By default False.  If True, nested dirs are created using umask of current user, instead of dirsMode.
    Given a seed name and a directory, hashes seed to a path nested under dir, creates any nested dirs in path that do not yet exist,
    and sets the perms of those subdirs created to 02777 by default.
    returns: absolute path to seed, including seed, e.g. /absolute/path/to/<dir>/a5/01/foo
    '''
    pathComponents = _getNestedComponents(seed, nesting)
    path = dir
    for component in pathComponents:
        path = os.path.join(path, component)
        if makeDirs and not os.path.exists(path):
            # attempt to avoid race conditions by double checking existence of path when an exception occurs.
            try:
                # mkdir automatically uses umask of the user.
                os.mkdir(path, dirsMode)
            except:
                if not os.path.exists(path):
                    raise
            if not dirsUseUmask:
                os.chmod(path, dirsMode)
    return os.path.abspath(path)


def getNestedPath(name, dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS):
    '''
    name: basename of file or dir to nest.  should not end with a slash.  E.g. not 'foo/'.
    Useful for checking if name exists, nested within dir, without creating any directories.
    does not create any directories, files, etc.  Simply hashes name to a path.
    returns: absolute path to nested path.
    '''
    basename = os.path.basename(name)
    if not basename:
        raise Exception('Error generating basename for %s.  Name arg should not end in a slash.')
    path = os.path.join(getNestedSeedDir(basename, dir=dir, nesting=nesting), basename)
    return path


def getNestedSeedDir(seed, dir=DEFAULT_TMP_DIR, nesting=DEFAULT_NESTED_LEVELS):
    '''
    returns: absolute path to nested seed dir
    '''
    components = _getNestedComponents(seed, nesting=nesting)
    return os.path.abspath(os.path.join(dir, *components))

    
def _getNestedComponents(seed, nesting=DEFAULT_NESTED_LEVELS):
    '''for a given seed, returns a list of path components the seed should be nested under.  e.g ['a1', '6f']'''
    # 40 char hexidecimal string, so maximum of 20 2-char nested dir levels.
    assert nesting >= 0
    assert nesting <= 20
    digest = hashlib.sha1(seed).hexdigest() # sha.new(seed).hexdigest()
    components = []
    for i in xrange(nesting):
        components.append(digest[(2*i):(2*i+2)])
    return components





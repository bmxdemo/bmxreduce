#! /usr/bin/python
#
#### /gpfs01/astro/packages/anaconda/default/bin/python
#
# Make simple and hardcode the directories. This assumes that the
# directories will not change.
#
# Need working directory - this will be where the lock file is placed
# need source and destination directories
#
# Call the rsync.py script as, for example
#    ./rsync.py -t /home/throwe/bmx -s data -d tgt
#
# This version runs on astroftp
#

import os
import os.path
import io
import sys
import subprocess
import time

workingDir = '/astro/u/bmx/bmxdaq'
lockFile = 'bmxrsync.lck'
sourceDir = 'bmxdaq/data'
destDir = '/astro/u/bmx/bmxdata/incoming'
rsync_py = '/astro/u/bmx/tgt/rsync.py'

maxLockAge = 3600

def touch(path):
    """
    Emulates the 'touch' command by creating the file at *path* if it does not
    exist.  If the file exist its modification time will be updated.
    """
    with io.open(path, 'ab'):
        os.utime(path, None)

# cd to topDir
try:
    os.chdir(workingDir)
except:
    print "Can not enter workingDir."
    sys.exit(1)

now = int(time.time())

if os.path.exists(lockFile):
    statinfo = os.stat(lockFile)
    if (now - statinfo.st_ctime) > maxLockAge:
        # notify someone.
        pass
    print "Lock file exists."
    sys.exit(1)

# If get here, create the lock file, run the rsync.py script and 
# then delete the lock file.
touch(lockFile)
try:
    # Need to do two hosts -m option to command
    #  -m bmxfer1 and -m bmxfer2
    # the "-r" option given below tells rsync to remove source files
    cmd = "%s -m bmxfer1 -t %s -s %s -d %s -r" % (rsync_py, workingDir, sourceDir, destDir,)
    #print "cmd: %s" % (cmd,)\
    #print "Running:",cmd
    #answer = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    answer = subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)
    #print "Answer:",answer
    # Now the second machine
    cmd = "%s -m bmxfer2 -t %s -s %s -d %s -r" % (rsync_py, workingDir, sourceDir, destDir,)
    #answer = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    answer = subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)
except:
    print "try error"
    pass

os.unlink(lockFile)


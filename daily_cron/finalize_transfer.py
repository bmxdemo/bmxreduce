#
# Script to walk the "incoming" directory from the BMX rsync and copy
# the files there to their final destination.
#
#  Files will be of the form YYMMDD_HHMM.data
#
# Want to go into
#
#  /gpfs01/astro/workarea/bmxdata/YYMMDD/YYMMDD_HHMM.data
#

import argparse
import os
import os.path

incomingDir = '/gpfs01/astro/workarea/bmxdata/incoming'
destDirTop = '/gpfs01/astro/workarea/bmxdata/raw'

parser = argparse.ArgumentParser(description='Process rynced files.')
parser.add_argument('-d','--dryrun', dest='dryrun', action='store_true',
                    help='Print out what would be done')
parser.add_argument('-v','--verbose', dest='verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()

# cd to incomingDir
try:
    os.chdir(incomingDir)
except:
    if args.verbose:
        print("cd failed.")
    sys.exit(1)

# The following assumes that the "incoming" dirstory has one level
for dirpath, dirnames, filenames in os.walk(incomingDir):
    for name in filenames:
        fullname=os.path.join(dirpath,name)
        # split name at the "_"
        dayDir = name.split('_')[0]
        try:
            int(dayDir)
        except ValueError:
            print("Do not recognize day ID",dayDir)
            continue
        monthDir=dayDir[:4]
        
        destDir = os.path.join(destDirTop,monthDir)
        rfiDir = os.path.join(destDirTop,monthDir,'rfi')
        ringDir = os.path.join(destDirTop,monthDir,'ring')

        # Verify destDir exists
        if  not os.path.exists(destDir):
            # Create it
            if args.dryrun:
                print("Would create directory %s" % (destDir,))
                print("Would create directory %s" % (rfiDir,))
                print("Would create directory %s" % (ringDir,))
            else:
                if args.verbose:
                    print("Creating %s" % (destDir))
                os.mkdir(destDir,0o755)
                os.mkdir(rfiDir,0o755)
                os.mkdir(ringDir,0o755)
        _,ext=os.path.splitext(name)
        destName = None
        if ext=='.data':
            destName = os.path.join(destDir,name)
        elif ext=='.rfi':
            destName = os.path.join(rfiDir,name)
        elif ext=='.ring':
            destName = os.path.join(ringDir,name)
        else:
            print("Unrecognized extension",ext)
            continue
        if args.dryrun:
            print("Would rename %s to %s" % (fullname, destName,))
        else:
            try:
                if args.verbose:
                    print("Renaming:",fullname,"->",destName)
                os.rename(fullname,destName)
            except:
                print("Cannot rename",fullname,"->",destName)


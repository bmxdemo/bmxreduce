import os
import string
import random
import numpy as np
from glob import glob
from datetime import datetime
import ntpath
import time
import subprocess



class farmit(object):

    def __init__(self, script, args=None, reqs=None, jobname=None, resubmit=False):
        """e.g. 
        f=farmit.farmit('test.py', jobname='test', args={'fmin':[1,2],
        'fmax':[5,6]}, reqs={'N':2})
        f.writejobfiles()
        f.runjobs()

        Values of argument dictonary must be arrays, even if only length 1
        """
        self.jobfilepath = os.getenv('HOME')+'/jobfiles/'

        if resubmit:
            #First arg is wildcard to existing job files. Just use these and do
            #nothing else
            self.jobfilenames = glob(script)
            return

        self.script = script
        self.args = args

        # Always forward X
        self.reqs = {'X':1}
        if reqs is not None:
            for k,val in enumerate(reqs):
                self.reqs[val]=reqs[val]

        if jobname is None:
            self.jobname = self.datestr()+'_'+self.randstring()
        else:
            self.jobname = jobname
        
        self.prepargs()

        self.getjobfilenames()
        self.runpath = os.getcwd()

        return

    def prepargs(self):
        """Make arguments all numpy arrays"""
        if self.args is not None:
            nvals = []
            for k,val in enumerate(self.args):
                x = self.args[val]
                if np.size(x) == 1:
                    x = np.reshape(np.array(x),(1))
                    self.args[val] = x
                    nvals.append(1)
                else:
                    x = np.array(x)
                    self.args[val] = x
                    nvals.append(x.size)
            self.njobs = np.max(np.array(nvals))
        else:
            self.njobs = 1

    def getjobfilenames(self):
        """Get job file names, sets self.jobfilename"""
        self.jobfilenames = []
        for k in range(self.njobs):
            fn = self.jobname + '_' + '{:04d}'.format(k) + '.job'
            self.jobfilenames.append(self.jobfilepath + fn)


    def randstring(self, size=4):
        """Generate random string of size size"""
        chars = string.ascii_uppercase + string.digits
        return ''.join(random.choice(chars) for _ in range(size))


    def datestr(self):
        """Get current date and time"""
        return datetime.now().strftime('%Y%m%d_%H%M%S')


    def getcmd(self, i=0):
        """Get command line command to issue. Args can come as an array of
        values, so use the ith value"""
        cmd = []
        cmd.append('cd ' + self.runpath)

        cmd0 = 'python '
        cmd0 += self.script
        if self.args is not None:
            for k,val in enumerate(self.args):
                cmd0 += ' -'
                cmd0 += val
                cmd0 += ' '
                if self.args[val].size > 1:
                    cmd0 += str(self.args[val][i])
                else:
                    cmd0 += str(self.args[val][0])

        cmd.append(cmd0)

        cmd.append('bash -c "if [ $? == "0" ]; then rm -f {0}; fi"'.format(self.jobfilenames[i]))

        return cmd


    def writejobfiles(self):
        """Write job files to $HOME/jobfiles/"""
        for k,val in enumerate(self.jobfilenames):

            f = open(val, 'w')
            f.write('command: |\n')
            cmd = self.getcmd(k)

            for j,cmd0 in enumerate(cmd):
                f.write('    ' + cmd0+'\n')

            if self.reqs is not None:
                for j,req in enumerate(self.reqs):
                    f.write(self.reqstring(req))
            f.close()

    def reqstring(self, req):
        return '{0}: {1}\n'.format(req, self.reqs[req])

    def runjobs(self, maxjobs=500):
        """Submit jobs"""

        if maxjobs is None:
            # Submit all jobs
            for k,val in enumerate(self.jobfilenames):
                self.submitjob(val)
        else:
            # Only submit a few at a time since wq cannot handle too many
            # submitted jobs, even if they are not running
            i = 0
            ntotal = len(self.jobfilenames)
            while True:
                njobs = self.getnjobs()
                nsubmit = maxjobs - njobs
                for k in range(nsubmit):
                    if i<ntotal:
                        self.submitjob(self.jobfilenames[i])
                        print('submitting job {0} of {1}'.format(i+1, ntotal))
                        i += 1
                    else:
                        return
                time.sleep(5)

    def submitjob(self, fn):
        """Submit a single job file"""
        cmd = 'wq sub -b {:s}'.format(fn)
        os.system(cmd)
    
    def waituntildone(self):
        """Wait until the jobs are done"""
        jfn = np.array(self.jobfilenames)
        for k,val in enumerate(jfn):
            jfn[k] = ntpath.basename(val)

        while True:
            ls = np.array(glob(self.jobfilepath+'*.job'))
            for k,val in enumerate(ls):
                ls[k] = ntpath.basename(val)
            njobs = len(np.intersect1d(jfn, ls))
            if njobs==0:
                break
            time.sleep(5)
        return

    def getnjobs(self):
        """Get number of running jobs"""
        res = subprocess.check_output('wq ls -u {:s}'.format(os.getlogin()), shell=True)
        ind1 = res.find('Jobs:')
        ind2 = res.find('Running:')
        njobs = int(res[(ind1+5):ind2])
        return njobs

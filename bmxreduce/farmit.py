
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

    def __init__(self, script, args=None, reqs=None, names=None, resubmit=False,
                 OMP_NUM_THREADS=None):
        """e.g. 
        f=farmit.farmit('test.py', name='test', args={'fmin':[1,2],
        'fmax':[5,6]}, reqs={'N':2})
        f.writejobfiles()
        f.runjobs()

        Values of argument dictonary must be arrays, even if only length 1
        """
        self.OMP_NUM_THREADS = OMP_NUM_THREADS
        self.jobfilepath = os.getenv('HOME')+'/jobfiles/'

        if not os.path.exists(self.jobfilepath):
            print ("Creating directory ",self.jobfilepath)
            os.makedirs(self.jobfilepath)

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

        self.prepargs()

        # Get jobnames
        if names is None:
            self.names = [self.datestr()+'_'+self.randstring() for k in range(self.njobs)]
        else:
            self.names = names

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
            fn = self.names[k] + '.job'
            self.jobfilenames.append(self.jobfilepath + fn)


    def randstring(self, size=6):
        """Generate random string of size size"""
        chars = string.ascii_uppercase + string.digits
        return ''.join(random.choice(chars) for _ in range(size))


    def datestr(self):
        """Get current date and time"""
        return datetime.now().strftime('%Y%m%d_%H%M%S')


    def getcmd(self, i=0, fname='logfile'):
        """Get command line command to issue. Args can come as an array of
        values, so use the ith value"""
        cmd = []
        
        if self.OMP_NUM_THREADS is not None:
            cmd.append('setenv OMP_NUM_THREADS {:d}'.format(self.OMP_NUM_THREADS))

        cmd.append('cd ' + self.runpath)

        cmd0 = 'python '
        cmd0 += self.script
        if self.args is not None:
            for k,val in enumerate(self.args):
                if len(val) == 1:
                    cmd0 += ' -'
                else:
                    cmd0 += ' --'
                cmd0 += val
                cmd0 += ' '
                if self.args[val].size > 1:
                    cmd0 += str(self.args[val][i]).replace(' ','').replace('[','').replace(']','')
                else:
                    cmd0 += str(self.args[val][0]).replace(' ','').replace('[','').replace(']','')
        cmd0 += ' >& {:s}.wqlog'.format(fname)
        cmd.append(cmd0)

        cmd.append('bash -c "if [ $? == "0" ]; then rm -f {0}; fi"'.format(self.jobfilenames[i]))

        return cmd


    def writejobfiles(self):
        """Write job files to $HOME/jobfiles/"""

        for k,val in enumerate(self.jobfilenames):

            f = open(val, 'w')
            f.write('command: |\n')
            cmd = self.getcmd(k, fname=val)

            for j,cmd0 in enumerate(cmd):
                f.write('    ' + cmd0+'\n')

            if self.reqs is not None:
                for j,req in enumerate(self.reqs):
                    f.write(self.reqstring(req))
            f.write('job_name: {0}\n'.format(self.names[k]))
            f.close()

    def reqstring(self, req):
        return '{0}: {1}\n'.format(req, self.reqs[req])

    def runjobs(self, maxjobs=512):
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
                #try:
                njobs = self.getnjobs()
                #except:
                #    # If error, don't submit more jobs and try again
                #    print('wq ls failed, trying again...')
                #    njobs = maxjobs
                nsubmit = maxjobs - njobs
                for k in range(nsubmit):
                    if i<ntotal:
                        self.submitjob(self.jobfilenames[i])
                        print('submitting job {0} of {1}'.format(i+1, ntotal))
                        time.sleep(0.2)
                        i += 1
                    else:
                        return
                time.sleep(5)

    def submitjob(self, fn):
        """Submit a single job file"""
        cmd = 'source ~/.bashrc; wq sub -b {:s} 2>&1'.format(fn)
        #cmd = 'nohup /astro/u/astrodat/local/bin/wq.exe sub {:s} 2>&1 >{:s}.wqlog &'.format(fn,fn)
        #cmd = 'nohup /astro/u/astrodat/local/bin/wq.exe sub {:s} 2>&1 > /dev/null &'.format(fn)
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
        res = subprocess.check_output('/astro/u/astrodat/local/bin/wq.exe ls -u {:s}'.format(os.getlogin()), shell=True)
        res=res.decode("utf-8") 
        ind1 = res.find('Jobs:')
        ind2 = res.find('Running:')
        njobs = int(res[(ind1+5):ind2])
        return njobs

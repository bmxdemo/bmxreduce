import numpy as np
import bmxdata
from datamanager import datamanager
from reduce_plot import genplots
import os


class reduce(object):

    def __init__(self, tag):
        """ Init with a tag string E.g.
        d = reduc('170928_1000')
        """

        # Initialize empty data manager so we can use its functions
        dm = datamanager()

        # Store tag string and get filename
        self.tag = tag
        self.rawfname = dm.getrawfname(tag)
        self.redfname = dm.getreducedfname(tag)
        self.d = bmxdata.BMXFile(self.rawfname)

        # Construct raw data time array (seconds from start)
        self.t = np.arange(self.d.nSamples)*self.d.deltaT
        
        # Now do reduction
        self.doreduce()


    def doreduce(self):
        """Do full reduction of a chunk of spectrometer data"""
        self.getind()
        self.maskdata()
        self.getcal()
        genplots(self)
        self.savedata()


    def getind(self):
        """Get array of start/stop indices for cal and data. Not fully tested
        for data with cal starting or stopping on first or last indices."""

        ind = self.d.data['lj_diode']
        
        # Huge kludge because indices don't synch up with data correctly 
        # (CDS 10/2/17)
        ind = np.roll(ind,1)
        ind[0] = ind[1]

        dind = ind - np.roll(ind,1)
        chind = np.where(dind != 0)[0]
        chind = chind[chind>0] # Change on first index is not true

        sc = [] # Start cal
        ec = [] # End cal
        sd = [] # Start data
        ed = [] # End data

        for k,val in enumerate(chind):
            if dind[val]==1:
                # Cal start index
                sc.append(val)
                ed.append(val-1)
                if len(sd)==0:
                    sd.append(0)
            else:
                # Cal end index
                ec.append(val-1)
                sd.append(val)
                if len(sc)==0:
                    sc.append(0)

        if len(ec)<len(sc):
            ec.append(self.d.nSamples-1)
        if len(ed)<len(sd):
            ed.append(self.d.nSamples-1)

        self.ind = {'sc':np.array(sc), 'ec':np.array(ec), 'sd':np.array(sd), 'ed':np.array(ed)}

        
    def getcal(self):
        """Get gain cal"""

        ncal = self.ind['sc'].size
        nf = self.d.freq[0].size

        self.Tcal = 2000 # Assume calibrator adds this power

        self.cal = np.zeros((self.d.nChan,ncal,nf)) # Mean of cal chunks
        self.dat = np.zeros((self.d.nChan,ncal,nf)) # Mean of surrounding data chunks
        self.calt = np.zeros(ncal)

        for k in range(ncal):
            # Cal start/stop indices
            s = self.ind['sc'][k]
            e = self.ind['ec'][k]

            # Data immediately after
            ind = np.where(self.ind['sd'] > e)[0]
            if len(ind) > 0:
                sd = self.ind['sd'][ind[0]]
                ed = self.ind['ed'][ind[0]]
            else:
                # Data immediately before
                ind = np.where(self.ind['ed'] < s)[0]
                sd = self.ind['sd'][ind[-1]]
                ed = self.ind['ed'][ind[-1]]
                        
            # Time is mean of cal times
            self.calt[k] = np.mean(self.t[s:(e+1)])

            for j in range(self.d.nChan):
                chn = self.getchname(j)
                self.cal[j,k,:] = np.mean(self.d.data[chn][s:(e+1),:],0)
                self.dat[j,k,:] = np.mean(self.d.data[chn][sd:(ed+1),:],0)


        # Get cal factor in ADU^2 / K
        self.g = (self.cal - self.dat)/self.Tcal

    def getchname(self,chanind):
        """Return channel name string given channel INDEX (i.e. starting from zero)"""
        return 'chan{:d}_0'.format(chanind+1)

    def maskdata(self):
        """Mask raw data based on 1st time deriv"""
        
        self.mask = np.ones((self.d.nChan,self.d.nSamples, self.d.freq[0].size))

        for j in range(self.d.nChan):
            chn = self.getchname(j)
            dx = self.d.data[chn] - np.roll(self.d.data[chn],1,axis=0)
            dx[0,:] = 0

            # Loop over frequencies
            for k in range(dx.shape[1]):
                v = dx[:,k]
                v = v-np.nanmean(v)
                
                # Get middle 80th percentile 
                p = np.percentile(v,[10,90])
                ind = (v>p[0]) & (v<p[1])
                std0 = np.nanstd(v[ind])

                # Mask 7 sigma outliers
                maskind = np.abs(v) > 7*std0

                # Unmark known calibrator on/off swings
                maskind[self.ind['sc']]   = 0
                maskind[self.ind['ec']+1] = 0

                # Expand by some number of samples on either side
                ker = np.ones(5)
                maskind = np.convolve(maskind.astype('float'), ker, 'same')
                maskind[maskind != 0] = 1
                maskind = maskind.astype('bool')

                self.mask[j, maskind, k] = 0

    def savedata(self):
        """Save reduced data after stripping raw data"""
        

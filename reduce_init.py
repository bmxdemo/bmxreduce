import numpy as np
import bmxdata
from datamanager import datamanager
from reduce_plot import genplots
import os
import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

# Initialize empty data manager so we can use its functions
dm = datamanager()

class reduce(object):

    def __init__(self, tag):
        """ Init with a tag string or data structure, e.g.
        r = reduce('170928_1000')
        or
        d = bmxdata.BMXFile('data/2017/170928_1000.data')
        r = reduce(d)
        """


        if type(tag) == bmxdata.BMXFile:
            # Raw data already loaded
            self.tag = dm.fname2tag(tag.fname)
            self.d = tag
            self.rawfname = dm.getrawfname(self.tag)
            self.redfname = dm.getreducedfname(self.tag)
        else:
            # Store tag string, get filename, and load data
            self.tag = tag
            self.rawfname = dm.getrawfname(self.tag)
            self.redfname = dm.getreducedfname(self.tag)
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
        self.calind = np.roll(ind,1)
        self.calind[0] = self.calind[1]

        dind = self.calind - np.roll(self.calind,1)
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

        # Calibrator ENR and any attenuation. Can be actually measured later
        self.ENR = 15 # dB
        self.calAtten = -10.0 # dB
        
        # Load cal port coupling
        x = dm.loadcsvbydate('S21_calport', self.tag)
        
        # Physical calibrator temperature. Should be measured temp as function
        # of time when we have it.
        self.Tcalphys = 290.0

        # Injected noise temperature. Interpolate measured cal port S21 onto
        # frequency grid
        f =  self.d.freq[0]
        fi = x[:,0]/1e6
        S21x = np.interp(f, fi, x[:,1])
        S21y = np.interp(f, fi, x[:,2])
        self.Tcalfast = np.zeros((self.d.nChan, self.d.nSamples, nf))

        # This is the excess temperature with the calibrator on
        self.Tcalfast[0, :, :] = self.Tcalphys*10**((S21x+self.calAtten)/10.)*(10**(self.ENR/10.)-1)
        self.Tcalfast[1, :, :] = self.Tcalphys*10**((S21y+self.calAtten)/10.)*(10**(self.ENR/10.)-1)
        

        ############################
        # Now take mean of cal and data

        self.cal = np.zeros((self.d.nChan,ncal,nf)) # Mean of cal chunks
        self.dat = np.zeros((self.d.nChan,ncal,nf)) # Mean of surrounding data chunks
        self.Tcal = np.zeros((self.d.nChan,ncal,nf)) # Mean of surrounding data chunks
        self.caltime = np.zeros(ncal)

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
            self.caltime[k] = np.mean(self.t[s:(e+1)])

            for j in range(self.d.nChan):
                chn = self.getchname(j)
                self.cal[j,k,:]  = np.nanmean(self.d.data[chn][s:(e+1),:],0)
                self.dat[j,k,:]  = np.nanmean(self.d.data[chn][sd:(ed+1),:],0)
                # Use mean of cal temperature during data, which has more samples 
                #and is less noisy. Really we probably want to smooth first.
                self.Tcal[j,k,:] = np.nanmean(self.Tcalfast[j,sd:(ed+1),:],0) 

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
                
                # Get middle 80th percentile of cal off data
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

    def getradec(self):
        """Get RA/Dec coordinates"""
        time = Time(self.d.data['mjd'], format='mjd')
        point = AltAz(alt=90*u.deg, az=0*u.deg, location=telescope_loc, obstime=time)
        sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
        self.ra = sky.ra
        self.dec = sky.dec


    def savedata(self):
        """Save reduced data after stripping raw data"""
        

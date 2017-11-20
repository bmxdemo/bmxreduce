import numpy as np
import bmxdata
from datamanager import datamanager
from reduce_plot import genplots
import os
import cPickle as cP
import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
import time

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

# Initialize empty data manager so we can use its functions
dm = datamanager()
dm.gettags()

class reduce(object):

    def __init__(self, tag):
        """ Init with a tag string E.g.
        r = reduce('170928_1000')
        """

        # Store tag string and get filename
        self.tag = tag
        self.rawfname = dm.getrawfname(tag)
        self.redfname = dm.getreducedfname(tag) 

        print('loading data...')
        t = time.time()
        self.d = bmxdata.BMXFile(self.rawfname)
        print('...took {:0.1f} sec'.format(time.time()-t))

        # Get frequency array
        self.f = self.d.freq[0]
        self.dtraw = self.d.deltaT[0]

        # Now concatenate on either side with data from adjacent tags if the
        # data are consecutive
        #self.paddata()


    def doreduce(self):
        """Do full reduction of a chunk of spectrometer data"""
        self.getind()

        # Plot unmasked raw data
        p = genplots(self)

        print('plotting raw waterfall...')
        t = time.time()
        p.plotrawwf(fext='_nomask')
        p.plotrawspec()
        print('...took {:0.1f} sec'.format(time.time()-t))

        print('masking data...')
        t = time.time()
        self.maskdata()
        print('...took {:0.1f} sec'.format(time.time()-t))

        print('plotting masked waterfall...')
        t = time.time()
        p.plotrawwf(fext='_mask')
        print('...took {:0.1f} sec'.format(time.time()-t))

        print('calibrating data...')
        t = time.time()
        self.getcal()
        self.applycal()
        print('...took {:0.1f} sec'.format(time.time()-t))

        print('plotting calibrated waterfall (log and lin)...')
        t = time.time()
        p.plotrawwf(fext='_cal')
        p.plotrawwf(fext='_callin', cscale='lin')
        print('...took {:0.1f} sec'.format(time.time()-t))

        print('downsampling and filtering data for plots...')
        t = time.time()
        self.downsample()
        self.meanfilter()
        self.svdfilter()
        print('...took {:0.1f} sec'.format(time.time()-t))

        print('plotting downsampled waterfall and power spectra...')
        t = time.time()
        p.plotcalwf()
        p.plotps()
        print('...took {:0.1f} sec'.format(time.time()-t))

        print('saving data...')
        self.savedata()
        print('...took {:0.1f} sec'.format(time.time()-t))

    def paddata(self):
        """Concatenate data before and after if it exists and is
        consecutive. Not sure if this currently works."""
        ind = np.where(dm.tags == self.tag)[0][0]
        
        # Store start/stop time of data
        self.mjdspan = [self.d.data['mjd'][[0,-1]]] 

        if ind >= 1:
            # Load tag before
            tag0 = dm.tags[ind-1]
            da = bmxdata.BMXFile(dm.getrawfname(tag0))
            self.d.data=np.hstack((self.d.data, da.data))

        self.d.nSamples=len(self.d.data)


    def getind(self):
        """Get array of start/stop indices for cal and data. Not fully tested
        for data with cal starting or stopping on first or last indices."""

        ind = self.d.data['lj_diode'].astype('bool')
        
        # Huge kludge because indices don't synch up with data correctly 
        # (CDS 10/2/17)
        self.calind = np.roll(ind,1)
        self.calind[0] = self.calind[1]

        dind = self.calind*1.0 - np.roll(self.calind,1)
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
        nf = self.f.size

        # Calibrator ENR and any attenuation. Can be actually measured later
        self.ENR = 15 # dB
        self.calAtten = 0.0 # dB
        
        # Load cal port coupling
        x = dm.loadcsvbydate('S21_calport', self.tag)
        
        # Physical calibrator temperature. Should be measured temp as function
        # of time when we have it.
        self.Tcalphys = 290.0

        # Injected noise temperature. Interpolate measured cal port S21 onto
        # frequency grid
        f =  self.f
        fi = x[:,0]/1e6
        #S21x = np.interp(f, fi, x[:,1])
        #S21y = np.interp(f, fi, x[:,2])
        S21x = -30.0
        S21y = -30.0
        self.Tcalfast = np.zeros((self.d.nChan, self.d.nSamples, nf))

        # This is the excess temperature with the calibrator on
        self.Tcalfast[0, :, :] = self.Tcalphys*10**((S21x+self.calAtten+self.ENR)/10.)
        self.Tcalfast[1, :, :] = self.Tcalphys*10**((S21y+self.calAtten+self.ENR)/10.)        

        ############################
        # Now take mean of cal and data

        self.cal = np.zeros((self.d.nChan,ncal,nf)) # Mean of cal chunks
        self.dat = np.zeros((self.d.nChan,ncal,nf)) # Mean of surrounding data chunks
        self.Tcal = np.zeros((self.d.nChan,ncal,nf)) # Mean of surrounding data chunks

        for k in range(ncal):
            # Cal start/stop indices
            s = self.ind['sc'][k]
            e = self.ind['ec'][k]

            # Data immediately before
            indb = np.where(self.ind['ed'] < s)[0]
            # Data immediately after
            inda = np.where(self.ind['sd'] > e)[0]

            if (len(inda) > 0) & (len(indb) > 0):
                sdb = self.ind['sd'][indb[-1]]
                edb = self.ind['ed'][indb[-1]]
                sda = self.ind['sd'][inda[0]]
                eda = self.ind['ed'][inda[0]]
            elif len(inda) > 0:
                # use after twice
                sda = self.ind['sd'][inda[0]]
                eda = self.ind['ed'][inda[0]]
                sdb=sda
                edb=eda
            elif len(indb) > 0:
                # use before twice
                sdb = self.ind['sd'][indb[-1]]
                edb = self.ind['ed'][indb[-1]]
                sda=sdb
                eda=edb
                        
            for j in range(self.d.nChan):
                chn = self.getchname(j)
                self.cal[j,k,:] = np.nanmean(self.d.data[chn][s:(e+1),:],0)
                self.dat[j,k,:] = 0.5 * (np.nanmean(self.d.data[chn][sda:(eda+1),:],0)+
                                         np.nanmean(self.d.data[chn][sdb:(edb+1),:],0))
                self.Tcal[j,k,:] = np.nanmean(self.Tcalfast[j,sdb:(eda+1),:],0)

        # Get cal factor in ADU^2 / K
        self.g = (self.cal - self.dat)/self.Tcal


    def applycal(self):
        """Apply calibration, for now mean of gain over time."""

        for j in range(self.d.nChan):

            # Mean of gain over time
            gmean = np.nanmean(self.g[j], 0)
            
            # ADU^2 -> K
            chn = self.getchname(j)
            self.d.data[chn] = self.d.data[chn]/gmean


    def getchname(self,chanind):
        """Return channel name string given channel INDEX (i.e. starting from zero)"""
        return 'chan{:d}_0'.format(chanind+1)


    def maskdata(self):
        """Mask raw data, data and cal data separately."""

        # Initialize mask array (1 for good, 0 for bad)
        mask = np.ones((self.d.nChan,self.d.nSamples, self.f.size)).astype('bool')

        # Create time array
        t0 = np.arange(self.d.nSamples)*self.d.deltaT

        # Load the frequency mask that applies to all times and all
        # channels. These frequencies will never make it through.
        mf = dm.loadcsvbydate('freqmask', self.tag)

        # Loop over channels
        for j in range(self.d.nChan):

            chn = self.getchname(j)
            x = self.d.data[chn]

            ###############
            # Mask Outliers
            # Loop over frequencies
            for k in range(x.shape[1]):

                # Initialize mask array
                maskind = np.zeros(x.shape[0]).astype('bool')

                # Do the masking step this many times
                nmask = 2
                for nm in range(nmask):
                
                    # Loop over cal on and cal off separately
                    for l in range(2):

                        if l==0:
                            ci = ~self.calind # Cal off data
                        else:
                            ci = self.calind # Cal on data

                        # Mask various outliers in v, dv/dt, etc.
                        for m in range(2):

                            if m==0:
                                # Get data itself
                                v = x[ci,k]
                            else:
                                # Get 1st deriv of data. There are gaps in the data for cal
                                # on/off but they're short so that's okay.
                                v = x[ci,k] - np.roll(x[ci,k], 1, axis=0)
                                v[0] = 0

                            # Times
                            t = t0[ci]

                            # Find where data is not already masked
                            gi = np.where(np.isfinite(v))
                            gt = t[gi]; gv = v[gi]
                            if len(gv) < 0.2*len(v):
                                # If more than 80% of the data is already NaN,
                                # mask the whole thing
                                maskind[ci] = True
                            else:
                                # Poly subtract data
                                p = np.polyfit(gt, gv, 4)
                                v = v - np.polyval(p, t)

                                # Get std of middle 80th percentile of values
                                p = np.nanpercentile(v,[10,90])
                                ind = (v>p[0]) & (v<p[1])
                                std0 = np.nanstd(v[ind])

                                # Mask 5 sigma outliers
                                maskind[ci] = (maskind[ci]) | (np.abs(v) > 5*std0)


                # Expand masked data by some number of samples on either side
                #ker = np.ones(2)
                #maskind = np.convolve(maskind.astype('float'), ker, 'same')
                #maskind[maskind != 0] = 1
                #maskind = maskind.astype('bool')

                # Mask first and last of cal ind
                maskind[self.ind['sc']] = True
                maskind[self.ind['ec']] = True
                maskind[self.ind['sd']] = True
                maskind[self.ind['ed']] = True

                mask[j, maskind, k] = False


            ###############################
            # Always mask these frequencies
            f  = self.f
            df = f[1] - f[0]
            fe = np.linspace(f[0]-df/2, f[-1]+df/2, f.size+1)
            felo = fe[0:-1]
            fehi = fe[1:]
            for k in range(mf.shape[0]):
                f0 = mf[k,0]
                f1 = mf[k,1]
                ind = np.where((fehi>=f0) & (felo<=f1))
                mask[j, :, ind] = False

            # Set masked data to NaN
            self.d.data[chn][~mask[j,:,:]] = np.nan

        ############
        # Store mask
        if hasattr(self,'mask'):
            self.mask = self.mask & mask
        else:
            self.mask = mask


    def getradec(self):
        """Get RA/Dec coordinates"""
        time = Time(self.d.data['mjd'], format='mjd')
        point = AltAz(alt=90*u.deg, az=0*u.deg, location=telescope_loc, obstime=time)
        sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
        self.ra = sky.ra
        self.dec = sky.dec


    def meanfilter(self):
        """Subtract mean over time from each frequency bin."""

        # Empty array
        self.data_mf = np.zeros(self.data.shape)

        # Mean filter
        for k,val in enumerate(self.data):
            self.data_mf[k] = val - np.nanmean(val,0)


    def svdfilter(self):
        """Singular value decompose data, zero first 10 modes, and convert
        back"""        
        self.data_svd = np.zeros(self.data.shape)

        # SVD filter
        for k,val in enumerate(self.data):
            # Pull out data
            v = self.data[k]*1.0

            # First mean filter
            v = v - np.nanmean(v, 0)

            # Normalize data to avoid RFI messing up SVD component separation
            norm = np.nanstd(v,0)**2

            # Prevent NaNs
            v[~np.isfinite(v)] = 0
            norm[~np.isfinite(norm)] = np.inf
            rat = v/norm
            rat[~np.isfinite(rat)]


            # SVD decomp
            U,s,V = np.linalg.svd(rat,full_matrices=True)

            # Now zero out first 10 components and convert back to data
            ss=np.zeros(v.shape)
            np.fill_diagonal(ss,s)
            ss[0:10,0:10] = 0
            vv = U.dot(ss.dot(V))

            # Now replace data
            norm[~np.isfinite(norm)] = np.nan
            self.data_svd[k] = vv* norm


    def downsample(self, dt=5.0):
        """Downsample data, dt in seconds."""
        
        self.dt = dt

        # MJD array has steps, make it not so
        mjd = self.d.data['mjd']
        mjd = np.linspace(mjd[0], mjd[-1], len(mjd))
        
        # Binned mjd edges
        mjdbin = np.arange(mjd[0], mjd[-1], dt/3600./24.)
        if mjdbin[-1] != mjd[-1]:
            np.append(mjdbin, mjd[-1])

        # Bin data
        nbin = len(mjdbin)-1
        self.data  = np.full((self.d.nChan, nbin, self.f.size), np.nan)
        self.var   = np.full((self.d.nChan, nbin, self.f.size), np.nan)
        self.nhits = np.full((self.d.nChan, nbin, self.f.size), np.nan)
        self.mjd   = np.full(nbin, np.nan)

        for k in range(nbin):
            ind = np.where( (mjd >= mjdbin[k]) & (mjd < mjdbin[k+1]) )[0]
            self.mjd[k] = np.nanmean(mjd[ind])

            for j in range(self.d.nChan):
                chn = self.getchname(j)
                dat = self.d.data[chn][ind, :]
                
                # Mask cal on data
                dat[self.calind[ind],:] = np.nan

                # Get stats
                mu = np.nanmean(dat, 0) # Mean
                var = np.nanvar(dat, 0) # Variance
                nhits = np.nansum(~np.isnan(dat), 0) # N hits
                
                # Require Nhits >= 50%
                ndat = dat.shape[0]
                badind = nhits < 0.5*ndat
                mu[badind] = np.nan
                var[badind] = np.inf

                # Store
                self.data[j, k, :]  = mu
                self.var[j, k, :]   = var
                self.nhits[j, k, :] = nhits


    def savedata(self):
        """Save reduced data after stripping raw data. Using numpy instead of
        cPickle because it seems to be much faster and will save time when
        coadding over many files, especially if doing so for many sim
        realizations."""
        np.savez(self.redfname, cal=self.cal, calAtten=self.calAtten,
                 dat=self.dat, data=self.data, dt=self.dt, ENR=self.ENR,
                 g=self.g, mjd=self.mjd, nhits=self.nhits, tag=self.tag,
                 Tcal=self.Tcal, Tcalphys=self.Tcalphys, f=self.f, dtraw=self.dtraw) 

import numpy as np
from glob import glob
import os
import astropy.units as u
from datetime import datetime, timedelta
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
import time
from datamanager import datamanager

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

class coaddbygroup(datamanager):
    
    def __init__(self, tags, mapdef, sn=None, fields=None):
        """Tag list over which to coadd, mapdef from mapmanager. If sn (serial
        number) is provided, we know this is a simulation and a list of fields
        must be provided (e.g. ['colore','gsync','hi4pi']). Ordering of fields
        does not matter, it gets sorted."""

        # Set up directories
        self.getdirs()

        # Save tags
        self.tags = tags
        self.m = mapdef

        # Do loading
        self.loadtags(sn, fields)
        self.getradec()
        self.windowdata()
        self.caldata()
        if sn is None:
            # Don't deglitch simulated noiseless data.
            self.deglitch()
        self.filterdata()

        return

    def loadtags(self, sn, fields):
        """Load tags, undo calibration, and concatenate"""
        for k,val in enumerate(self.tags):

            # Load
            if sn is None:
                fn = self.getreducedfname(val)
            else:
                fn = self.getreducedsimfname(val, sn, fields)
            print('loading {0}'.format(fn))
            x0 = np.load(fn)

            # Convert to dictionary so we can manipulate it
            x = {}
            for fld in x0.keys():
                x[fld] = x0[fld]

            # Undo gain calibration
            nchan = x['data'].shape[0]
            for l in range(nchan):
                x['data'][l] = x['data'][l] * np.nanmean(x['g'][l],0)

            # Concatenate
            for fld in ['data','g','mjd','nhits']:            
                if k==0:
                    setattr(self,fld,x[fld])
                else:
                    setattr(self,fld,np.hstack((getattr(self,fld),x[fld])))
        
        self.nchan = self.data.shape[0]
        self.f = x['f']

        return

    def getradec(self):
        """Get ra/dec"""
        time = Time(self.mjd, format='mjd')
        point = AltAz(alt=90*u.deg, az=0*u.deg, location=telescope_loc, obstime=time)
        sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
        self.ra = sky.ra.value
        self.dec = sky.dec.value
        return

    def windowdata(self):
        """Cut down data to map region"""
        indf = (self.f >= self.m['frange'][0]) & (self.f <= self.m['frange'][1])
        indt = (self.ra >= self.m['rarange'][0]) & (self.ra <= self.m['rarange'][1])
        indch = np.arange(self.nchan)

        self.data = self.data[np.ix_(indch,indt,indf)]
        self.g = self.g[:,:,indf]
        self.nhits = self.nhits[np.ix_(indch,indt,indf)]
        self.f = self.f[indf]
        self.mjd = self.mjd[indt]
        self.ra = self.ra[indt]
        self.dec = self.dec[indt]

    def caldata(self):
        """Apply calibration"""
        self.T = self.data
        for k in range(self.nchan):
            self.T[k] = self.T[k] / np.nanmean(self.g[k],0)

    def deglitch(self):
        """Remove outliers over time"""
        for k in range(self.nchan):
            v = self.data[k]
            for j in range(self.f.size):
                vv = v[:,j]*1.0

                if all(~np.isfinite(vv)):
                    continue
                
                # Only  mask on dv/dt for Galactic HI
                if self.f[j] < 1418.0:
                    vv = vv - np.nanmean(vv)
                    p = np.percentile(vv[np.isfinite(vv)],90)
                    sig = np.nanstd(vv[vv<=p])
                    maskind1 = vv>4*sig
                else:
                    maskind1 = np.zeros_like(vv).astype('bool')

                dvv = vv - np.roll(vv,1)
                if all(~np.isfinite(dvv)):
                    continue
                p = np.percentile(dvv[np.isfinite(dvv)],90)
                sig = np.nanstd(dvv[dvv<=p])
                maskind2 = np.abs(dvv)>4*sig

                # Join v and dv/dt mask index arrays
                maskind = (maskind1) | (maskind2)

                # Expand mask index
                nind = 10
                x = np.zeros_like(maskind)
                x[maskind] = 1.0
                x = np.convolve(x, np.ones(nind/2),'same')
                x[x>0] = 1.0
                maskind = x.astype('bool')

                # Apply mask
                v[maskind,j] = np.nan

    def filterdata(self):
        """Filter data"""
        
        # Compute template as median spectrum over time. Fit this to each
        # spectrum individually plus a polynomial
        self.mod = np.zeros_like(self.data)

        for k in range(self.nchan):
            v = self.data[k]
            # MEDIAN IS BAD!!! NON-LINEAR!!! WARNING!!! Need to improve glitch
            # cutting and RFI excision so we can use nanmean.
            #temp = np.nanmean(v,0)
            temp = np.nanmedian(v,0)
            
            for j in range(self.mjd.size):
                # Set up least sqares regression
                b = v[j]
                if all(~np.isfinite(b)):
                    continue
                a = np.vstack((temp,np.ones_like(temp),self.f,self.f**2,self.f**3))
                ind = np.isfinite(b)
                x = np.linalg.lstsq(a[:,ind].T,b[ind])
                self.mod[k,j,:] = (a.T).dot(x[0])

        return

import numpy as np
from glob import glob
import os
import astropy.units as u
from datetime import datetime, timedelta
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
import time
from datamanager import datamanager
from time import time
from sklearn.linear_model import Ridge
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

class coaddbygroup(datamanager):
    
    def __init__(self, tags, mapdef, sn=None, fields=None):
        """Tag list over which to coadd, mapdef from mapmanager. If sn (serial
        number) is provided, we know this is a simulation and a list of fields
        must be provided (e.g. ['colore','gsync','hi4pi']). Ordering of fields
        does not matter, it gets sorted."""

        # Set up directories
        self.getdirs()
        self.sn = sn

        # Save tags
        self.tags = tags
        self.m = mapdef

        # Do loading
        self.loadtags(fields)
        self.getradec()

        return

    def reduce(self, dodeglitch=True, dofilter=True, docpm=True):
        """Reduction steps"""
        self.windowdata()
        self.caldata()
        if self.sn is None:
            # Don't deglitch simulated noiseless data.
            self.deglitch(dodeglitch)
        self.filterdata(dofilter)
        self.binra()
        self.cpm(docpm)
        self.interpmod()

        return

    def loadtags(self, fields):
        """Load tags, undo calibration, and concatenate"""
        for k,val in enumerate(self.tags):

            # Load
            if self.sn is None:
                fn = self.getreducedfname(val)
            else:
                fn = self.getreducedsimfname(val, self.sn, fields)
            print('loading {0}'.format(fn))
            x0 = np.load(fn)

            # Convert to dictionary so we can manipulate it
            x = {}
            for fld in x0.keys():
                x[fld] = x0[fld]

            # fld gain calibration
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
        for k in range(self.nchan):
            self.data[k] = self.data[k] / np.nanmean(self.g[k],0)

    def deglitch(self, dodeglitch=True):
        """Remove outliers over time"""

        if ~dodeglitch:
            return

        for k in range(self.nchan):
            v = self.data[k]
            nh = self.nhits[k]
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
                nh[maskind,j] = 0

    def filterdata(self, dofilter=True):
        """Filter data"""
        
        # Compute template as median spectrum over time. Fit this to each
        # spectrum individually plus a polynomial
        self.mod = np.zeros_like(self.data)
        
        if ~dofilter:
            return

        for k in range(self.nchan):
            v = self.data[k]
            # MEDIAN IS BAD!!! NON-LINEAR!!! WARNING!!! Need to improve glitch
            # cutting and RFI excision so we can use nanmean.
            #temp = np.nanmean(v,0)
            temp = np.nanmedian(v,0)
            
            # Deglitch templates
            dtdx = temp - np.roll(temp,1)
            p = np.nanpercentile(dtdx,99)
            ind = np.abs(dtdx) > p
            temp[ind] = np.interp(self.f[ind], self.f[~ind], temp[~ind])

            for j in range(self.mjd.size):
                # Set up least sqares regression
                b = v[j]*1.0
                
                # Don't fit galactic HI
                b[(self.f>1420.3) & (self.f<1420.8)] = np.nan

                if all(~np.isfinite(b)):
                    continue

                ff = self.f - self.f.mean()
                ff = ff/np.max(ff)
                a = np.vstack((temp,np.ones_like(temp),ff,ff**2,ff**3))
                ind = np.isfinite(b) & np.isfinite(a).all(0)
                if all(~ind):
                    continue
                x = np.linalg.lstsq(a[:,ind].T,b[ind])
                self.mod[k,j,:] = (a.T).dot(x[0])
        return

    def cpm(self, docpm=True):
        """CPM model of data"""
        
        # Initialize
        self.modcpm = np.zeros_like(self.data)

        if ~docpm:
            return

        rr=Ridge(alpha=1, normalize=False, fit_intercept=False)
        print('ridge regression alpha = {:0.2E}'.format(rr.alpha))

        # Time axis in minutes
        t = (self.mjd-self.mjd[0])*24*60

        # Don't use any data within dt (minutes) and df (MHz) of target datum
        # for regression 
        dt = 30.0
        df = 5.0

        fi = np.arange(self.f.size)

        for k in range(self.nchan):

            v = self.data[k] - self.mod[k]
            
            for i,t0 in enumerate(t):

                doindt = np.where(np.abs(t-t0)>dt)[0]
                print(i)

                for j,f0 in enumerate(self.f):

                    s = time()

                    doindf = np.where(np.abs(self.f - f0)>df)[0]

                    # Construct regressor
                    #X = v[doindt][:,doindf]
                    
                    # Data to fit
                    #y = v[doindt][:,j]

                    # Equivalent to above but faster?
                    X0 = v.take(doindt,0)
                    X = X0.take(doindf,1)
                    y = X0.take(j,1)

                    # Only fit if data exists
                    ind = np.isfinite(y)
                    y = y[ind]
                    X = X[ind,:]

                    # Don't allow NaNs in regressor

                    X[~np.isfinite(X)] = 0

                    # If there's nothing to fit, skip
                    if len(y) == 0:
                        continue

                    # Do regression
                    #b = np.linalg.lstsq(X, y)[0]
                    rr.fit(X, y); b=rr.coef_

                    # Get prediction
                    X = v[i,doindf]
                    X[~np.isfinite(X)] = 0
                    self.modcpm[k,i,j] = X.dot(b)

                    # Timing
                    e = time()
                    #print('regress {:f}'.format(e-s))


    def interpmod(self):
        """Interpolate model over 1420 line"""

        ind = (self.f>1420.2) & (self.f<1422.0)
        for k in range(self.nchan):
            mod = self.modcpm[k]
            fi = interp2d(self.f[~ind], self.mjd, mod[:,~ind], kind='linear')
            self.modcpm[k][:,ind] = fi(self.f[ind], self.mjd)

    def binra(self):
        """Bin in RA"""
        
        databin = np.zeros((self.data.shape[0], len(self.m['ra']), self.data.shape[2]))
        nhitsbin = np.zeros((self.data.shape[0], len(self.m['ra']), self.data.shape[2]))
        modbin = np.zeros((self.data.shape[0], len(self.m['ra']), self.data.shape[2]))
        mjdbin = np.zeros(len(self.m['ra']))

        for k in range(len(self.m['ra'])):
            ralo = self.m['rabe'][k]
            rahi = self.m['rabe'][k+1]
        
            ind = np.where( (self.ra>=ralo) & (self.ra<rahi))[0]

            mjdbin[k] = np.mean(self.mjd[ind])

            for j in range(self.data.shape[0]):
                databin[j, k, :] = np.nanmean(self.data[j][ind, :],0)
                if hasattr(self,'mod'):
                    modbin[j, k, :] = np.nanmean(self.mod[j][ind, :],0)
                if hasattr(self,'nhits'):
                    nhitsbin[j,k,:] = np.nansum(self.nhits[j][ind, :],0)

        self.mod = modbin
        self.data = databin
        self.nhits = nhitsbin
        self.mjd = mjdbin
        self.ra = self.m['ra']
        
        if len(self.dec)==0:
            self.dec = []
        else:
            self.dec = np.ones_like(self.ra)*self.dec[0]

    def save(self):
        """Save"""
        fdir = 'maps/bmx/{:s}/'.format(self.m['mapdefn'])
        if  not os.path.isdir(fdir):
            os.mkdir(fdir)
        fn = '{:s}/{:s}_map.npz'.format(fdir,self.tags[0][0:6])

        np.savez(fn, data=self.data, dataroot=self.dataroot, dec=self.dec,
                 f=self.f, g=self.g, m=self.m, mjd=self.mjd, mod=self.mod,
                 modcpm=self.modcpm, nchan=self.nchan, nhits=self.nhits, ra=self.ra,
                 reducedroot=self.reducedroot,
                 reducedsimroot=self.reducedsimroot, tags=self.tags)

    def load(self, tag):
        """Load"""

        x = np.load(tag) 
        self.data = x['data']
        self.dataroot = x['dataroot']
        self.dec = x['dec']
        self.f = x['f']
        self.g = x['g']
        self.m = x['m'].item()
        self.mjd = x['mjd']
        self.mod = x['mod']
        self.modcpm = x['modcpm']
        self.nchan = x['nchan']
        self.nhits = x['nhits']
        self.ra = x['ra']
        self.reducedroot = x['reducedroot']
        self.reducedsimroot = x['reducedsimroot']
        self.tags = x['tags']



class coaddbyday(coaddbygroup):

    def __init__(self, fn, dodeglitch=True):
        """Initialize with list of filenames:
        e.g. 
        fn = sort(glob('maps/bmx/*.npz'))
        coaddbyday(fn)
        """
        
        for k,val in enumerate(fn):

            print('coadding {:s}'.format(val))

            self.load(val)
            self.data = self.data - self.mod - self.modcpm
            if dodeglitch:
                self.deglitch()
            self.getw()
            self.docoadd()
            self.finalize()


    def deglitch(self):
        """2nd round deglitching"""

        for k in range(self.nchan):

            v = self.data[k]*1.0

            for j in range(3):
                
                if j<2:
                    # Deriv
                    vv = v - np.roll(v,1,j)
                else:
                    vv = v
                p = np.nanpercentile(vv,90)
                sig = np.nanstd(vv[np.abs(vv)<p])
                ind = np.abs(vv) > 7*sig
                
                # Don't cut galactic 21-cm
                ind[:, (self.f>1420.2) & (self.f<1420.8)] = False

                v[ind] = np.nan

                #plt.clf()
                #plt.imshow(v);plt.colorbar()
                #plt.draw()
                #raw_input("Press Enter to continue...")


            # Cut any row/column with >50% NaNs
            t = 0.5
            fracgoodx = np.sum(np.isfinite(v),0)*1.0 / v.shape[0]
            fracgoody = np.sum(np.isfinite(v),1)*1.0 / v.shape[1]
            v[fracgoody<t,:] = np.nan
            v[:,fracgoodx<t] = np.nan

            # Put back
            self.data[k] = v
            
    def getw(self):
        """Get weight as simply 1/var where var is computed over frequency at
        each time"""
        
        self.var0 = np.zeros_like(self.data)
        self.w0 = np.zeros_like(self.data)

        nf = self.data.shape[2]

        for k in range(self.nchan):
            v = self.data[k]
            self.var0[k,:,:] = np.repeat(np.nanvar(v, 1)[:,np.newaxis],nf,1)
            self.w0[k,:,:] = np.repeat(1/np.nanvar(v, 1)[:,np.newaxis],nf,1)
            
        self.w0[~np.isfinite(self.data)] = 0


    def docoadd(self):
        """Do the coaddition"""

        if not hasattr(self,'w'):
            # First loop, initialize
            self.w = np.zeros_like(self.data)
            self.wvar = np.zeros_like(self.data)
            self.wv = np.zeros_like(self.data)
            

        self.var0[~np.isfinite(self.var0)] = 0
        self.data[~np.isfinite(self.data)] = 0

        self.w += self.w0
        self.wvar += self.w0 * self.var0
        self.wv += self.w0 * self.data

    def finalize(self):
        """Finish"""
        
        delattr(self, 'w0')
        delattr(self, 'var0')
        delattr(self, 'data')

        self.T = self.wv / self.w
        self.var = self.wvar / self.w

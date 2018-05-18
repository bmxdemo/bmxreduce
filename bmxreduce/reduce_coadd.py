import numpy as np
from glob import glob
import os
import astropy.units as u
from datetime import datetime, timedelta
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
import time
from mapmanager import mapmanager
from time import time
from sklearn.linear_model import Ridge
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import sys

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

class coaddbygroup(mapmanager):
    
    def __init__(self, tags, mapdef, sn=None, fields=None):
        """Tag list over which to coadd, mapdef from mapmanager. If sn (serial
        number) is provided, we know this is a simulation and a list of fields
        must be provided (e.g. ['colore','gsync','hi4pi']). Ordering of fields
        does not matter, it gets sorted."""

        # Set up directories
        self.getdirs()
        self.sn = sn
        self.fields = fields

        # Save tags
        self.tags = tags
        self.getmapdefn(mapdef)

        # Do loading
        self.loadtags()
        self.getradec()

        # Horrible horrible kludge to deal with data with Xpol (vertical) on
        # channel 2 and Ypol (horizontal) on channel 1. Other times must be cut
        # in cuttimes.csv.
        if np.int(self.tags[0][0:6]) <= 180206:
            self.reversechannels()

        return

    def reversechannels(self):
        """Reverse channel order"""
        self.data = self.data[::-1]
        self.g = self.g[::-1]
        self.nhits = self.nhits[::-1]

    def reduce(self, dodeglitch=True, dofilter=True, docpm=True, cpmdtmin=30.0,
               cpmdtmax=np.inf, cpmdf=10.0, cpmalpha=1e-4):
        """Reduction steps"""
        self.windowdata()
        self.caldata()
        if self.sn == 'real':
            # Only deglitch real data
            self.deglitch(dodeglitch)
        self.filterdata(dofilter)
        self.binra()
        self.cpm(docpm, cpmdtmin=cpmdtmin, cpmdtmax=cpmdtmax, cpmdf=cpmdf, cpmalpha=cpmalpha)
        self.interpmod()
        if self.sn == 'real':
            self.getcutstats()
        return

    def loadtags(self):
        """Load tags, undo calibration, and concatenate"""
        for k,val in enumerate(self.tags):

            # Load
            if self.sn == 'real':
                fn = self.getreducedfname(val)
            else:
                fn = self.getreducedsimfname(val, self.sn, self.fields)
            print('loading {0}'.format(fn))
            sys.stdout.flush()
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
            
            # Simopts
            if (k==0) & (self.sn != 'real'):
                self.simopts = x['simopts'].item()

        self.nchan = self.data.shape[0]
        self.f = x['f']

        return

    def getradec(self):
        """Get ra/dec"""

        print('calculating ra/dec...')
        sys.stdout.flush()

        time = Time(self.mjd, format='mjd')
        point = AltAz(alt=90*u.deg, az=0*u.deg, location=telescope_loc, obstime=time)
        sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
        self.ra = sky.ra.value
        self.dec = sky.dec.value
        return

    def windowdata(self):
        """Cut down data to map region"""
        
        print('windowing data...')
        sys.stdout.flush()

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
        print('applying calibration')
        sys.stdout.flush()
        for k in range(self.nchan):
            self.data[k] = self.data[k] / np.nanmean(self.g[k],0)

    def deglitch(self, dodeglitch=True):
        """Remove outliers over time"""

        if ~dodeglitch:
            return

        print('deglitching...')
        sys.stdout.flush()

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
        
        if not dofilter:
            return

        print('poly + template flitering data...')

        for k in range(self.nchan):
            v = self.data[k]
            # MEDIAN IS BAD!!! NON-LINEAR!!! WARNING!!! Need to improve glitch
            # cutting and RFI excision so we can use nanmean.
            temp = np.nanmean(v,0)
            #temp = np.nanmedian(v,0)
            
            for j in range(self.mjd.size):
                # Set up least sqares regression
                b = v[j]*1.0
                
                # Don't fit galactic HI
                b[(self.f>1420.3) & (self.f<1420.8)] = np.nan

                if all(~np.isfinite(b)):
                    continue

                ff = self.f - self.f.mean()
                ff = ff/np.max(ff)
                #a = np.vstack((temp,np.ones_like(temp),ff,ff**2,ff**3))
                a = np.vstack((np.ones_like(temp),ff,ff**2,ff**3))
                ind = np.isfinite(b) & np.isfinite(a).all(0)
                if all(~ind):
                    continue
                x = np.linalg.lstsq(a[:,ind].T,b[ind])
                self.mod[k,j,:] = (a.T).dot(x[0])
                
            temp = np.nanmean(self.data[k]-self.mod[k], 0)
            self.mod[k] += np.repeat(temp[np.newaxis,:], self.mjd.size, 0)
            
        return

    def cpm(self, docpm=True, cpmdtmin=30.0, cpmdtmax=np.inf, cpmdf=10.0, cpmalpha=1e-4):
        """CPM model of data.
        Only use data separated from target datum by more than dt[0] and less
        than dt[1] (minutes). Don't use any data within df (MHz) of target
        datum frequency. cpmalpha is alpha param for ridge regression.
        """
        
        # Initialize
        self.cpmdtmin = cpmdtmin
        self.cpmdtmax = cpmdtmax
        self.cpmdf = cpmdf
        self.cpmalpha = cpmalpha
        self.modcpm = np.zeros_like(self.data)

        if not docpm:
            return


        rr=Ridge(alpha=self.cpmalpha, normalize=True, fit_intercept=False)
        print('CPM ridge regression alpha = {:0.2E}'.format(rr.alpha))
        sys.stdout.flush()

        # Time axis in minutes
        t = (self.mjd-self.mjd[0])*24*60

        for k in range(self.nchan):

            v = self.data[k] - self.mod[k]
            
            for i,t0 in enumerate(t):
                
                absdt = np.abs(t-t0)
                doindt = np.where( (absdt>=self.cpmdtmin) & (absdt<=self.cpmdtmax))[0]

                print('{:d} of {:d} time steps'.format(i,len(t)))
                sys.stdout.flush()

                dt = [] # timing
                for j,f0 in enumerate(self.f):

                    s = time()

                    doindf = np.where(np.abs(self.f - f0)>self.cpmdf)[0]

                    # Construct regressor
                    X = v[doindt][:,doindf]
                    
                    # Data to fit
                    y = v[doindt][:,j]

                    # Autoregressive
                    #Xauto = np.zeros((len(doindt),len(doindt)))
                    #for jj,val in enumerate(doindt):
                    #    Xauto[][s:e] = v[s:e, j]

                    # Concatenate
                    #X = np.hstack((X,Xauto))

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
                    dt.append(e-s)
                    
                print('mean time per freq. is {:f} s'.format(np.array(dt).mean()))
                sys.stdout.flush()

    def interpmod(self):
        """Interpolate model over 1420 line"""

        print('interpolating CPM model over 1420 line...')
        sys.stdout.flush()

        ind = (self.f>1419.5) & (self.f<1421.5)

        for k in range(self.nchan):
            mod = self.modcpm[k]
            fi = interp2d(self.f[~ind], self.mjd, mod[:,~ind], kind='linear')
            y = fi(self.f[ind], self.mjd)
            self.modcpm[k][:,ind] = y

    def binra(self):
        """Bin in RA"""
        
        print('binning in RA...')
        sys.stdout.flush()

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

    def getmapfname(self):
        """Get map filename"""

        fdir = 'maps'
        if self.sn is None:
            sn = 'real'
        else:
            sn = self.sn
        fld = self.m['mapdefn']
        day = self.tags[0][0:6]

    def getcutstats(self):
        """Get cut statistics"""

        print('calculating cut statistics...')
        sys.stdout.flush()

        self.cutstat = []

        for k in range(self.nchan):

            cs = {}

            # 5th, 50th, and 95th percentile of median T data
            x = np.nanmedian(self.data[k],0)
            cs['medT05'] = np.nanpercentile(x, 5)
            cs['medT50'] = np.nanpercentile(x, 50)
            cs['medT95'] = np.nanpercentile(x, 95)

            # 5th, 50th, and 95th percentile of all T data
            x = self.data[k]
            cs['T05'] = np.nanpercentile(x, 5)
            cs['T50'] = np.nanpercentile(x, 50)
            cs['T95'] = np.nanpercentile(x, 95)
            
            # Max median gain (below 99th percentile)
            # Min gain in last 100 freq bins
            x = np.nanmedian(self.g[k],0)
            cs['maxmedg'] = np.nanmax(x[x<np.nanpercentile(x,99)])
            cs['minmedg'] = np.nanmin(x[-100:])

            # 5th, 50th, and 95th percentile of all g
            x = self.g[k]
            cs['g05'] = np.nanpercentile(x, 5)
            cs['g50'] = np.nanpercentile(x, 50)
            cs['g95'] = np.nanpercentile(x, 95)
            
            # std of the middle 95 percent of data
            x = self.data[k]
            ind = (x>=np.nanpercentile(x,2.5)) & (x<=np.nanpercentile(x,97.5))
            cs['std'] = np.nanstd(x[ind])

            self.cutstat.append(cs)

    def save(self):
        """Save"""
        
        print('saving data')
        sys.stdout.flush()

        fdir = 'maps/bmx/{:s}/{:s}/'.format(self.sn, self.m['mapdefn'])
        if not os.path.isdir(fdir):
            os.makedirs(fdir)

        fn = '{:s}/{:s}_map.npz'.format(fdir,self.tags[0][0:6])
        np.savez(fn, data=self.data, dataroot=self.dataroot, dec=self.dec,
                 f=self.f, m=self.m, mjd=self.mjd, mod=self.mod,
                 modcpm=self.modcpm, nchan=self.nchan, nhits=self.nhits, ra=self.ra,
                 reducedroot=self.reducedroot,
                 reducedsimroot=self.reducedsimroot, tags=self.tags)

        if self.sn == 'real':
            fn = '{:s}/{:s}_cutstat.npz'.format(fdir,self.tags[0][0:6])
            np.savez(fn, cutstat=self.cutstat)


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


    def load(self, fn):
        """Load"""

        x = np.load(fn) 
        y = dict(x)
        x.close()

        for k,val in enumerate(y.keys):
            setattr(self,val,y[val])

        if hasattr(self, simopts):
            self.simopts = self.simopts.item()

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

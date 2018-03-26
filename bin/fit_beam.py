#!/usr/bin/env python
import sys
sys.path+=['.','..','bmxreduce']
import bmxreduce as br
import numpy as np
from numpy import sin,cos
import matplotlib.pyplot as plt
import glob
from scipy.optimize import minimize
##
## Ok, this first try is now intentionally hackish, we'll do a better
## job later.
##

## Let's load data, extract some frequencies and then see where we stand

fmin,fmax = 1215.,1240.
#fmin,fmax = 1300.,1400

flist=sorted(glob.glob('data/reduced/1803/18030[4,5]*.npz'))

mjd=[]
sig=[]

for fn in flist[:24]:
    print fn
    da=np.load(fn)
    cmjd=da['mjd']
    wm=np.where((np.isfinite(cmjd)))
    wf=np.where((da['f']>fmin)&(da['f']<fmax))
    dcut=(da['data'][0])[:,wf]
    dcut=np.nanmean(dcut[:,0,:],axis=1)
    mjd.append(cmjd)
    sig.append(dcut)

mjd=np.hstack(mjd)
sig=np.hstack(sig)

## get ephemerides
e=br.sat_ephem()
sats=e.get_tracks(mjd,max_zen_distance=20*np.pi/180)

if True:

    for prn, ndx, dla, dlng in sats:
        cls=0.3/np.sqrt(dla**2+dlng**2).min()
        print cls
        plt.subplot(2,1,1)
        plt.plot(mjd[ndx],3.*np.ones(len(mjd[ndx])),label="#"+str(prn),lw=min(cls,10))
        plt.subplot(2,1,2)
        plt.plot(mjd[ndx],3.*np.ones(len(mjd[ndx])),label="#"+str(prn),lw=min(cls,10))


    print mjd.shape, sig.shape
    plt.subplot(2,1,1)
    plt.plot(mjd,sig-85.)
    plt.ylabel('sig-85')
    plt.legend()

    plt.subplot(2,1,2)
    plt.plot(mjd,sig-85.)
    plt.ylabel('sig-85')
    plt.semilogy()
    plt.xlabel('time [mjd]')
    plt.show()


## now, let's model


meanx=0.0
meany=0.0
sigmaxx=3**2
sigmayy=sigmaxx
sigmaxy=0.0
ofs=85.

deg2rad=1/180.*np.pi

def msig(x):
    meanx,meany,sxx,syy,sxy,ofs=x[:6]
    amps=x[6:6+len(sats)]
    meanx*=deg2rad
    meany*=deg2rad
    sxx*=deg2rad**2
    sxy*=deg2rad**2
    syy*=deg2rad**2
    msig=ofs*np.ones(len(mjd))
    
    for i,(prn, ndx, dla, dlng) in enumerate(sats):
        dy=dla-meany
        dx=dlng-meanx
        if (sxy**2>sxx*syy*0.99):
            sxy=sxx*syy*0.99
        #print sxx,syy,sxy, meanx,meany
        beam=np.exp(-0.5*((dx*dx*syy+dy*dy*sxx+2*dx*dy*(-sxy))/(sxx*syy-sxy**2)))
        #if (np.isnan(beam).any()):
        #    continue
        msig[ndx]+=amps[i]*beam
    return msig


def chi2(x):
    mod=msig(x)
    if (np.isnan(mod).any()):
        return 1e30
    chi2=np.nansum((sig-mod)**2)
    print chi2
    return chi2


xg=[meanx,meany,sigmaxx,sigmayy,sigmaxy,ofs]+list(300*np.ones(len(sats)))

mod=msig(xg)
plt.plot(mjd,sig,label='data')
plt.plot(mjd,mod,label='first trial')
#xopt=minimize(chi2,xg,method='powell')
xopt=minimize(chi2,xg)
mod=msig(xopt.x)
plt.plot(mjd,mod,label='optimized')
plt.legend()

meanx,meany,sxx,syy,sxy,ofs=xopt.x[:6]
print meanx,meany
print np.sqrt(sxx), np.sqrt(syy), sxy/np.sqrt(sxx*syy)

plt.show()

        










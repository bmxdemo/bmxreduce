import matplotlib.pyplot as plt
import numpy as np
import reduce_init
from astropy.time import Time

plt.ioff()

def nanhist(x, **kwargs):
    plt.hist(x[np.isfinite(x)],**kwargs)

class genplots():
    
    def __init__(self, r, fext=''):
        """Takes reduce_init.reduce object as input, generates reduc plots"""

        # Reduced data
        self.r = r
        
        # Choose this (must be float)
        self.dpi = 80.0

        # Get useful info
        self.fmin = self.r.f[0]
        self.fmax = self.r.f[-1]
        self.tmin = 0
        self.tmax = self.r.d.deltaT[0]*self.r.d.nSamples / 60. # minutes

        self.ts = Time(self.r.d.data['mjd'][0],  format='mjd')
        self.te = Time(self.r.d.data['mjd'][-1], format='mjd')

        # File naming
        self.filebase = self.r.tag

        # Close all open plots
        plt.close('all')

        return

    def plotrawwf(self, dochan=None, fext='', cscale='log'):
        """Plot raw waterfall plot. Takes chan index."""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)
            v = self.r.d.data[chn]

            # Size of data
            sz = np.array(v.shape)

            # Let's make a figure of a standard size but with hugely increased
            # dpi. This will hopefully keep labels looking readable when zoomed
            # out. 
            padx = 0.25
            pady = 0.25
            figsz = (1+np.array([padx,pady]))*sz/self.dpi
            fac = 20.0
            fig = plt.figure(figsize=figsz/fac, dpi=self.dpi*fac)

            if cscale=='log':
                plt.imshow(np.log10(v.T), extent=(self.tmin,self.tmax,self.fmax,self.fmin))
                titlab = 'log10[ADU^2]'
            else:
                plt.imshow(v.T, extent=(self.tmin,self.tmax,self.fmax,self.fmin))
                plt.clim(0,200)
                titlab = 'ADU^2'

            plt.grid('on')
            plt.xlabel('t (minutes)');
            plt.ylabel('freq (MHz)');
            plt.title('raw {:s} adc ({:s}) -- {:s} - {:s} (UTC)'.format(chn,titlab,self.ts.iso,self.te.iso))
            plt.colorbar(pad=0)

            fname = self.filebase + '_'+chn+'_wfraw'+fext+'.jpg'
            plt.savefig(fname, dpi=self.dpi*fac, bbox_inches='tight')
            plt.close(fig)

    def plotrawspec(self, dochan=None):
        """Plot raw time average spectrum"""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)
            v = self.r.d.data[chn]

            # Get median spectrum
            soff = np.nanmedian(v[~self.r.calind,:],0)
            son = np.nanmedian(v[self.r.calind,:],0)
            
            # Plot
            fig = plt.figure(figsize=(7,7))
            
            plt.subplot(2,1,1)
            plt.semilogy(self.r.f, soff, label='cal off')
            plt.semilogy(self.r.f, son, label='cal on')
            plt.xlabel('f (MHz)')
            plt.ylabel('ADU^2')
            plt.grid('on')
            plt.title('Median raw spectrum, {:s}, {:s}'.format(self.r.tag,chn))
            plt.legend()

            plt.subplot(2,1,2)
            plt.plot(self.r.f, son-soff, label='on - off')
            plt.xlabel('f (MHz)')
            plt.grid('on')
            yl=plt.ylim()
            plt.ylim((0,yl[1]))
            plt.legend()

            fname = self.filebase + '_'+chn+'_specraw.png'
            plt.savefig(fname, bbox_inches='tight')
            plt.close(fig)


    def plotcalwf(self, dochan=None):
        """Plot calibrated, downsampled waterfall plot."""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)

            for k in range(3):
                if k==0:
                    v = self.r.data[chan].T
                    titlab = 'T (K)'
                    cl = (0,200)
                    fext = 'data'
                if k==1:
                    v = np.log10(self.r.g[chan].T)
                    titlab = 'log10(gain) (ADU^2/K)'
                    cl = None
                    fext = 'gain'
                if k==2:
                    v = (self.r.nhits[chan]/self.r.var[chan]).T
                    titlab = 'weight = nhits/variance (1/K^2)'
                    cl = None
                    fext = 'weight'

                if k==0:
                    # Size of data
                    sz = np.array(v.shape)
                    # Let's make a figure of a standard size but with hugely increased
                    # dpi. This will hopefully keep labels looking readable when zoomed
                    # out. 
                    padx = 0.25
                    pady = 0.25
                    figsz = (1+np.array([padx,pady]))*sz/self.dpi
                    fac = 2.0

                # Plot waterfall
                fig = plt.figure(figsize=figsz/fac, dpi=self.dpi*fac)
                plt.imshow(v.T, extent=(self.fmin,self.fmax,self.tmax,self.tmin),aspect='auto')
                if cl is not None:
                    plt.clim(*cl)

                plt.grid('on')
                plt.ylabel('t (minutes)');
                plt.xlabel('freq (MHz)');
                plt.title('{:s} {:s} -- {:s} - {:s} (UTC)'.format(chn,titlab,self.ts.iso,self.te.iso))
                plt.colorbar(pad=0)

                fname = self.filebase + '_'+chn+'_wfcal_'+fext+'.png'
                plt.savefig(fname, dpi=self.dpi*fac, bbox_inches='tight')
                plt.close(fig)

                # Plot median
                fig = plt.figure(figsize=(7,5))
                plt.plot(self.r.f, np.nanmedian(v,1))
                plt.xlabel('freq (MHz)');
                plt.ylabel(titlab)
                plt.title('{:s} median {:s}'.format(chn,titlab,self.r.tag))
                if cl is not None:
                    plt.ylim(*cl)

                fname = self.filebase + '_'+chn+'_medcal_'+fext+'.png'
                plt.savefig(fname, bbox_inches='tight')
                plt.close(fig)
                

    def plotvariance(self, dochan=None):
        """Plot variance and expected variance from radiometer equation"""
        
        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)

            # Actual variance
            v = self.r.var[chan]*1.0
            v[~np.isfinite(v)]=np.nan
            v = np.nanmean(v,0)

            # Radiometer equation
            df = (self.r.f[1]-self.r.f[0])*1e6 # In Hz
            dt = self.r.dtraw # In sec
            T = np.nanmean(self.r.data[chan], 0)
            dT = T/np.sqrt(dt*df)

            # Plot
            fig = plt.figure(figsize=(7,7))

            plt.subplot(2,1,1)
            plt.plot(self.r.f, np.sqrt(v), label='std(T)')
            plt.plot(self.r.f, dT, label='T/sqrt(df*dt)')
            plt.xlabel('freq (MHz)');
            plt.ylabel('Kelvin')
            plt.title('{:s} {:s}'.format(chn,self.r.tag))
            plt.legend()
            plt.ylim(0,10)

            plt.subplot(2,1,2)
            plt.plot(self.r.f, np.sqrt(v)/dT, label='ratio (std(T)/[T/sqrt(df*dt)])')
            plt.xlabel('freq (MHz)');
            plt.ylabel('ratio')
            plt.legend()
            plt.ylim(0,3)
            plt.plot([self.r.f[0],self.r.f[-1]],[1,1],':k')

            fname = self.filebase + '_'+chn+'_variance.png'
            plt.savefig(fname, bbox_inches='tight')
            plt.close(fig)
            


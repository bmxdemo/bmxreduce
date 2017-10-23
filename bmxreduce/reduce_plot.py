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
        self.fmin = self.r.d.freq[0][0]
        self.fmax = self.r.d.freq[0][-1]
        self.tmin = 0
        self.tmax = self.r.d.deltaT[0]*self.r.d.nSamples / 60. # minutes

        self.ts = Time(self.r.d.data['mjd'][0],  format='mjd')
        self.te = Time(self.r.d.data['mjd'][-1], format='mjd')

        # File naming
        self.filebase = self.r.tag

        # Close all open plots
        plt.close('all')

        return

    def plotrawdata(self, dochan=None, fext=''):
        """Plot raw data, maksed and unmasked. Takes chan index."""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif ~iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)
            v = self.r.d.data[chn]

            # Size of data
            sz = np.array(v.shape)

            # Let's make a figure of a standard size but with hugely increased
            # dpi. This will hopefully keep labels looking readable when zoomed
            # out. 
            padx = 0.2
            pady = 0.25
            figsz = (1+np.array([padx,pady]))*sz/self.dpi
            fac = 20.0
            fig = plt.figure(figsize=figsz/fac, dpi=self.dpi*fac)

            ax = plt.Axes(fig, [padx/2 - padx/4, pady/2+.02, 1-padx+padx/4, 1-pady])
            fig.add_axes(ax)
            imax = ax.imshow(np.log10(v.T), extent=(self.tmin,self.tmax,self.fmax,self.fmin))
            ax.grid('on')
            imax.axes.set_xlabel('t (minutes)');
            imax.axes.set_ylabel('freq (MHz)');
            imax.axes.set_title('raw adc (log10[ADU^2]) -- {:s} - {:s} (UTC)'.format(self.ts.iso,self.te.iso))

            cpad = 0.01
            cbaxes = fig.add_axes([1-padx/2+cpad, pady/2, padx/2-cpad-padx/4-.02, 1-pady]) 
            cax = fig.colorbar(imax, cax=cbaxes)

            fname = self.filebase + '_'+chn+'_wfraw'+fext+'.png'

            plt.savefig(fname, dpi=self.dpi*fac)
            plt.close(fig)

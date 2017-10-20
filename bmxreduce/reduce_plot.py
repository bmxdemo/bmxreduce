import matplotlib.pyplot as plt
import numpy as np
import reduce_init

plt.ioff()

def nanhist(x, **kwargs):
    plt.hist(x[np.isfinite(x)],**kwargs)

class genplots():
    
    def __init__(self, r):
        """Takes reduce_init.reduce object as input, generates reduc plots"""

        # Reduced data
        self.r = r
        
        # Choose this (must be float)
        self.dpi = 80.0

        # Close all open plots
        plt.close('all')

        return

    def plotrawdata(self, chan):
        """Plot raw data, maksed and unmasked. Takes chan index."""
        
        chn = self.r.getchname(chan)
        v = self.r.d.data[chn]
        
        # Size of data
        sz = np.array(v.shape)

        # Let's make a figure of the right size with maybe some padding all around
        padx = 0.025
        pady = 0.1
        figsz = (1+np.array([padx,pady]))*sz/self.dpi
        fig = plt.figure(figsize=figsz, dpi=self.dpi)
        
        ax = plt.Axes(fig, [padx/2, pady/2, 1-padx, 1-pady])
        fig.add_axes(ax)
        ax.imshow(v.T)
        plt.savefig('test.png', dpi=self.dpi)

        plt.close(fig)

from glob import glob
import numpy as np
import matplotlib.pyplot as plt

def getcuts(mapdir, gettags=False):
    """Plot all cuts contained in directory mapdir,
    i.e. plotcuts('maps/bmx/real/SDSSlowz/')"""

    # Get cutstat files
    fn = np.sort(glob(mapdir + '/*_cutstat.npz'))

    tags = []
    for k,val in enumerate(fn):

        if gettags:
            # Get first tag, which is annoyingly saved only in the map file
            x = np.load(val.replace('cutstat','map'))
            tags.append(x['tags'][0])
            x.close()

        # Load cuts
        x = np.load(val)
        c = x['cutstat']
        x.close()

        # Number of channels
        nchan = len(c)


        if k==0:
            cnames = list(c[0].keys())
            cc = c
            for n in range(nchan):
                for cn in cnames:
                    cc[n][cn] = np.array([cc[n][cn]])
        else:
            for n in range(nchan):
                for cn in cnames:
                    cc[n][cn] = np.append(cc[n][cn], c[n][cn])

    if gettags:
        tags = np.array(tags)
        return cc, tags
    else:
        return cc

def getclim():
    """Get cut limits"""

    cl = {}
    cl['std']         = [0, np.inf]
    cl['medT05']      = [35., 80.]
    cl['medT50']      = [40., 77.]
    cl['medT95']      = [45., 90.]
    cl['T05']         = [35., 80.]
    cl['T50']         = [42., 75.]
    cl['T95']         = [40., 100.]
    cl['g05']         = [5e10, 3e11]
    cl['g50']         = [1e11, 5e11]
    cl['g95']         = [2e11, 8e11]
    cl['minmedg']     = [5e10, 3e11]
    cl['maxmedg']     = [2e11, 8e11]



    return cl

def plotcuts(c, clim=None):
    """Plot cuts"""

    cnames = list(c[0].keys())
    cnames.sort()
    cnames = cnames[0:6] + cnames[7:10] + [cnames[6]] + cnames[10:]
    nchan = len(c)

    plt.close(1)
    plt.figure(1,figsize=(12,10))

    plt.close(2)
    plt.figure(2, figsize=(12,10))

    # Initialize mask
    ndays = len(c[0][cnames[0]])*1.0
    mask = np.ones((nchan,ndays), dtype=bool)

    for k,cn in enumerate(cnames):

        for n in range(nchan):
            plt.figure(1)
            plt.subplot(4,3,k+1)

            y = c[n][cn]
            plt.plot(y,'.', label='chan{:d}'.format(n))
            plt.title(cn)
            if clim is not None:
                lo = clim[cn][0]
                hi = clim[cn][1]
                plt.ylim(lo/5.,hi*5)
                xl = plt.xlim()
                plt.plot(xl, [lo,lo],':k')
                plt.plot(xl, [hi,hi],':k')

                # Update mask
                mask[n] = (mask[n]) & (y>=lo) & (y<=hi)

            if k==0:
                plt.legend()                

            # Print some info about cumulative cuts
            if clim is not None:
                remfrac = len(np.where(mask[n])[0]) / ndays
                st = 'chan{:d} pass frac = {:0.1f}%'.format(n, remfrac*100)
                plt.text(0.03, 0.98-n*0.08, st, transform=plt.gca().transAxes, fontsize=10,
                         verticalalignment='top', horizontalalignment='left')


            plt.figure(2)
            plt.subplot(4,3,k+1)
            if clim is None:
                plt.hist(c[n][cn], label='chan{:d}'.format(n), alpha=0.5)
            else:
                plt.hist(c[n][cn], label='chan{:d}'.format(n), range=(lo/5.,hi*5),
                         alpha=0.5, bins=50)
                yl = plt.ylim()
                plt.plot([lo,lo],yl,':k')
                plt.plot([hi,hi],yl,':k')
            if k==0:
                plt.legend()

            plt.title(cn)

    return mask

def getmaskval(cs, cl):
    """Return True if data should be retained, False if it should be cut"""
    maskval = True
    for cn in list(cs.keys()):
        lo = cl[cn][0]
        hi = cl[cn][1]
        if (cs[cn]<lo) | (cs[cn]>hi):
            maskval = False
            break
    return maskval

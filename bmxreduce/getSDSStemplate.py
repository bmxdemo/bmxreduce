import numpy as np

def getSDSStemplate(c):
    """Takes as input coaddbyday object"""

    rasdss,decsdss,zsdss,zerr,u,g,dum,i,zz = np.loadtxt('personal/SDSS_galaxies_z1.csv',delimiter=',',comments='#',unpack=True)
    fsdss = 1420.406 / (1+zsdss)

    fwhm = 4.0
    ddec = 0.0
    dra = -1.0

    BNLlat = c.dec.mean() + ddec
    w = np.exp(-(decsdss-BNLlat)**2 / (2*(fwhm/2.355)**2))

    # Just saves some time in the histogram step
    dec0 = BNLlat - 5*fwhm/2
    dec1 = BNLlat + 5*fwhm/2
    ind = (decsdss>dec0) & (decsdss<dec1)

    rabe = c.m.item()['rabe']
    df = c.f[1]-c.f[0]
    fbe = np.hstack((c.f-df/2,c.f[-1]+df/2))

    N,xe,ye = np.histogram2d(rasdss[ind] + dra, fsdss[ind], bins=[rabe,fbe],
                             normed=True, weights=w[ind])

    N = N.T

    # Convolve in RA with beam
    xc = (xe[0:-1] + xe[1:])/2
    xc = (xc - xc.mean())*np.cos(BNLlat*np.pi/180)
    sigma = fwhm / 2.355
    beam = np.exp(-xc**2/(2*sigma**2))
    beam = beam/np.sum(beam)
    for k,val in enumerate(N):
        N[k] = np.convolve(N[k],beam,'same')


    # Do a little smoothing in frequency
    kerdf = 0.5 # MHz
    kersz = 3*kerdf / df
    kerx = np.linspace(-1,1,kersz)*df
    ker = np.exp(-kerx/(2*kerdf**2))
    ker = ker/np.sum(ker)
    for k in range(N.shape[1]):
        N[:,k] = np.convolve(N[:,k],ker,'same')


    N = N.T


    return N


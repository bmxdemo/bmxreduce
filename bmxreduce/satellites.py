#
# satellite predictor
#
import numpy as np
import glob
from orbit_predictor.sources import EtcTLESource
from orbit_predictor.locations import Location
from orbit_predictor.predictors import Position
from .telescope import BMXLatLon
from astropy.time import Time

## stolen from Will Tyndal

def llh_to_altaz(loc1, loc2, radian=False):
    '''Return the altaz looking from loc2'''
    loc1_llh = loc1.position_llh
    loc2_llh = loc2.position_llh
    loc1_xyz = loc1.position_ecef
    loc2_xyz = loc2.position_ecef

    if radian:
        coeff = 1
    else:
        coeff = 180 / np.pi

    dx = loc1_llh[1] - loc2_llh[1]
    dy = loc1_llh[0] - loc2_llh[0]
    az = np.arctan(np.float64(dx) / np.float64(dy)) * coeff
    if dy < 0:
        az += np.pi * coeff
    if az > np.pi * coeff:
        az -= 2 * np.pi * coeff

    # Earth ellipsoid parameters
    a = 6378.1370
    b = 6356.752314
    # earth radius
    n1 = np.sqrt(loc1_xyz[0]**2 + loc1_xyz[1]**2 + loc1_xyz[2]**2)
    n2 = np.sqrt(loc2_xyz[0]**2 + loc2_xyz[1]**2 + loc2_xyz[2]**2)
    dist = np.sqrt((loc1_xyz[0] - loc2_xyz[0])**2 + (loc1_xyz[1] - loc2_xyz[1])**2 + (loc1_xyz[2] - loc2_xyz[2])**2)\
    # cosA = (b^2 + c^2 - a^2) / 2bc
    cosalt_center = (n2**2 + dist**2 - n1**2) / (2 * n2 * dist)
    #print(cosalt, n1+loc1_llh[2], n2+loc2_llh[2], dist)
    alt_center = np.pi - np.arccos(cosalt_center)

    lat2_center = np.arctan(loc2_xyz[2] / np.sqrt(loc2_xyz[0]**2 + loc2_xyz[1]**2))
    lat2_corr = loc2_llh[0] / 180 * np.pi - lat2_center

    alt = (0.5*np.pi - (alt_center + lat2_corr)) * coeff
    #print(alt, lat2_corr)

    return alt, az




class Satellites:
    def __init__ (self,mjds, logfn):
        """ Pass list of mjds and we predict satellites """
        self.log=logfn
        self.mjds=mjds
        self.almanac='/astro/u/bmx/bmxreduce/data/almanac/TLE/'
        lat,lon=BMXLatLon()
        self.loc = Location("BNL", latitude_deg=lat, longitude_deg=lon, elevation_m=79)
        ## find central mjd and the closest TLE file
        mean_mjd=mjds.mean()
        self.tledate=self.find_tle_date(mean_mjd)
        ## now we need to convert mjds to datetime

    def find_tle_date(self,target_mjd):
        datelist=[]
        for fname in glob.glob(self.almanac+'/*/*.tle'):
            datelist.append(fname[-10:-4])
        datelist=sorted(set(datelist))
        bestdif=1e10
        self.log("Found %i TLE files"%len(datelist))
        for date in datelist:
            timestr='20%s-%s-%sT07:30:00'%(date[0:2],date[2:4],date[4:6]) ##interpolate between daylight and not
            mjd=Time(timestr, format='isot', scale='utc').mjd
            dt=mjd-target_mjd
            if abs(dt)<bestdif:
                bestdif=abs(dt)
                bestdate=date
        self.log("Found closest TLE date: %s"%bestdate)
        if (bestdif>10):
            self.log("Warning, best TLE more than 10 days away")

        return bestdate

        
    def get_predictions(self):
        outlist=[]

        dtimes=Time(self.mjds,format='mjd').to_datetime()
        ## let's make another set that is every 5 mins. If in range, we'll recalculate in full
        dtimes_sparse=Time(np.arange(self.mjds[0],self.mjds[-1],5/(60*24)) ,format='mjd').to_datetime()

        for satype in 'GPS,GAL,GLO,BEI'.split(','):
            ## first load the sources
            fn=self.almanac+'20%s/%s%s.tle'%(self.tledate[:2],satype,self.tledate)
            self.log("Loading %s ..."%fn)
            satlist=[x[:-1] for x in open(fn).readlines()[::3]]
            self.log("Found %i satelites. "%(len(satlist)))
            ## this requires multi-tle-support branch of https://github.com/jamlamberti/orbit-predictor
            ## hopefully to be fixed soon
            sources=EtcTLESource(fn)
            for sa in satlist:
                sanam=sa.strip().replace(' ','_') ## clean name
                pred=sources.get_predictor(sa)
                ## now something like this
                alt,az=np.array([llh_to_altaz(pred.get_position(when_utc=utc),self.loc)
                                 for utc in dtimes_sparse]).T
                if (np.any(alt>70)):
                    alt,az=np.array([llh_to_altaz(pred.get_position(when_utc=utc),self.loc)
                                 for utc in dtimes]).T
                    altmax=alt.max()
                    self.log ("Found transiting satellite: ",sanam, 'max alt (%2.0f deg)'%altmax)
                    outlist.append((altmax,sanam,alt,az))
        return outlist
                



    


        

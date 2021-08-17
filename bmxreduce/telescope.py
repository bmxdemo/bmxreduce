#
#
# Utility function for transformations, etc.
import os, pickle
import numpy as np
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
import astropy.units as u


def BMXLatLon():
    # latitude=40.878180, longitude_deg=-72.856640
    #lat=40.87792
    #lon=-72.85852
    #   fixed numbers based on POC email from Tue Sep 1st, 2020 
    lat = 40.869951
    lon =-72.866072
    return lat,lon



def location():
    lat,lon=BMXLatLon()
    return EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=0*u.m)


def time2coord(time):
    point=AltAz(alt=90*u.deg, az=0*u.deg, location=location(), obstime=time)
    sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
    gal = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='galactic'))
    ra = sky.ra.rad 
    dec = sky.dec.rad 
    gall = gal.l.rad
    galb = gal.b.rad
    return ra,dec,gall,galb

def times2coords(times,usecache=True):
    cafn="cache/times2coords.pickle"
    cad={}
    ncount=0

    if usecache and os.path.isfile(cafn):
        cad=pickle.load(open(cafn,'rb'))
    for t in times:
        ts=str(t)
        if not ts in cad:
            cad[ts]=time2coord(t)
            ncount+=1
    if usecache and ncount>0:
        print ("Calculated %i new coord conversions for a total of %i."%(ncount,len(cad)))
        pickle.dump(cad,open(cafn,'wb'))
    res=np.array([cad[str(t)] for t in times])
    ra,dec,gall,galb=res.T
    return ra,dec,gall, galb
    
        
def mjd2coords(mjd):
    """Get RA/Dec coordinates"""
    times = Time(mjd, format='mjd')
    return times2coords(times)




    

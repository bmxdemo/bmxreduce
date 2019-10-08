#
#
# Utility function for transformations, etc.
#
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
import astropy.units as u


def location():
    return EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)


def Time2coords(times):
    points = [AltAz(alt=90*u.deg, az=0*u.deg, location=location(), obstime=time) for time in times]
    sky = [point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs')) for point in points]
    gal = [point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='galactic')) for point in points]
    ra = [s.ra.rad for s in sky]
    dec = [s.dec.rad for s in sky]
    gall = [g.l.rad for g in gal]
    galb = [g.b.rad for g in gal]
    return ra,dec,gall, galb
    
        
def mjd2coords(mjd):
    """Get RA/Dec coordinates"""
    times = Time(mjd, format='mjd')
    return Time2radec(times)




    

import numpy as np
from glob import glob
import os
import astropy.units as u
from datetime import datetime, timedelta
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
import time
from datamanager import datamanager

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

class mapmanager(datamanager):
    
    def __init__(self):

        # Set up directories
        self.getdirs()

        return

    def getmapdefn(self, type='SDSSlowz'):
        """Get map type"""

        m = {}

        if type == 'SDSSlowz':
            m['rarange'] = np.array([100.0,260.0])
            m['dra'] = 1.0
            m['frange'] = np.array([1290.0, 1510.0])
            m['df'] = None

        m['rabe'] = np.arange(m['rarange'][0],m['rarange'][1]+m['dra'],m['dra'])
        m['ra'] = (m['rabe'][0:-1] + m['rabe'][1:])/2.
        if m['df'] is not None:
            m['fbe'] = np.arange(m['frange'][0],m['frange'][1]+m['df'],m['df'])
            m['f'] = (m['fbe'][0:-1] + m['fbe'][1:])/2.
        else:
            m['fbe'] = None
            m['f'] = None

        self.m = m

    def getmaptags(self):
        """Return a list of tag lists, where the sub-lists are temporally
        contiguous tags containing data in the map field."""
        # Get tags
        self.gettags(reduced=True)
        
        # Get start/stop times
        self.start = []
        self.stop = []
        for k,tag in enumerate(self.tags):
            # Get start/stop datetime
            self.start.append(self.tag2datetime(self.tags[k]))
            # Assume it's an hour
            self.stop.append(self.start[k] + timedelta(hours=1))

        # Get ra/dec of start/stop times
        self.rastart,decstart = self.dt2radec(self.start)
        self.rastop,decstop  = self.dt2radec(self.stop)

        # Return tags with data w/in range
        ind = (self.rastart <= self.m['rabe'][-1]) & (self.rastop >= self.m['rabe'][0])
        self.tags = self.tags[ind]
        self.start = np.array(self.start)[ind]
        self.stop = np.array(self.stop)[ind]
        self.rastart = self.rastart[ind]
        self.rastop = self.rastop[ind]

        # Now split into temporally contiguous groups
        xinner = []
        xouter = []
        tinner = []
        for k,tag in enumerate(self.tags):

            if k==0:
                xinner.append(tag)
                continue
            
            # This logic will break if the map field spans RA=0. Defer fixing
            # until then.
            if self.rastart[k] >= self.rastop[k-1]:
                xinner.append(tag)
            else:
                xouter.append(xinner)
                xinner = [tag]
    
        self.tagnest = xouter

    def dt2radec(self, dt):
        """Datetime to RA/Dec"""
        time = Time(dt, format='datetime', scale='utc')
        point = AltAz(alt=90*u.deg, az=0*u.deg, location=telescope_loc, obstime=time)
        sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
        return sky.ra.value, sky.dec.value


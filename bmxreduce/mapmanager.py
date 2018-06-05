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

        if type == 'galcross':
            m['rarange'] = np.array([40,180])
            m['dra'] = 1.0
            m['frange'] = np.array([1290.0, 1510.0])
            m['df'] = None

        if type == 'SDSS':
            m['rarange'] = np.array([100.0,260.0])
            m['dra'] = 1.0
            m['frange'] = np.array([1100.0, 1510.0])
            m['df'] = None

        if type == 'SDSSlowz':
            m['rarange'] = np.array([100.0,260.0])
            m['dra'] = 1.0
            m['frange'] = np.array([1290.0, 1510.0])
            m['df'] = None

        if type == 'SDSSlowzlowra':
            m['rarange'] = np.array([100.0,180])
            m['dra'] = 1.0
            m['frange'] = np.array([1290.0, 1510.0])
            m['df'] = None

        if type == 'SDSSlowzra120':
            m['rarange'] = np.array([120.0,130])
            m['dra'] = 1.0
            m['frange'] = np.array([1290.0, 1510.0])
            m['df'] = None

        if type == 'GPS':
            m['rarange'] = np.array([55.0,65.0])
            m['dra'] = 1.0
            m['frange'] = np.array([1110.0, 1510.0])
            m['df'] = None


        m['rabe'] = np.arange(m['rarange'][0],m['rarange'][1]+m['dra'],m['dra'])
        m['ra'] = (m['rabe'][0:-1] + m['rabe'][1:])/2.
        if m['df'] is not None:
            m['fbe'] = np.arange(m['frange'][0],m['frange'][1]+m['df'],m['df'])
            m['f'] = (m['fbe'][0:-1] + m['fbe'][1:])/2.
        else:
            m['fbe'] = None
            m['f'] = None

        m['mapdefn'] = type

        self.m = m

    def getmaptags(self, hasmap='all', sn='real', canonical=False):
        """Return a list of tag lists, where the sub-lists are temporally
        contiguous tags containing data in the map field. If canonical=True,
        only return tags that fall within the canonical sim date range.

        hasmap = 'all', True, False
        sn is map serial number, relevant only if hasmap != 'all'
        """
        # Get tags
        self.gettags(reduced=True, applycuts=True)
        
        # Trim down to canonical
        if canonical:
            year,month,day,hr,min = self.parsetags(self.tags)
            ind = np.where((year==18) & (month==4) & (day>=1) & (day<=2))[0]
            self.tags = self.tags[ind]

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
    
        xouter.append(xinner)
        self.tagnest = xouter
        self.tagnest = np.array(self.tagnest)

        # Only keep longest tag list if canonical
        if canonical:
            taglen = np.array([len(x) for x in self.tagnest])
            ind = np.where(taglen==np.max(taglen))[0]
            self.tagnest = self.tagnest[ind]

        # Get rid of tags with no map
        if hasmap != 'all':
            keepind = []
            for k,tags in enumerate(self.tagnest):
                fn = self.getmapfname(sn, tags)
                keepind.append(os.path.isfile(fn))
            keepind = np.array(keepind)

            if hasmap:
                self.tagnest = self.tagnest[keepind]
            else:
                self.tagnest = self.tagnest[~keepind]

    def dt2radec(self, dt):
        """Datetime to RA/Dec"""
        time = Time(dt, format='datetime', scale='utc')
        point = AltAz(alt=90*u.deg, az=0*u.deg, location=telescope_loc, obstime=time)
        sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
        return sky.ra.value, sky.dec.value

    def getmapfname(self, sn, tags):
        """Get maap filename"""
        fdir = 'maps/bmx/{:s}/{:s}/'.format(sn, self.m['mapdefn'])
        if not os.path.isdir(fdir):
            os.makedirs(fdir)
        fn = '{:s}/{:s}_map.npz'.format(fdir, tags[0][0:6])
        return fn

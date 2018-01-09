import numpy as np
from glob import glob
import os
from datetime import datetime, timedelta

class datamanager(object):
    
    def __init__(self, dataroot=None):
        # Set up dirs
        self.getdirs(dataroot=dataroot)
        return

    def getdirs(self, dataroot=None):
        """Set up directories. Put here instead of __init__ so other classes can
        inherit this"""
        if dataroot is None:
            if os.environ.has_key('BMXDATA'):
                dataroot=os.environ['BMXDATA']
            else:
                dataroot='data/raw'  ## make sure data softlinks to some actual position
        self.dataroot=dataroot
        if os.environ.has_key('BMXREDUCED'):
            self.reducedroot=os.environ['BMXREDUCED']
        else:
            self.reducedroot='data/reduced'
        return

    
    def gettags(self, new=False, reduced=False):
        """Get a list of all data tags with spectrometer files and in "good"
        time periods
        
        if new = True, only return data tags that have no reduced data
        associated with them.

        if reduced = True, only return tags with reduced data.

        Returns list of tag strings"""

        # Get list of directories with raw data, which let's say are directories
        # matching pattern 20??
        dirs = glob(self.dataroot+'/[0-9][0-9][0-9][0-9]')

        # Initialize tag list
        tags = []

        # Loop over dirs
        for k,val in enumerate(dirs):
            fn = glob(os.path.join(val,'*.data'))
            fn = [self.fname2tag(x) for x in fn]
            for j in fn:
                tags.append(j)

        # Load times to cut, convert to datetime objects
        x = np.loadtxt('auxdata/cuttimes.csv', delimiter=',', comments='#', dtype='string')
        if len(x)>0:

            sc = [datetime.strptime(k.strip(), '%Y-%m-%d:%H:%M:%S') for k in x[:,0] ]
            ec = [datetime.strptime(k.strip(), '%Y-%m-%d:%H:%M:%S') for k in x[:,1] ]
    
            # Convert tags to datetime objects. Right way to do this would be to
            # load data and get actual start/top times down to the hour, optionally
            # saving a separate csv file to avoid repeating. We will instead be
            # cheesy and get start times from the filename, assume they occur
            # exactly on the minute, and assume all data files last one hour.
            st = [self.tag2datetime(k) for k in tags]
            et = [k + timedelta(hours=1) for k in st]

            # Turn into numpy arrays
            sc = np.array(sc)
            ec = np.array(ec)
            st = np.array(st)
            et = np.array(et)
            tags = np.array(tags)

            # Now go through and get rid of tags that fall within the cut times
            for k in range(sc.size):
                ind1 = (st>=sc[k]) & (st<ec[k])
                ind2 = (et>=sc[k]) & (et<ec[k])
                ind = ind1 | ind2
                tags = tags[~ind]
                et = et[~ind]
                st = st[~ind]

        if new:
            # Only return tags lacking a reduced file
            ind = np.array([os.path.isfile(self.getreducedfname(k)) for k in tags])
            tags = tags[~ind]

        if reduced:
            # Only return tags with a reduced file
            ind = np.array([os.path.isfile(self.getreducedfname(k)) for k in tags])
            tags = tags[ind]
            
        self.tags = np.sort(tags)
        self.fnames = [self.getrawfname(k) for k in self.tags]


    def tag2datetime(self, tag):
        """Turn tag string into datetime object"""
        return datetime.strptime(tag.strip(), '%y%m%d_%H%M')

    def fname2tag(self, fn):
        """Turn filename into tag string. Assumes path is
        "dir/totallyarbitrary/YYMMDD_HHMM+whatever.anything" and returns the
        "YYMMDD_HHMM" portion"""
        return os.path.split(fn)[1][0:11]

    def parsetag(self, tag):
        """Parse a tag, returns year, month, day, hr, min
        e.g. yr,mn,d,h,min = parsetag('170901_1000')
        """
        yr = tag[0:2]
        month = tag[2:4]
        day = tag[4:6]
        hr = tag[7:9]
        minute = tag[9:11]

        return yr, month, day, hr, minute


    def getrawfname(self, tag):
        """Get filename from tag"""
        taginfo = self.parsetag(tag)
        return os.path.join(self.dataroot,taginfo[0]+taginfo[1], tag+'.data')


    def getreducedfname(self, tag):
        """Get filename from tag"""
        taginfo = self.parsetag(tag)
        return os.path.join(self.reducedroot,taginfo[0]+taginfo[1], tag+'_reduced.data.npz')

    def loadcsvbydate(self, fname, tag):
        """Get dated csv file closest in time and before tag matching pattern:
        auxdata/fname_YYYYMMDD.csv"""
        fn = np.array(glob('auxdata/{:s}_*.csv'.format(fname)))
        fndate = np.array([np.float(x[-10:-4]) for x in fn])
        tagdate = np.float(tag[0:6])
        ind = np.where(fndate < tagdate)[0]

        fn = fn[ind]
        fndate = fndate[ind]
        
        ind = np.where(fndate == np.max(fndate))[0][0]
        
        print('loading {:s}'.format(fn[ind]))
        x = np.loadtxt(fn[ind], delimiter=',', comments='#')

        return x


    def getradec(self):
        """Get RA/Dec coordinates"""
        time = Time(self.d.data['mjd'], format='mjd')
        point = AltAz(alt=90*u.deg, az=0*u.deg, location=telescope_loc, obstime=time)
        sky = point.transform_to(SkyCoord(0*u.deg, 0*u.deg, frame='icrs'))
        self.ra = sky.ra
        self.dec = sky.dec

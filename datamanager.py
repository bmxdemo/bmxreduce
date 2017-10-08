import numpy as np
from glob import glob
import os
from datetime import datetime, timedelta

class datamanager(object):
    
    def __init__(self, tag=None, dataroot=None):
        # Do nothing
        if dataroot is None:
            if os.environ.has_key('BMXDATA'):
                dataroot=os.env['BMXDATA']
            else:
                dataroot='/gpfs01/astro/workarea/bmxdata'
        self.dataroot=dataroot
        return


    def gettags(self, new=False):
        """Get a list of all data tags with spectrometer files and in "good"
        time periods
        
        if new = True, only return data tags that have no reduced data
        associated with them.

        Returns list of tag strings"""

        # Get list of directories with raw data, which let's say are directories
        # matching pattern 20??
        dirs = glob(self.dataroot+'/raw/20[0-9][0-9]')

        # Initialize tag list
        tags = []

        # Loop over dirs
        for k,val in enumerate(dirs):
            fn = glob(os.path.join(val,'*.data'))
            fn = [self.fn2tag(x) for x in fn]
            for j in fn:
                tags.append(j)
        
        # Load times to cut, convert to datetime objects
        x = np.loadtxt('auxdata/cuttimes.csv', delimiter=',', comments='#', dtype='string')
        sc = [datetime.strptime(k.strip(), '%Y-%m-%d:%H:%M:%S') for k in x[:,0] ]
        ec = [datetime.strptime(k.strip(), '%Y-%m-%d:%H:%M:%S') for k in x[:,1] ]
        
        # Convert tags to datetime objects. Right way to do this would be to
        # load data and get actual start/top times down to the hour, optionally
        # saving a separate csv file to avoid repeating. We will instead be
        # cheesy and get start times from the filename, assume they occur
        # exactly on the minute, and assume all data files last one hour.
        st = [self.tag2datetime(k) for k in tags]
        et = [k + timedelta(hours=1) for k in st]

        # Now go through and get rid of tags that fall within the cut times
        sc = np.array(sc)
        ec = np.array(ec)
        st = np.array(st)
        et = np.array(et)
        tags = np.array(tags)

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

        self.tags = np.sort(tags)
        self.fnames = [self.getrawfname(k) for k in self.tags]


    def tag2datetime(self, tag):
        """Turn tag string into datetime object"""
        return datetime.strptime(tag.strip(), '%y%m%d_%H%M')

    def fn2tag(self, fn):
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
        return os.path.join(self.dataroot,'raw','20'+taginfo[0], tag+'.data')


    def getreducedfname(self, tag):
        """Get filename from tag"""
        taginfo = self.parsetag(tag)
        return os.path.join(self.dataroot,'reduced','20'+taginfo[0], tag+'_reduced.data')

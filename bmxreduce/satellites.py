#
# Based on code by Fergus Noble <fergus@swift-nav.com>
# Copyright (C) 2011-2014 Swift Navigation Inc.
#
#

from __future__ import print_function, division
import numpy as np
from numpy import sin,cos
import os, time
from datetime import datetime, timedelta
from bmx_const import BMX_LAT, BMX_LONG


NAV_GM = 3.986005e14
NAV_OMEGAE_DOT = 7.2921151467e-005
GPS_L1_HZ = 1.57542e9
NAV_C = 299792458.0

RANDOM_SUNDAY = 1283644784



def longlat2ECEF(dlat,dlng):
    """ Convert lng lat in degreess to ECEF coordinate system """
    lat,lng=dlat/180*np.pi, dlng/180*np.pi
    #Following : http://mathforum.org/library/drmath/view/51832.html
    R=6378137
    f= 1./298.257224
    C= 1/np.sqrt(cos(lat)**2+(1-f)**2*sin(lat)**2)
    S=(1-f)**2*C
    location=R*C*cos(lat)*cos(lng), R*C*cos(lat)*sin(lng), R*S*sin(lat)
    return location

## various time conversions utils
def mjd2gps(mjd_):
    mjd=mjd_+16./(24.*60.*60.)
    dayofweek=(mjd+3.0)%7
    timeofweek = dayofweek*24*60*60
    gpsweek=int((mjd+3.0)/7)-6321
    return gpsweek,timeofweek
                  
def mjd2yearday(mjd):
    ti=mjd2time(mjd)
    ts=time.gmtime(ti)
    return ts.tm_year, ts.tm_yday
    
def time2mjd(time):
    return time/86400.0  + 40587.0

def mjd2time(mjd):
    return (mjd-40587.0)*86400.0


def time_of_week(time):
    return (time - RANDOM_SUNDAY) % (7 * 24 * 60 * 60)


# main class

class sat_ephem(object):
    
    def __init__(self, location=None):
        """ if location is None assume BMX position """
        
        # Set up dirs, improve this eventually
        self.almanacroot = 'data/almanac'
        self.gpsalmanac = self.almanacroot+'/gps'
        self.location = longlat2ECEF (BMX_LAT, BMX_LONG)
        self.nlocation = self.location / np.sqrt(np.dot(self.location, self.location))
        #lat,lng=BMX_LAT/180*np.pi, BMX_LONG/180*np.pi
        #self.nlocation = (cos(lat)*cos(lng), cos(lat)*sin(lng), sin(lat))
        self.phivec = np.cross ([0,0,1],self.nlocation)
        self.phivec /= np.sqrt(np.dot(self.phivec,self.phivec))
        self.thetavec = np.cross (self.nlocation,self.phivec)
        self.thetavec /= np.sqrt(np.dot(self.thetavec,self.thetavec))
        

        self.debug=False
        
    def get_tracks (self, mjd, max_zen_distance=20*np.pi/180):
        """ For a given list of MJDs produce "tracks" of satellites transiting bmx at less than max_zen_distance
            All units in radians.
            It assumed MJDs are monotonically increasing.
            Returns list of tuples. Each tupple is
            (satid, indices, ddec, dra)
        """
        ## first chop things into individual MJDs
        imjd=mjd.astype("int")
        mjdl=sorted(set(imjd))
        ## now get starting indices
        edgendx=[next(i for i in xrange(len(imjd)) if imjd[i] == m) for m in mjdl]
        edgendx.append(len(mjd))
        MXSat=100
        ndx=[np.zeros(0,int) for x in range(MXSat)] ## max 100 satellites
        dlatl=[np.zeros(0,float) for x in range(MXSat)]
        dlngl=[np.zeros(0,float) for x in range(MXSat)]
        
        for ci,imjd in enumerate(mjdl):
            mjds=mjd[edgendx[ci]:edgendx[ci+1]]
            ## Load almanac
            self.load_almanac(imjd)
            for si,s in enumerate(self.sats):
                sat_pos=s.ECEF(mjds)
                ## correct for our location
                # if s.prn==26:
                #     print (sat_pos, self.location)


                sat_pos -= self.location

                ### vectors on the unit sphere
                sat_posr = np.sqrt((sat_pos**2).sum(axis=1))
                ### differences with normal vectorr
                dsphere_pos = np.array(sat_pos)
                dsphere_pos[:,0] /= sat_posr
                dsphere_pos[:,1] /= sat_posr
                dsphere_pos[:,2] /= sat_posr
                ## find vectors that are aligned enough
                w=np.where(np.dot(dsphere_pos,self.nlocation)>np.cos(max_zen_distance))[0]
                # if s.prn==26:
                #     print (sat_pos,'Y')
                #     print (sat_posr,'Y')
                #     print (dsphere_pos,'Y')
                #     #print (dsphere_pos**2).sum(axis=1)
                #     print (np.arccos(np.dot(dsphere_pos,self.nlocation)),'XX')
                if (len(w>0)):
                    dsphere_pos = dsphere_pos[w]
                    dsphere_pos-= self.nlocation
                    ## project onto phi and theta unit vectors, correct for projection angle
                    dlat = np.arcsin(np.dot(dsphere_pos, self.thetavec))
                    dlng = np.arcsin(np.dot(dsphere_pos, self.phivec))
                    ndx[si]=np.concatenate((ndx[si],w))
                    dlatl[si]=np.concatenate((dlatl[si],dlat))
                    dlngl[si]=np.concatenate((dlngl[si],dlng))
                # if s.prn==26:
                #     print (dlat, dlng,np.sqrt(dlat**2+dlng**2))
                    #stop()
        ## now split into individual tracks
        res=[]
        for si in range(MXSat):
            if len(ndx[si])>0:
                nd=ndx[si]
                dlat=dlatl[si]
                dlng=dlngl[si]
                li=0
                while li<len(nd):
                    hi=next((i for i in xrange(li+1,len(nd)) if nd[i] != nd[i-1]+1),len(nd) )
                    res.append((self.sats[si].prn, nd[li:hi],dlat[li:hi], dlng[li:hi]))
                    li=hi
        return res

    def load_almanac(self,mjd):
        """ Loads almanac for particular MJD. Sateliates are reloaded for this 
            particular mjd """
        
        yr,dy=mjd2yearday(mjd)
        fn=self.gpsalmanac+"/almanac_%4d_%03d.alm"%(yr,dy)
        print (fn)
        blocks = []
        yuma=open(fn).readlines()
        for (i, line) in enumerate(yuma):
            if line[:3] == "ID:":
                blocks += [yuma[i:i + 13]]
        self.sats = map(lambda bl: Sat(bl), blocks)
                
            
class Sat:
    def __init__(self, yuma_block):
        fields = map(lambda x: x[25:], yuma_block)
        self.prn = int(fields[0])
        self.healthy = (int(fields[1]) == 0)
        self.ecc = float(fields[2])
        self.toa = float(fields[3])
        self.inc = float(fields[4])
        self.rora = float(fields[5])
        self.a = float(fields[6])**2
        self.raaw = float(fields[7])
        self.argp = float(fields[8])
        self.ma = float(fields[9])
        self.af0 = float(fields[10])
        self.af1 = float(fields[11])
        self.week = int(fields[12])

    def ECEF(self, mjd):
        tow = time_of_week(mjd2time(mjd))
                                   
        tdiff = tow - self.toa
        # if self.prn==26:
        #     print (mjd,tow,tdiff)
        tdiff[tdiff > 302400.0] -= 604800.0
        tdiff[tdiff < -302400.0] += 604800.0

        ma_dot = np.sqrt(NAV_GM / self.a**3)
        ma = self.ma + ma_dot * tdiff

        # Iteratively solve for the Eccentric Anomaly (from Keith Alter and David Johnston)
        ea = ma  # Starting value for E
        ea_old = ea + 1

        while (any(abs(ea - ea_old) > 1.0E-14)):
            ea_old = ea
            tempd1 = 1.0 - self.ecc * np.cos(ea_old)
            ea = ea + (ma - ea_old + self.ecc * np.sin(ea_old)) / tempd1
        ea_dot = ma_dot / tempd1

        # Begin calc for True Analomay and Argument of Latitude
        tempd2 = np.sqrt(1.0 - self.ecc * self.ecc)
        # [rad] Argument of Latitude = True Anomaly + Argument of Perigee
        al = np.arctan2(tempd2 * np.sin(ea), np.cos(ea) - self.ecc) + self.argp
        al_dot = tempd2 * ea_dot / tempd1

        # Calculate corrected radius based on argument of latitude
        r = self.a * tempd1
        r_dot = self.a * self.ecc * np.sin(ea) * ea_dot

        # Calculate inclination based on argument of latitude
        inc = self.inc

        # Calculate position and velocity in orbital plane
        x = r * np.cos(al)
        y = r * np.sin(al)
        x_dot = r_dot * np.cos(al) - y * al_dot
        y_dot = r_dot * np.sin(al) + x * al_dot

        # Corrected longitude of ascending node
        om_dot = self.rora - NAV_OMEGAE_DOT
        om = self.raaw + tdiff * om_dot - NAV_OMEGAE_DOT * self.toa

        sat_pos = np.zeros((len(mjd),3))

        # Compute the satellite's position in Earth-Centered Earth-Fixed coordiates
        sat_pos[:,0] = x * np.cos(om) - y * np.cos(inc) * np.sin(om)
        sat_pos[:,1] = x * np.sin(om) + y * np.cos(inc) * np.cos(om)
        sat_pos[:,2] = y * np.sin(inc)

        return sat_pos

                     
    def packed(self):
        import struct
        return struct.pack("<ddddddddddHBBB", self.ecc, self.toa, self.inc,
                           self.rora, self.a, self.raaw, self.argp, self.ma,
                           self.af0, self.af1, self.week, self.prn,
                           self.healthy, 1)



        
        

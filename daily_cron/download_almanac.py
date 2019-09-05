#!/usr/bin/env python
import sys, datetime, os
n = datetime.datetime.now()
d = '/astro/u/bmx/bmxdata/almanac/TLE/'+str(n.year)
if not os.path.exists(d):
    os.mkdir(d)
assert(os.path.exists(d))

todownload = [ 
('http://celestrak.com/NORAD/elements/gps-ops.txt','GPS'), ## GPS
('http://celestrak.com/NORAD/elements/glo-ops.txt','GLO'), ## GLONASS
('http://celestrak.com/NORAD/elements/galileo.txt','GAL'), ## GALILEO
('http://celestrak.com/NORAD/elements/beidou.txt','BEI')  ##BEIDOU
]
for link,id in todownload:
    filename="%s%02d%02d%02d.tle"%(id,n.year-2000,n.month,n.day)
    print (filename)
    target=os.path.join(d,filename)
    os.system ('wget -q -O %s %s'%(target,link))



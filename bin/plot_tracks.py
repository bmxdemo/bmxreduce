#!/usr/bin/env python
import sys
sys.path+=['.','..','bmxreduce']
import bmxreduce as br
import numpy as np
from numpy import sin,cos
import matplotlib.pyplot as plt

e=br.sat_ephem()
mjdstart=58189
mjds=np.linspace(mjdstart,mjdstart+2.99, 24*60*60)
for prn, ndx, dla, dlng in e.get_tracks(mjds,max_zen_distance=10*np.pi/180):
    dla*=180/np.pi
    dlng*=180/np.pi
    plt.plot(dlng,dla,'-',label=prn)

for r in range(1,11):
    phi=np.linspace(0,2*np.pi,1000)
    if (r%2):
        plt.plot(r*sin(phi), r*cos(phi),"k:")
    else:
        plt.plot(r*sin(phi), r*cos(phi),"k--")

plt.plot(0,0,'ro',markersize=6)
plt.xlim(-12,12)
plt.ylim(-12,12)
plt.axes().set_aspect('equal')
plt.legend()
plt.show()


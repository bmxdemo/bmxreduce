import mapmanager
import numpy as np
import farmit

cpmdtmin = 30.0
cpmdtmax = np.inf
cpmdf = 10.0
cpmalpha = 1e-4

sn = 'real'
mapdef = 'SDSSlowz'

m = mapmanager.mapmanager()
m.getmapdefn(mapdef)
m.getmaptags()

# Turn tags into string
tags = [str(k) for k in m.tagnest]

# Farm
f = farmit.farmit('driver/coaddbygroup.py', 
                  args={'cpmdtmin':[cpmdtmin],
                        'cpmdtmax':[cpmdtmax],
                        'cpmdf':[cpmdf],
                        'cpmalpha':[cpmalpha],
                        'sn':[sn],
                        'mapdef':[mapdef],
                        'tags':tags})

f.writejobfiles()


    

#!/usr/bin/env python
import sys,os
reduce_home=os.path.abspath(os.path.dirname(sys.argv[0])+"/..")
sys.path.append(reduce_home)
import bmxreduce as br
dm = br.datamanager()
dm.setTags(new=True)
dm.setPassages()


#f = br.farmit.farmit('bin/reduce_batch.py', args={'t':dm.tags},
#                  names=['BMX_'+tag for tag in dm.tags],
#                  reqs={'N':2,'X':0,'priority':'low','mode':'bycore1'})
#f.writejobfiles()
#f.runjobs(maxjobs=500)


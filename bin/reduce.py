#!/usr/bin/env python
import sys
sys.path+=['.','..','bmxreduce']
import bmxreduce as br
import farmit

h = br.genhtml()
h.genindex()
h.gentagindex()
h.gentagpages()

dm = br.datamanager()
dm.gettags(new=True)
f = farmit.farmit('bin/reduce_batch.py', args={'t':dm.tags}, reqs={'N':4})
f.writejobfiles()
f.runjobs(maxjobs=500)


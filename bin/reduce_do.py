#!/usr/bin/env python
#
# This is a low level driver, not to be used standalone
#
#
import sys,os
reduce_home=os.path.abspath(os.path.dirname(sys.argv[0])+"/..")
sys.path.append(reduce_home)
import bmxreduce as br

typ,nam,tags=sys.argv[1:4]
r=br.reduce(typ,nam,tags.split(','))
r.go()


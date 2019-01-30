#!/usr/bin/env python
import sys
sys.path+=['.','..','bmxreduce']
import bmxreduce as br
import mapmanager
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", dest="tags", type="str", default='')
parser.add_option("-m", dest="mdef", type="str", default='')
(o, args) = parser.parse_args()

m = mapmanager.mapmanager()
m.getmapdefn(o.mdef)

tags = o.tags.replace("'",'')
tags = tags.split(',')

r = br.reduce_coadd.coaddbygroup(tags, m.m)
r.reduce(dodeglitch=False, dofilter=False, docpm=False)
r.save()


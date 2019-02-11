#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import sys
sys.path+=['.','..','bmxreduce','../bmxdaq/py']
import bmxreduce as br
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", dest="tag", type="str", default='')
(o, args) = parser.parse_args()

r = br.reduce(o.tag)
r.doreduce()


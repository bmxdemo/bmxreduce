#!/usr/bin/env python
#
# This is a low level driver, not to be used standalone
#
#
import sys,os
reduce_home=os.path.abspath(os.path.dirname(sys.argv[0])+"/..")
sys.path.append(reduce_home)
import bmxreduce as br
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-x", dest="typ", type="str", default='pas')
parser.add_option("-t", dest="tags", type="str", default='')
parser.add_option("-n", dest="name", type="str", default='')
parser.add_option("-y", dest="extraname", type="str", default='')
(o, args) = parser.parse_args()



r=br.reduce(o.typ,o.name,o.tags.split(','),extraname=o.extraname)
r.go()


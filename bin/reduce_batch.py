#
# deprecated
#
#!/usr/bin/env python
stop()
import matplotlib as mpl
mpl.use('Agg')
import bmxreduce as br
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", dest="tag", type="str", default='')
(o, args) = parser.parse_args()

r = br.reduce(o.tag)
r.doreduce()


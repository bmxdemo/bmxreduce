import reduce_init
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", dest="tag", type="str", default='')
(o, args) = parser.parse_args()

r=reduce_init.reduce(o.tag)
r.doreduce()


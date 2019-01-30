import argparse
import numpy as np
import reduce_coadd

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--tags", dest="tags", required=True,
                    help="str(taglist) where taglist is list of strings", type=str)
parser.add_argument("--mapdef", dest="mapdef", default='SDSSlowz',
                    help="mapdefinition string", type=str)
parser.add_argument("--sn", dest="sn", default='real',
                    help="serial number string, i.e. 'real' or '00001'", type=str)
parser.add_argument("--fields", dest="fields", default=None,
                    help="'+' separated field list for use with sims, e.g. colore+gsync", type=str)
parser.add_argument("--cpmdf", dest="cpmdf", default=10.0,
                    help="do not use data within cpmdf (MHz) of target datum in CPM filter", type=float)
parser.add_argument("--cpmdtmin", dest="cpmdtmin", default=30.0,
                    help="do not use data within cpmdtmin (minutes) of target datum in CPM filter", type=float)
parser.add_argument("--cpmdtmax", dest="cpmdtmax", default=np.inf,
                    help="do not use data farther than cpmdtmin (minutes) from target datum in CPM filter", type=float)
parser.add_argument("--cpmalpha", dest="cpmalpha", default=1e-4,
                    help="alpha param in CPM ridge regression", type=float)

o = parser.parse_args()

# Tag string to list
tags = o.tags.replace('[','').replace(']','').replace(' ','').replace("'",'').split(',')

r = reduce_coadd.coaddbygroup(tags, o.mapdef, sn=o.sn, fields=o.fields)
r.reduce(cpmdtmin=o.cpmdtmin, cpmdtmax=o.cpmdtmax, cpmdf=o.cpmdf, cpmalpha=o.cpmalpha)
r.save()


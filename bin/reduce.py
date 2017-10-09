#!/usr/bin/env python
import sys
sys.path+=['.','..']
import bmxreduce as br
import argparse
# argument parser
parser = argparse.ArgumentParser(
    description="Reduce BMX data.\n")
parser.add_argument('tags', metavar='tags', type=str, nargs='+',
                    help='tags to reduce, e.g. 170928_1000')
o = parser.parse_args()

# set up reduction objet
br.reduce(o.tags)
# and reduce...
br.doreduce()



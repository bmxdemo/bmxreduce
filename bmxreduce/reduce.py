#
# Main reduction drivers
#

import numpy as np
import bmxdata
from .datamanager import datamanager
import os


class reduce:
    def __init__(self,typ, name, tags, from_stage=0, to_stage=1):
        self.typ=typ
        assert (typ in ['pas','fra'])
        self.name=name
        self.tags=tags
        self.dm=datamanager()
        self.cstage = self.dm.reduce_stage (typ,name)
        self.root= self.dm.reduce_stage_root(typ,name)
        self.from_stage = from_stage
        self.to_stage = to_stage

        if not os.path.isdir(self.root):
            os.mkdir(self.root)
            self.set_stage(0)

    def set_stage(self,i):
        open(self.root+"/stage",'w').write(str(i)+"\n")

    def go(self):
        stage=self.from_stage
        while stage<self.to_stage:
            self.set_stage(stage)
            stage+=1
            if (stage==1):
                self.stage_one()

    def stage_one(self):
        pass
        

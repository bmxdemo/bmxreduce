#
# Main reduction drivers
#

import numpy as np
import bmxdata 
from .datamanager import datamanager
import os
import datetime

class reduce:
    def __init__(self,typ, name, tags, from_stage=0, to_stage=1):
        self.typ=typ
        assert (typ in ['pas','fra'])
        self.name=name
        self.tags=tags
        self.dm=datamanager()
        self.root= self.dm.reduce_root(typ,name)
        self.cstage = self.dm.reduce_stage (typ,name)
        self.load_meta()
        if typ=="pas": self.add_meta("passage")
        if typ=="fra": self.add_meta("fragment")
        self.from_stage = from_stage
        self.to_stage = to_stage
        self.set_logto(None)
        
        if not os.path.isdir(self.root):
            os.mkdir(self.root)
            self.set_stage(0)
            with open(self.root+"/tags",'w') as ft:
                for t in tags:
                    ft.write(t+"\n")
        self.log ("Initialized reduce object. Root = %s "%self.root)

    def set_logto(self,logtag):
        if (self.logfile is not None):
            self.logfile.close()
        self.logtag=logtag
        if logtag is not None:
            self.logfile = open(os.path.join(self.root,self.logtag+".log"),'w')
    
    def log(self,*m):
        if self.logfile is not None:
            print (*m, file=self.logfile)
        if self.logtag is None:
            print (*m)
        else:
            print ("%s:"%self.logtag,*m)
            
    def set_stage(self,i):
        open(self.root+"/stage",'w').write(str(i)+"\n")

    def go(self):
        stage=self.from_stage
        while stage<self.to_stage:
            self.set_stage(stage)
            stage+=1
            if (stage==1):
                ok=self.stage_one():

            self.set_logto(None)
            if not ok:
                print ("Stage failed.")
                return
            

    def load_meta(self):
        tagfn=self.root+"/tags")
        if os.path.isfile(tafn):
            self.meta=set([tag.replace("\n","") for tag in open(tagn).readlines()])
        else:
            self.meta=set())

    def add_meta(self, meta):
        if meta in self.meta:
            return
        self.meta.add(meta)
        with f as open(self.root+"/meta",'w'):
            for m in sorted(self.meta):
                f.write(m+"\n")

    def timenow_string(self):
        return datetime.datetime.now.strftime("%Y-%m-%d %H:%M"))

    def stage_one(self):
        self.logto("stage1")
        self.log("Starting Stage1 reduction ",self.timenow_string)
        for ext in 'D1','D2':
            data=[]
            for fn in [self.dm.getrawfname(tag) for tag self.tags]:
                if ext=='D2':
                    fn.replace('D1','D2')
                self.log("Loading %s..."%fn)
                data.append(bmxdata.BMXFile(tag))
            ncuts=data[0].ncuts
            nchan=data[0].nchan
            assert (nchan==4)
            ## first check for continuity of mjd
            mjd=data[0].data['mjd']
            deltamjd=(mjd[1:]-mjd[:-1]).mean()
            for i in range(len(data)-1):
                if data[i+1].data['mjd'][0]-data[i].data['mjd'][-1]>deltamjd*1.1:
                    self.log("Tags %s %s discontinous in MJD.",self.tags[i],self.tags[i+1])
                    self.log("Quitting!")
                    return False
            mjd=self.vstack([d.data['mjd'] for d in data])
            if ext=='D1':
                mjd_d1=mjd
                ##
                ## Now caclulate ra/dec so that we can get cut into our passage
                ##
                -----------work here
            else:
                ## we now have both, let's combine
                if (mjd_d1[0]-mjd[0])>deltamjd*0.5: ## should be aligned at least by half.
                    self.log("Computers not aligned in MJD.")
                    self.log("Quitting!")
                
                
                
        
            
        
        

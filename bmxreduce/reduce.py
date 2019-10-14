#
# Main reduction drivers
#

import numpy as np
import os
import datetime

import bmxdata 
import fitsio
from .datamanager import datamanager
from . import telescope

from scipy.interpolate import interp1d

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
        self.log ("\n- - - - - - ")
        self.log ("Initialized reduce object. Root = %s "%self.root)

    def set_logto(self,logtag):
        if not hasattr(self,'logfile'):
            self.logfile=None
        
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
                ok=self.stage_one()

            self.set_logto(None)
            if not ok:
                print ("Stage failed.")
                return
            

    def load_meta(self):
        tagfn=self.root+"/tags"
        if os.path.isfile(tagfn):
            self.meta=set([tag.replace("\n","") for tag in open(tagn).readlines()])
        else:
            self.meta=set()

    def add_meta(self, meta):
        if meta in self.meta:
            return
        self.meta.add(meta)
        with open(self.root+"/meta",'w') as f:
            for m in sorted(self.meta):
                f.write(m+"\n")

    def timenow_string(self):
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    def stage_one(self):
        ## debug
        ##D##self.tags=self.tags[:2]+self.tags[-2:]

        self.set_logto("stage1")
        self.log("Starting Stage1 reduction on %s"%self.timenow_string())
        for ext in 'D1','D2':
            data=[]
            for fn in [self.dm.getrawfname(tag) for tag in self.tags]:
                if ext=='D2':
                    fn=fn.replace('D1','D2')
                self.log("Loading %s..."%fn)
                data.append(bmxdata.BMXFile(fn))

            ncuts=data[0].ncuts
            assert (data[0].nChan==2)
            assert (data[0].nCards==2)
            ## first check for continuity of mjd
            mjd=data[0].data['mjd']
            deltamjd=(mjd[1:]-mjd[:-1]).mean()
            for i in range(len(data)-1):
                if data[i+1].data['mjd'][0]-data[i].data['mjd'][-1]>deltamjd*1.1:
                    self.log("Tags %s %s discontinous in MJD."%(self.tags[i],self.tags[i+1]))
                    self.log("Quitting!")
                    return False
            mjd=np.hstack([d.data['mjd'] for d in data])
            if ext=='D1':
                mjd_d1=mjd
                ##
                ## Now caclulate ra/dec so that we can get cut into our passage
                ##
                self.log('Calculating celestical coords...')
                mjdlist=np.arange(mjd[0],mjd[-1]+0.005,0.005)## every 7 mins or so
                ra,dec,gall,galb=telescope.mjd2coords(mjdlist)
                ## set up interpolator
                raint=interp1d(mjdlist,ra)
                decint=interp1d(mjdlist,dec)
                gallint=interp1d(mjdlist,gall)
                galbint=interp1d(mjdlist,galb)
                radeg=raint(mjd)*180/np.pi
                if self.typ=='pas':
                    istart=np.where(radeg>270.0)[0][0]-1 ## we should be safely away from corner cases
                    iend=np.where(radeg<330.0)[0][-1]+1
                else:
                    istart=0
                    iend=len(mjd)
                print ("Cutting %i %i samples at beginning / end..."%(istart,len(mjd)-iend))

            else:
                ## we now have both, let's combine
                if (mjd_d1[0]-mjd[0])>deltamjd*0.2: ## should be aligned at least by 20%
                    self.log("Computers not aligned in MJD.")
                    self.log("Quitting!")
                mjd=0.5*(mjd+mjd_d1)
                mjd=mjd[istart:iend]
                ra=raint(mjd)
                dec=decint(mjd)
                gall=gallint(mjd)
                galb=galbint(mjd)
                self.log("Writing MJD file...")
                fitsio.write(self.root+"/mjd.fits", np.array(mjd,dtype=[('mjd','f4')]))
                self.log("Writing coordinate files")
                fitsio.write(self.root+"/coords.fits",
                             np.rec.fromarrays([ra,dec,gall,galb],
                            dtype=[('ra','f4'),('dec','f4'),('lgal','f4'),('bgal','f4')]))
            
                
                
                
        
            
        
        

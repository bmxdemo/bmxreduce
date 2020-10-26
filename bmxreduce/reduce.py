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
from .satellites import Satellites
from scipy.interpolate import interp1d

class reduce:
    def __init__(self,typ, name, tags, extraname="", from_stage=0, to_stage=1):
        print ("\n- - - - - - \nStarting ... ")
        self.typ=typ
        assert (typ in ['pas','fra'])
        self.name=name
        self.tags=tags
        self.dm=datamanager(extraname=extraname)
        self.root= self.dm.reduce_root(typ,name)
        if not os.path.isdir(self.root):
            os.mkdir(self.root)
            self.set_stage(0)
            with open(self.root+"/tags",'w') as ft:
                for t in tags:
                    ft.write(t+"\n")

        self.cstage = self.dm.reduce_stage (typ,name)
        self.load_meta()
        if typ=="pas": self.add_meta("passage")
        if typ=="fra": self.add_meta("fragment")
        self.from_stage = from_stage
        self.to_stage = to_stage
        self.set_logto(None)
        self.logtag_extra=None
        self.log ("Initialized reduce object. Root = %s "%self.root)

    def set_logto(self,logtag):
        if not hasattr(self,'logfile'):
            self.logfile=None
        
        if (self.logfile is not None):
            self.logfile.close()
            self.logfile=None
        self.logtag=logtag
        if logtag is not None:
            self.logfile = open(os.path.join(self.root,self.logtag+".log"),'w')

    def set_logtag_extra(self,logtag_extra):
        self.logtag_extra=logtag_extra
    
    def log(self,*m):
        if self.logfile is not None:
            if (self.logtag_extra is not None):
                self.logfile.write("(%s) "%self.logtag_extra)
            print (*m, file=self.logfile)
        if self.logtag is None:
            print (*m)
        else:
            if self.logtag_extra is not None:
                print ("%s (%s):"%(self.logtag,self.logtag_extra),*m)
            else:
                print ("%s:"%self.logtag,*m)
            
    def set_stage(self,i):
        open(self.root+"/stage",'w').write(str(i)+"\n")

    def go(self):
        stage=self.from_stage
        self.set_stage(stage)
        while stage<self.to_stage:
            stage+=1
            if (stage==1):
                ok=self.stage_one()

            self.set_logto(None)
            if not ok:
                self.log ("Stage failed.")
                break
            else:
                self.log("Stage complete.")
                self.set_stage(stage)
        self.log("All done.")

    def load_meta(self):
        metafn=self.root+"/meta"
        if os.path.isfile(metafn):
            self.meta=set([meta.replace("\n","") for meta in open(metafn).readlines()])
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
        debug=False
        if debug:
            self.tags=self.tags[:1]+self.tags[-1:]

        self.set_logto("stage1")
        self.log("Starting Stage1 reduction on %s"%self.timenow_string())
        cutdirs={}
        for ext in 'D1','D2':
            data=[]
            for fn in [self.dm.getrawfname(tag) for tag in self.tags]:
                if ext=='D2':
                    fn=fn.replace('D1','D2')
                self.log("Loading %s ..."%fn)
                data.append(bmxdata.BMXFile(fn))

            if (hasattr(data[0],'wires')):
                self.log("Writing wires file ...")
                with open(self.root+"/wires",'wb' if ext=='D1' else 'ab') as f:
                    for i in range(1,5):
                        f.write (data[0].wires[i]+b"\n")

            ncuts=data[0].ncuts
            for cut in range(ncuts):
                if cut not in cutdirs:
                    cdir=self.root+"/cut%i"%cut
                    cutdirs[cut]=cdir
                    if not os.path.isdir(cdir):
                        self.log("Creating %s ..."%cdir)
                        os.mkdir(cdir)
            assert (data[0].nChan==2)
            assert (data[0].nCards==2)
            ## first check for continuity of mjd
            mjd=data[0].data['mjd']
            deltamjd=data[0].deltaT/(24*3600)
            #for i in range(len(data)-1):
            #    if data[i+1].data['mjd'][0]-data[i].data['mjd'][-1]>deltamjd*1.1:
            #        self.log("Tags %s %s discontinous in MJD."%(self.tags[i],self.tags[i+1]))
            #        if not debug: ## during debug we work with discontinous files
            #            self.log("Quitting!")
            #            return False
            mjd=np.hstack([d.data['mjd'] for d in data])
            if  data[0].version>=10:
                temp_fgpa = np.vstack([d.data['temp_fgpa'] for d in data])
                temp_adc = np.vstack([d.data['temp_adc'] for d in data])
                temp_frontend = np.vstack([d.data['temp_frontend'] for d in data])
                have_temps = True
            else:
                have_temps = False
            if ext=='D1':
                mjd_d1=mjd
                if have_temps:
                    temp_fgpa_D1 = temp_fgpa
                    temp_adc_D1 = temp_adc
                    temp_frontend_D1 =temp_frontend
                ##
                ## Now caclulate ra/dec so that we can get cut into our passage
                ##
                self.log('Calculating celestical coords...')
                mjdlist=np.arange(mjd[0]-1e-5,mjd[-1]+0.005,0.005)## every 7 mins or so
                ra,dec,gall,galb=telescope.mjd2coords(mjdlist)
                ## set up interpolator
                ## we now need to make sure ra is always rising
                while (np.any( (ra[1:]-ra[:-1])<0)):
                    i=np.argmax((ra[1:]-ra[:-1])<0)
                    ra[:i+1]-=2*np.pi
                raint=interp1d(mjdlist,ra)
                decint=interp1d(mjdlist,dec)
                gallint=interp1d(mjdlist,gall)
                galbint=interp1d(mjdlist,galb)
                radeg=raint(mjd)*180/np.pi
                if self.typ=='pas':
                    istart=np.where(radeg>270.0-360)[0][0]-1 ## we should be safely away from corner cases
                    iend=np.where(radeg<330.0)[0][-1]+1
                else:
                    istart=0
                    iend=len(mjd)-1 ## cut one sample at the end to make sure D2 is OK
                print ("Cutting %i %i samples at beginning / end..."%(istart,len(mjd)-iend))
                print ("Number of samples: ",iend-istart)
            else:
                ## we now have both, let's combine
                if (mjd_d1[0]-mjd[0])>deltamjd*0.2: ## should be aligned at least by 20%
                    self.log("Computers not aligned in MJD.")
                    self.log("Quitting!")
                mjd_d1=mjd_d1[istart:iend]
                mjd_d2=mjd[istart:iend]
                mjd=0.5*(mjd_d1+mjd_d2)
                ## Let's fix MJD, to be really continous
                mjdfix=np.arange(len(mjd))*deltamjd
                mjdfix+=mjd.mean()-mjdfix.mean()
                #print (deltamjd, mjdfix-mjd)
                # up to 0.5s discrepancy
                if np.any(np.abs(mjdfix-mjd)>max(deltamjd,0.5/(3600*24))):
                    self.log("Warning: Tag perhaps discontinous in MJD.")
                    #if not debug: ## during debug we work with discontinous files
                    #    self.log("Quitting!")
                    #    return False
                mjd=mjdfix
                ra=raint(mjd)
                dec=decint(mjd)
                gall=gallint(mjd)
                galb=galbint(mjd)
                self.log("Writing MJD file...",mjd.shape,mjd_d1.shape,mjd_d2.shape)
                fitsio.write(self.root+"/mjd.fits", np.rec.fromarrays([mjd,mjd_d1,mjd_d2],
                            dtype=[('mjd','f8'),('mjd_d1','f8'),('mjd_d2','f8')]),
                            clobber=True)
                self.log("Writing coordinate file...")
                fitsio.write(self.root+"/coords.fits",
                             np.rec.fromarrays([ra,dec,gall,galb],
                             dtype=[('ra','f4'),('dec','f4'),('lgal','f4'),('bgal','f4')]),
                             clobber=True)
                ## now also write the labjack
                labjack=np.hstack([d.data['lj_diode'] for d in data])[istart:iend]
                if labjack.sum()==0:
                    self.log("No labjack diode...")
                    self.add_meta("no_diode")
                else:
                    outfn=self.root+"/diode.fits"
                    self.log("Writing %s ... "%outfn)
                    fitsio.write(outfn,labjack,clobber=True)
                if have_temps:
                    self.add_meta("have_temperatures")
                    self.log("Writing temperature file")
                    fitsio.write(self.root+"/temperatures.fits",
                                 np.rec.fromarrays(
                                     [np.hstack((temp_fgpa_D1[istart:iend,:],temp_fgpa[istart:iend,:])),
                                      np.hstack((temp_adc_D1[istart:iend,:],temp_adc[istart:iend,:])),
                                      np.hstack((temp_frontend_D1[istart:iend,:],temp_frontend[istart:iend,:]))],
                                     dtype=[('fgpa','4f4'),('adc','4f4'),('frontend','4f4')]),
                             clobber=True)
                    
                ## Now, set up satellites
                self.set_logtag_extra("sats")
                sats=Satellites(mjdfix, logfn=self.log)
                satpred=sats.get_predictions()
                self.set_logtag_extra(None)
                self.log("Writing satellites file...")
                clobber=True
                for altmax,name,alt,az in sorted(satpred,reverse=True):
                    fitsio.write(self.root+"/satellites.fits",
                                 np.rec.fromarrays([alt,az],
                                 dtype=[('alt','f4'),('az','f4')]),
                                 clobber=clobber,header={'SAT_ID':name, 'MAX_ALT:':altmax})
                    clobber=False
                    
            ## Now let's dump our guys
            for cut in range(ncuts):
                nfreq=data[0].nP[cut]
                if (ext=='D1'):
                    outfn=cutdirs[cut]+"/freq.fits"
                    self.log("Writing %s ... "%outfn)
                    fitsio.write(outfn,data[0].freq[cut],clobber=True)
                    if (cut==1) and abs(data[0].freq[cut].mean()-1420.0)<5:
                        self.log("Galactic 21cm cut present ...")
                        self.add_meta("galactic21_cut")
                for ch1 in range(1,5): ## hardcoded, yes, but can become more flexible later
                    for ch2 in range(ch1,5):
                        if (ch1==ch2):
                            name="chan%i_%i"%(ch1,cut)
                            dataR=np.vstack([d.data[name] for d in data])[istart:iend]
                            outfn=cutdirs[cut]+'/auto_%i.fits'%(ch1+4*(ext=='D2'))
                            self.log("Writing %s ... "%outfn)
                            fitsio.write(outfn,dataR,clobber=True)
                        else:
                            nameR='chan%ix%iR_%i'%(ch1,ch2,cut)
                            dataR=np.vstack([d.data[nameR] for d in data])[istart:iend]
                            nameI='chan%ix%iI_%i'%(ch1,ch2,cut)
                            dataI=np.vstack([d.data[nameI] for d in data])[istart:iend]
                            outfn=cutdirs[cut]+'/cross_%i%i.fits'%(ch1+4*(ext=='D2'),ch2+4*(ext=='D2'))
                            self.log("Writing %s ..."%outfn)
                            fitsio.write(outfn,dataR,clobber=True)
                            fitsio.write(outfn,dataI)
        #done
        self.log("Complete.")
        return True

                            
                            
                
                
        
            
        
        

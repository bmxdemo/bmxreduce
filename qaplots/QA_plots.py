#!/cvmfs/astro.sdcc.bnl.gov/SL73/packages/anaconda3/2019.07/bin/python
# coding: utf-8

import numpy as np
from numpy.fft import rfft,irfft, rfftfreq
import scipy
import scipy.interpolate as scint
import astropy
from astropy.io import fits
from astropy.table import Table

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy.fft import rfft2,rfft
import scipy.linalg as la
dT=2**25/1.1e9*32
from matplotlib.colors import LogNorm
import os
import os.path
import sys

import argparse

parser = argparse.ArgumentParser(description='BMX QA plot generator.')
parser.add_argument('-d',"--dataset", dest='dataset', help="Dataset in format YYMMDD_HHMM")
parser.add_argument('-f',"--force", dest='force', action='store_true', help="Overwite existing files")
parser.add_argument('-v',"--verbose", dest='verbose', action='store_true', help="Display messages")

args = parser.parse_args()

QA_out_top_dir = '/gpfs02/astro/www/bmx/qaplots'
reduced_pas_dir = '/gpfs02/astro/workarea/bmxdata/reduced/pas'

#dataset = args.dataset

def getWork():
    ''' The original plan will not work since structure of reduced/pas
        directories changed over time. The earliest that works with this
        code is 191014_2100. So need to select directories after this one and
        ssume that all after will work.
    Get list of "1x..." and "2x..." directories in the dataset directory.
        Get list of "1x..." and "2x..." directories in the qa web directory.
        Remove all of the web directory entries from the dataset list and
        return the shortened list. Lese are th datbasets that need to be
        processed.
        '''

    first_valid = '191014_2100'
    ds_list = list()
    qa_list = list()
    with os.scandir(reduced_pas_dir) as ds_it:
        for ds_entry in ds_it:
            if ds_entry.name < first_valid:
                continue
            if (ds_entry.name.startswith('1') or ds_entry.name.startswith('2')) and ds_entry.is_dir():
                ds_list.append(ds_entry.name)
    with os.scandir(QA_out_top_dir) as qa_it:
        for qa_entry in qa_it:
            if (qa_entry.name.startswith('1') or qa_entry.name.startswith('2')) and qa_entry.is_dir():
                qa_list.append(qa_entry.name)

    if args.verbose:
        print("ds_list len:", len(ds_list))
        print("qa_list len:", len(qa_list))
    
    # remove qa_list from ds_list
    for qa_dir in qa_list:
        if qa_dir in ds_list:
            ds_list.remove(qa_dir)

    return ds_list

    
def makePlots(dataset):
    dir = os.path.join(reduced_pas_dir, dataset)
    QA_out_dir = os.path.join(QA_out_top_dir, dataset)
    if os.path.exists(QA_out_dir):
        if not args.force:
            if args.verbose:
                print("Output exists for", dataset, "aborting.")
            sys.exit(1)
    else:
        os.mkdir(QA_out_dir)
        if args.verbose:
            print("Created", QA_out_dir)
    
    wires=[x.replace('\n','') for x in open(dir+'/wires').readlines()]
    hdu_freq = fits.open(dir+'/cut1/freq.fits')
    freq = hdu_freq[0].data
    hdu_coords = fits.open(dir+'/coords.fits')
    coords = hdu_coords[1].data
    ra=coords['ra']

    if os.path.isfile(dir+'/temperatures.fits'):
        hdu_temp = fits.open(dir+'/temperatures.fits')
        tempe = hdu_temp[1].data
    else:
        tempe = None

    if os.path.isfile(dir+'/diode.fits'):
        hdu_diode = fits.open(dir+'/diode.fits')
        diode = hdu_diode[0].data
    else:
        print ("No diode")
        diode = None

    if diode is not None:
        i=0
        N=len(diode)
        donoff=[]
        while (diode[i]>0):
            i+=1
        while True:
            while (diode[i]==0):
                i+=1
                if i==N: break
            if i==N: break
            st=i
            while (diode[i]>0):
                i+=1
                if i==N: break
            if i==N: break
            donoff.append((st,i))
        donoff=donoff[1:-1] ## skip first and last one

        def process_diode(da,diode):
            di=[]
            for i,j in donoff:
                h=(j-i)
                a=i-h
                b=j+h
                diff=da[i:j].mean(axis=0)-0.5*(da[a:i].mean(axis=0)+da[j:b].mean(axis=0))
                di.append(diff)
                da[i:i+h//2]=da[i-h//2:i]
                da[i+h//2:j]=da[j:j+(j-i-h//2)]
            di=np.array(di)
            afreq=di.mean(axis=0)
            atime=di.mean(axis=1)
            return da,afreq, atime

    def measure_rfi(vec,returnbad=False):
        vecout=np.copy(vec)
        bad=[]
        for i in range(2,len(vec)-2):
            #print (i,(vec[i]+vec[i+1])/(vec[i-1]+vec[i+2]))
            if i in bad:
                continue
            if (freq[i]>1420.3)and(freq[i]<1420.5):
                continue ## actual peak of 21cm, hard to tell
            if (vec[i]/vec[i-1])>1.005:
                ## jump up
            
                if (vec[i]/vec[i+1])>1.005:
                    bad.append(i)
                    vecout[i]=0.5*(vec[i-1]+vec[i+1])
                elif (vec[i]/vec[i+2])>1.005:
                    bad.append(i)
                    bad.append(i+1)
                    vecout[i]=0.5*(vec[i-1]+vec[i+2])
                    vecout[i+1]=0.5*(vec[i-1]+vec[i+2])
        if returnbad:
            return bad,vecout
        rfi=(vec-vecout).sum()/vecout.sum()
        nrfi=len(bad)
        snr=vecout.max()/vecout.min()
        return rfi,nrfi,snr


    rdata=[]
    for i in range(1,9):
        hdu_da = fits.open(dir+'/cut1/auto_%i.fits'%i)
        #da=fitsio.read(dir+'/cut1/auto_%i.fits'%i)
        da = hdu_da[0].data
        if diode is not None:
            da,afreq,atime=process_diode(da,diode)
        else:
            afreq,atime = None, None
        rdata.append((da,afreq,atime))

    for i,(da,_,_) in enumerate(rdata):
        fig = plt.figure(figsize=(40,6))
        plt.subplot(1,3,1)
        da=da[:da.shape[0]//128*128,:]
        da=da.reshape(-1,128,da.shape[1]).mean(axis=1)
        plt.imshow(da.T,aspect='auto',norm=LogNorm(), extent=(ra[0],ra[-1],freq[0],freq[-1]),origin='lower' )
        plt.colorbar()
        plt.ylabel(wires[i])
        plt.subplot(1,3,2)
        v=da.mean(axis=0)
        plt.plot(freq,v)
        b,_=measure_rfi(v,returnbad=True)
        plt.plot(freq[b],v[b],'ro')
        plt.subplot(1,3,3)
        plt.plot(da.mean(axis=1))
        plt.tight_layout()
        # Write files to web space
        pngfile = 'da_%i.png' % (i,)
        pngout = os.path.join(QA_out_dir, pngfile)
        fig.savefig(pngout, format='png')
        plt.close()
        if args.verbose:
            print("Wrote", pngfile)

    if diode is not None:
        fig = plt.figure(figsize=(20,5))
        for i,(_,afreq,_) in enumerate(rdata):
            plt.plot(freq,afreq,label=wires[i])
        plt.legend()

        # Write file to web space
        pngfile = 'afreq.png'
        pngout = os.path.join(QA_out_dir, pngfile)
        fig.savefig(pngout, format='png')
        plt.close()
        if args.verbose:
            print("Wrote", pngout)

        fig = plt.figure(figsize=(20,5))
        for i,(_,_,atime) in enumerate(rdata):
            plt.plot(atime,label=wires[i])
        plt.semilogy()
        plt.ylim(1e9,2e10)
        # Write file to web space
        pngfile = 'atime.png'
        pngout = os.path.join(QA_out_dir, pngfile)
        fig.savefig(pngout, format='png')
        plt.close()
        if args.verbose:
            print("Wrote", pngout)
    else:
        if tempe is not None:
            # box average tempe
            def boxavg(data):
                avg = data[:len(data)//100*100]
                avg = avg.reshape((-1,10)).mean(axis=1)
                return avg
            davg= boxavg(ra)
            fig = plt.figure(figsize=(20,5))
            plt.plot(davg,boxavg(tempe['fgpa'][:,0]),label='FPGA1')
            plt.plot(davg,boxavg(tempe['fgpa'][:,1]),label='FPGA2')
            plt.plot(davg,boxavg(tempe['adc'][:,0]),label='ADC1')
            plt.plot(davg,boxavg(tempe['adc'][:,1]),label='ADC2')
            plt.plot(davg,boxavg(tempe['frontend'][:,0]),label='FRONT1')
            plt.plot(davg,boxavg(tempe['frontend'][:,1]),label='FRONT2')
            plt.legend()
            plt.xlabel('ra')
            plt.ylabel('T [K]')
            pngfile = 'afreq.png'
            pngout = os.path.join(QA_out_dir, pngfile)
            fig.savefig(pngout, format='png')
            fig = plt.figure(figsize=(20,5))
            plt.plot(davg,boxavg(tempe['fgpa'][:,2]),label='FPGA3')
            plt.plot(davg,boxavg(tempe['fgpa'][:,3]),label='FPGA4')
            plt.plot(davg,boxavg(tempe['adc'][:,2]),label='ADC3')
            plt.plot(davg,boxavg(tempe['adc'][:,3]),label='ADC4')
            plt.plot(davg,boxavg(tempe['frontend'][:,2]),label='FRONT3')
            plt.plot(davg,boxavg(tempe['frontend'][:,3]),label='FRONT4')
            plt.legend()
            plt.xlabel('ra')
            plt.ylabel('T [K]')
            pngfile = 'atime.png'
            pngout = os.path.join(QA_out_dir, pngfile)
            fig.savefig(pngout, format='png')

        else:
            fig = plt.figure(figsize=(20,5))
            plt.text(0.3,0.5, "No diode data.",fontsize=20)
            pngfile = 'atime.png'
            pngout = os.path.join(QA_out_dir, pngfile)
            fig.savefig(pngout, format='png')
            pngfile = 'afreq.png'
            pngout = os.path.join(QA_out_dir, pngfile)
            fig.savefig(pngout, format='png')
        


    # Table data
    tableFile = 'table.html'
    tableFull = os.path.join(QA_out_dir, tableFile)
    tableout = open(tableFull, 'w')

    for i,(da,_,_) in enumerate(rdata):
        v=da.mean(axis=0)
        rfi,nrfi,snr=measure_rfi(v)
        #print("Receiver %s : Ampl: %f, RFI badness: %f, RFI #ch:%i,  MW SNR: %f"%(wires[i],v.min()/1e11,rfi*1e3,nrfi,snr))
        # write table data
        tableout.write("    <tr>\n")
        outStr = '      <td class="py-0">%s</td>\n' % (wires[i],)
        tableout.write(outStr)
        outStr = '      <td class="py-0">%f</td>\n' % (v.min()/1e11,)
        tableout.write(outStr)
        outStr = '      <td class="py-0">%f</td>\n' % (rfi*1e3,)
        tableout.write(outStr)
        outStr = '      <td class="py-0">%i</td>\n' % (nrfi,)
        tableout.write(outStr)
        outStr = '      <td class="py-0">%f</td>\n' % (snr,)
        tableout.write(outStr)
        tableout.write("    </tr>\n")

    # Close table
    tableout.close()

    def gettau (data,df):
        N=len(data)
        adata = np.hstack((data*np.hamming(N),np.zeros(8*N)))
        ft = np.abs(np.fft.fft(adata)**2)
        tau = np.fft.fftfreq(len(adata),df)[np.argmax(ft)]*1000
        return tau

    ##  tau data
    cygi1 = np.where(ra>-1.03)[0][0]
    cygi2 = np.where(ra>-1.03+2*np.pi)[0][0]
    print (cygi1,cygi2)
    df = freq[1]-freq[0]
    print (df)
    hdu_da = fits.open(dir+'/cut1/cross_24.fits')
    data = hdu_da[0].data+1j*hdu_da[1].data
    tau12_1 = gettau(data[cygi1, :],df)
    tau12_2 = gettau(data[cygi2, :],df)
    hdu_da = fits.open(dir+'/cut1/cross_68.fits')
    data = hdu_da[0].data+1j*hdu_da[1].data
    tau34_1 = gettau(data[cygi1, :],df)
    tau34_2 = gettau(data[cygi2, :],df)

    tauFile = 'tau.html'
    tauFull = os.path.join(QA_out_dir, tauFile)
    tauout = open(tauFull, 'w')
    if abs(tau12_1)>100 or abs(tau12_2)>100:
        tauout.write ('<font color="red"> <b>Digitizer 1/2 offset: </b> %3.2f / %3.2f </font> <br>'
                      %(tau12_1,tau12_2))
    else:
        tauout.write ("<b>Digitizer 1/2 offset: </b> %3.2f / %3.2f </font> <br>"%(tau12_1,tau12_2))
    if abs(tau34_1)>100 or abs(tau34_2)>100:
        tauout.write ('<font color="red"> <b>Digitizer 3/4 offset: </b> %3.2f / %3.2f </font> <br>'
                      %(tau34_1,tau34_2))
    else:
        tauout.write ("<b>Digitizer 3/4 offset: </b> %3.2f / %3.2f </font> <br>"%(tau34_1,tau34_2))
    tauout.close()

    print ("Tau diff:",tau12_1,tau12_2,tau34_1,tau34_2)

    if args.verbose:
        print("Wrote table file", tableFull)

if __name__ == "__main__":
    dss = list()
    if args.dataset:
        dss.append(args.dataset)
    else:
        dss = getWork()

    if args.verbose:
        print("dss len:", len(dss))
    for ds in dss:
        if args.verbose:
            print("Processing", ds)
        makePlots(ds)
        
    sys.exit(0)


from __future__ import print_function, division
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['image.interpolation']='none' # Turn off interpolation for image display
mpl.rcParams['image.aspect']='auto' # Auto aspect ratio by default
import matplotlib.pyplot as plt
import numpy as np
import reduce_init
from astropy.time import Time
import os
import datamanager
import gc

plt.rc('legend',**{'fontsize':10}) # Legend fonts are too big
plt.ioff()

def nanhist(x, **kwargs):
    plt.hist(x[np.isfinite(x)],**kwargs)

class genplots():
    
    def __init__(self, r, fext=''):
        """Takes reduce_init.reduce object as input, generates reduc plots"""

        # Reduced data
        self.r = r
        
        # Choose this (must be float)
        self.dpi = 80.0

        # Get useful info
        self.fmin = self.r.f[0]
        self.fmax = self.r.f[-1]
        self.tmin = 0
        self.tmax = self.r.d.deltaT[0]*self.r.d.nSamples / 60. # minutes

        self.ts = Time(self.r.d.data['mjd'][0],  format='mjd')
        self.te = Time(self.r.d.data['mjd'][-1], format='mjd')

        # File naming
        self.filebase = self.r.tag
        self.plotdir = 'browser/plots'

        # Close all open plots
        plt.close('all')
        gc.collect()

        return

    def plotrawwf(self, dochan=None, fext='', cscale='log'):
        """Plot raw waterfall plot. Takes chan index."""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)
            v = self.r.d.data[chn]

            # Size of data
            sz = np.array(v.shape)

            # Let's make a figure of a standard size but with hugely increased
            # dpi. This will hopefully keep labels looking readable when zoomed
            # out. 
            padx = 0.25
            pady = 0.25
            figsz = (1+np.array([padx,pady]))*sz/self.dpi
            fac = 20.0
            fig = plt.figure(figsize=figsz/fac, dpi=self.dpi*fac)

            if cscale=='log':
                plt.imshow(np.log10(v.T), extent=(self.tmin,self.tmax,self.fmax,self.fmin))
                titlab = 'log10[ADU^2]'
            else:
                plt.imshow(v.T, extent=(self.tmin,self.tmax,self.fmax,self.fmin))
                plt.clim(0,200)
                titlab = 'ADU^2'

            if 'cal' in fext:
                titlab.replace('ADU','K')

            plt.grid('on')
            plt.xlabel('t (minutes)');
            plt.ylabel('freq (MHz)');
            plt.title('raw {:s} adc ({:s}) -- {:s} - {:s} (UTC)'.format(chn,titlab,self.ts.iso,self.te.iso))
            plt.colorbar(pad=0)

            fname       = self.filebase + '_'+chn+'_wfraw'+fext+'.jpg'
            fname_thumb = self.filebase + '_'+chn+'_wfraw'+fext+'_thumbnail.jpg'
            plt.savefig(os.path.join(self.plotdir,fname), dpi=self.dpi*fac, bbox_inches='tight')
            plt.savefig(os.path.join(self.plotdir,fname_thumb), dpi=self.dpi, bbox_inches='tight')
            plt.close(fig)
            gc.collect()

    def plotrawspec(self, dochan=None):
        """Plot raw time average spectrum"""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)
            v = self.r.d.data[chn]

            # Get median spectrum
            soff = np.nanmedian(v[~self.r.calind,:],0)
            son = np.nanmedian(v[self.r.calind,:],0)
            
            # Plot
            fig = plt.figure(figsize=(7,7))
            
            plt.subplot(2,1,1)
            plt.semilogy(self.r.f, soff, label='cal off')
            plt.semilogy(self.r.f, son, label='cal on')
            plt.xlabel('f (MHz)')
            plt.ylabel('ADU^2')
            plt.grid('on')
            plt.title('Median raw spectrum, {:s}, {:s}'.format(self.r.tag,chn))
            plt.legend()

            plt.subplot(2,1,2)
            plt.plot(self.r.f, son-soff, label='on - off')
            plt.xlabel('f (MHz)')
            plt.grid('on')
            yl=plt.ylim()
            plt.ylim((0,yl[1]))
            plt.legend()

            fname = self.filebase + '_'+chn+'_specraw.png'
            plt.savefig(os.path.join(self.plotdir,fname), bbox_inches='tight')
            plt.close(fig)
            gc.collect()


    def plotcalwf(self, dochan=None):
        """Plot calibrated, downsampled waterfall plot."""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)

            for k in range(5):
                if k==0:
                    v = self.r.data[chan].T
                    titlab = 'T (K)'
                    cl = (0,200)
                    fext = 'data'
                if k==1:
                    v = self.r.data_mf[chan].T
                    titlab = 'T (K)'
                    cl = (-1,1)
                    fext = 'data_mf'
                if k==2:
                    v = self.r.data_svd[chan].T
                    titlab = 'T (K)'
                    cl = (-1,1)
                    fext = 'data_svd'
                if k==3:
                    v = np.log10(self.r.g[chan].T)
                    titlab = 'log10(gain) (ADU^2/K)'
                    cl = None
                    fext = 'gain'
                if k==4:
                    v = (self.r.nhits[chan]/self.r.var[chan]).T
                    titlab = 'weight = nhits/variance (1/K^2)'
                    cl = None
                    fext = 'weight'

                if k==0:
                    # Size of data
                    sz = np.array(v.shape)
                    # Let's make a figure of a standard size but with hugely increased
                    # dpi. This will hopefully keep labels looking readable when zoomed
                    # out. 
                    padx = 0.25
                    pady = 0.25
                    figsz = (1+np.array([padx,pady]))*sz/self.dpi
                    fac = 2.0

                # Plot waterfall
                fig = plt.figure(figsize=figsz/fac, dpi=self.dpi*fac)
                plt.imshow(v.T, extent=(self.fmin,self.fmax,self.tmax,self.tmin),aspect='auto')
                if cl is not None:
                    plt.clim(*cl)

                plt.grid('on')
                plt.ylabel('t (minutes)');
                plt.xlabel('freq (MHz)');
                plt.title('{:s} {:s} -- {:s} - {:s} (UTC)'.format(chn,titlab,self.ts.iso,self.te.iso))
                plt.colorbar(pad=0)

                fname       = self.filebase + '_'+chn+'_wfcal_'+fext+'.png'
                fname_thumb = self.filebase + '_'+chn+'_wfcal_'+fext+'_thumbnail.png' 
                plt.savefig(os.path.join(self.plotdir,fname), dpi=self.dpi*fac, bbox_inches='tight')
                plt.savefig(os.path.join(self.plotdir,fname_thumb), dpi=self.dpi, bbox_inches='tight')

                plt.close(fig)
                gc.collect()

                # Plot median
                fig = plt.figure(figsize=(7,5))
                plt.plot(self.r.f, np.nanmedian(v,1))
                plt.xlabel('freq (MHz)');
                plt.ylabel(titlab)
                plt.title('{:s} median {:s}'.format(chn,titlab,self.r.tag))
                if cl is not None:
                    plt.ylim(*cl)

                fname = self.filebase + '_'+chn+'_medcal_'+fext+'.png'
                plt.savefig(os.path.join(self.plotdir,fname), bbox_inches='tight')

                plt.close(fig)
                gc.collect()


    def plotvariance(self, dochan=None, fld_in=['data_mf','data_svd']):
        """Plot variance and expected variance from radiometer equation"""
        
        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)

            for fld in fld_in:
                # Actual variance
                v = getattr(self.r,fld)[chan]*1.0
                v[~np.isfinite(v)]=np.nan
                v = np.nanvar(v,0)

                # Radiometer equation
                df = (self.r.f[1]-self.r.f[0])*1e6 # In Hz
                dt = self.r.dtraw * self.r.nhits[chan].mean(0) # In sec
                T = np.nanmean(self.r.data[chan], 0)
                dT = T/np.sqrt(dt*df)

                # Plot
                fig = plt.figure(figsize=(7,7))

                plt.subplot(2,1,1)
                plt.plot(self.r.f, np.sqrt(v), label='std(T)')
                plt.plot(self.r.f, dT, label='T/sqrt(df*dt)')
                plt.xlabel('freq (MHz)');
                plt.ylabel('Kelvin')
                plt.title('{:s}, {:s}, {:s}'.format(chn,fld,self.r.tag))
                plt.legend()
                plt.ylim(0,10)

                plt.subplot(2,1,2)
                plt.plot(self.r.f, np.sqrt(v)/dT, label='ratio (std(T)/[T/sqrt(df*dt)])')
                plt.xlabel('freq (MHz)');
                plt.ylabel('ratio')
                plt.legend()
                plt.ylim(0,10)
                plt.plot([self.r.f[0],self.r.f[-1]],[1,1],':k')

                fname = self.filebase + '_'+chn+'_variance_'+fld+'.png'
                plt.savefig(os.path.join(self.plotdir,fname), bbox_inches='tight')
                plt.close(fig)
                gc.collect()


    def plotps(self, dochan=None, fld_in=['data_mf','data_svd']):
        """Plot power spectrum over time for each frequency"""

        if dochan is None:
            dochan = range(self.r.d.nChan)
        elif not np.iterable(dochan):
            dochan = [dochan]

        for chan in dochan:

            chn = self.r.getchname(chan)

            for fld in fld_in:
                # Do rfft
                v = getattr(self.r,fld)[chan]*1.0
                v[~np.isfinite(v)] = 0.0
                ft = np.fft.rfft(v,norm='ortho',axis=0)
                p = np.real(ft*np.conj(ft))
                f = np.linspace(0,1/(2*self.r.dt),p.shape[0])

                # Plot waterfall
                sz = np.array(v.T.shape)
                padx = 0.25
                pady = 0.25
                figsz = (1+np.array([padx,pady]))*sz/self.dpi
                fac = 2.0

                fig = plt.figure(figsize=figsz/fac, dpi=self.dpi*fac)
                plt.imshow(np.log10(p), extent=(self.fmin,self.fmax,f[-1],f[0]),aspect='auto')
                plt.clim(-10,5)

                plt.grid('on')
                plt.ylabel('freq (Hz)');
                plt.xlabel('freq (MHz)');
                plt.title('Power spectrum (K^2), {:s}, {:s}'.format(chn,self.r.tag))
                plt.colorbar(pad=0)

                fname       = self.filebase + '_'+chn+'_pswf_'+fld+'.png'
                fname_thumb = self.filebase + '_'+chn+'_pswf_'+fld+'_thumbnail.png'
                plt.savefig(os.path.join(self.plotdir,fname_thumb), dpi=self.dpi, bbox_inches='tight')
                plt.savefig(os.path.join(self.plotdir,fname), dpi=self.dpi*fac, bbox_inches='tight')
                plt.close(fig)
                gc.collect()


                # Plot linear
                fig = plt.figure(figsize=(7,5))
                # Plot power spectrum of these frequency bins
                ind = [100, 500, 1000, 2000, 3000]
                # Plot, omit DC bin
                for k in ind:
                    if np.all(p[1:,k] <= 0):
                        continue
                    plt.loglog(f[1:], p[1:,k],label='{:0.1f}'.format(self.r.f[k]))
                plt.xlabel('freq (Hz)');
                plt.ylabel('Power (K^2)')
                plt.title('Power spectrum, {:s}, {:s}'.format(chn, self.r.tag))
                plt.ylim(1e-7,1e5)
                plt.grid('on')
                plt.legend(loc='lower left')

                fname = self.filebase + '_'+chn+'_psloglog_'+fld+'.png'
                plt.savefig(os.path.join(self.plotdir,fname), bbox_inches='tight')
                plt.close(fig)
                gc.collect()



class genhtml():
    
    def __init__(self):
        """Generate reduc plot html"""

        # Get all tags with raw data that are not cut
        dm = datamanager.datamanager()
        dm.gettags()
        # Sort tags with latest one first
        self.tags = dm.tags[::-1]
        self.htmldir = 'browser'

    def maybe_create_dirs(self):
        if not os.path.exists(self.htmldir):
            print ("Creating directory structure under",self.htmldir)
            os.makedirs(self.htmldir)
            os.makedirs(os.path.join(self.htmldir,'pages'))
            os.makedirs(os.path.join(self.htmldir,'plots'))
            
            
    def genindex(self):
        """Generate index page"""

        # Last tag
        tag = self.tags[0]

        self.maybe_create_dirs()
        f = open(os.path.join(self.htmldir,'index.html'),'w')

        f.write('<html>\n')

        f.write('<head><title>BMX data browser</title></head>\n')

        f.write('<SCRIPT LANGUAGE="JavaScript">\n')
        f.write('<!--\n')
        f.write('type="medcal_data";\n')
        f.write('chan="chan1_0";\n')

        f.write('tag=\'\';\n')
        f.write('url=\'\';\n')
        f.write('prefix="../plots/";\n')

        f.write('function plupdate(){\n')
        f.write(' url_base = prefix + "/";\n')
        f.write(' if (type.includes("wfraw")) {\n')
        f.write('   ext = ".jpg";\n')
        f.write(' } else {\n')
        f.write('   ext = ".png";\n')
        f.write(' }\n')
        f.write('\n')
        f.write(' if (type.includes("wf")) {\n')
        f.write('   tagpage.document["reduc"].style.width="100%";\n')
        f.write('   linkurl = url_base + tag + "_" + chan + "_" + type + ext;\n')
        f.write('   url     = url_base + tag + "_" + chan + "_" + type + "_thumbnail"+ext;\n')
        f.write(' } else {\n')
        f.write('   tagpage.document["reduc"].style.width="auto";\n')
        f.write('   url     = url_base + tag + "_" + chan + "_" + type + ext;\n')
        f.write('   linkurl = url_base + tag + "_" + chan + "_" + type + ext;\n')
        f.write(' }\n')
        f.write('\n')
        f.write(' tagpage.document["reduc"].src=url;\n')
        f.write(' tagpage.document.getElementById("reduc_link").href=linkurl;\n')
        f.write('}\n')
        f.write('\n')
        f.write('function set_type(plot_type){\n')
        f.write('  type=plot_type;\n')
        f.write('  plupdate();\n')
        f.write('}\n')
        f.write('\n')
        f.write('function set_chan(chan_num){\n')
        f.write('  chan=chan_num;\n')
        f.write('  plupdate();\n')
        f.write('}\n')
        f.write('\n')
        f.write('function set_tag(tagstr){\n')
        f.write('  tag=tagstr;\n')
        f.write('  plupdate();\n')
        f.write('}\n')
        f.write('\n')
        f.write('//-->\n')
        f.write('</SCRIPT>\n')
        f.write('\n')
        f.write('<frameset noresize="noresize" cols="220,*">\n')
        f.write('<frame src="tag_index.html#left_top">\n')
        f.write('<frame src="pages/{:s}.html" name="tagpage">\n'.format(tag))
        f.write('</frameset>\n')
        f.write('\n')
        f.write('\n')
        f.write('</html>\n')

        f.close()


    def gentagindex(self):
        """Generate tag index page"""

        self.maybe_create_dirs()
        f = open(os.path.join(self.htmldir,'tag_index.html'),'w')

        f.write('<html>\n')
        f.write('<body bgcolor="#d0d0d0">\n')
        f.write('<a name="left_top">\n')
        f.write('<pre>\n')

        tag0 = 'dummy'

        for tag in self.tags:

            if tag[0:6] != tag0[0:6]:
                # New day, new table
                if tag0 != 'dummy':
                    f.write('</table>\n')
                f.write('<font color="green">20{:s}</font>\n'.format(tag[0:6]))
                f.write('<table>\n')
                k = 0

            if np.mod(k,6) == 0:
                f.write('<tr>\n')
            f.write('<td><a href="pages/{:s}.html" target="tagpage">{:s}</a>\n'.format(tag,tag[7:]))
            k += 1
            tag0 = tag

    def gentagpages(self):
        """Generate individual tag pages"""

        ntags = len(self.tags)

        self.maybe_create_dirs()
        for k,tag in enumerate(self.tags):

            f = open(os.path.join(self.htmldir,'pages','{:s}.html'.format(tag)),'w')

            f.write('<html>\n')
            f.write('\n')
            f.write('<head>\n')
            f.write('\n')
            f.write('</head>\n')
            f.write('<body>\n')
            f.write('\n')
            f.write('<center><h2>BMX tag browser</h2></center>\n')
            if (k+1) < ntags:
                f.write('    <center>[<a href="{:s}.html">-tag</a>] \n'.format(self.tags[k+1]))
            else:
                f.write('    <center>[-tag] \n')
            f.write('      <b>{:s}</b>\n'.format(tag))
            if (k-1) >=0:
                f.write('      [<a href="{:s}.html">+tag</a>] \n'.format(self.tags[k-1]))
            else:
                f.write('      [+tag] \n'.format(self.tags[k-1]))
            f.write('    </center>\n')
            f.write('<hr>\n')
            f.write('<SCRIPT LANGUAGE="JavaScript">\n')
            f.write('<!--\n')
            f.write('function set_type(tt){\n')
            f.write('  parent.set_type(tt);\n')
            f.write('}\n')
            f.write('function set_chan(cc){\n')
            f.write('  parent.set_chan(cc);\n')
            f.write('}\n')
            f.write('//-->\n')
            f.write('</SCRIPT>\n')
            f.write('\n')
            f.write('<center>\n')
            f.write('\n')
            f.write('<table border="0" cellspacing="1" cellpadding="3">\n')
            f.write('<tr><th>Median spectrum: </th>\n')
            f.write('<td><a href="javascript:set_type(\'specraw\');">raw</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'medcal_data\');">calibrated</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'medcal_data_mf\');">mean filtered</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'medcal_data_svd\');">SVD filtered</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'medcal_gain\');">gain</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'medcal_weight\');">weight</a></td>\n')
            f.write('</tr>\n')
            f.write('</table>\n')
            f.write('\n')
            f.write('<table border="0" cellspacing="1" cellpadding="3">\n')
            f.write('<tr><th>Calibrated, downsampled waterfall plots: </th>\n')
            f.write('<td><a href="javascript:set_type(\'wfcal_data\');">data</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'wfcal_data_mf\');">mean filter</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'wfcal_data_svd\');">mean+SVD filter</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'wfcal_gain\');">gain</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'wfcal_weight\');">weight</a></td>\n')
            f.write('</tr>\n')
            f.write('</table>\n')
            f.write('\n')
            f.write('<table border="0" cellspacing="1" cellpadding="3">\n')
            f.write('<tr><th>Undownsampled waterfall plots: </th>\n')
            f.write('<td><a href="javascript:set_type(\'wfraw_nomask\');">raw</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'wfraw_mask\');">after masking</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'wfraw_cal\');">after cal (log)</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'wfraw_callin\');">after cal (lin)</a></td>\n')
            f.write('</tr>\n')
            f.write('</table>\n')
            f.write('\n')
            f.write('<table border="0" cellspacing="1" cellpadding="3">\n')
            f.write('<tr><th>Power spectra: </th>\n')
            f.write('<td><a href="javascript:set_type(\'pswf_data_mf\');">WF</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'pswf_data_svd\');">WF (after SVD)</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'psloglog_data_mf\');">loglog</a></td>\n')
            f.write('<td><a href="javascript:set_type(\'psloglog_data_svd\');">loglog (after SVD)</a></td>\n')
            f.write('</tr>\n')
            f.write('</table>\n')
            f.write('\n')
            f.write('<table border="0" cellspacing="1" cellpadding="3">\n')
            f.write('<tr><th>channel: </th><td><a href="javascript:set_chan(\'chan1_0\');">1</a></td>\n')
            f.write('<td><a href="javascript:set_chan(\'chan2_0\');">2</a></td>\n')
            f.write('</tr>\n')
            f.write('</table>\n')
            f.write('\n')
            f.write('<a href="../plots/{:s}_chan1_0_medcal_data.png" id="reduc_link">\n'.format(tag))
            f.write('<img src="../plots/{:s}_chan1_0_medcal_data.png" width=""  name="reduc">\n'.format(tag))
            f.write('</a>\n')
            f.write('</center>\n')
            f.write('\n')
            f.write('<SCRIPT LANGUAGE="JavaScript">\n')
            f.write('<!--\n')
            f.write('parent.tag=\'{:s}\';\n'.format(tag))
            f.write('parent.plupdate(); \n')  
            f.write('//-->\n')
            f.write('</SCRIPT>\n')
            f.write('\n')
            f.write('\n')
            f.write('</body>\n')
            f.write('\n')
            f.write('</html>\n')

            f.close()


            
            
        
                 
        

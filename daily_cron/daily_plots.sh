#!/bin/bash
#source ~astrodat/setup/setup-wq.sh
#alias wq='wq.exe'
source /astro/u/bmx/.bashrc
export BMXDATA=/gpfs01/astro/workarea/bmxdata/raw
export BMXREDUCED=/gpfs01/astro/workarea/bmxdata/reduced
export PATH=/gpfs01/astro/packages/anaconda/default/bin:$PATH
export LOGFN=/direct/astro+u/bmx/bmxreduce/daily_cron/logs/`date +%y%m%d`.log
export LOGFNE=/direct/astro+u/bmx/bmxreduce/daily_cron/logs/`date +%y%m%d`.err
echo "---------- daily_plots.sh: `date` ----------"  >>$LOGFN 2>>$LOGFNE
echo "----------------- calling QA plots --------------------" >>$LOGFN 2>>$LOGFNE
/astro/u/bmx/bmxreduce/qaplots/QA_plots.py -v  >>$LOGFN 2>>$LOGFNE
echo "------------------- DONE -----------------------" >>$LOGFN 2>>$LOGFNE



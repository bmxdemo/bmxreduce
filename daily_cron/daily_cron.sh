#!/bin/bash
source ~astrodat/setup/setup-wq.sh
alias wq='wq.exe'
export BMXDATA=/gpfs01/astro/workarea/bmxdata/raw
export BMXREDUCED=/gpfs01/astro/workarea/bmxdata/reduced
export PATH=/gpfs01/astro/packages/anaconda/default/bin:$PATH
export LOGFN=/direct/astro+u/bmx/daily_cron/logs/`date +%y%m%d`.log
export LOGFNE=/direct/astro+u/bmx/daily_cron/logs/`date +%y%m%d`.err
echo "-----------------finalize_transfer------------" >$LOGFN 2>$LOGFNE
/direct/astro+u/bmx/bmxreduce/daily_cron/finalize_transfer.py -v >> $LOGFN 2>$LOGFNE
echo "-----------------BMXREDUCE--------------------" >>$LOGFN 2>$LOGFNE
cd /direct/astro+u/bmx/bmxreduce 
bin/reduce.py >>$LOGFN 2>$LOGFNE
echo "-------------------DONE-----------------------" >>$LOGFN 2>$LOGFNE


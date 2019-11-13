#!/bin/bash
#source ~astrodat/setup/setup-wq.sh
#alias wq='wq.exe'
source /astro/u/bmx/.bashrc
export BMXDATA=/gpfs01/astro/workarea/bmxdata/raw
export BMXREDUCED=/gpfs01/astro/workarea/bmxdata/reduced
export PATH=/gpfs01/astro/packages/anaconda/default/bin:$PATH
export LOGFN=/direct/astro+u/bmx/bmxreduce/daily_cron/logs/`date +%y%m%d`.log
export LOGFNE=/direct/astro+u/bmx/bmxreduce/daily_cron/logs/`date +%y%m%d`.err
echo "-----------------finalize_transfer------------" >$LOGFN 2>$LOGFNE
/direct/astro+u/bmx/bmxreduce/daily_cron/finalize_transfer.py -v >> $LOGFN 2>>$LOGFNE
echo "------------------downloading almanac ------------">$LOGFN 2>$LOGFNE
/direct/astro+u/bmx/bmxreduce/daily_cron/download_almanac.py -v >> $LOGFN 2>>$LOGFNE
echo "-----------------BMXREDUCE--------------------" >>$LOGFN 2>>$LOGFNE
#echo "            --- disabled --- " >> >>$LOGFN 2>>$LOGFNE
which python >>$LOGFN 2>>$LOGFNE
#echo $PATH  >>$LOGFN 2>>$LOGFNE
#export PYTHONPATH=/astro/u/bmx/bmxdaq/py:/astro/u/bmx/bmxdaq:/astro/u/bmx/bmxreduce/bmxreduce:/astro/u/bmx/bmxreduce:/astro/u/astrodat/local/lib/python2.7/site-packages:$PYTHONPATH
#export PATH=/astro/u/astrodat/local/wq/binary/bin:$PATH
cd /direct/astro+u/bmx/bmxreduce 
echo "-----------------calling reduce--------------------" >>$LOGFN 2>>$LOGFNE
bin/reduce.py  >>$LOGFN 2>>$LOGFNE
echo "-------------------DONE-----------------------" >>$LOGFN 2>>$LOGFNE


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
dayName=`date +%a | tr '[:lower:]' '[:upper:]'`
if [ $dayName = "TUE" ]
then
wget https://datacenter.iers.org/data/9/finals2000A.all -O //gpfs02/astro/workarea/bmxdata/almanac/finals2000A.all >> $LOGFN 2>>$LOGFNE
fi
echo "-----------------BMXREDUCE--------------------" >>$LOGFN 2>>$LOGFNE
#echo "            --- disabled --- " >> >>$LOGFN 2>>$LOGFNE
declare -x ANACONDA_DIR="/gpfs01/astro/packages/anaconda/default"
declare -x CONDA_EXE="/gpfs01/astro/packages/anaconda3/bin/conda"
declare -x CONDA_PREFIX="/astro/u/bmx/.conda/envs/bmxpy"
declare -x CONDA_PYTHON_EXE="/gpfs01/astro/packages/anaconda3/bin/python"
declare -x LOADEDMODULES="afs/SL64:anaconda/gpfs:astrodat/1.0:wq/binary"
declare -x PATH="/astro/u/bmx/.conda/envs/bmxpy/bin:/gpfs01/astro/packages/anaconda3/condabin:/gpfs01/astro/packages/anaconda3/bin:/astro/u/bmx/local/bin:/astro/u/astrodat/local/bin:/astro/u/astrodat/local/wq/binary/bin:/gpfs01/astro/packages/anaconda/default/bin:/astro/u/bmx/bin:/opt/astro/SL64/bin:/usr/afsws/bin:/usr/bin:/etc:/usr/sbin:/usr/ucb:/astro/u/bmx/bin:/usr/bin/X11:/sbin:."
which python >>$LOGFN 2>>$LOGFNE
#echo $PATH  >>$LOGFN 2>>$LOGFNE
#export PYTHONPATH=/astro/u/bmx/bmxdaq/py:/astro/u/bmx/bmxdaq:/astro/u/bmx/bmxreduce/bmxreduce:/astro/u/bmx/bmxreduce:/astro/u/astrodat/local/lib/python2.7/site-packages:$PYTHONPATH
#export PATH=/astro/u/astrodat/local/wq/binary/bin:$PATH
cd /direct/astro+u/bmx/bmxreduce 
echo "-----------------calling reduce--------------------" >>$LOGFN 2>>$LOGFNE
bin/reduce.py cron  >>$LOGFN 2>>$LOGFNE
echo "-----------------calling QA plots--------------------" >>$LOGFN 2>>$LOGFNE
/astro/u/bmx/bmxreduce/qaplots/QA_plots.py -v  >>$LOGFN 2>>$LOGFNE
echo "-------------------DONE-----------------------" >>$LOGFN 2>>$LOGFNE


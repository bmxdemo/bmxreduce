#!/bin/bash
/usr/bin/ssh bmxfer1 cat bmxdaq/log/daq.log | /usr/bin/tac > /gpfs02/astro/www/bmx/daqlog/bmx1.log
/usr/bin/ssh bmxfer2 cat bmxdaq/log/daq.log | /usr/bin/tac > /gpfs02/astro/www/bmx/daqlog/bmx2.log
source /astro/u/bmx/.bashrc
python  /astro/u/bmx/bmxreduce/astroftp_cron/check_log.py >  /astro/u/bmx/bmxreduce/astroftp_cron/check_log.log 2>/astro/u/bmx/bmxreduce/astroftp_cron/check_log.err

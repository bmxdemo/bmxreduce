#!/bin/bash
/usr/bin/ssh bmxfer1 cat bmxdaq/log/daq.log | /usr/bin/tac > /gpfs02/astro/www/bmx/daqlog/bmx1.log
/usr/bin/ssh bmxfer2 cat bmxdaq/log/daq.log | /usr/bin/tac > /gpfs02/astro/www/bmx/daqlog/bmx2.log

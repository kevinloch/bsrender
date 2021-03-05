#!/bin/sh
# 
# This script downloads the full set of pre-generated GEDR3 data files for bsrender.
# You do not have to use 'gaia-edr3-extract.sh' and 'mkgalaxy' if you have downloaded these
# files.
#
# Note: the lower parallax quality files are very large.  Total size of all files is 43GB.
# You do not need all of these files, only down to the parallax quality level you intend
# to use in bsrender.cfg.  pq100 through pq010 is a total of 3.2GB.
#
datahost="https://kevinloch.com"
#
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-external.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq100.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq050.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq030.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq020.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq010.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq005.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq003.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq002.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq001.dat
wget --no-check-certificate ${datahost}/bsrender/sample_data/0.9.0/galaxy-pq000.dat

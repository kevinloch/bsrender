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
datadir="/bsrender/sample_data/0.9.0-dev-42"
#
wget ${datahost}${datadir}/README.txt
wget ${datahost}${datadir}/external.csv
wget ${datahost}${datadir}/galaxy-external.dat
wget ${datahost}${datadir}/galaxy-pq100.dat
wget ${datahost}${datadir}/galaxy-pq050.dat
wget ${datahost}${datadir}/galaxy-pq030.dat
wget ${datahost}${datadir}/galaxy-pq020.dat
wget ${datahost}${datadir}/galaxy-pq010.dat
wget ${datahost}${datadir}/galaxy-pq005.dat
wget ${datahost}${datadir}/galaxy-pq003.dat
wget ${datahost}${datadir}/galaxy-pq002.dat
wget ${datahost}${datadir}/galaxy-pq001.dat
wget ${datahost}${datadir}/galaxy-pq000.dat

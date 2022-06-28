#!/bin/sh
# 
# This script downloads the full set of pre-generated GDR3 data files for bsrender.
# You do not have to use 'gaia-dr3-extract.sh' and 'mkgalaxy' if you have downloaded these
# files. These sample files are in little-endian byte-order.
#
# Note: the lower parallax quality files are very large. Total size of all files is 47GB.
# You do not need all of these files, only down to the parallax quality level you intend
# to use in bsrender.cfg. pq100 through pq010 is a total of 3.2GB.
#
dlprog="wget"
datahost="https://bsrender.io"
datadir="/sample_data/1.0-dev"
#
${dlprog} ${datahost}${datadir}/README.txt
${dlprog} ${datahost}${datadir}/external.csv
${dlprog} ${datahost}${datadir}/galaxy-external-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq100-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq050-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq030-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq020-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq010-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq005-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq003-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq002-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq001-le.bsr
${dlprog} ${datahost}${datadir}/galaxy-gdr3-pq000-le.bsr

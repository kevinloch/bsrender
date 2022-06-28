#!/bin/sh
#
# Billion Star 3D Rendering for Gaia DR3
# Kevin M. Loch
#
# Data extractor for the ESA Gaia DR3 star dataset
#
# This script expects the compressed Gaia DR3 source fiiles to be in ./gaia_source
# It will extract the necessary fields into extracted.csv for processing by mkgalaxy
#

extract=gaia-dr3-extracted.csv

:>${extract}

for f in `ls gaia_source`
do
  echo "Extracting from ${f}"
  gzcat gaia_source/${f} | grep -v '#' | cut -f 3,6,8,10,12,36,38,39,67,72,77,116,132,141,147 -d "," >> ${extract}
done

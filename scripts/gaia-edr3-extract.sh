#!/bin/sh
#
# Billion Star 3D Rendering Engine Proof of Concept
# Kevin M. Loch
#
# Data extractor for the ESA Gaia EDR3 star dataset
#
# This script expects the compressed Gaia EDR3 source fiiles to be in ./gaia_source
# It will extract the necessary fields into extracted.csv for processing by mkgalaxy
#

extract=gaia-edr3-extracted.csv

:>${extract}

for f in `ls gaia_source`
do
  echo "Extracting from ${f}"
  gzcat gaia_source/${f} | cut -f 4,6,8,10,12,36,38,39,67,72,77,99 -d "," >> ${extract}
done

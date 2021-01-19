#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>

void printUsage() {

printf("bsrender version %s\n", BSR_VERSION);
printf("\n\
NAME\n\
     bsrender -- Render PNG image from 3D star database\n\
\n\
SYNOPSIS\n\
     bsrender [-c filename] [-d path] [-h]\n\
 \n\
OPTIONS:\n\
     -c filename\n\
          Set configuration file name (default: ./bsrender.cfg)\n\
\n\
     -d path\n\
          Set path to data files (default: ./galaxydata)\n\
\n\
     -h\n\
          Show help\n\
\n\
DESCRIPTON\n\
 bsrender (Billion Star Rendering Engine) is designed to handle billions of stars such as those in the ESA's Gaia EDR3 dataset.\n\
 \n");
}

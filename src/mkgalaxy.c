//
// Billion Star 3D Rendering Engine
// Kevin M. Loch
//
// 3D rendering engine for the ESA Gaia DR3 star dataset

/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2021, Kevin Loch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//
// pre-processor for the ESA Gaia DR3 star dataset
//
// This program creates binary data files for use by the rendering engine
// Data is sorted into separate files by "Gaia paralax quality" which is the value of the 'parallax_over_error' field in the GDR3 data set
// galaxy-pq030 for example has stars with a minimum parallax quality of 30.0, but less than the next highest tier (50.0)
//

//#define DEBUG

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "bandpass-ratio.h"
#include "util.h"

void printUsage() {
  printf("mkgalaxy version %s\n", BSR_VERSION);
  printf("\n\
NAME\n\
     mkgalaxy -- create binary data files for use with bsrender\n\
\n\
SYNOPSIS\n\
     mkgalaxy [-b] [-w] [-d] [-p] [-c] [-n] [-m] [-l] [-g] [-h]\n\
 \n\
OPTIONS:\n\
     -b\n\
          Use bp/G and rp/G bandpass ratios to estimate apparent temperature (default)\n\
\n\
     -w\n\
          Use Gaia DR3 color wavenumber field to estimate apparent temperature\n\
\n\
     -d\n\
          Use Gaia DR3 'gspphot_distance' field instead of 'parallax' when available\n\
\n\
     -p\n\
          Only use Gaia DR3 'parallax' for star distance (default)\n\
\n\
     -c\n\
          Enable Lindegren et al parallax calibration\n\
\n\
     -n\n\
          Do not use Lindegren et al parallax calibration (default)\n\
\n\
     -m\n\
          Maximum star distance from Earth. Stars further away will be fixed to this value. Set to zero to disable (default is 50,000 parsecs)\n\
\n\
     -l\n\
          Force output to little-endian format (default is to match this platform)\n\
\n\
     -g\n\
          Force output big-endian format (default is to match this platform)\n\
\n\
     -h\n\
          Show help\n\
\n\
DESCRIPTON\n\
 mkgalaxy processes extracted fields from ESA's Gaia DR3 dataset for use with bsrender. Uses the output from 'gaia-dr3-extract.sh' in the bsrender package\n\
 \n");
}

int setDefaults(mkg_config_t *mkg_config) {
  int little_endian;

  mkg_config->use_bandpass_ratios=1;
  mkg_config->use_gspphot_distance=0;
  mkg_config->calibrate_parallax=0;
  mkg_config->enable_maximum_distance=1;
  mkg_config->maximum_distance=50000.0;

  little_endian=littleEndianTest();
  if (little_endian == 1) {
    mkg_config->output_little_endian=1;
  } else {
    mkg_config->output_little_endian=0;
  }

  return(0);
}

int processCmdArgs(mkg_config_t *mkg_config, int argc, char **argv) {
  int i;
  char *option_start;
  char tmpstr[32];
  int option_length;

  if (argc == 1) {
    return(0);
  } else {
    for (i=1; i <= (argc - 1); i++) {
      if (argv[i][1] == 'b') {
        // use bandpass ratios for apparent temperature
        mkg_config->use_bandpass_ratios=1;
      } else if (argv[i][1] == 'w') {
        // use Gaia DR3 wavenumber for apparent temperature
        mkg_config->use_bandpass_ratios=0;
      } else if (argv[i][1] == 'd') {
        // use DR3 'gspphot_distance' instead of 'parallax' when available 
        mkg_config->use_gspphot_distance=1;
      } else if (argv[i][1] == 'p') {
        // only use DR3 'parallax'
        mkg_config->use_gspphot_distance=0;
      } else if (argv[i][1] == 'c') {
        // enable Lindegren et al parallax calibration mode
        mkg_config->calibrate_parallax=1;
      } else if (argv[i][1] == 'n') {
        // do not enable Lindegren et al parallax calibration mode
        mkg_config->calibrate_parallax=0;
      } else if (argv[i][1] == 'm') {
        // maximum distance (parsecs)
        if (argv[i][2] != 0) {
          // option concatenated onto switch
          option_start=argv[i];
          option_length=strnlen(option_start + (size_t)2, 31);
          if (option_length > 31) {
            option_length=31;
          }
          strncpy(tmpstr, (option_start + (size_t)2), option_length);
          tmpstr[option_length]=0;
          mkg_config->maximum_distance=strtod(tmpstr, NULL);
          if (mkg_config->maximum_distance <=0.0) {
            mkg_config->enable_maximum_distance=0;
            mkg_config->maximum_distance=0.0;
          } else {
            mkg_config->enable_maximum_distance=1;
          }
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          option_length=strnlen(option_start, 31);
          if (option_length > 31) {
            option_length=31;
          }
          strncpy(tmpstr, option_start, option_length);
          tmpstr[option_length]=0;
          mkg_config->maximum_distance=strtod(tmpstr, NULL);
          if (mkg_config->maximum_distance <=0.0) {
            mkg_config->enable_maximum_distance=0;
          } else {
            mkg_config->enable_maximum_distance=1;
          }
        } // end if no space
      } else if (argv[i][1] == 'l') {
        // force output to littl-endian
        mkg_config->output_little_endian=1;
      } else if (argv[i][1] == 'g') {
        // force output to big-endian
        mkg_config->output_little_endian=0;
      } else if (argv[i][1] == 'h') {
        // print help
        printUsage();
        exit(0);
      } // end which option
    } // end for argc
  } // end if any options
  return(0);
}

double calibrateParallax(double *parallax, int astrometric_params_solved, double G, double neff, double ecl_lat) {
  //
  // from Lindegren et. al., Gaia Early Data Release 3 Parallax bias versus magnitude, colour, and position, Astronomy & Astrophysics manuscript no. DR3-Parallaxes December 8, 2020
  //
  int j;
  int k;
  double Z;
  double q[5][3];
  double c[5];
  double b[3];
  double new_parallax;
  double low; // interpolate low ratio
  double high; // interpolate high ratio
  // note: G=magnitude, neff=nu_eff or pseudocolor

  //
  // initialize q
  //
  for (j=0; j <= 4; j++) {
    for (k=0; k <=2; k++) {
      q[j][k]=0.0;
    }
  }
  if (astrometric_params_solved == 31) {
    // Z5 coefficients for five-parameter solution based on magnitude G
    if (G < 6.0) {
      q[0][0]= -26.98;
      q[0][1]= -9.62;
      q[0][2]= +27.40;
      q[1][0]= -25.1;
      q[1][1]=  +15.7;
      q[2][0]= -1257.0;
    } else if (G < 10.8) {
      high=(G - 6.0) / 4.8;
      low=1.0 - high;
      q[0][0]= (low * -26.98) + (high * -27.23);
      q[0][1]= (low * -9.62)  + (high * -3.07);
      q[0][2]= (low * +27.40) + (high * +23.04);
      q[1][0]= (low * -25.1)  + (high * +35.3);
      q[1][1]= high * +15.7;
      q[2][0]= -1257.0;
    } else if (G < 11.2) {
      high=(G - 10.8) / 0.4;
      low=1.0 - high;
      q[0][0]= (low * -27.23) + (high * -30.33);
      q[0][1]= (low * -3.07)  + (high * -9.23);
      q[0][2]= (low * +23.04) + (high * +9.08);
      q[1][0]= (low * +35.3)  + (high * -88.4);
      q[1][1]= (low * +15.7)  + (high * -11.8);
      q[2][0]= -1257.0;
   } else if (G < 11.8) {
      high=(G - 11.2) / 0.6;
      low=1.0 - high;
      q[0][0]= (low * -30.33) + (high * -33.54);
      q[0][1]= (low * -9.23)  + (high * -10.08);
      q[0][2]= (low * +9.08)  + (high * +13.28);
      q[1][0]= (low * -88.4)  + (high * -126.7);
      q[1][1]= (low * -11.8)  + (high * +11.6);
      q[2][0]= -1257.0;
   } else if (G < 12.2) {
      high=(G - 11.8) / 0.6;
      low=1.0 - high;
      q[0][0]= (low * -33.54) + (high * -13.65);
      q[0][1]= (low * -10.08) + (high * -0.07);
      q[0][2]= (low * +13.28) + (high * +9.35);
      q[1][0]= (low * -126.7) + (high * -111.4);
      q[1][1]= (low * +11.6)  + (high * +40.6);
      q[2][0]= -1257.0;
   } else if (G < 12.9) {
      high=(G - 12.2) / 0.7;
      low=1.0 - high;
      q[0][0]= (low * -13.65) + (high * -19.53);
      q[0][1]= (low * -0.07)  + (high * -1.64);
      q[0][2]= (low * +9.35)  + (high * +15.86);
      q[1][0]= (low * -111.4) + (high * -66.8);
      q[1][1]= (low * +40.6)  + (high * +20.6);
      q[2][0]= -1257.0;
   } else if (G < 13.1) {
      high=(G - 12.9) / 0.2;
      low=1.0 - high;
      q[0][0]= (low * -19.53) + (high * -37.99);
      q[0][1]= (low * -1.64)  + (high * +2.63);
      q[0][2]= (low * +15.86) + (high * +16.14);
      q[1][0]= (low * -66.8)  + (high * -5.7);
      q[1][1]= (low * +20.6)  + (high * +14.0);
      q[2][0]= -1257.0;
      q[3][0]= high * +107.9;
      q[4][0]= high * +104.3;
   } else if (G < 15.9) {
      high=(G - 13.1) / 2.8;
      low=1.0 - high;
      q[0][0]= (low * -37.99)  + (high * -38.33);
      q[0][1]= (low * +2.63)   + (high * +5.61);
      q[0][2]= (low * +16.14)  + (high * +15.42);
      q[1][0]= low * -5.7;
      q[1][1]= (low * +14.0)   + (high * +18.7);
      q[2][0]= (low * -1257.0) + (high * -1189.0);
      q[3][0]= (low * +107.9)  + (high * +243.8);
      q[4][0]= (low * +104.3)  + (high * +155.2);
   } else if (G < 16.1) {
      high=(G - 15.9) / 0.2;
      low=1.0 - high;
      q[0][0]= (low * -38.33)  + (high * -31.05);
      q[0][1]= (low * +5.61)   + (high * +2.83);
      q[0][2]= (low * +15.42)  + (high * +8.59);
      q[1][1]= (low * +18.7)   + (high * +15.5);
      q[2][0]= (low * -1189.0) + (high * -1404.0);
      q[3][0]= (low * +243.8)  + (high * +105.5);
      q[4][0]= (low * +155.2)  + (high * +170.7);
   } else if (G < 17.5) {
      high=(G - 16.1) / 1.4;
      low=1.0 - high;
      q[0][0]= (low * -31.05)  + (high * -29.18);
      q[0][1]= (low * +2.83)   + (high * -0.09);
      q[0][2]= (low * +8.59)   + (high * +2.41);
      q[1][1]= (low * +15.5)   + (high * +24.5);
      q[2][0]= (low * -1404.0) + (high * -1165.0);
      q[3][0]= (low * +105.5)  + (high * +189.7);
      q[4][0]= (low * +170.7)  + (high * +325.0);
   } else if (G < 19.0) {
      high=(G - 17.5) / 1.5;
      low=1.0 - high;
      q[0][0]= (low * -29.18) + (high * -18.40);
      q[0][1]= (low * -0.09)  + (high * +5.98);
      q[0][2]= (low * +2.41)  + (high * -6.46);
      q[1][1]= (low * +24.5)  + (high * +5.5);
      q[2][0]= low * -1165.0;
      q[3][0]= low * +189.7;
      q[4][0]= (low * +325.0) + (high * +276.6);
   } else if (G < 20.0) {
      high=(G - 19.0) / 1.0;
      low=1.0 - high;
      q[0][0]= (low * -18.40) + (high * -12.65);
      q[0][1]= (low * +5.98)  + (high * -4.57);
      q[0][2]= (low * -6.46)  + (high * -7.46);
      q[1][1]= (low * +5.5)   + (high * +97.9);
      q[4][0]= low * +276.6;
   } else if (G < 21.0) {
      high=(G - 20.0) / 1.0;
      low=1.0 - high;
      q[0][0]= (low * -12.65) + (high * -18.22);
      q[0][1]= (low * -4.57)  + (high * -15.24);
      q[0][2]= (low * -7.46)  + (high * -18.54);
      q[1][1]= (low * +97.9)  + (high * +128.2);
   } else {
      q[0][0]= -18.22;
      q[0][1]= -15.24;
      q[0][2]= -18.54;
      q[1][1]= +128.2;
   }
  } else if (astrometric_params_solved == 95) {
    // Z6 coefficients for six-parameter solution based on magnitude G
    if (G < 6.0) {
      q[0][0]= -27.85;
      q[0][1]= -7.78;
      q[0][2]= +27.47;
      q[1][0]= -32.1;
      q[1][1]= +14.4;
      q[1][2]= +9.5;
      q[2][0]= -67.0;
    } else if (G < 10.8) {
      high=(G - 6.0) / 4.8;
      low=1.0 - high;
      q[0][0]= (low * -27.85) + (high * -28.91);
      q[0][1]= (low * -7.78)  + (high * -3.57);
      q[0][2]= (low * +27.47) + (high * +22.92);
      q[1][0]= (low * -32.1)  + (high * +7.7);
      q[1][1]= (low * +14.4)  + (high * +12.6);
      q[1][2]= (low * +9.5)   + (high * +1.6);
      q[2][0]= (low * -67.0)  + (high * -572.0);
    } else if (G < 11.2) {
      high=(G - 10.8) / 0.4;
      low=1.0 - high;
      q[0][0]= (low * -28.91) + (high * -26.72);
      q[0][1]= (low * -3.57)  + (high * -8.74);
      q[0][2]= (low * +22.92) + (high * +9.36);
      q[1][0]= (low * +7.7)   + (high * -30.3);
      q[1][1]= (low * +12.6)  + (high * +5.6);
      q[1][2]= (low * +1.6)   + (high * +17.2);
      q[2][0]= (low * -572.0) + (high * -1104.0);
    } else if (G < 11.8) {
      high=(G - 11.2) / 0.6;
      low=1.0 - high;
      q[0][0]= (low * -26.72)  + (high * -29.04);
      q[0][1]= (low * -8.74)   + (high * -9.69);
      q[0][2]= (low * +9.36)   + (high * +13.63);
      q[1][0]= (low * -30.3)   + (high * -49.4);
      q[1][1]= (low * +5.6)    + (high * +36.3);
      q[1][2]= (low * +17.2)   + (high * +17.7);
      q[2][0]= (low * -1104.0) + (high * -1129.0);
    } else if (G < 12.2) {
      high=(G - 11.8) / 0.4;
      low=1.0 - high;
      q[0][0]= (low * -29.04)  + (high * -12.39);
      q[0][1]= (low * -9.69)   + (high * -2.16);
      q[0][2]= (low * +13.63)  + (high * +10.23);
      q[1][0]= (low * -49.4)   + (high * -92.6);
      q[1][1]= (low * +36.3)   + (high * +19.8);
      q[1][2]= (low * +17.7)   + (high * +27.6);
      q[2][0]= (low * -1129.0) + (high * -365.0);
    } else if (G < 12.9) {
      high=(G - 12.2) / 0.7;
      low=1.0 - high;
      q[0][0]= (low * -12.39) + (high * -18.99);
      q[0][1]= (low * -2.16)  + (high * -1.93);
      q[0][2]= (low * +10.23) + (high * +15.90);
      q[1][0]= (low * -92.6)  + (high * -57.2);
      q[1][1]= (low * +19.8)  + (high * -8.0);
      q[1][2]= (low * +27.6)  + (high * +19.9);
      q[2][0]= (low * -365.0) + (high * -554.0);
    } else if (G < 13.1) {
      high=(G - 12.9) / 0.2;
      low=1.0 - high;
      q[0][0]= (low * -18.99) + (high * -38.29);
      q[0][1]= (low * -1.93)  + (high * +2.59);
      q[0][2]= (low * +15.90) + (high * +16.20);
      q[1][0]= (low * -57.2)  + (high * -10.5);
      q[1][1]= (low * -8.0)   + (high * +1.4);
      q[1][2]= (low * +19.9)  + (high * +0.4);
      q[2][0]= (low * -554.0) + (high * -960.0);
    } else if (G < 15.9) {
      high=(G - 13.1) / 2.8;
      low=1.0 - high;
      q[0][0]= (low * -38.29) + (high * -36.83);
      q[0][1]= (low * +2.59)  + (high * +4.20);
      q[0][2]= (low * +16.20) + (high * +15.76);
      q[1][0]= (low * -10.5)  + (high * +22.3);
      q[1][1]= (low * +1.4)   + (high * +11.1);
      q[1][2]= (low * +0.4)   + (high * +10.0);
      q[2][0]= (low * -960.0) + (high * -1367.0);
    } else if (G < 16.1) {
      high=(G - 15.9) / 0.2;
      low=1.0 - high;
      q[0][0]= (low * -36.83)  + (high * -28.37);
      q[0][1]= (low * +4.20)   + (high * +1.99);
      q[0][2]= (low * +15.76)  + (high * +9.28);
      q[1][0]= (low * +22.3)   + (high * +50.4);
      q[1][1]= (low * +11.1)   + (high * +17.2);
      q[1][2]= (low * +10.0)   + (high * +13.7);
      q[2][0]= (low * -1367.0) + (high * -1351.0);
    } else if (G < 17.5) {
      high=(G - 16.1) / 1.4;
      low=1.0 - high;
      q[0][0]= (low * -28.37) + (high * -24.68);
      q[0][1]= (low * +1.99)  + (high * -1.37);
      q[0][2]= (low * +9.28)  + (high * +3.52);
      q[1][0]= (low * +50.4)  + (high * +86.8);
      q[1][1]= (low * +17.2)  + (high * +19.8);
      q[1][2]= (low * +13.7)  + (high * +21.3);
      q[2][0]=(low * -1351.0) + (high * -1380.0);
    } else if (G < 19.0) {
      high=(G - 17.5) / 1.5;
      low=1.0 - high;
      q[0][0]= (low * -24.68)  + (high * -15.3);
      q[0][1]= (low * -1.37)   + (high * +4.01);
      q[0][2]= (low * +3.52)   + (high * -6.03);
      q[1][0]= (low * +86.8)   + (high * +29.2);
      q[1][1]= (low * +19.8)   + (high * +14.1);
      q[1][2]= (low * +21.3)   + (high * +0.4);
      q[2][0]= (low * -1380.0) + (high * -563.0);
    } else if (G < 20.0) {
      high=(G - 19.0) / 1.0;
      low=1.0 - high;
      q[0][0]= (low * -15.32) + (high * -13.73);
      q[0][1]= (low * +4.01)  + (high * -10.92);
      q[0][2]= (low * -6.03)  + (high * -8.3);
      q[1][0]= (low * +29.2)  + (high * -74.4);
      q[1][1]= (low * +14.1)  + (high * +196.4);
      q[1][2]= (low * +0.4)   + (high * -42.0);
      q[2][0]= (low * -563.0) + (high * +536.0);
    } else if (G < 21.0) {
      high=(G - 20.0) / 1.0;
      low=1.0 - high;
      q[0][0]= (low * -13.73) + (high * -29.53);
      q[0][1]= (low * -10.92) + (high * -20.34);
      q[0][2]= (low * -8.3)   + (high * -18.74);
      q[1][0]= (low * -74.4)  + (high * -39.5);
      q[1][1]= (low * +196.4) + (high * +326.8);
      q[1][2]= (low * -42.0)  + (high * -262.3);
      q[2][0]= (low * +536.0) + (high * +1598.0);
    } else {
      q[0][0]= -29.53;
      q[0][1]= -20.34;
      q[0][2]= -18.74;
      q[1][0]= -39.5;
      q[1][1]= +326.8;
      q[1][2]= -262.3;
      q[2][0]= +1598.0;
    }
  } else {
    return(1); // we should not get here but just in case
  }

  //
  // initialize c
  //
  c[0]=1.0;
  if (neff <= 1.24) {
    c[1]=-0.24;
  } else if (neff <= 1.72) {
    c[1]=neff - 1.48;
  } else {
    c[1]=0.24;
  }
  if (neff <= 1.24) {
    c[2]=pow(0.24, 3.0);
  } else if (neff <= 1.48) {
    c[2]=pow((1.48 - neff), 3.0);
  } else {
    c[2]=0.0;
  }
  if (neff <= 1.24) {
    c[3]=neff - 1.24;
  } else {
    c[3]=0.0;
  }
  if (neff <= 1.72) {
    c[4]=0.0;
  } else {
    c[4]=neff - 1.72;
  }

  //
  // initialize b
  //
  b[0]=1.0;
  b[1]=sin(ecl_lat * M_PI / 180.0);
  b[2]=pow(sin((ecl_lat * M_PI / 180.0) - (1.0 / 3.0)), 2.0);

  //
  // apply calibration
  //
  Z=0;
  for (j=0; j <= 4; j++) {
    for (k=0; k <= 2; k++) {
      Z+=q[j][k] * c[j] * b[k];
    }
  }
  Z=Z / 1000.0; // Z is in uas, parallax is in miliarcseconds
  new_parallax = *parallax - Z;
  *parallax=new_parallax;

  return(new_parallax); // necessary for -Ofast to not skip function
}

int main(int argc, char **argv) {
  struct timespec overall_starttime;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  double overall_elapsed_time;
  FILE *input_file;
  char file_name[256];
  FILE *output_file_pq000;
  FILE *output_file_pq001;
  FILE *output_file_pq002;
  FILE *output_file_pq003;
  FILE *output_file_pq005;
  FILE *output_file_pq010;
  FILE *output_file_pq020;
  FILE *output_file_pq030;
  FILE *output_file_pq050;
  FILE *output_file_pq100;
  char *input_line_p;
  char input_line[256];
  char *field_start;
  char *field_end;
  size_t field_length;
  char tmpstr[32];
  double color_wavenumber;
  uint64_t input_count;
  uint64_t pq000_count;
  uint64_t pq001_count;
  uint64_t pq002_count; 
  uint64_t pq003_count;
  uint64_t pq005_count;
  uint64_t pq010_count;
  uint64_t pq020_count;
  uint64_t pq030_count;
  uint64_t pq050_count;
  uint64_t pq100_count;
  uint64_t total_output_count;
  uint64_t discard_no_flux_count;
  uint64_t discard_parms_count;
  uint64_t discard_parallax_count;
  uint64_t distance_override_negative_count;
  uint64_t distance_override_toohigh_count;
  uint64_t temperature_from_bp_G_count;
  uint64_t temperature_from_rp_G_count;
  uint64_t temperature_from_bp_rp_count;
  uint64_t temperature_from_nu_eff_count;
  uint64_t temperature_from_pseudocolor_count;
  uint64_t min_temp_count;
  uint64_t max_temp_count;
  uint64_t unreddened_min_temp_count;
  uint64_t unreddened_max_temp_count;
  uint64_t gspphot_distance_count;
  uint64_t undimmed_count;
  uint64_t unreddened_count;
  const double flux_to_vega=5.3095E-11; // approximate conversion factor for phot_g_mean_flux to intensity relative to Vega
  float linear_1pc_intensity;
  float linear_1pc_intensity_undimmed;
  uint64_t color_temperature;
  uint64_t color_temperature_unreddened;
  double bp_over_G_ref[32768];
  double rp_over_G_ref[32768];
  double bp_over_rp_ref[32768];
  int i;
  double rp_over_G;
  double bp_over_G;
  double bp_over_rp;
  int bestmatch_temperature;
  int all_invalid;
  int bp_over_G_invalid;
  int rp_over_G_invalid;
  int bp_over_rp_invalid;
  mkg_config_t mkg_config;
  double icrs_x;
  double icrs_y;
  double icrs_z;
  char star_record[BSR_STAR_RECORD_SIZE];
  size_t star_record_size=(size_t)BSR_STAR_RECORD_SIZE;
  char file_header[BSR_FILE_HEADER_SIZE];
  size_t file_header_size;
  int little_endian;
  int same_endian; // 1 == output endianness is same as arch

  // fields imported from GDR3
  uint64_t source_id;
  double ra;
  double dec;
  double parallax;
  double parallax_over_error;
  int astrometric_params_solved;
  double nu_eff_used_in_astrometry;
  double pseudocolor;
  double phot_G_mean_flux;
  double phot_bp_mean_flux;
  double phot_rp_mean_flux;
  double ecl_lat;
  double teff_gspphot;
  double gspphot_distance;
  double ag_gspphot;

  // temp working vars
  double distance;
  double ra_rad;
  double dec_rad;
  double magnitude;
  double linear_intensity; // intensity relative to vega
  double linear_intensity_undimmed; // intensity relative to vega with extinction removed (according to ag_gspphot)
  uint64_t tmp64;
  uint32_t tmp32;
  uint16_t tmp16;

  //
  // initialize timers
  //
  clock_gettime(CLOCK_REALTIME, &overall_starttime);
  clock_gettime(CLOCK_REALTIME, &starttime);

  //
  // set default options
  //
  setDefaults(&mkg_config);

  //
  // proces command line options
  processCmdArgs(&mkg_config, argc, argv);

  //
  // print version and options
  //
  printf("mkgalaxy version %s\n", BSR_VERSION);
  if (mkg_config.use_bandpass_ratios == 1) {
    printf("Star temperatures determined by r/G and bp/G ratios when available\n");
  } else {
    printf("Star temperatures determined by 'nu_eff_used_in_astrometery' or 'pseudocolor'\n");
  }
  if (mkg_config.use_gspphot_distance == 1) {
    printf("Star distance derived from DR3 'gspphot_distance' when available instead of 'parallax'\n");
  } else {
    printf("Star distance derived from DR3 'parallax' only\n");
  }
  if (mkg_config.calibrate_parallax == 1) {
    printf("Lindegren et. al. parallax calibration enabled\n");
  } else {
    printf("Lindegrenn et. al. parallax calibration disabled\n");
  }
  if (mkg_config.enable_maximum_distance == 1) {
    printf("Maximum distance of %.1e parsecs will be enforced\n", mkg_config.maximum_distance);
  } else {
    printf("Maximum distance enforcement disabled\n");
  }
  if (mkg_config.output_little_endian == 1) {
    printf("Output data files will be in little-endian format\n");
  } else {
    printf("Output data files will be in big-endian format\n");
  }

  //
  // init bandpass ratio tables
  //
  initBandpassRatioTables(rp_over_G_ref, bp_over_G_ref, bp_over_rp_ref);

  //
  // attempt to open input file
  //
  printf("Opening input file gaia-dr3-extracted.csv\n");
  fflush(stdout);
  input_file=fopen("gaia-dr3-extracted.csv", "rb");
  if (input_file == NULL) {
    printf("Error: could not open gaia-dr3-extracted.csv\n");
    fflush(stdout);
    return(1);
  }

  //
  // attempt to open ouptut files
  //

  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq000-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq000-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq000=fopen(file_name, "wb");
  if (output_file_pq000 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq001-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq001-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq001=fopen(file_name, "wb");
  if (output_file_pq001 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq002-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq002-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq002=fopen(file_name, "wb");
  if (output_file_pq002 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq003-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq003-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq003=fopen(file_name, "wb");
  if (output_file_pq003 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq005-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq005-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq005=fopen(file_name, "wb");
  if (output_file_pq005 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq010-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq010-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq010=fopen(file_name, "wb");
  if (output_file_pq010 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq020-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq020-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq020=fopen(file_name, "wb");
  if (output_file_pq020 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq030-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq030-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq030=fopen(file_name, "wb");
  if (output_file_pq030 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq050-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq050-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq050=fopen(file_name, "wb");
  if (output_file_pq050 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-pq100-%s.%s", BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-pq100-%s.%s", BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file_pq100=fopen(file_name, "wb");
  if (output_file_pq100 == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }

  //
  // check endianness
  //
  little_endian=littleEndianTest();
  if ((mkg_config.output_little_endian ^ little_endian) == 0) {
    same_endian=1;
  } else {
    same_endian=0;
  }

  //
  // write file headers
  //
  if (mkg_config.output_little_endian == 1) {
    snprintf(file_header, BSR_FILE_HEADER_SIZE, "%s, mkgalaxy version: %s, use_bandpass_ratios: %d, use_gspphot_distance: %d, calibrate_parallax: %d, enable_maximum_distance: %d, maximum_distance: %.1e\n", BSR_MAGIC_NUMBER_LE, BSR_VERSION, mkg_config.use_bandpass_ratios, mkg_config.use_gspphot_distance, mkg_config.calibrate_parallax, mkg_config.enable_maximum_distance, mkg_config.maximum_distance);
  } else {
    snprintf(file_header, BSR_FILE_HEADER_SIZE, "%s, mkgalaxy version: %s, use_bandpass_ratios: %d, use_gspphot_distance: %d, calibrate_parallax: %d, enable_maximum_distance: %d, maximum_distance: %.1e\n", BSR_MAGIC_NUMBER_BE, BSR_VERSION, mkg_config.use_bandpass_ratios, mkg_config.use_gspphot_distance, mkg_config.calibrate_parallax, mkg_config.enable_maximum_distance, mkg_config.maximum_distance);
  } 
  // pad the rest of file_header with zeros
  file_header_size=strnlen(file_header, (BSR_FILE_HEADER_SIZE - 1));
  for (i=(int)file_header_size; i < BSR_FILE_HEADER_SIZE; i++) {
    file_header[i]=0;
  }
  printf("Writing file headers\n");
  fflush(stdout);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq100);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq050);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq030);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq020);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq010);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq005);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq003);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq002);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq001);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file_pq000);

  //
  // read and process each line of input file
  //
  input_count=0;
  discard_no_flux_count=0;
  discard_parms_count=0;
  discard_parallax_count=0;
  distance_override_negative_count=0;
  distance_override_toohigh_count=0;
  pq000_count=0;
  pq001_count=0; 
  pq002_count=0; 
  pq003_count=0; 
  pq005_count=0; 
  pq010_count=0; 
  pq020_count=0;
  pq030_count=0;
  pq050_count=0;
  pq100_count=0;
  total_output_count=0;
  temperature_from_bp_G_count=0;
  temperature_from_rp_G_count=0;
  temperature_from_bp_rp_count=0;
  temperature_from_nu_eff_count=0;
  temperature_from_pseudocolor_count=0;
  min_temp_count=0;
  max_temp_count=0;
  unreddened_min_temp_count=0;
  unreddened_max_temp_count=0;
  gspphot_distance_count=0;
  undimmed_count=0;
  unreddened_count=0;
  input_line_p=fgets(input_line, 256, input_file);
  printf("Beginning input file processing\n");
  fflush(stdout);
  while (input_line_p != NULL) {

#ifdef DEBUG
    printf("debug, input_line: %s", input_line_p);
    fflush(stdout);
#endif

    //
    // source_id,ra,dec,parallax,parallax_over_error,astrometric_params_solved,nu_eff_used_in_astrometry,pseudocolour,phot_g_mean_flux,phot_bp_mean_flux,phot_rp_mean_flux,ecl_lat,teff_gspphot,distance_gspphot,ag_gspphot
    //
    if (input_line[0] != 'r') { // skip csv header lines
      input_count++;
      field_start=input_line;

      //
      // source_id
      //
      field_end=strchr(field_start, ','); 
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      source_id=strtoull(tmpstr, NULL, 10);
      field_start=(field_end+1);

      //
      // ra
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      ra=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // dec
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      dec=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // parallax
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        parallax=0.0;
      } else {
        parallax=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // parallax_over_error
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        parallax_over_error=0.0;
      } else {
        parallax_over_error=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // astrometric_params_solved
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      astrometric_params_solved=strtol(tmpstr, NULL, 10);
      field_start=(field_end+1);

      //
      // nu_eff_used_in_astrometry
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        nu_eff_used_in_astrometry=0.0;
      } else {
        nu_eff_used_in_astrometry=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // pseudocolor
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        pseudocolor=0.0;
      } else {
        pseudocolor=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // phot_G_mean_flux
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        phot_G_mean_flux=0.0;
      } else {
        phot_G_mean_flux=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // phot_bp_mean_flux
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        phot_bp_mean_flux=0.0;
      } else {
        phot_bp_mean_flux=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // phot_rp_mean_flux
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        phot_rp_mean_flux=0.0;
      } else {
        phot_rp_mean_flux=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // ecl_lat
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      ecl_lat=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // teff_gspphot
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        teff_gspphot=0.0;
      } else {
        teff_gspphot=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // gspphot_distance
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        gspphot_distance=0.0;
      } else {
        gspphot_distance=strtod(tmpstr, NULL);
      }
      field_start=(field_end+1);

      //
      // ag_gspphot
      //
      field_end=strchr(field_start, '\0'); // special processing for last field
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      if (tmpstr[0] == 'n') { // null record
        ag_gspphot=0.0;
      } else {
        ag_gspphot=strtod(tmpstr, NULL);
      }

      //
      // only continue if record has parallax (5 parms solved or 6 parms solved) and phot_G_mean_flux > 0
      //
      if (((astrometric_params_solved == 31) || (astrometric_params_solved == 95)) && (phot_G_mean_flux > 0.0)) {
        //
        // transform flux to linear intensity relative to vega (do this before parallax calibration)
        //
        linear_intensity=phot_G_mean_flux * flux_to_vega;
        magnitude=-2.5 * log10(linear_intensity); // some stars have blank phot_G_mean_magnitude so we derive from the more reliable flux column

        //
        // select correct color wavenumber variable (do this before parallax calibration)
        //
        if (astrometric_params_solved == 31) {
          color_wavenumber=nu_eff_used_in_astrometry;
        } else if (astrometric_params_solved == 95) {
          color_wavenumber=pseudocolor;
        }

        //
        // optionally calibrate parallax according to Lindegren et. al
        //
        if (mkg_config.calibrate_parallax == 1) {
          calibrateParallax(&parallax, astrometric_params_solved, magnitude, color_wavenumber, ecl_lat);
        }

        //
        // only continue if gspphot_distance (if used) or parallax is valid
        //
        if ((mkg_config.enable_maximum_distance == 1) || ((mkg_config.use_gspphot_distance == 1) && (gspphot_distance > 0.0)) || (parallax > 0.0)) {
          //
          // select desired distance source
          //
          if ((mkg_config.use_gspphot_distance == 1) && (gspphot_distance > 0.0)) {
            distance=gspphot_distance;
            gspphot_distance_count++;
          } else if (parallax <= 0.0) {
            // should only get here if enable_maximum_distance == 1 
            distance=mkg_config.maximum_distance;
            distance_override_negative_count++;
          } else {
            //distance=M_PI / (648000.0 * tan(M_PI * parallax / 648000000)); // in parsecs, full calculation
            distance=1000.0 / parallax; // in parsecs, approximation ignoring small angle tan()
          }

          //
          // optionally override maximum distance
          //
          if ((mkg_config.enable_maximum_distance == 1) && (distance > mkg_config.maximum_distance)) {
            distance=mkg_config.maximum_distance;
            distance_override_toohigh_count++;
          }

          //
          // transform spherical icrs to euclidian icrs
          //
          ra_rad=ra * M_PI / 180.0;
          dec_rad=dec * M_PI / 180.0;
          icrs_x=distance * cos(dec_rad) * cos(ra_rad);
          icrs_y=distance * cos(dec_rad) * sin(ra_rad);
          icrs_z=distance * sin(dec_rad);

          //
          // convert linear intensity to intensity at 1pc and undimmed intensity at 1pc (use ag_gspphot available);
          //
          linear_1pc_intensity=(float)(linear_intensity * pow(distance, 2.0));
          if (ag_gspphot > 0.0) {
            linear_intensity_undimmed=pow(100.0, (-(magnitude - ag_gspphot) / 5.0));
            linear_1pc_intensity_undimmed=(float)(linear_intensity_undimmed * pow(distance, 2.0));
            undimmed_count++;
          } else {
            linear_1pc_intensity_undimmed=linear_1pc_intensity;
          }

    //
    // determine star apparent color temperature
    // try to use rp/G and bp/G ratios to determine apparent temperature if this mode is enabled and rp and bp have flux
    //
          //
          // get ratios
          //
          rp_over_G=phot_rp_mean_flux / phot_G_mean_flux;
          bp_over_G=phot_bp_mean_flux / phot_G_mean_flux;
          bp_over_rp=(bp_over_G / rp_over_G);

          //
          // check for invalid parameters outside of allowable ranges
          // these ranges were determined by looking at the values of the given ratios
          // from 500K-32767K and adjusting limits emperically to minimize the
          // number of stars hitting max low / max high temperature (indicating a likely invalid temperature match).
          //
          bp_over_G_invalid=0;
          rp_over_G_invalid=0;
          bp_over_rp_invalid=0;
          all_invalid=0;
          if ((mkg_config.use_bandpass_ratios != 1) || ((phot_bp_mean_flux <= 0.0) && (phot_rp_mean_flux <= 0.0))) {
             all_invalid=1;
          }
          if (phot_bp_mean_flux <= 0.0) {
            bp_over_G_invalid=1;
            bp_over_rp_invalid=1;
          }
          if (phot_rp_mean_flux <= 0.0) {
            rp_over_G_invalid=1;
            bp_over_rp_invalid=1;
          }
          if ((bp_over_rp < 3.06E-6) || (bp_over_rp > 4.02888)) {
            bp_over_rp_invalid=1;
          }
          if ((bp_over_G < 8.16E-6) || (bp_over_G > 0.933898)) {
            bp_over_G_invalid=1;
          }
          if ((rp_over_G < 0.231800) || (rp_over_G > 2.664)) {
            rp_over_G_invalid=1;
          }

          //
          // select best match for apparent temperature
          //
          if ((all_invalid == 0) && ((phot_bp_mean_flux > 0.0) || (phot_rp_mean_flux > 0.0))) {
            //
            // we have at least one parameter to match against Planck spectrum reference
            //
            bestmatch_temperature=0;
            if (bp_over_rp_invalid == 0) {
              // bp over rp generally gives the best results
              for (i=500; ((i < 32768) && (bestmatch_temperature == 0)); i++) {
                if (bp_over_rp_ref[i] > bp_over_rp) {
                  bestmatch_temperature=i;
                }
              }
              if (bestmatch_temperature == 0) {
                bestmatch_temperature=32767;
              }
              temperature_from_bp_rp_count++;
            } else if (bp_over_G_invalid == 0) {
              // next best is bp over G
              for (i=500; ((i < 32768) && (bestmatch_temperature == 0)); i++) {
                if (bp_over_G_ref[i] > bp_over_G) {
                  bestmatch_temperature=i;
                }
              }   
              if (bestmatch_temperature == 0) {
                bestmatch_temperature=32767;
              }
              temperature_from_bp_G_count++;
            } else if (rp_over_G_invalid == 0) {
              // rp over G is least likely to give accurate apparent temp due to background infrared
              for (i=500; ((i < 32768) && (bestmatch_temperature == 0)); i++) {
                if (rp_over_G_ref[i] < rp_over_G) {
                  bestmatch_temperature=i;
                }
              }
              if (bestmatch_temperature == 0) {
                bestmatch_temperature=32767;
              }
              temperature_from_rp_G_count++;
            // below we handle cases where none of the three parameters were within valid ranges but we have data so we set to min/max value
            } else if ((phot_bp_mean_flux > 0.0) && (phot_rp_mean_flux > 0.0)) {
              // bp/rp out of range but has value so set to max
              bestmatch_temperature=32767;
              temperature_from_bp_rp_count++;
            } else if (phot_bp_mean_flux > 0.0) {
              // blue out of range but has value so set to max
              bestmatch_temperature=32767;
              temperature_from_bp_G_count++;
            } else if (phot_rp_mean_flux > 0.0) {
              // red out of range but has value so set to min
              bestmatch_temperature=500;
              temperature_from_rp_G_count++;
            } else {
              // catchall, should never get here but just in case set to default temp
              bestmatch_temperature=4300;
            }
            // range checks 
            if (bestmatch_temperature < 500) {
              bestmatch_temperature=500;
            } else if (bestmatch_temperature > 32767) {
              bestmatch_temperature=32767;
            }
            if (bestmatch_temperature == 500) {
              min_temp_count++;
            } else if (bestmatch_temperature == 32767) {
              max_temp_count++;
            }
            color_temperature=(uint64_t)bestmatch_temperature;
          } else {
            //
            // no parameter with flux value, transofrm color_wavenumber to integer Kelvin blackbody temperature and clip extraneous values
            //
            color_temperature=(uint64_t)((2897.771955 * color_wavenumber) + 0.5); // wein's displacement to convert to Kelvin
            if (color_temperature < 500ul) {
              color_temperature=500ul;
              min_temp_count++;
            } else if (color_temperature > 32767ul) {
              color_temperature=32767ul;
              max_temp_count++;
            }
            if (astrometric_params_solved == 31) {
              temperature_from_nu_eff_count++;
            } else if (astrometric_params_solved == 95) {
              temperature_from_pseudocolor_count++;
            }
          } // end if some parameter is valid

          //
          // set unreddened color temperature from teff_gspphot
          //
          if (teff_gspphot > 0.0) {
            color_temperature_unreddened=(uint64_t)(teff_gspphot + 0.5);
            if (color_temperature_unreddened < 500ul) {
              color_temperature_unreddened=500ul;
              unreddened_min_temp_count++;
            } else if (color_temperature_unreddened > 32767ul) {
              color_temperature_unreddened=32767ul;
              unreddened_max_temp_count++;
            }
            unreddened_count++;
          } else {
            color_temperature_unreddened=color_temperature;
          }

#ifdef DEBUG
          printf("debug, source_id: %ull, parallax_over_error: %.4e, distance: %.4e, icrs_x: %.4e, icrs_y: %.4e, icrs_z: %.4e, linear_1pc_intensity: %.4e, linear_1pc_intensity_undimmed: %.4e, color_temperature: %d, color_temperature_unreddened: %d\n", source_id, parallax_over_error, distance, icrs_x, icrs_y, icrs_z, linear_1pc_intensity, linear_1pc_intensity_undimmed, color_temperature, color_temperature_unreddened);
          fflush(stdout);
#endif

          // Binary data files have a fixed-length 256-bit ascii header (including the file identifier in the first 11 bytes), followed
          // by a variable number of 33-byte star records.
          //
          // Each star record includes a 64-bit unsigned integer for Gaia DR3 'source_id', three 40-bit truncaed doubles for x,y,z,
          // a 24-bit truncated float for linear_1pc_intensity, a 24-bit truncated float for linear_1pc_intensity_undimmed,
          // a 16-bit unsigned int for color_temperature, and a 16-bit unsigned int for color_temperature_unreddened.
          // these are packed into a 33 byte star record with each field encoded in the selected byte order.
          //
          // +---------------+---------+---------+---------+-----+-----+---+---+
          // |   source_id   |    x    |    y    |    z    | li  |li-u | c |c-u|
          // +---------------+---------+---------+---------+-----+-----+---+---+
          // |      8        |    5    |    5    |    5    |  3  |  3  | 2 | 2 |
          //                               bytes
          //
          // mkgalaxy and mkexternal have options to generate either little-endian or big-endian files but the bye-order of the file
          // must match the architecture it is used on with bsrender. This is to avoid unnessary operations in the performance
          // critical inner-loop of processStars() which iterates over potentially billions of star records.
          //

          //
          // pack star record fields into 33-byte star_record
          //
          if (same_endian == 1) {
            //
            // output byte order is same as this arch
            //

            // source_id
            tmp64=0;
            tmp64=source_id;
            star_record[0]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[1]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[2]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[3]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[4]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[5]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[6]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[7]=*(char *)&tmp64;

            // icrs_x
            tmp64=0;
            tmp64=*(uint64_t *)&icrs_x;
            if (little_endian == 1) {
              tmp64 >>= 24; // skip 24 lsb if source is little-endian
            }
            star_record[8]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[9]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[10]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[11]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[12]=*(char *)&tmp64;

            // icrs_y
            tmp64=0;
            tmp64=*(uint64_t *)&icrs_y;
            if (little_endian == 1) {
              tmp64 >>= 24; // skip 24 lsb if source is little-endian
            }
            star_record[13]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[14]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[15]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[16]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[17]=*(char *)&tmp64;

            // icrs_z
            tmp64=0;
            tmp64=*(uint64_t *)&icrs_z;
            if (little_endian == 1) {
              tmp64 >>= 24; // skip 24 lsb if source is little-endian
            }
            star_record[18]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[19]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[20]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[21]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[22]=*(char *)&tmp64;

            // linear_1pc_intensity
            tmp32=0;
            tmp32=*(uint32_t *)&linear_1pc_intensity;
            if (little_endian == 1) {
              tmp32 >>= 8; // skip 8 lsb if source is little-endian
            }
            star_record[23]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[24]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[25]=*(char *)&tmp32;

            // linear_1pc_intensity_undimmed
            tmp32=0;
            tmp32=*(uint32_t *)&linear_1pc_intensity_undimmed;
            if (little_endian == 1) {
              tmp32 >>= 8; // skip 8 lsb if source is little-endian
            }
            star_record[26]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[27]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[28]=*(char *)&tmp32;

            // color_temperature
            tmp16=0;
            tmp16=*(uint16_t *)&color_temperature;
            star_record[29]=*(char *)&tmp16;
            tmp16 >>= 8;
            star_record[30]=*(char *)&tmp16;

            // color_temperature_unreddened
            tmp16=0;
            tmp16=*(uint16_t *)&color_temperature_unreddened;
            star_record[31]=*(char *)&tmp16;
            tmp16 >>= 8;
            star_record[32]=*(char *)&tmp16;
          } else {
            //
            // output byte order is opposite this arch
            //

            // source_id
            tmp64=0;
            tmp64=source_id;
            star_record[7]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[6]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[5]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[4]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[3]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[2]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[1]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[0]=*(char *)&tmp64;

            // icrs_x
            tmp64=0;
            tmp64=*(uint64_t *)&icrs_x;
            if (little_endian == 1) {
              tmp64 >>= 24; // skip 24 lsb if source is little-endian
            }
            star_record[12]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[11]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[10]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[9]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[8]=*(char *)&tmp64;

            // icrs_y
            tmp64=0;
            tmp64=*(uint64_t *)&icrs_y;
            if (little_endian == 1) {
              tmp64 >>= 24; // skip 24 lsb if source is little-endian
            }
            star_record[17]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[16]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[15]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[14]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[13]=*(char *)&tmp64;

            // icrs_z
            tmp64=0;
            tmp64=*(uint64_t *)&icrs_z;
            if (little_endian == 1) {
              tmp64 >>= 24; // skip 24 lsb if source is little-endian
            }
            star_record[22]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[21]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[20]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[19]=*(char *)&tmp64;
            tmp64 >>= 8;
            star_record[18]=*(char *)&tmp64;

            // linear_1pc_intensity
            tmp32=0;
            tmp32=*(uint32_t *)&linear_1pc_intensity;
            if (little_endian == 1) {
              tmp32 >>= 8; // skip 8 lsb if source is little-endian
            }
            star_record[25]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[24]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[23]=*(char *)&tmp32;

            // linear_1pc_intensity_undimmed
            tmp32=0;
            tmp32=*(uint32_t *)&linear_1pc_intensity_undimmed;
            if (little_endian == 1) {
              tmp32 >>= 8; // skip 8 lsb if source is little-endian
            }
            star_record[28]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[27]=*(char *)&tmp32;
            tmp32 >>= 8;
            star_record[26]=*(char *)&tmp32;

            // color_temperature
            tmp16=0;
            tmp16=*(uint16_t *)&color_temperature;
            star_record[30]=*(char *)&tmp16;
            tmp16 >>= 8;
            star_record[29]=*(char *)&tmp16;

            // color_temperature_unreddend
            tmp16=0;
            tmp16=*(uint16_t *)&color_temperature_unreddened;
            star_record[32]=*(char *)&tmp16;
            tmp16 >>= 8;
            star_record[31]=*(char *)&tmp16;
          }

          //
          // output star_record to correct output file
          //
          if (parallax_over_error >= 100.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq100);
            pq100_count++;
            total_output_count++;
          } else if (parallax_over_error >= 50.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq050);
            pq050_count++;
            total_output_count++;
          } else if (parallax_over_error >= 30.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq030);
            pq030_count++;
            total_output_count++;
          } else if (parallax_over_error >= 20.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq020);
            pq020_count++;
            total_output_count++;
          } else if (parallax_over_error >= 10.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq010);
            pq010_count++;
            total_output_count++;
          } else if (parallax_over_error >= 5.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq005);
            pq005_count++;
            total_output_count++;
          } else if (parallax_over_error >= 3.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq003);
            pq003_count++;
            total_output_count++;
          } else if (parallax_over_error >= 2.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq002);
            pq002_count++;
            total_output_count++;
          } else if (parallax_over_error >= 1.0) {
            fwrite(star_record, star_record_size, 1, output_file_pq001);
            pq001_count++;
            total_output_count++;
          } else {
            fwrite(star_record, star_record_size, 1, output_file_pq000);
            pq000_count++;
            total_output_count++;
          }

        } else {
          discard_parallax_count++;
#ifdef DEBUG
          printf("debug, ignoring star due to missing or negative parallax\n");
          fflush(stdout);
#endif
        } // end ignore if zero or negative parallax after corrections
      } else {
        if (phot_G_mean_flux <= 0.0) {
          discard_no_flux_count++;
#ifdef DEBUG
          printf("debug, ignoring star due to no G-band flux\n");
          fflush(stdout);
#endif
        } else {
          discard_parms_count++;
#ifdef DEBUG
          printf("debug, ignoring star due to insufficient astrometric params\n");
          fflush(stdout);
#endif
        }
      } // end ignore if not parms=5 or 6 or no G flux
#ifdef DEBUG
    } else {
      printf("debug, ignoring csv header line\n");
      fflush(stdout);
#endif
    } // end ignore csv header lines

    //
    // periodic status
    //
    if ((input_count > 0) && ((input_count % 1000000) == 0)) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      overall_elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(overall_starttime.tv_sec - 1500000000) + ((double)overall_starttime.tv_nsec) / 1.0E9);
      printf("------\nInput records: %9lu\n", input_count);
      printf("\nOutput by parallax quality\n");
      printf("  pq000: %lu\n", pq000_count);
      printf("  pq001: %lu\n", pq001_count);
      printf("  pq002: %lu\n", pq002_count);
      printf("  pq003: %lu\n", pq003_count);
      printf("  pq005: %lu\n", pq005_count);
      printf("  pq010: %lu\n", pq010_count);
      printf("  pq020: %lu\n", pq020_count);
      printf("  pq030: %lu\n", pq030_count);
      printf("  pq050: %lu\n", pq050_count);
      printf("  pq100: %lu\n", pq100_count);
      printf("  Total: %lu\n", total_output_count);
      printf("\nDiscards\n");
      printf("  2-parameter solution (no parallax): %lu\n", discard_parms_count);
      printf("  5 or 6-parameter solution but no G-band flux: %lu\n", discard_no_flux_count);
      if (mkg_config.enable_maximum_distance == 1) {
        printf("  Negative parallax: (maximum distance override enabled)\n");
      } else {
        printf("  Negative parallax: %lu\n", discard_parallax_count);
      }
      printf("\nValues derived from Gaia DR3 GSPPhot fields\n");
      if (mkg_config.use_gspphot_distance == 0) {
        printf("  Distance: (disabled)\n");
      } else {
        printf("  Distance: %lu\n", gspphot_distance_count);
      }
      printf("  Undimmed intensity: %lu\n", undimmed_count);
      printf("  Unreddened temperature: %lu\n", unreddened_count);
      if (mkg_config.enable_maximum_distance == 1) {
        printf("\nDistance override (max=%.1e parsecs)\n", mkg_config.maximum_distance);
        printf("  Negative parallax: %lu\n", distance_override_negative_count);
        printf("  Distance too high: %lu\n", distance_override_toohigh_count);
      }
      printf("\nApparent temperature derived from\n");
      printf("  bp/rp: %lu\n", temperature_from_bp_rp_count);
      printf("  bp/G: %lu\n", temperature_from_bp_G_count);
      printf("  rp/G: %lu\n", temperature_from_rp_G_count);
      printf("  nu_eff_used_in_astrometry: %lu\n", temperature_from_nu_eff_count);
      printf("  pseudocolor: %lu\n", temperature_from_pseudocolor_count);
      printf("\nTemperature min/max override\n");
      printf("  Apparent temperature < 500K: %lu\n", min_temp_count);
      printf("  Apparent temperature > 32767K: %lu\n", max_temp_count);
      printf("  Unreddened temperature < 500K: %lu\n", unreddened_min_temp_count);
      printf("  Unreddened temperature > 32767K: %lu\n", unreddened_max_temp_count);
      printf("\nIncremental time: %.3fs, total time: %.3fs\n", elapsed_time, overall_elapsed_time);
      fflush(stdout);
      clock_gettime(CLOCK_REALTIME, &starttime);
    }

    input_line_p=fgets(input_line, 256, input_file);
  }  // end while read file
      
  //
  // print final status 
  //
  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
  overall_elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(overall_starttime.tv_sec - 1500000000) + ((double)overall_starttime.tv_nsec) / 1.0E9);
  printf("------\nInput records: %9lu\n", input_count);
  printf("\nOutput by parallax quality\n");
  printf("  pq000: %lu\n", pq000_count);
  printf("  pq001: %lu\n", pq001_count);
  printf("  pq002: %lu\n", pq002_count);
  printf("  pq003: %lu\n", pq003_count);
  printf("  pq005: %lu\n", pq005_count);
  printf("  pq010: %lu\n", pq010_count);
  printf("  pq020: %lu\n", pq020_count);
  printf("  pq030: %lu\n", pq030_count);
  printf("  pq050: %lu\n", pq050_count);
  printf("  pq100: %lu\n", pq100_count);
  printf("  Total: %lu\n", total_output_count);
  printf("\nDiscards\n");
  printf("  2-parameter solution (no parallax): %lu\n", discard_parms_count);
  printf("  5 or 6-parameter solution but no G-band flux: %lu\n", discard_no_flux_count);
  if (mkg_config.enable_maximum_distance == 1) {
    printf("  Negative parallax: (maximum distance override enabled)\n");
  } else {
    printf("  Negative parallax: %lu\n", discard_parallax_count);
  }
  printf("\nValues derived from Gaia DR3 GSPPhot fields\n");
  if (mkg_config.use_gspphot_distance == 0) {
    printf("  Distance: (disabled)\n");
  } else {
    printf("  Distance: %lu\n", gspphot_distance_count);
  }
  printf("  Undimmed intensity: %lu\n", undimmed_count);
  printf("  Unreddened temperature: %lu\n", unreddened_count);
  if (mkg_config.enable_maximum_distance == 1) {
    printf("\nDistance override (max=%.1e parsecs)\n", mkg_config.maximum_distance);
    printf("  Negative parallax: %lu\n", distance_override_negative_count);
    printf("  Distance too high: %lu\n", distance_override_toohigh_count);
  }
  printf("\nApparent temperature derived from\n");
  printf("  bp/rp: %lu\n", temperature_from_bp_rp_count);
  printf("  bp/G: %lu\n", temperature_from_bp_G_count);
  printf("  rp/G: %lu\n", temperature_from_rp_G_count);
  printf("  nu_eff_used_in_astrometry: %lu\n", temperature_from_nu_eff_count);
  printf("  pseudocolor: %lu\n", temperature_from_pseudocolor_count);
  printf("\nTemperature min/max override\n");
  printf("  Apparent temperature < 500K: %lu\n", min_temp_count);
  printf("  Apparent temperature > 32767K: %lu\n", max_temp_count);
  printf("  Unreddened temperature < 500K: %lu\n", unreddened_min_temp_count);
  printf("  Unreddened temperature > 32767K: %lu\n", unreddened_max_temp_count);
  printf("\nIncremental time: %.3fs, total time: %.3fs\n", elapsed_time, overall_elapsed_time);
  fflush(stdout);

  //
  // clean up
  //
  fclose(input_file);
  fclose(output_file_pq000);
  fclose(output_file_pq001);
  fclose(output_file_pq002);
  fclose(output_file_pq003);
  fclose(output_file_pq005);
  fclose(output_file_pq010);
  fclose(output_file_pq020);
  fclose(output_file_pq030);
  fclose(output_file_pq050);
  fclose(output_file_pq100);

  return(0);
}

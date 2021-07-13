//
// Billion Star 3D Rendering Engine
// Kevin M. Loch
//
// 3D rendering engine for the ESA Gaia EDR3 star dataset

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

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <math.h>
#include "Gaia-passbands.h"

int initBandpassRatioTables(double *rp_over_G_ref, double *bp_over_G_ref, double *bp_over_rp_ref) {
  //
  // generate Grp/G, and Ggp/G ratios for a given blackbody temperature
  // This is used in mkgalaxy to determine star temperature by finding best fit to rp/G and bp/G ratios
  //
  int i;
  const double kb=1.380649E-23;
  const double h=6.62607015E-34;
  const double c=299792458.0;
  double temp;
  int wavelength_start;
  int wavelength_end;
  int wavelength;
  double specific_intensity;
  double G_intensity;
  double rp_intensity;
  double bp_intensity;

  //
  // Range of Gaia EDR3 passband data file is 320-1100nm
  //
  wavelength_start=320;
  wavelength_end=1100;

  //
  // calculate rp, g, bp values for each integer Kelvin temp from 0 - 32767K
  //
  for (i=0; i < 32768; i++) {
    temp=(double)i;

    //
    // scan over Gaia passband range and assign intensity chunks to appropriate channels
    //
    rp_intensity=0.0;
    G_intensity=0.0;
    bp_intensity=0.0;
    for (wavelength=wavelength_start; wavelength <= wavelength_end; wavelength++) {
      // omit unnecessary constants since we normalize relative to G-band after integrating
      specific_intensity=1.0 / (pow((wavelength * 1.0E-9), 5.0) * (exp(h * c  / (wavelength * 1.0E-9 * kb * temp)) - 1));
      rp_intensity+=(getGaiaTransmissivityRp(wavelength) * specific_intensity);
      G_intensity+=(getGaiaTransmissivityG(wavelength) * specific_intensity);
      bp_intensity+=(getGaiaTransmissivityBp(wavelength) * specific_intensity);
    } // end for wavelength
      
    //
    // normalize intensity values by setting G-band=1.0;
    //
    if (G_intensity != 0.0) {
      rp_intensity=(rp_intensity / G_intensity);
      bp_intensity=(bp_intensity / G_intensity);
    } else {
      rp_intensity=0.0;
      bp_intensity=0.0;
    }

    //
    // store in bandpass ratio arrays
    //
    rp_over_G_ref[i]=rp_intensity;
    bp_over_G_ref[i]=bp_intensity;
    bp_over_rp_ref[i]=bp_intensity / rp_intensity; // G cancells out

/*
printf("temp: %d, bp/rp: %.6e, bp/G: %.6e, rp/G: %.6e\n", i, bp_over_rp_ref[i], bp_over_G_ref[i], rp_over_G_ref[i]);
fflush(stdout);
*/

  } // end for i

  return(0);
}

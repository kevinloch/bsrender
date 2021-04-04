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

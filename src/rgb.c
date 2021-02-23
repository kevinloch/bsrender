#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <math.h>

int initRGBTables(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
//
// generate RGB color values for a given blackbody temperature
//

  int i;
  const double kb=1.380649E-23;
  const double h=6.62607015E-34;
  const double c=299792458.0;
  double temp;
  double red_wb_factor;
  double green_wb_factor;
  double blue_wb_factor;
  double color_max;
  double color_min;
  double color_mid;

  const double Gaia_Gband_long_limit=1050.0;
  const double Gaia_Gband_short_limit=330.0;
  double wavelength_start;
  double wavelength_end;
  double wavelength;
  double specific_intensity;
  double Gband_intensity;
  double red_intensity;
  double green_intensity;
  double blue_intensity;
  const int wavelength_increments=200;
  double wavelength_increment;

  //
  // determine wavelength scan range and increment
  //
  wavelength_start=Gaia_Gband_long_limit;
  if (bsr_config->red_filter_long_limit > wavelength_start) {
    wavelength_start=bsr_config->red_filter_long_limit;
  }
  if (bsr_config->green_filter_long_limit > wavelength_start) {
    wavelength_start=bsr_config->green_filter_long_limit;
  }
  if (bsr_config->blue_filter_long_limit > wavelength_start) {
    wavelength_start=bsr_config->blue_filter_long_limit;
  }
  wavelength_end=Gaia_Gband_short_limit;
  if (bsr_config->red_filter_short_limit < wavelength_end) {
    wavelength_end=bsr_config->red_filter_short_limit;
  }
  if (bsr_config->green_filter_short_limit < wavelength_end) {
    wavelength_end=bsr_config->green_filter_short_limit;
  }
  if (bsr_config->blue_filter_short_limit < wavelength_end) {
    wavelength_end=bsr_config->blue_filter_short_limit;
  }
  wavelength_increment = (wavelength_start - wavelength_end) / (double)wavelength_increments;

  if (bsr_config->cgi_mode == 0) {
    printf(" wavelength_start: %.3e, wavelength_end: %.3e, wavelength_increment: %.3e ", wavelength_start, wavelength_end, wavelength_increment);
    fflush(stdout);
  }

  //
  // calculate white balance factors
  //
  if (bsr_config->camera_wb_enable == 1) {
    temp=bsr_config->camera_wb_temp;

    //
    // scan over wavelength range and assign intensity chunks to appropriate channels
    //
    Gband_intensity=0.0;
    red_intensity=0.0;
    green_intensity=0.0;
    blue_intensity=0.0;
    for (wavelength=wavelength_start; wavelength >= wavelength_end; wavelength-=wavelength_increment) {
      // omit unnecessary constants since we normalize relative to Gband after integrating
      specific_intensity=1.0 / (pow((wavelength * 1.0E-9), 5.0) * (exp(h * c  / (wavelength * 1.0E-9 * kb * temp)) - 1));
      if ((wavelength <= Gaia_Gband_long_limit) && (wavelength >= Gaia_Gband_short_limit)) {
        Gband_intensity+=specific_intensity;
      }
      if ((wavelength <= bsr_config->red_filter_long_limit) && (wavelength >= bsr_config->red_filter_short_limit)) {
        red_intensity+=specific_intensity;
      }
      if ((wavelength <= bsr_config->green_filter_long_limit) && (wavelength >= bsr_config->green_filter_short_limit)) {
        green_intensity+=specific_intensity;
      }
      if ((wavelength <= bsr_config->blue_filter_long_limit) && (wavelength >= bsr_config->blue_filter_short_limit)) {
        blue_intensity+=specific_intensity;
      }
    }

    //
    // normalize intensity values by comparing to G-band intensity because rgb values are mltiplied by G-band flux during rendering
    //
    if (Gband_intensity != 0.0) {
      red_intensity = red_intensity / Gband_intensity;
      green_intensity = green_intensity / Gband_intensity;
      blue_intensity = blue_intensity / Gband_intensity;
    }
    red_wb_factor=green_intensity / red_intensity;
    green_wb_factor=1.0;
    blue_wb_factor=green_intensity / blue_intensity;
  } else {
    // wb disabled
    red_wb_factor=1.0;
    green_wb_factor=1.0;
    blue_wb_factor=1.0;
  } // end if wb enabled

  //
  // calculate rgb values for each integer Kelvin temp from 0 - 32767K
  //
  for (i=0; i < 32768; i++) {
    temp=(double)i;

    //
    // scan over wavelength range and assign intensity chunks to appropriate channels
    //
    Gband_intensity=0.0;
    red_intensity=0.0;
    green_intensity=0.0;
    blue_intensity=0.0;
    for (wavelength=wavelength_start; wavelength >= wavelength_end; wavelength-=wavelength_increment) {
      // omit unnecessary constants since we normalize relative to Gband after integrating
      specific_intensity=1.0 / (pow((wavelength * 1.0E-9), 5.0) * (exp(h * c  / (wavelength * 1.0E-9 * kb * temp)) - 1));
      if ((wavelength <= Gaia_Gband_long_limit) && (wavelength >= Gaia_Gband_short_limit)) {
        Gband_intensity+=specific_intensity;
      }
      if ((wavelength <= bsr_config->red_filter_long_limit) && (wavelength >= bsr_config->red_filter_short_limit)) {
        red_intensity+=specific_intensity;
      }
      if ((wavelength <= bsr_config->green_filter_long_limit) && (wavelength >= bsr_config->green_filter_short_limit)) {
        green_intensity+=specific_intensity;
      }
      if ((wavelength <= bsr_config->blue_filter_long_limit) && (wavelength >= bsr_config->blue_filter_short_limit)) {
        blue_intensity+=specific_intensity;
      }
    }
      
    //
    // normalize intensity values by comparing to G-band intensity because rgb values are mltiplied by G-band flux during rendering
    //
    if (Gband_intensity != 0.0) {
      red_intensity = red_wb_factor * red_intensity / Gband_intensity;
      green_intensity = green_wb_factor * green_intensity / Gband_intensity;
      blue_intensity = blue_wb_factor * blue_intensity / Gband_intensity;
    }

    //
    // apply camera color saturation adjustment
    //
    color_max=-1.0;
    if (red_intensity > color_max) {
      color_max=red_intensity;
    } 
    if (green_intensity > color_max) {
      color_max=green_intensity;
    } 
    if (blue_intensity > color_max) {
      color_max=blue_intensity;
    } 
    color_min=1.0E99;
    if (red_intensity < color_min) {
      color_min=red_intensity;
    }
    if (green_intensity < color_min) {
      color_min=green_intensity;
    }
    if (blue_intensity < color_min) {
      color_min=blue_intensity;
    } 
    color_mid=(color_max + color_min) / 2.0;
    red_intensity=color_mid + (bsr_config->camera_color_saturation * (red_intensity - color_mid));
    if (red_intensity < 0) {
      red_intensity=0;
    }
    green_intensity=color_mid + (bsr_config->camera_color_saturation * (green_intensity - color_mid));
    if (green_intensity < 0) {
      green_intensity=0;
    }
    blue_intensity=color_mid + (bsr_config->camera_color_saturation * (blue_intensity - color_mid));
    if (blue_intensity < 0) {
      blue_intensity=0;
    }

    //
    // store in rgb arrays
    //
    bsr_state->rgb_red[i]=red_intensity;
    bsr_state->rgb_green[i]=green_intensity;
    bsr_state->rgb_blue[i]=blue_intensity;
  } // end for i

  return(0);
}

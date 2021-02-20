#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <math.h>

int initRGBTables(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
//
// generate RGB color values for a given blackbody temperature
//

  int i;
  double red_freq;
  double green_freq;
  double blue_freq;
  double red_intensity;
  double green_intensity;
  double blue_intensity;
  const double kb=1.380649E-23;
  const double h=6.62607015E-34;
  const double c=299792458.0;
  double temp;
  double red_wb_factor;
  double green_wb_factor;
  double blue_wb_factor;
  double normalization_factor;
  double color_max;
  double color_min;
  double color_mid;

  //
  // convert reference wavelengths to frequency
  //
  red_freq=c / bsr_config->red_center;
  green_freq=c / bsr_config->green_center;
  blue_freq=c / bsr_config->blue_center;

  //
  // calculate white balance factors
  //
  red_intensity=  (2.0 * h * pow(red_freq,   3.0) / pow(c, 2.0)) / (exp(h * red_freq   / (kb * bsr_config->camera_wb_temp)) - 1);
  green_intensity=(2.0 * h * pow(green_freq, 3.0) / pow(c, 2.0)) / (exp(h * green_freq / (kb * bsr_config->camera_wb_temp)) - 1);
  blue_intensity= (2.0 * h * pow(blue_freq,  3.0) / pow(c, 2.0)) / (exp(h * blue_freq  / (kb * bsr_config->camera_wb_temp)) - 1);
  green_wb_factor=1.0 / green_intensity;
  red_wb_factor=1.0 / red_intensity;
  blue_wb_factor=1.0 / blue_intensity;
  
  //
  // calculate rgb values for each integer Kelving temp from 0 - 32767K
  //
  for (i=0; i < 32768; i++) {
    temp=(double)i;
    red_intensity=red_wb_factor     * (2.0 * h * pow(red_freq,   3.0) / pow(c, 2.0)) / (exp(h * red_freq   / (kb * temp)) - 1);
    green_intensity=green_wb_factor * (2.0 * h * pow(green_freq, 3.0) / pow(c, 2.0)) / (exp(h * green_freq / (kb * temp)) - 1);
    blue_intensity= blue_wb_factor  * (2.0 * h * pow(blue_freq,  3.0) / pow(c, 2.0)) / (exp(h * blue_freq  / (kb * temp)) - 1);

    //
    // calculate normalization factor so r+g+b=1
    //
    if (red_intensity == 0) {
      red_intensity=1.0;
    }
    normalization_factor=1.0 / (red_intensity + green_intensity + blue_intensity);
    red_intensity=red_intensity * normalization_factor;
    green_intensity=green_intensity * normalization_factor;
    blue_intensity=blue_intensity * normalization_factor;

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
    color_min=2.0;
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
    // re-normalize and store in rgb arrays
    //
    normalization_factor=1.0 / (red_intensity + green_intensity + blue_intensity);
    bsr_state->rgb_red[i]=red_intensity * normalization_factor;
    bsr_state->rgb_green[i]=green_intensity * normalization_factor;
    bsr_state->rgb_blue[i]=blue_intensity * normalization_factor;
  } // end for i

  return(0);
}

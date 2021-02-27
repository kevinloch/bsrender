#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <math.h>
#include "Bessel.h"

double makeAiryMap(double *Airymap, int max_r, int half_oversampling, double pixel_scaling_factor, double I0) {
  double *Airymap_p;
  int max_xy;
  int pixel_x_index;
  int pixel_y_index;
  double pixel_x;
  double pixel_y;
  double pixel_r;
  int oversampling;
  int oversample_x_index;
  int oversample_y_index;
  double oversample_x;
  double oversample_y;
  double oversample_r;
  double Bessel_x;
  int Bessel_x_index;
  double Airymap_sum;

  //
  // set shorthand variables
  //
  max_xy = (max_r * 2) + 1;
  oversampling=(half_oversampling * 2) +1;

  //
  // generate Airy disk map
  //
  Airymap_p=Airymap;
  Airymap_sum=0.0;
  for (pixel_y_index=0; pixel_y_index < max_xy; pixel_y_index++) {
    for (pixel_x_index=0; pixel_x_index < max_xy; pixel_x_index++) {
      pixel_x=(double)(pixel_x_index - max_r);
      pixel_y=(double)(pixel_y_index - max_r);
      pixel_r=sqrt((pixel_x * pixel_x) + (pixel_y * pixel_y));
      if ((pixel_r <= (double)max_r) && ((pixel_r * pixel_scaling_factor) <= 2370)) {
        *Airymap_p=0.0;
        for (oversample_y_index=0; oversample_y_index < oversampling; oversample_y_index++) {
          for (oversample_x_index=0; oversample_x_index < oversampling; oversample_x_index++) {
            oversample_x=(pixel_x + (((double)oversample_x_index - (double)half_oversampling) / (double)oversampling));
            oversample_y=(pixel_y + (((double)oversample_y_index - (double)half_oversampling) / (double)oversampling));
            oversample_r=sqrt((oversample_x * oversample_x) + (oversample_y * oversample_y));
            Bessel_x=oversample_r * pixel_scaling_factor;
            Bessel_x_index=(int)((Bessel_x * 10) + 0.5);
            if ((oversample_r == 0) || (Bessel_x_index == 0)) {
              *Airymap_p+=I0;
            } else if (Bessel_x_index > 23710) {
              *Airymap_p=0.0; // ignore if beyond range of Bessel.h (too many orders of diffraction)
              oversample_x_index=(oversampling + 1);
              oversample_y_index=(oversampling + 1);
            } else {
              *Airymap_p+=(I0 * pow(2.0 * Bessel_J1[Bessel_x_index] / ((double)Bessel_x_index / 10.0), 2.0));
            }
          } // end for oversample_x
        } // end for oversample_y
      } else {
        *Airymap_p=0.0; // ignore if outside max radius
      } // end if pixel_r < max_r
      if (*Airymap_p > 0.0) {
        Airymap_sum+=*Airymap_p;
      }
      Airymap_p++;
    } // end for x
  } // end for y

  return(Airymap_sum);
}

int initAiryMaps(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  double pixel_scaling_factor_red;
  double pixel_scaling_factor_green;
  double pixel_scaling_factor_blue;
  int half_oversampling_red;
  int half_oversampling_green;
  int half_oversampling_blue;
  int oversampling_red;
  int oversampling_green;
  int oversampling_blue;
  double I0_red;
  double I0_green;
  double I0_blue;
  double Airymap_sum;
  double red_center;
  double green_center;
  double blue_center;

  //
  // calculate center wavelengths for each color channel
  //
  red_center=bsr_config->red_filter_short_limit + ((bsr_config->red_filter_long_limit - bsr_config->red_filter_short_limit) / 2.0);
  green_center=bsr_config->green_filter_short_limit + ((bsr_config->green_filter_long_limit - bsr_config->green_filter_short_limit) / 2.0);
  blue_center=bsr_config->blue_filter_short_limit + ((bsr_config->blue_filter_long_limit - bsr_config->blue_filter_short_limit) / 2.0);

  //
  // calculate Airy disk scaling factors (pixels) for each color.  Green is defined in config with Airy_disk_first_null
  //
  pixel_scaling_factor_green=3.8317 / (double)bsr_config->Airy_disk_first_null;
  pixel_scaling_factor_red=pixel_scaling_factor_green * green_center / red_center;
  pixel_scaling_factor_blue=pixel_scaling_factor_green * green_center / blue_center;
  
  // calculate pixel oversampling for each color to make full use of 10x Bessel function resolution from Bessel.h
  // and minimum of 11x11
  //
  half_oversampling_red=(int)((pixel_scaling_factor_red * 10.0) + 0.5);
  if (half_oversampling_red < 5) {
    half_oversampling_red=5;
  }
  oversampling_red=(half_oversampling_red * 2) +1;
  half_oversampling_green=(int)((pixel_scaling_factor_green * 10.0) + 0.5);
  if (half_oversampling_green < 5) {
    half_oversampling_green=5;
  }
  oversampling_green=(half_oversampling_green * 2) +1;
  half_oversampling_blue=(int)((pixel_scaling_factor_blue * 10.0) + 0.5);
  if (half_oversampling_blue < 5) {
    half_oversampling_blue=5;
  }
  oversampling_blue=(half_oversampling_blue * 2) +1;

  //
  // calculate center intensity calibration for each color
  //
  I0_green=1.16823 / pow((double)(bsr_config->Airy_disk_first_null * oversampling_green), 2.0);
  I0_red=1.16823 * pow(green_center, 2.0) / (pow(red_center, 2.0) * pow((double)(bsr_config->Airy_disk_first_null * oversampling_red), 2.0));
  I0_blue=1.16823 * pow(green_center, 2.0) / (pow(blue_center, 2.0)  *pow((double)(bsr_config->Airy_disk_first_null * oversampling_blue), 2.0));

  //
  // generate Airy disk map for each color
  //
  Airymap_sum=makeAiryMap(bsr_state->Airymap_red, bsr_config->Airy_disk_max_extent, half_oversampling_red, pixel_scaling_factor_red, I0_red);
  if (bsr_config->cgi_mode != 1) {
    printf(" Airymap_sum_red: %.6e,", Airymap_sum);
    fflush(stdout);
  }
  Airymap_sum=makeAiryMap(bsr_state->Airymap_green, bsr_config->Airy_disk_max_extent, half_oversampling_green, pixel_scaling_factor_green, I0_green);
  if (bsr_config->cgi_mode != 1) {
    printf(" Airymap_sum_green: %.6e,", Airymap_sum);
    fflush(stdout);
  }
  Airymap_sum=makeAiryMap(bsr_state->Airymap_blue, bsr_config->Airy_disk_max_extent, half_oversampling_blue, pixel_scaling_factor_blue, I0_blue);
  if (bsr_config->cgi_mode != 1) {
    printf(" Airymap_sum_blue: %.6e", Airymap_sum);
    fflush(stdout);
  }

  return 0;
}

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "util.h"
#include "Lanczos.h"
#include "Gaussian-blur.h"

int postProcess(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  int i;
  int current_image_x;
  int current_image_y;
  pixel_composition_t *current_image_p;
  double inv_camera_pixel_limit;
  double pixel_r;
  double pixel_g;
  double pixel_b;

  //
  // shortcuts
  //
  inv_camera_pixel_limit = 1.0 / bsr_config->camera_pixel_limit;

  //
  // get current image pointer
  //
  current_image_p=bsr_state->current_image_buf;

  //
  // display status message if not in cgi mode
  //
  if (bsr_config->cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Applying camera gamma and intensity limit...");
    fflush(stdout);
  }

  //
  // apply cmaera_gamma and intensity limiting
  //
  current_image_x=0;
  current_image_y=0;
  for (i=0; i < (bsr_state->current_image_res_x * bsr_state->current_image_res_y); i++) {

    //
    // set new row pointer if we have reached end of row
    //
    if (current_image_x == bsr_state->current_image_res_x) {
      current_image_x=0;
      current_image_y++;
    }

    //
    // convert pixel values to output range ~0-1.0 with camera sensitivity reference level = 1.0
    //
    pixel_r=current_image_p->r * inv_camera_pixel_limit;
    pixel_g=current_image_p->g * inv_camera_pixel_limit;
    pixel_b=current_image_p->b * inv_camera_pixel_limit;

    //
    // optionally apply camera gamma setting
    //
    if (bsr_config->camera_gamma != 1.0) { // this is expensive so only if not 1.0
      pixel_r=pow(pixel_r, bsr_config->camera_gamma);
      pixel_g=pow(pixel_g, bsr_config->camera_gamma);
      pixel_b=pow(pixel_b, bsr_config->camera_gamma);
    }

    //
    // limit pixel intensity to range [0..1]
    //
    if (bsr_config->camera_pixel_limit_mode == 0) {
      limitIntensity(&pixel_r, &pixel_g, &pixel_b);
    } else if (bsr_config->camera_pixel_limit_mode == 1) {
      limitIntensityPreserveColor(&pixel_r, &pixel_g, &pixel_b);
    }

    //
    // copy back to current image buf
    //
    current_image_p->r=pixel_r;
    current_image_p->g=pixel_g;
    current_image_p->b=pixel_b;

    current_image_x++;
    current_image_p++;
  } // end for i

  //
  // output execution time if not in cgi mode
  //
  if (bsr_config->cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.4fs)\n", elapsed_time);
    fflush(stdout);
  }

  //
  // optionally blur image
  //
  if (bsr_config->Gaussian_blur_radius > 0) {
    GaussianBlur(bsr_config, bsr_state); 
  }

  //
  // optionally resize image
  //
  if (bsr_config->output_scaling_factor != 1.0) {
    resizeLanczos(bsr_config, bsr_state);
  }

  return(0);
}

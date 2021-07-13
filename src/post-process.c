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
#include <time.h>
#include "util.h"
#include "Lanczos.h"
#include "Gaussian-blur.h"
#include "overlay.h"

int postProcess(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  long long image_offset;
  int current_image_x;
  int current_image_y;
  pixel_composition_t *current_image_p;
  double inv_camera_pixel_limit;
  double pixel_r;
  double pixel_g;
  double pixel_b;

  //
  // main thread: camera_gamma, intensity limiting
  //
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
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
    for (image_offset=0; image_offset < ((long long)bsr_state->current_image_res_x * (long long)bsr_state->current_image_res_y); image_offset++) {

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
  } // end if main thread

  //
  // all threads: optionally blur image
  //
  if (bsr_config->Gaussian_blur_radius > 0.0) {
    GaussianBlur(bsr_config, bsr_state); 
  }

  //
  // all threads: optionally resize image
  //
  if (bsr_config->output_scaling_factor != 1.0) {
    resizeLanczos(bsr_config, bsr_state);
  }

  //
  // main thread: optionally draw overlays
  //
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    if (bsr_config->draw_crosshairs == 1) {
      drawCrossHairs(bsr_config, bsr_state);
    }
    if (bsr_config->draw_grid_lines == 1) {
      drawGridLines(bsr_config, bsr_state);
    }
  } // end if main thread

  return(0);
}

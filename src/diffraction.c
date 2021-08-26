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
#include "Bessel.h"
#include "util.h"

int makeAiryMap(bsr_state_t *bsr_state, double *Airymap, int max_r, int half_oversampling, float pixel_scaling_factor, float I0, float obs_ratio) {
  double *Airymap_p;
  int max_xy;
  int map_offset;
  int map_index_x;
  int map_index_y;
  float pixel_x;
  float pixel_y;
  float pixel_r;
  int oversampling;
  int oversample_x_index;
  int oversample_y_index;
  float oversample_x;
  float oversample_y;
  float oversample_r;
  float Bessel_x;
  int Bessel_x_index;
  int lines_per_thread;
  float obs_I0_factor=0.0f;
  int obs_x_index=0;

  //
  // set shorthand variables
  //
  max_xy=max_r + 1;
  oversampling=(half_oversampling * 2) +1;
  lines_per_thread=(int)ceilf(((float)max_xy / ((float)bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }
  if (obs_ratio > 0.0) {
    obs_I0_factor=1.0f / powf((1.0f - (obs_ratio * obs_ratio)), 2.0f);
  }

  //
  // generate Airy disk map
  //
  map_index_x=0;
  map_index_y=(bsr_state->perthread->my_thread_id * lines_per_thread);
  Airymap_p=Airymap + (max_xy * map_index_y);
  for (map_offset=0; ((map_offset < (max_xy * lines_per_thread)) && (map_index_y < max_xy)); map_offset++) {
    pixel_x=(float)map_index_x;
    pixel_y=(float)map_index_y;
    pixel_r=sqrtf((pixel_x * pixel_x) + (pixel_y * pixel_y));
    if ((pixel_r <= (float)max_r) && ((pixel_r * pixel_scaling_factor) < 12800)) {
      *Airymap_p=0.0;
      for (oversample_y_index=0; oversample_y_index < oversampling; oversample_y_index++) {
        for (oversample_x_index=0; oversample_x_index < oversampling; oversample_x_index++) {
          oversample_x=(pixel_x + (((float)oversample_x_index - (float)half_oversampling) / (float)oversampling));
          oversample_y=(pixel_y + (((float)oversample_y_index - (float)half_oversampling) / (float)oversampling));
          oversample_r=sqrtf((oversample_x * oversample_x) + (oversample_y * oversample_y));
          Bessel_x=oversample_r * pixel_scaling_factor;
          Bessel_x_index=(int)((Bessel_x * 10.0f) + 0.5);
          if ((oversample_r == 0.0f) || (Bessel_x_index == 0)) {
            *Airymap_p+=(double)I0;
          } else if (Bessel_x_index >= 128000) {
            *Airymap_p=0.0; // ignore if beyond range of Bessel.h (too many orders of diffraction)
            // skip to next map pixel
            oversample_x_index=(oversampling + 1);
            oversample_y_index=(oversampling + 1);
          } else if (obs_ratio > 0.0) {
            obs_x_index=(int)((obs_ratio * Bessel_x * 10.0f) + 0.5);
            *Airymap_p+=(double)(I0 * obs_I0_factor * powf((((2.0f * Bessel_J1[Bessel_x_index]) - (2.0f * obs_ratio * Bessel_J1[obs_x_index])) / Bessel_x), 2.0f));
          } else {
            *Airymap_p+=(double)(I0 * powf((2.0f * Bessel_J1[Bessel_x_index] / Bessel_x), 2.0f));
          }
        } // end for oversample_x
      } // end for oversample_y
    } else {
      *Airymap_p=0.0; // ignore if outside max radius
    } // end if pixel_r < max_r

    map_index_x++;
    if (map_index_x == max_xy) {
      map_index_x=0;
      map_index_y++;
    }
    Airymap_p++;
  } // end for pixel_index

  return(0);
}

int initAiryMaps(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  float pixel_scaling_factor_red;
  float pixel_scaling_factor_green;
  float pixel_scaling_factor_blue;
  int half_oversampling_red;
  int half_oversampling_green;
  int half_oversampling_blue;
  int oversampling_red;
  int oversampling_green;
  int oversampling_blue;
  float I0_red;
  float I0_green;
  float I0_blue;
  float red_center;
  float green_center;
  float blue_center;
  int i;
  float obs_ratio;

  //
  // main thread: display status message if not in cgi mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing Airy disk maps...");
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_AIRY_MAP_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_AIRY_MAP_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: calculate center wavelengths for each color channel
  //
  red_center=(float)bsr_config->red_filter_short_limit + (((float)bsr_config->red_filter_long_limit - (float)bsr_config->red_filter_short_limit) / 2.0f);
  green_center=(float)bsr_config->green_filter_short_limit + (((float)bsr_config->green_filter_long_limit - (float)bsr_config->green_filter_short_limit) / 2.0f);
  blue_center=(float)bsr_config->blue_filter_short_limit + (((float)bsr_config->blue_filter_long_limit - (float)bsr_config->blue_filter_short_limit) / 2.0f);

  //
  // all threads: calculate Airy disk scaling factors (pixels) for each color.  Green is defined in config with Airy_disk_first_null
  //
  pixel_scaling_factor_green=3.8317f / (float)bsr_config->Airy_disk_first_null;
  pixel_scaling_factor_red=pixel_scaling_factor_green * green_center / red_center;
  pixel_scaling_factor_blue=pixel_scaling_factor_green * green_center / blue_center;

  // all threads: calculate pixel oversampling for each color to make full use of 10x Bessel function resolution from Bessel.h
  // and minimum of 11x11
  //
  half_oversampling_red=(int)((pixel_scaling_factor_red * 10.0f) + 0.5f);
  if (half_oversampling_red < 5) {
    half_oversampling_red=5;
  }
  oversampling_red=(half_oversampling_red * 2) + 1;
  half_oversampling_green=(int)((pixel_scaling_factor_green * 10.0f) + 0.5f);
  if (half_oversampling_green < 5) {
    half_oversampling_green=5;
  }
  oversampling_green=(half_oversampling_green * 2) +1;
  half_oversampling_blue=(int)((pixel_scaling_factor_blue * 10.0f) + 0.5f);
  if (half_oversampling_blue < 5) {
    half_oversampling_blue=5;
  }
  oversampling_blue=(half_oversampling_blue * 2) +1;

  //
  // all threads: set obs_ratio if central obstruction configured
  //
  if (bsr_config->Airy_disk_obstruction > 0.0) {
    if (bsr_config->Airy_disk_obstruction > 0.99) {
      obs_ratio=0.99;
    } else {
      obs_ratio=bsr_config->Airy_disk_obstruction;
    }
  } else {
    obs_ratio=0.0;
  }

  //
  // all threads: calculate center intensity calibration for each color
  //
  I0_green=1.16823f / powf(((float)bsr_config->Airy_disk_first_null * oversampling_green), 2.0f);
  I0_red=1.16823f * powf(green_center, 2.0f) / (powf(red_center, 2.0f) * powf(((float)bsr_config->Airy_disk_first_null * oversampling_red), 2.0f));
  I0_blue=1.16823f * powf(green_center, 2.0f) / (powf(blue_center, 2.0f)  * powf(((float)bsr_config->Airy_disk_first_null * oversampling_blue), 2.0f));

  //
  // all threads: generate Airy disk map for each color
  //
  makeAiryMap(bsr_state, bsr_state->Airymap_red, bsr_config->Airy_disk_max_extent, half_oversampling_red, pixel_scaling_factor_red, I0_red, obs_ratio);
  makeAiryMap(bsr_state, bsr_state->Airymap_green, bsr_config->Airy_disk_max_extent, half_oversampling_green, pixel_scaling_factor_green, I0_green, obs_ratio);
  makeAiryMap(bsr_state, bsr_state->Airymap_blue, bsr_config->Airy_disk_max_extent, half_oversampling_blue, pixel_scaling_factor_blue, I0_blue, obs_ratio);

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_AIRY_MAP_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_AIRY_MAP_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_AIRY_MAP_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_AIRY_MAP_CONTINUE;
    }
  } // end if not main thread

  //
  // main thread: output execution time if not in cgi mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1)) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  return 0;
}

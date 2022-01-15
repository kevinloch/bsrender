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

int makeAiryMap(bsr_state_t *bsr_state, double *Airymap, int max_extent, int half_oversampling, double pixel_scaling_factor, double I0, double obs_ratio) {
  double *Airymap_p;
  int Airymap_max_width;
  int map_offset;
  int map_index_x;
  int map_index_y;
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
  int lines_per_thread;
  double obs_I0_factor=0.0;
  int obs_x_index=0;
  double oversampling_factor;
  double oversample;

  //
  // set shorthand variables
  //
  Airymap_max_width=max_extent + 1;
  oversampling=(half_oversampling * 2) +1;
  lines_per_thread=(int)ceil(((double)Airymap_max_width / ((double)bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }
  if (obs_ratio > 0.0) {
    obs_I0_factor=1.0 / pow((1.0 - (obs_ratio * obs_ratio)), 2.0);
  }
  oversampling_factor=1.0 / oversampling;

  //
  // generate Airy disk map
  //
  map_index_x=0;
  map_index_y=(bsr_state->perthread->my_thread_id * lines_per_thread);
  Airymap_p=Airymap + (Airymap_max_width * map_index_y);
  for (map_offset=0; ((map_offset < (Airymap_max_width * lines_per_thread)) && (map_index_y < Airymap_max_width)); map_offset++) {
    pixel_x=(double)map_index_x;
    pixel_y=(double)map_index_y;
    pixel_r=sqrt((pixel_x * pixel_x) + (pixel_y * pixel_y));
    if ((pixel_r <= (double)max_extent) && ((pixel_r * pixel_scaling_factor) < 12800)) {
      *Airymap_p=0.0;
      oversample=0.0;
      for (oversample_y_index=0; oversample_y_index < oversampling; oversample_y_index++) {
        for (oversample_x_index=0; oversample_x_index < oversampling; oversample_x_index++) {
          oversample_x=(pixel_x + (((double)oversample_x_index - (double)half_oversampling) * oversampling_factor));
          oversample_y=(pixel_y + (((double)oversample_y_index - (double)half_oversampling) * oversampling_factor));
          oversample_r=sqrt((oversample_x * oversample_x) + (oversample_y * oversample_y));
          Bessel_x=oversample_r * pixel_scaling_factor;
          Bessel_x_index=(int)((Bessel_x * 10.0) + 0.5);
          if ((oversample_r == 0.0) || (Bessel_x_index == 0)) {
            oversample=I0;
          } else if (Bessel_x_index >= 128000) {
            oversample=0.0; // ignore if beyond range of Bessel.h (too many orders of diffraction)
            // skip to next map pixel
            oversample_x_index=(oversampling + 1);
            oversample_y_index=(oversampling + 1);
          } else if (obs_ratio > 0.0) {
            obs_x_index=(int)((obs_ratio * Bessel_x * 10.0) + 0.5);
            oversample=(I0 * obs_I0_factor * pow((((20.0 * Bessel_J1[Bessel_x_index]) - (20.0 * obs_ratio * Bessel_J1[obs_x_index])) / (double)Bessel_x_index), 2.0));
          } else {
            oversample=(I0 * pow((20.0 * Bessel_J1[Bessel_x_index] / (double)Bessel_x_index), 2.0));
          }
          *Airymap_p+=oversample;
        } // end for oversample_x
      } // end for oversample_y
    } else {
      *Airymap_p=0.0; // ignore if outside max radius
    } // end if pixel_r < max_extent

    map_index_x++;
    if (map_index_x == Airymap_max_width) {
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
  double red_center;
  double green_center;
  double blue_center;
  int i;
  double obs_ratio;
  const double I0_calibration=1.1675;

  //
  // main thread: display status message if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
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
  red_center=(double)bsr_config->red_filter_short_limit + (((double)bsr_config->red_filter_long_limit - (double)bsr_config->red_filter_short_limit) / 2.0);
  green_center=(double)bsr_config->green_filter_short_limit + (((double)bsr_config->green_filter_long_limit - (double)bsr_config->green_filter_short_limit) / 2.0);
  blue_center=(double)bsr_config->blue_filter_short_limit + (((double)bsr_config->blue_filter_long_limit - (double)bsr_config->blue_filter_short_limit) / 2.0);

  //
  // all threads: calculate Airy disk scaling factors (pixels) for each color. Green is defined in config with Airy_disk_first_null
  //
  pixel_scaling_factor_green=3.8317 / bsr_config->Airy_disk_first_null;
  pixel_scaling_factor_red=pixel_scaling_factor_green * green_center / red_center;
  pixel_scaling_factor_blue=pixel_scaling_factor_green * green_center / blue_center;

  // all threads: calculate pixel oversampling for each color to make full use of 10x Bessel function resolution from Bessel.h
  // and minimum of 11x11
  //
  half_oversampling_red=(int)((pixel_scaling_factor_red * 10.0) + 0.5);
  if (half_oversampling_red < 5) {
    half_oversampling_red=5;
  }
  oversampling_red=(half_oversampling_red * 2) + 1;
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
  I0_green=I0_calibration / pow((bsr_config->Airy_disk_first_null * oversampling_green), 2.0);
  I0_red=I0_calibration * pow(green_center, 2.0) / (pow(red_center, 2.0) * pow((bsr_config->Airy_disk_first_null * oversampling_red), 2.0));
  I0_blue=I0_calibration * pow(green_center, 2.0) / (pow(blue_center, 2.0)  * pow((bsr_config->Airy_disk_first_null * oversampling_blue), 2.0));

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
  // main thread: output execution time if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  return 0;
}

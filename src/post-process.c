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
  uint64_t image_offset;
  int current_image_x;
  int current_image_y;
  pixel_composition_t *current_image_p;
  double inv_camera_pixel_limit;
  double pixel_r;
  double pixel_g;
  double pixel_b;
  int current_image_res_x;
  int current_image_res_y;
  int lines_per_thread;
  int i;

  //
  // main thread: display status message if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Applying camera gamma and intensity limit...");
    fflush(stdout);
  }

  //
  // all threads: get current image resolution and lines per thread
  //
  current_image_res_x=bsr_state->current_image_res_x;
  current_image_res_y=bsr_state->current_image_res_y;
  lines_per_thread=(int)ceil(((double)current_image_res_y / (double)(bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }
  inv_camera_pixel_limit = 1.0 / bsr_state->camera_pixel_limit;

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_POST_PROCESS_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_POST_PROCESS_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: normalize pixels to 1.0 reference, and apply cmaera_gamma
  //
  current_image_x=0;
  current_image_y=bsr_state->perthread->my_thread_id * lines_per_thread;
  current_image_p=bsr_state->current_image_buf + ((uint64_t)current_image_res_x * (uint64_t)current_image_y);
  for (image_offset=0; ((image_offset < ((uint64_t)bsr_state->current_image_res_x * (uint64_t)lines_per_thread)) && (current_image_y < current_image_res_y)); image_offset++) {
    //
    // normalize pixel values to camera saturation reference level = 1.0
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
    // copy back to current image buf
    //
    current_image_p->r=pixel_r;
    current_image_p->g=pixel_g;
    current_image_p->b=pixel_b;

    // if end of this line, move to next line
    current_image_x++;
    if (current_image_x == bsr_state->current_image_res_x) {
      current_image_x=0;
      current_image_y++;
    }
    current_image_p++;
  } // end for i

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_POST_PROCESS_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_POST_PROCESS_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_POST_PROCESS_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_POST_PROCESS_CONTINUE;
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

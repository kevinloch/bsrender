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
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "util.h"

int resizeLanczos(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  pixel_composition_t *current_image_p;
  pixel_composition_t *image_resize_p;
  uint64_t current_image_offset;
  int resize_res_x;
  int resize_res_y;
  double source_w;
  double half_source_w;
  double source_x_center;
  double source_y_center;
  int source_x;
  int source_y;
  int resize_x;
  int resize_y;
  double L_kernel;
  double L_x_r;
  double L_x_g;
  double L_x_b;
  double L_y_r;
  double L_y_g;
  double L_y_b;
  double L_distance_x;
  double L_distance_y;
  uint64_t resize_i;
  int i_x;
  int i_y;
  int Lanczos_order;
  int current_image_res_x;
  int current_image_res_y;
  int lines_per_thread;
  int i;
  int current_image_x;
  int current_image_y;
  uint64_t image_offset;

  //
  // all threads: get current image resolution and calculate resize resolution
  //
  current_image_res_x=bsr_state->current_image_res_x;
  current_image_res_y=bsr_state->current_image_res_y;
  resize_res_x=bsr_state->resize_res_x;
  resize_res_y=bsr_state->resize_res_y;
  source_w=1.0 / bsr_config->output_scaling_factor;
  half_source_w=source_w / 2.0;

  //
  // main thread: display status message if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Resizing image from %dx%d to %dx%d...", current_image_res_x, current_image_res_y, resize_res_x, resize_res_y);
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_LANCZOS_PREP_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_LANCZOS_PREP_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: convert to log scale to reduce clipping artifacts. This will be undone
  // at the end of resizeLanczos()
  //
  current_image_x=0;
  lines_per_thread=(int)ceil(((double)current_image_res_y / (double)(bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }
  current_image_y=bsr_state->perthread->my_thread_id * lines_per_thread;
  current_image_p=bsr_state->current_image_buf + ((uint64_t)current_image_res_x * (uint64_t)current_image_y);
  for (image_offset=0; ((image_offset < ((uint64_t)bsr_state->current_image_res_x * (uint64_t)lines_per_thread)) && (current_image_y < current_image_res_y)); image_offset++) {
    L_x_r=current_image_p->r;
    L_x_g=current_image_p->g;
    L_x_b=current_image_p->b;
    L_x_r=log(BSR_RESIZE_LOG_OFFSET + L_x_r); 
    L_x_g=log(BSR_RESIZE_LOG_OFFSET + L_x_g);
    L_x_b=log(BSR_RESIZE_LOG_OFFSET + L_x_b);
    current_image_p->r=L_x_r;
    current_image_p->g=L_x_g;
    current_image_p->b=L_x_b;

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
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_LANCZOS_PREP_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_LANCZOS_RESAMPLE_BEGIN);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_LANCZOS_PREP_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_LANCZOS_RESAMPLE_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: copy rendered image to resize buffer using Lanczos interpolation
  //
  if (bsr_config->Lanczos_order < 2) {
    Lanczos_order=2;
  } else if (bsr_config->Lanczos_order > 10) {
    Lanczos_order=10;
  } else {
    Lanczos_order=bsr_config->Lanczos_order;
  }
  resize_x=0;
  lines_per_thread=(int)ceil(((double)resize_res_y / (double)(bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }
  resize_y=bsr_state->perthread->my_thread_id * lines_per_thread;
  image_resize_p=bsr_state->image_resize_buf + ((uint64_t)resize_res_x * (uint64_t)resize_y);
  for (resize_i=0; ((resize_i < ((uint64_t)resize_res_x * (uint64_t)lines_per_thread)) && (resize_y < resize_res_y)); resize_i++) {
    source_x_center=((double)resize_x * source_w) + half_source_w - 0.5;
    source_y_center=((double)resize_y * source_w) + half_source_w - 0.5;
    L_y_r=0.0;
    L_y_g=0.0;
    L_y_b=0.0;
    for (i_y=((int)source_y_center - Lanczos_order + 1); i_y <= ((int)source_y_center + Lanczos_order); i_y++) {
      source_y=i_y;

      //
      // L_x
      //
      L_x_r=0.0;
      L_x_g=0.0;
      L_x_b=0.0;
      for (i_x=((int)source_x_center - Lanczos_order + 1); i_x <= ((int)source_x_center + Lanczos_order); i_x++) {
        source_x=i_x;
        if ((source_x >= 0) && (source_x < current_image_res_x) && (source_y >= 0) && (source_y < current_image_res_y)) {
          current_image_offset=((uint64_t)source_y * (uint64_t)current_image_res_x) + (uint64_t)source_x;
          current_image_p=bsr_state->current_image_buf + current_image_offset;
          L_distance_x=source_x_center - (double)i_x;
          if (L_distance_x == 0.0) {
            L_x_r+=current_image_p->r;
            L_x_g+=current_image_p->g;
            L_x_b+=current_image_p->b;
          } else if ((L_distance_x >= -(double)Lanczos_order) && (L_distance_x <= (double)Lanczos_order)) {
            L_kernel=Lanczos_order * sin(M_PI * L_distance_x) * sin(M_PI * L_distance_x / (double)Lanczos_order) / (M_PI * M_PI * L_distance_x * L_distance_x);
            L_x_r+=(current_image_p->r * L_kernel);
            L_x_g+=(current_image_p->g * L_kernel);
            L_x_b+=(current_image_p->b * L_kernel);
          } else {
            // L_x == 0
          } // end if L_distance_x
        } // end if within composition buffer bounds
      } // end for i_x

      //
      // L_y
      //
      L_distance_y=source_y_center - (double)i_y;
      if (L_distance_y == 0.0) {
        L_y_r+=L_x_r;
        L_y_g+=L_x_g;
        L_y_b+=L_x_b;
      } else if ((L_distance_y >= -(double)Lanczos_order) && (L_distance_y <= (double)Lanczos_order)) {
        L_kernel=Lanczos_order * sin(M_PI * L_distance_y) * sin(M_PI * L_distance_y / (double)Lanczos_order) / (M_PI * M_PI * L_distance_y * L_distance_y);
        L_y_r+=(L_x_r * L_kernel);
        L_y_g+=(L_x_g * L_kernel);
        L_y_b+=(L_x_b * L_kernel);
      } else {
        // L_y = 0
      } // end if L_distance_y
    } // end for i_y

    //
    // undo log scaling
    //
    L_y_r=exp(L_y_r) - BSR_RESIZE_LOG_OFFSET;
    L_y_g=exp(L_y_g) - BSR_RESIZE_LOG_OFFSET;
    L_y_b=exp(L_y_b) - BSR_RESIZE_LOG_OFFSET;

    //
    // handle negative clipping, which is common
    //
    if (L_y_r < 0.0) {
      L_y_r=0.0;
    }
    if (L_y_g < 0.0) {
      L_y_g=0.0;
    }
    if (L_y_b < 0.0) {
      L_y_b=0.0;
    }

    //
    // copy to output buffer
    //
    image_resize_p->r=L_y_r;
    image_resize_p->g=L_y_g;
    image_resize_p->b=L_y_b;

    // if end of this line, move to next line
    resize_x++;
    if (resize_x == resize_res_x) {
      resize_x=0;
      resize_y++;
    }
    image_resize_p++;
  }

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_LANCZOS_RESAMPLE_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_LANCZOS_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_LANCZOS_RESAMPLE_COMPLETE);
    //
    // main thread: update current_image_buf pointer
    //
    if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
      bsr_state->current_image_buf=bsr_state->image_resize_buf;
      bsr_state->current_image_res_x=resize_res_x;
      bsr_state->current_image_res_y=resize_res_y;
    }
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_LANCZOS_CONTINUE;
    }
  } // end if not main thread

  //
  // main thread: output execution time if not in CGI mode
  //
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    if ((bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.3fs)\n", elapsed_time);
      fflush(stdout);
    }
  } // end if main thread

  return(0);
}

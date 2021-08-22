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
#include <unistd.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "util.h"

int GaussianBlur(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  int i;
  double radius;
  int lines_per_thread;
  int sample_size;
  int half_sample_size;
  int blur_res_x;
  int blur_res_y;
  long long blur_i;
  int blur_x;
  int blur_y;
  int source_x;
  int source_y;
  int kernel_i;
  int current_image_res_x;
  int current_image_res_y;
  long long current_image_offset;
  double *G_kernel_array;
  double *G_kernel_p;
  double G_kernel_sum;
  double G;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  pixel_composition_t *current_image_p;
  pixel_composition_t *image_blur_p;
  double G_r;
  double G_g;
  double G_b;

  //
  // all threads: determine Gaussian kernel size
  //
  radius=bsr_config->Gaussian_blur_radius;
  sample_size=((int)ceil(radius) * 6) + 1;
  half_sample_size=((int)ceil(radius) * 3) + 1;

  //
  // main thread: display status message if not in cgi mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Applying Gaussian blur with radius %.3e...", radius);
    fflush(stdout);
  }

  //
  // all threads: allocate memory for 1D Gaussian kernel
  //
  G_kernel_array=(double *)malloc(sample_size * sizeof(double));
  if (G_kernel_array == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate memory for Gaussian blur kernel\n");
      fflush(stdout);
    }
    return(1);
  }

  //
  // all threads: generate 1D Gaussian kernel
  //
  G_kernel_sum=0.0;
  G_kernel_p=G_kernel_array;
  for (kernel_i=-half_sample_size + 1; kernel_i < half_sample_size; kernel_i++) {
    G=exp(-(double)kernel_i * (double)kernel_i / (2.0 * radius * radius)) / sqrt(2.0 * M_PI * radius * radius);
    G_kernel_sum+=G;
    *G_kernel_p=G;
    G_kernel_p++;
  } // end for kernel

  //
  // all threads: normalize Gaussian kernel
  //
  G_kernel_p=G_kernel_array;
  for (kernel_i=-half_sample_size + 1; kernel_i < half_sample_size; kernel_i++) {
    *G_kernel_p/=G_kernel_sum;
    G_kernel_p++;
  } // end for kernel

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_GAUSSIAN_BLUR_HORIZONTAL_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_GAUSSIAN_BLUR_HORIZONTAL_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: get current image pointer and resolution
  //
  current_image_p=bsr_state->current_image_buf;
  current_image_res_x=bsr_state->current_image_res_x;
  current_image_res_y=bsr_state->current_image_res_y;
  blur_res_x=current_image_res_x;
  blur_res_y=current_image_res_y;
  lines_per_thread=(int)ceil(((double)current_image_res_y / (double)(bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }

  //
  // all threads: apply Gaussian 1D kernel to each pixel horizontally and put output in blur buffer
  //
  blur_x=0;
  blur_y=(bsr_state->perthread->my_thread_id * lines_per_thread);
  image_blur_p=bsr_state->image_blur_buf + ((long long)blur_res_x * (long long)blur_y);
  for (blur_i=0; ((blur_i < ((long long)blur_res_x * (long long)lines_per_thread)) && (blur_y < blur_res_y)); blur_i++) {

    //
    // apply Gaussian kernel to this pixel horizontally
    //
    G_r=0.0;
    G_g=0.0;
    G_b=0.0;
    G_kernel_p=G_kernel_array;
    for (kernel_i=-half_sample_size + 1; kernel_i < half_sample_size; kernel_i++) {
      source_x=blur_x + kernel_i;
      if ((source_x >= 0) && (source_x < current_image_res_x)) {
        current_image_offset=((long long)blur_y * (long long)blur_res_x) + (long long)source_x;
        current_image_p=bsr_state->current_image_buf + current_image_offset;
        G_r+=(current_image_p->r * *G_kernel_p);
        G_g+=(current_image_p->g * *G_kernel_p);
        G_b+=(current_image_p->b * *G_kernel_p);
      } // end if within current image bounds
      G_kernel_p++;
    } // end for kernel

    //
    // copy blurred pixel to blur buffer
    //
    image_blur_p->r=G_r;
    image_blur_p->g=G_g;
    image_blur_p->b=G_b;

    blur_x++;
    if (blur_x == blur_res_x) {
      blur_x=0;
      blur_y++;
    }
    image_blur_p++;
  } // end for blur_i

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_GAUSSIAN_BLUR_HORIZONTAL_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_GAUSSIAN_BLUR_VERTICAL_BEGIN);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_GAUSSIAN_BLUR_HORIZONTAL_COMPLETE);
    //
    // ready to continue, set all worker thread status to begin vertical
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_GAUSSIAN_BLUR_VERTICAL_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: apply Gaussian 1D kernel to each pixel vertically and put output back in 'current_image_buffer'
  // note in this step we use image_blur_buf as source and current_iamge_buf as dest so some variable names will be backwards
  //
  blur_x=0;
  blur_y=(bsr_state->perthread->my_thread_id * lines_per_thread);
  image_blur_p=bsr_state->current_image_buf + ((long long)blur_res_x * (long long)blur_y);
  for (blur_i=0; ((blur_i < ((long long)blur_res_x * (long long)lines_per_thread)) && (blur_y < blur_res_y)); blur_i++) {

    //
    // apply Gaussian kernel to this pixel vertically
    //
    G_r=0.0;
    G_g=0.0;
    G_b=0.0;
    G_kernel_p=G_kernel_array;
    for (kernel_i=-half_sample_size + 1; kernel_i < half_sample_size; kernel_i++) {
      source_y=blur_y + kernel_i;
      if ((source_y >= 0) && (source_y < current_image_res_y)) {
        current_image_offset=((long long)source_y * (long long)blur_res_x) + (long long)blur_x;
        current_image_p=bsr_state->image_blur_buf + current_image_offset;
        G_r+=(current_image_p->r * *G_kernel_p);
        G_g+=(current_image_p->g * *G_kernel_p);
        G_b+=(current_image_p->b * *G_kernel_p);
      } // end if within current image bounds
      G_kernel_p++;
    } // end for kernel

    //
    // copy blurred pixel to blur buffer
    //
    image_blur_p->r=G_r;
    image_blur_p->g=G_g;
    image_blur_p->b=G_b;

    blur_x++;
    if (blur_x == blur_res_x) {
      blur_x=0;
      blur_y++;
    }
    image_blur_p++;
  } // end for blur_i

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_GAUSSIAN_BLUR_VERTICAL_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_GAUSSIAN_BLUR_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_GAUSSIAN_BLUR_VERTICAL_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_GAUSSIAN_BLUR_CONTINUE;
    }
  } // end if not main thread

  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    //
    // main thread: display execution time if not in cgi mode
    //
    if (bsr_config->cgi_mode != 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.4fs)\n", elapsed_time);
      fflush(stdout);
    }
  } // end if main thread

  return(0);
}

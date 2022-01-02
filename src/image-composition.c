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
#include "Gaia-passbands.h"

int initImageCompositionBuffer(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  long long image_offset;
  int current_image_x;
  int current_image_y;
  pixel_composition_t *current_image_p;
  int current_image_res_x;
  int current_image_res_y;
  int lines_per_thread;
  int i;
  int skyglow_temp;
  double skyglow_intensity;
  double background_red;
  double background_green;
  double background_blue;

  //
  // all threads: get current image resolution and lines per thread
  //
  current_image_res_x=bsr_state->current_image_res_x;
  current_image_res_y=bsr_state->current_image_res_y;
  lines_per_thread=(int)ceil(((float)current_image_res_y / (float)(bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }

  //
  // main thread: display status message if not in cgi mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing image composition buffer %dx%d...", current_image_res_x, current_image_res_y);
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_INIT_IMAGECOMP_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_INIT_IMAGECOMP_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: set background value
  //
  if (bsr_config->skyglow_enable == 1) {
    // set skyglow rgb values
    // note: rgb lookup table values are already adjusted for Gaia Gband transmissivity, so we must uncorrect for that
    skyglow_temp=(int)(bsr_config->skyglow_temp + 0.5);
    skyglow_intensity=Gaia_Gband_scalar * pow(100.0, (-bsr_config->skyglow_per_pixel_mag / 5.0));
    background_red=skyglow_intensity * bsr_state->rgb_red[skyglow_temp];
    background_green=skyglow_intensity * bsr_state->rgb_green[skyglow_temp];
    background_blue=skyglow_intensity * bsr_state->rgb_blue[skyglow_temp];
  } else {
    // no skyglow set to zero
    background_red=0.0;
    background_green=0.0;
    background_blue=0.0;
  }

  //
  // all threads: initialize image composition buffer
  //
  current_image_x=0;
  current_image_y=(bsr_state->perthread->my_thread_id * lines_per_thread);
  current_image_p=bsr_state->current_image_buf + ((long long)current_image_res_x * (long long)current_image_y);
  for (image_offset=0; ((image_offset < ((long long)bsr_state->current_image_res_x * (long long)lines_per_thread)) && (current_image_y < current_image_res_y)); image_offset++) {
    //
    // set pixel rgb to background value (0.0 if now skyglow)
    //
    current_image_p->r=background_red;
    current_image_p->g=background_green;
    current_image_p->b=background_blue;
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
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_INIT_IMAGECOMP_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_INIT_IMAGECOMP_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_INIT_IMAGECOMP_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_INIT_IMAGECOMP_CONTINUE;
    }
  } // end if not main thread

  //
  // main thread: output execution time if not in cgi mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1)) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.4fs)\n", elapsed_time);
    fflush(stdout);
  }

  return(0);
}

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

int quantize(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  long long image_offset;
  pixel_composition_t *current_image_p;
  png_byte *image_output_p;
  const double one_over_2dot4=1.0 / 2.4;
  double pixel_r;
  double pixel_g;
  double pixel_b;
  int lines_per_thread;
  int i;
  int output_res_x;
  int output_res_y;
  int output_x;
  int output_y;
  int bpp; // bytes per pixel

  //
  // main thread: display status message if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    if (bsr_config->bits_per_color == 8) {
      printf("Converting to 8 bits per color...");
    } else if (bsr_config->bits_per_color == 16) {
      printf("Converting to 16 bits per color...");
    }
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_QUANTIZE_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_QUANTIZE_BEGIN;
    }
  } // end if not main thread

  //
  // all threads: get current image resolution and lines per thread
  //
  output_res_x=bsr_state->current_image_res_x;
  output_res_y=bsr_state->current_image_res_y;
  lines_per_thread=(int)ceil(((float)bsr_state->current_image_res_y / (float)(bsr_state->num_worker_threads + 1)));
  if (lines_per_thread < 1) {
    lines_per_thread=1;
  }

  //
  // all threads: quantize current_image_buf to 8 or 16 bits per pixel and store in image_output_buf
  // also update row_pointers for libpng
  //
  if (bsr_config->bits_per_color == 8) {
    bpp=3;
  } else {
    bpp=6;
  }
  output_x=0;
  output_y=(bsr_state->perthread->my_thread_id * lines_per_thread);
  image_output_p=bsr_state->image_output_buf + ((long long)output_res_x * (long long)output_y * (long long)bpp);
  current_image_p=bsr_state->current_image_buf + ((long long)output_res_x * (long long)output_y);
  if (output_y < output_res_y) {
    bsr_state->row_pointers[output_y]=image_output_p;
  }
  for (image_offset=0; ((image_offset < ((long long)output_res_x * (long long)lines_per_thread)) && (output_y < output_res_y)); image_offset++) {
    //
    // copy pixel data from current_image_buf
    //
    pixel_r=current_image_p->r;
    pixel_g=current_image_p->g;
    pixel_b=current_image_p->b;

    //
    // limit pixel intensity to range [0..1]
    //
    if (bsr_config->camera_pixel_limit_mode == 0) {
      limitIntensity(&pixel_r, &pixel_g, &pixel_b);
    } else if (bsr_config->camera_pixel_limit_mode == 1) {
      limitIntensityPreserveColor(&pixel_r, &pixel_g, &pixel_b);
    }

    //
    // apply encoding gamma for color spaces that use it
    //
    if ((bsr_config->icc_profile == 1) || (bsr_config->icc_profile == 2)) { // sRGB and Display-P3
      if (pixel_r <= 0.0031308) {
        pixel_r=pixel_r * 12.92;
      } else {
        pixel_r=(1.055 * pow(pixel_r, one_over_2dot4) - 0.055);
      }
      if (pixel_g <= 0.0031308) {
        pixel_g=pixel_g * 12.92;
      } else {
        pixel_g=(1.055 * pow(pixel_g, one_over_2dot4) - 0.055);
      }
      if (pixel_b <= 0.0031308) {
        pixel_b=pixel_b * 12.92;
      } else {
        pixel_b=(1.055 * pow(pixel_b, one_over_2dot4) - 0.055);
      }
    } else if (bsr_config->icc_profile == 3) { // Rec. 2020
      if (pixel_r < 0.018053968510807) {
        pixel_r=pixel_r * 4.5;
      } else {
        pixel_r=(1.09929682680944 * pow(pixel_r, 0.45) - 0.09929682680944);
      }
      if (pixel_g < 0.018053968510807) {
        pixel_g=pixel_g * 4.5;
      } else {
        pixel_g=(1.09929682680944 * pow(pixel_g, 0.45) - 0.09929682680944);
      }
      if (pixel_b < 0.018053968510807) {
        pixel_b=pixel_b * 4.5;
      } else {
        pixel_b=(1.09929682680944 * pow(pixel_b, 0.45) - 0.09929682680944);
      }
    }

    //
    // convert r,g,b to 8 or 16 bit values and copy to output buffer
    //
    if (bsr_config->bits_per_color == 8) {
      *image_output_p=(unsigned char)((pixel_r * 255.0) + 0.5);
      image_output_p++;
      *image_output_p=(unsigned char)((pixel_g * 255.0) + 0.5);
      image_output_p++;
      *image_output_p=(unsigned char)((pixel_b * 255.0) + 0.5);
      image_output_p++;
    } else if (bsr_config->bits_per_color == 16) {
      png_save_uint_16(image_output_p, (uint16_t)((pixel_r * 65535.0) + 0.5));
      image_output_p+=2;
      png_save_uint_16(image_output_p, (uint16_t)((pixel_g * 65535.0) + 0.5));
      image_output_p+=2;
      png_save_uint_16(image_output_p, (uint16_t)((pixel_b * 65535.0) + 0.5));
      image_output_p+=2;
    }

    //
    // set new row pointer if we have reached end of row
    //
    output_x++;
    if (output_x == output_res_x) {
      output_x=0;
      output_y++;
      if (((image_offset + (long long)1) < ((long long)output_res_x * (long long)lines_per_thread)) && (output_y < output_res_y)) {
        bsr_state->row_pointers[output_y]=image_output_p;
      }
    }
    current_image_p++;
  } // end for i

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_QUANTIZE_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_QUANTIZE_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_QUANTIZE_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_QUANTIZE_CONTINUE;
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

  return(0);
}

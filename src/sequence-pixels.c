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

int sequencePixels(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  //
  // This function takes pixel data from the current_image_buf after image generation and post processing
  // and converts it into the specific unsigned char sequence required by an image output format.
  // Image formats use a variety of number formats (integer/floating-point), encoding gamma, bit depth, color
  // channel order, and endianness. The result is stored in image_output_buf and row_pointers.
  //

  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  uint64_t image_offset;
  pixel_composition_t *current_image_p;
  unsigned char *image_output_p;
  unsigned char *image_output_R_p=NULL; // used by EXR since it groups same-channel pixel data together
  unsigned char *image_output_G_p=NULL; // used by EXR since it groups same-channel pixel data together
  unsigned char *image_output_B_p=NULL; // used by EXR since it groups same-channel pixel data together
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
  int bytes_per_pixel=0;
  int bytes_per_color=0;

  //
  // main thread: display status message if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    if (bsr_config->image_number_format == 0) {
      if (bsr_config->bits_per_color == 8) { 
        printf("Converting to 8-bit unsigned integer per color...");
      } else if (bsr_config->bits_per_color == 16) {
        printf("Converting to 16-bit unsigned integer per color...");
      } else if (bsr_config->bits_per_color == 32) {
        printf("Converting to 32-bit unsigned integer per color...");
      }
    } else if (bsr_config->image_number_format == 1) {
      if (bsr_config->bits_per_color == 16) {
        printf("Converting to 16-bit floating-point per color...");
      } else if (bsr_config->bits_per_color == 32) {
        printf("Converting to 32-bit floating-point per color...");
      }
    }
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_SEQUENCE_PIXELS_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_SEQUENCE_PIXELS_BEGIN;
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
  // all threads: convert current_image_buf to unsigned char byte sequence and store
  // in image_output_buf. Also update row_pointers
  //
  if (bsr_config->bits_per_color == 8) {
    bytes_per_color=1;
    bytes_per_pixel=3;
  } else if (bsr_config->bits_per_color == 16) {
    bytes_per_color=2;
    bytes_per_pixel=6;
  } else if (bsr_config->bits_per_color == 32) {
    bytes_per_color=4;
    bytes_per_pixel=12;
  }
  output_x=0;
  output_y=bsr_state->perthread->my_thread_id * lines_per_thread;
  image_output_p=bsr_state->image_output_buf + ((uint64_t)output_res_x * (uint64_t)output_y * (uint64_t)bytes_per_pixel);
  if (bsr_config->image_format == 1) {
    // EXR groups same channel pixel data together
    image_output_B_p=image_output_p;
    image_output_G_p=image_output_p + ((uint64_t)bytes_per_color * (uint64_t)output_res_x);
    image_output_R_p=image_output_p + (2ll * (uint64_t)bytes_per_color * (uint64_t)output_res_x);
  }
  current_image_p=bsr_state->current_image_buf + ((uint64_t)output_res_x * (uint64_t)output_y);
  if (output_y < output_res_y) {
    bsr_state->row_pointers[output_y]=image_output_p;
  }
  for (image_offset=0; ((image_offset < ((uint64_t)output_res_x * (uint64_t)lines_per_thread)) && (output_y < output_res_y)); image_offset++) {
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
      limitIntensity(bsr_config, &pixel_r, &pixel_g, &pixel_b);
    } else if (bsr_config->camera_pixel_limit_mode == 1) {
      limitIntensityPreserveColor(bsr_config, &pixel_r, &pixel_g, &pixel_b);
    }

    //
    // apply encoding gamma for image formats and color spaces that use it
    //
    if (bsr_config->image_format != 1) { // EXR does not use encoding gamma
      if ((bsr_config->icc_profile == -1) || (bsr_config->icc_profile == 1) || (bsr_config->icc_profile == 2)) { // PNG default, sRGB and Display-P3
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
      } else if ((bsr_config->icc_profile == 3) || (bsr_config->icc_profile == 4)\
              || (bsr_config->icc_profile == 5) || (bsr_config->icc_profile == 6)) { // Rec. 2020, Rec. 601 NTSC, Rec. 601 PAL, Rec. 709
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
      } else if (bsr_config->icc_profile == 7) { // flat 2.0 gamma
        pixel_r=pow(pixel_r, 0.5);
        pixel_g=pow(pixel_g, 0.5);
        pixel_b=pow(pixel_b, 0.5);
      } // end if icc_profile
    } // end if image_format

    //
    // convert r,g,b to output byte sequence and store in output buffer
    //
    if (bsr_config->image_format == 0) {
      // PNG format. Channels are RGB RGB RGB order and stored big-endian
      if (bsr_config->bits_per_color == 8) {
        *image_output_p=(unsigned char)((pixel_r * 255.0) + 0.5);
        image_output_p+=bytes_per_color;
        *image_output_p=(unsigned char)((pixel_g * 255.0) + 0.5);
        image_output_p+=bytes_per_color;
        *image_output_p=(unsigned char)((pixel_b * 255.0) + 0.5);
        image_output_p+=bytes_per_color;
      } else if (bsr_config->bits_per_color == 16) {
        storeU16BE(image_output_p,  (uint16_t)((pixel_r * 65535.0) + 0.5));
        image_output_p+=bytes_per_color;
        storeU16BE(image_output_p,  (uint16_t)((pixel_g * 65535.0) + 0.5));
        image_output_p+=bytes_per_color;
        storeU16BE(image_output_p,  (uint16_t)((pixel_b * 65535.0) + 0.5));
        image_output_p+=bytes_per_color;
      } // end if bits_per_color
    } else if (bsr_config->image_format == 1) {
      // EXR format. Channels are in BBB GGG RRR order and stored little-endian
      if (bsr_config->image_number_format == 1) {
        // floating-point
        if (bsr_config->bits_per_color == 16) {
          storeHalfLE(image_output_B_p, (float)pixel_b);
          image_output_B_p+=bytes_per_color;
          storeHalfLE(image_output_G_p, (float)pixel_g);
          image_output_G_p+=bytes_per_color;
          storeHalfLE(image_output_R_p, (float)pixel_r);
          image_output_R_p+=bytes_per_color;
          image_output_p+=bytes_per_pixel; // to keep beginning of row in sync
        } else if (bsr_config->bits_per_color == 32) {
          storeFloatLE(image_output_B_p, (float)pixel_b);
          image_output_B_p+=bytes_per_color;
          storeFloatLE(image_output_G_p, (float)pixel_g);
          image_output_G_p+=bytes_per_color;
          storeFloatLE(image_output_R_p, (float)pixel_r);
          image_output_R_p+=bytes_per_color;
          image_output_p+=bytes_per_pixel; // to keep beginning of row in sync
        } // end if bits_per_color
      } else if (bsr_config->image_number_format == 0) {
        // unsigned integer 
        storeU32LE(image_output_B_p, (uint32_t)((pixel_b * 4294967295.0) + 0.5));
        image_output_B_p+=bytes_per_color;
        storeU32LE(image_output_G_p, (uint32_t)((pixel_g * 4294967295.0) + 0.5));
        image_output_G_p+=bytes_per_color;
        storeU32LE(image_output_R_p, (uint32_t)((pixel_r * 4294967295.0) + 0.5));
        image_output_R_p+=bytes_per_color;
        image_output_p+=bytes_per_pixel; // to keep beginning of row in sync
      } // end if image_number_format
    } // end if image_format

    //
    // set new row pointer if we have reached end of row
    //
    output_x++;
    if (output_x == output_res_x) {
      output_x=0;
      output_y++;
      if (((image_offset + (uint64_t)1) < ((uint64_t)output_res_x * (uint64_t)lines_per_thread)) && (output_y < output_res_y)) {
        bsr_state->row_pointers[output_y]=image_output_p;
      }
      if (bsr_config->image_format == 1) {
        // update EXR output pointers
        image_output_B_p=image_output_p;
        image_output_G_p=image_output_p + ((uint64_t)bytes_per_color * (uint64_t)output_res_x);
        image_output_R_p=image_output_p + (2ll * (uint64_t)bytes_per_color * (uint64_t)output_res_x);
      } // end if image_format
    } //end if output_x

    current_image_p++;
  } // end for i

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_SEQUENCE_PIXELS_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_SEQUENCE_PIXELS_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_SEQUENCE_PIXELS_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_SEQUENCE_PIXELS_CONTINUE;
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

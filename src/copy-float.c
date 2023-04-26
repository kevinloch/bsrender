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

int storeHalfLE(unsigned char *dest, float src) {
  unsigned char *src_p;
  unsigned char *dest_p;
  unsigned char tmp8;

  //
  // convert binary32 to binary16 and store little-endian
  //
#ifdef BSR_LITTLE_ENDIAN_COMPILE
  // start at least significant byte of src and dest
  src_p=(unsigned char *)&src;
  dest_p=dest;

  src_p++; // move to second source byte
  tmp8=*src_p;
  tmp8 >>= 5; // 3 msb of second source byte becomes 3 lsb of binary16 significand
  *dest_p=tmp8;

  src_p++; // move to third source byte
  tmp8=*src_p;
  tmp8 <<= 3; // 5 lsb of third source byte becomes next 5 sb of binary16 significand
  *dest_p |= tmp8;
  dest_p++; // move to second dest byte
  tmp8=*src_p;
  tmp8 &= 0x7f; // suppress msb = least significant binary32 exponent bit
  tmp8 >>= 5; // bits 6,7 of third source byte are 2 msb of binary16 significand
  *dest_p=tmp8;

  src_p++; // move to fourth source byte
  tmp8=*src_p;
  tmp8 &= 0xfc; // 6 msb of fourth source byte has sign and all 5 binary16 exponent bits
  *dest_p |= tmp8;
#elif defined BSR_BIG_ENDIAN_COMPILE
  // tbd
#endif

  return(0);
}

int storeFloatLE(unsigned char *dest, float src) {
  unsigned char *src_p;
  unsigned char *dest_p;

  dest_p=dest;
#ifdef BSR_LITTLE_ENDIAN_COMPILE
  src_p=(unsigned char *)&src;
  *dest_p=*src_p;
  dest_p++;
  src_p++;
  *dest_p=*src_p;
  dest_p++;
  src_p++;
  *dest_p=*src_p;
  dest_p++;
  src_p++;
  *dest_p=*src_p;
#elif defined BSR_BIG_ENDIAN_COMPILE
  src_p=(unsigned char *)&src;
  src_p+=3;
  *dest_p=*src_p;
  dest_p++;
  src_p--;
  *dest_p=*src_p;
  dest_p++;
  src_p--;
  *dest_p=*src_p;
  dest_p++;
  src_p--;
  *dest_p=*src_p;
#endif

  return(0);
}

int copyFloat(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  long long image_offset;
  pixel_composition_t *current_image_p;
  unsigned char *image_output_p;
  double pixel_r;
  double pixel_g;
  double pixel_b;
  int lines_per_thread;
  int i;
  int output_res_x;
  int output_res_y;
  int output_x;
  int output_y;
  int bytes_per_pixel;

  //
  // main thread: display status message if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    if (bsr_config->bits_per_color == 32) {
      printf("Converting to 32-bit floating-point per color...");
    } else { // default 16 bits per color
      printf("Converting to 16-bit floating-point per color...");
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
  // all threads: copy current_image_buf to image_output_buf
  // also update row_pointers
  //
  if (bsr_config->bits_per_color == 32) {
    bytes_per_pixel=12;
  } else { // default 16 bits per color
    bytes_per_pixel=6;
  }
  output_x=0;
  output_y=bsr_state->perthread->my_thread_id * lines_per_thread;
  image_output_p=bsr_state->image_output_buf + ((long long)output_res_x * (long long)output_y * (long long)bytes_per_pixel);
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
    // limit pixel intensity to range [0..1] (> 0.0 for floating point)
    //
    if (bsr_config->camera_pixel_limit_mode == 0) {
      limitIntensity(bsr_config, &pixel_r, &pixel_g, &pixel_b);
    } else if (bsr_config->camera_pixel_limit_mode == 1) {
      limitIntensityPreserveColor(bsr_config, &pixel_r, &pixel_g, &pixel_b);
    }

    //
    // convert r,g,b to 32, or 16 bit floating-point values and copy to output buffer
    //
    // EXR is currently the only supported image format that uses floating-point. Channels are in BGR order and stored little-endian
    if (bsr_config->bits_per_color == 32) {
      storeFloatLE(image_output_p, (float)pixel_b);
      image_output_p+=4;
      storeFloatLE(image_output_p, (float)pixel_g);
      image_output_p+=4;
      storeFloatLE(image_output_p, (float)pixel_r);
      image_output_p+=4;
    } else { // default 16 bit half precision
      storeHalfLE(image_output_p, (float)pixel_b);
      image_output_p+=2;
      storeHalfLE(image_output_p, (float)pixel_g);
      image_output_p+=2;
      storeHalfLE(image_output_p, (float)pixel_r);
      image_output_p+=2;
    } // end if image_format

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

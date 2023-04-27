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
#include <string.h>
#include "util.h"
#include "cgi.h"
#include "icc-profiles.h"
#include "bsr-exr.h"
#include "sequence-pixels.h"

int storeI32LE(unsigned char *dest, int32_t src) {
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

  // return number of bytes stored
  return(4);
}

int storeU8(unsigned char *dest, unsigned char src) {
  // seems silly but helps keep outputEXRHeader clean

  *dest=src;

  // return number of bytes stored
  return(1);
}

int storeStr32(unsigned char *dest, char *src) {
  int size;

  snprintf((char *)dest, 32, "%s", src);

  // return number of bytes stored
  size=strlen(src) + 1;
  if (size > 32) {
    return(32);
  } else {
    return(size);
  }
}

int outputEXRHeader(bsr_config_t *bsr_config, FILE *output_file) {
  unsigned char header[4096];
  unsigned char *header_p;
  int pixel_type=EXR_PIXEL_FLOAT;

  // pixel_type shortcut
  if (bsr_config->image_number_format == 0) {
    if (bsr_config->bits_per_color == 32) {
      pixel_type=EXR_PIXEL_UINT;
    }
  } else if (bsr_config->image_number_format == 1) {
    if (bsr_config->bits_per_color == 16) {
      pixel_type=EXR_PIXEL_HALF;
    } else if (bsr_config->bits_per_color == 32) {
      pixel_type=EXR_PIXEL_FLOAT;
    }
  }

  header_p=header;

  // magic number
  header_p+=storeU8(header_p, 0x76);
  header_p+=storeU8(header_p, 0x2f);
  header_p+=storeU8(header_p, 0x31);
  header_p+=storeU8(header_p, 0x01);

  // version
  header_p+=storeU8(header_p, 0x02); // EXR version 2.0
  header_p+=storeU8(header_p, 0);    // version flags are all zero (no tiles, no long attribute names, no deep data, not multi-part)
  header_p+=storeU8(header_p, 0);    // reserved
  header_p+=storeU8(header_p, 0);    // reserved

  // channels
  header_p+=storeStr32(header_p, "channels");  // attribute name
  header_p+=storeStr32(header_p, "chlist");    // attribute type
  header_p+=storeI32LE(header_p, (int32_t)55); // attribute value length in bytes
    // Blue
    header_p+=storeStr32(header_p, "B");                       // channel name
    header_p+=storeI32LE(header_p, (int32_t)pixel_type);       // half, float, uint32
    header_p+=storeU8(header_p, EXR_PERCEPTUALLY_LOGARITHMIC); // perceptual treatment (RGB=logarithmic, XYZ=linear)
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeI32LE(header_p, (int32_t)1);                // x_sampling=1 for RGB (4:4:4)
    header_p+=storeI32LE(header_p, (int32_t)1);                // y_sampling=1 for RGB (4:4:4)
    // Green
    header_p+=storeStr32(header_p, "G");                       // channel name
    header_p+=storeI32LE(header_p, (int32_t)pixel_type);       // half, float, uint32
    header_p+=storeU8(header_p, EXR_PERCEPTUALLY_LOGARITHMIC); // perceptual treatment (RGB=logarithmic, XYZ=linear)
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeI32LE(header_p, (int32_t)1);                // x_sampling=1 for RGB (4:4:4)
    header_p+=storeI32LE(header_p, (int32_t)1);                // y_sampling=1 for RGB (4:4:4)
    // Red
    header_p+=storeStr32(header_p, "R");                       // channel name
    header_p+=storeI32LE(header_p, (int32_t)pixel_type);       // half, float, uint32
    header_p+=storeU8(header_p, EXR_PERCEPTUALLY_LOGARITHMIC); // perceptual treatment (RGB=logarithmic, XYZ=linear)
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeU8(header_p, 0);                            // reserved
    header_p+=storeI32LE(header_p, (int32_t)1);                // x_sampling=1 for RGB (4:4:4)
    header_p+=storeI32LE(header_p, (int32_t)1);                // y_sampling=1 for RGB (4:4:4)
    // null byte indicating no more channels
    header_p+=storeU8(header_p, 0);

  // compression
  header_p+=storeStr32(header_p, "compression");    // attribute name
  header_p+=storeStr32(header_p, "compression");    // attribute type
  header_p+=storeI32LE(header_p, (int32_t)1);       // attribute value length in bytes
  //header_p+=storeU8(header_p, EXR_COMPRESSION_NONE); // for testing
  //header_p+=storeU8(header_p, EXR_COMPRESSION_ZIPS); // deflate algorithm, 16 lines at a time
  header_p+=storeU8(header_p, EXR_COMPRESSION_ZIP); // deflate algorithm, one line at a time

  // data window
  header_p+=storeStr32(header_p, "dataWindow");                          // attribute name
  header_p+=storeStr32(header_p, "box2i");                               // attribute type
  header_p+=storeI32LE(header_p, (int32_t)16);                           // attribute value length in bytes
  header_p+=storeI32LE(header_p, (int32_t)0);                            // xmin
  header_p+=storeI32LE(header_p, (int32_t)0);                            // ymin
  header_p+=storeI32LE(header_p, (int32_t)(bsr_config->camera_res_x-1)); // xmax
  header_p+=storeI32LE(header_p, (int32_t)(bsr_config->camera_res_y-1)); // ymax

  // display window
  header_p+=storeStr32(header_p, "displayWindow");                       // attribute name
  header_p+=storeStr32(header_p, "box2i");                               // attribute type
  header_p+=storeI32LE(header_p, (int32_t)16);                           // attribute value length in bytes
  header_p+=storeI32LE(header_p, (int32_t)0);                            // xmin
  header_p+=storeI32LE(header_p, (int32_t)0);                            // ymin
  header_p+=storeI32LE(header_p, (int32_t)(bsr_config->camera_res_x-1)); // xmax
  header_p+=storeI32LE(header_p, (int32_t)(bsr_config->camera_res_y-1)); // ymax

  // line order
  header_p+=storeStr32(header_p, "lineOrder");             // attribute name
  header_p+=storeStr32(header_p, "lineOrder");             // attribute type
  header_p+=storeI32LE(header_p, (int32_t)1);              // attribute value length in bytes
  header_p+=storeU8(header_p, EXR_LINEORDER_INCREASING_Y); // line order

  // pixel aspect ratio
  header_p+=storeStr32(header_p, "pixelAspectRatio"); // attribute name
  header_p+=storeStr32(header_p, "float");            // attribute type
  header_p+=storeI32LE(header_p, (int32_t)4);         // attribute value length in bytes
  header_p+=storeFloatLE(header_p, 1.0f);             // pixel aspect ratio

  // screen window center
  header_p+=storeStr32(header_p, "screenWindowCenter"); // attribute name
  header_p+=storeStr32(header_p, "v2f");                // attribute type
  header_p+=storeI32LE(header_p, (int32_t)8);           // attribute value length in bytes
  header_p+=storeFloatLE(header_p, 0.0f);               // window center x?
  header_p+=storeFloatLE(header_p, 0.0f);               // window center y?

  // screen window width
  header_p+=storeStr32(header_p, "screenWindowWidth"); // attribute name
  header_p+=storeStr32(header_p, "float");             // attribute type
  header_p+=storeI32LE(header_p, (int32_t)8);          // attribute value length in bytes
  header_p+=storeFloatLE(header_p, 1.0f);              // window width


//
// incomplete - work in progress
//



  // output header
  if (bsr_config->cgi_mode == 1) {
    fwrite(header, (header_p - header), 1, stdout);
  } else {
    fwrite(header, (header_p - header), 1, output_file);
  }

  return(0);
}

int outputEXR(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  FILE *output_file=NULL;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  int i;

  //
  // main thread: display status update if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Writing %s...", bsr_config->output_file_name);
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_IMAGE_OUTPUT_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_IMAGE_OUTPUT_BEGIN;
    }
  } // end if not main thread

  //
  // main thread: open output file if not CGI mode
  //
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    output_file=fopen(bsr_config->output_file_name, "wb");
    if (output_file == NULL) {
      printf("Error: could not open %s for writing\n", bsr_config->output_file_name);
      fflush(stdout);
      exit(1);
    }
  }




//
// incomplete - work in progress
//





  //
  // main thread: construct and output EXR file header
  //
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    outputEXRHeader(bsr_config, output_file);
  } // end if main thread

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_IMAGE_OUTPUT_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_IMAGE_OUTPUT_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_IMAGE_OUTPUT_COMPLETE);
    //
    // ready to continue, set all worker thread status to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_IMAGE_OUTPUT_CONTINUE;
    }
  } // end if not main thread

  //
  // main thread: display status message and close output file if not CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1)) {
    if (bsr_config->print_status == 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.3fs)\n", elapsed_time);
      fflush(stdout);
    }

    // clean up
    fclose(output_file);
  }

  return(0);
}

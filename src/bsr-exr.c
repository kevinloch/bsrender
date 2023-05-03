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
#include <time.h>
#include "util.h"
#include "cgi.h"
#include "icc-profiles.h"
#include "bsr-exr.h"
#include "sequence-pixels.h"

int outputEXRHeader(bsr_config_t *bsr_config, bsr_state_t *bsr_state, FILE *output_file) {
  unsigned char header[4096];
  unsigned char *header_p;
  int pixel_type=EXR_PIXEL_FLOAT;
  chromaticities_t chromaticities=Rec709_c;
  int header_size;

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
  header_p+=storeU32LE(header_p, (uint32_t)0x1312F76);

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
  header_p+=storeStr32(header_p, "compression");       // attribute name
  header_p+=storeStr32(header_p, "compression");       // attribute type
  header_p+=storeI32LE(header_p, (int32_t)1);          // attribute value length in bytes
  header_p+=storeU8(header_p, EXR_COMPRESSION_NONE); // for testing
  //header_p+=storeU8(header_p, EXR_COMPRESSION_ZIPS); // deflate algorithm, one line at a time
  //header_p+=storeU8(header_p, EXR_COMPRESSION_ZIP);    // deflate algorithm, 16 lines at a time

  // data window
  header_p+=storeStr32(header_p, "dataWindow");                                // attribute name
  header_p+=storeStr32(header_p, "box2i");                                     // attribute type
  header_p+=storeI32LE(header_p, (int32_t)16);                                 // attribute value length in bytes
  header_p+=storeI32LE(header_p, (int32_t)0);                                  // xmin
  header_p+=storeI32LE(header_p, (int32_t)0);                                  // ymin
  header_p+=storeI32LE(header_p, (int32_t)(bsr_state->current_image_res_x-1)); // xmax
  header_p+=storeI32LE(header_p, (int32_t)(bsr_state->current_image_res_y-1)); // ymax

  // display window
  header_p+=storeStr32(header_p, "displayWindow");                             // attribute name
  header_p+=storeStr32(header_p, "box2i");                                     // attribute type
  header_p+=storeI32LE(header_p, (int32_t)16);                                 // attribute value length in bytes
  header_p+=storeI32LE(header_p, (int32_t)0);                                  // xmin
  header_p+=storeI32LE(header_p, (int32_t)0);                                  // ymin
  header_p+=storeI32LE(header_p, (int32_t)(bsr_state->current_image_res_x-1)); // xmax
  header_p+=storeI32LE(header_p, (int32_t)(bsr_state->current_image_res_y-1)); // ymax

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
  header_p+=storeI32LE(header_p, (int32_t)4);          // attribute value length in bytes
  header_p+=storeFloatLE(header_p, 1.0f);              // window width

  //
  // EXR does not support ICC profiles but it does have a standard optional attribute
  // for color space information. Default is no chromaticity header
  //
  if ((bsr_config->icc_profile >= 1) && (bsr_config->icc_profile <= 6)) {
    if (bsr_config->icc_profile == 1) {
      // sRGB
      chromaticities=sRGB_c;
    } else if (bsr_config->icc_profile == 2) {
      // Display-P3
      chromaticities=DisplayP3_c;
    } else if (bsr_config->icc_profile == 3) {
      // Rec. 2020
      chromaticities=Rec2020_c;
    } else if (bsr_config->icc_profile == 4) {
      // Rec. 601 NTSC
      chromaticities=Rec601NTSC_c;
    } else if (bsr_config->icc_profile == 5) {
      // Rec. 601 PAL
      chromaticities=Rec601PAL_c;
    } else if (bsr_config->icc_profile == 6) {
      // Rec. 709
      chromaticities=Rec709_c;
    }

    // chromaticities
    header_p+=storeStr32(header_p, "chromaticities");        // attribute name
    header_p+=storeStr32(header_p, "chromaticities");        // attribute type
    header_p+=storeI32LE(header_p, (int32_t)32);             // attribute value length in bytes
    header_p+=storeFloatLE(header_p, chromaticities.redX);   // redX
    header_p+=storeFloatLE(header_p, chromaticities.redY);   // redY
    header_p+=storeFloatLE(header_p, chromaticities.greenX); // greenX
    header_p+=storeFloatLE(header_p, chromaticities.greenY); // greenY
    header_p+=storeFloatLE(header_p, chromaticities.blueX);  // blueX
    header_p+=storeFloatLE(header_p, chromaticities.blueY);  // blueY
    header_p+=storeFloatLE(header_p, chromaticities.whiteX); // whiteX
    header_p+=storeFloatLE(header_p, chromaticities.whiteY); // whiteY
  }

  // null byte to signal end of attributes
  header_p+=storeU8(header_p, 0);

  // output header
  header_size=(int)(header_p - header);
  if (bsr_config->cgi_mode == 1) {
    fwrite(header, header_size, 1, stdout);
  } else {
    fwrite(header, header_size, 1, output_file);
  }

  // return number of bytes in header
  return(header_size);
}

int outputEXROffsetTable(bsr_config_t *bsr_config, bsr_state_t *bsr_state, FILE *output_file, int header_size) {
  int output_res_x;
  int output_res_y;
  uint64_t offset_table_size;
  uint64_t bytes_per_line=0;
  uint64_t chunk_start;
  uint64_t offset;
  unsigned char offset_buf[8];
  const int chunk_header_size=8; // 8 bytes: y coordinate + pixel_data_size
  int output_y;
 
  output_res_x=bsr_state->current_image_res_x;
  output_res_y=bsr_state->current_image_res_y;

  // bits per line
  if (bsr_config->bits_per_color == 8) {
    bytes_per_line=(uint64_t)(chunk_header_size + (3 * output_res_x));
  } else if (bsr_config->bits_per_color == 16) {
    bytes_per_line=(uint64_t)(chunk_header_size + (6 * output_res_x));
  } else if (bsr_config->bits_per_color == 32) {
    bytes_per_line=(uint64_t)(chunk_header_size + (12 * output_res_x));
  }

  // output entry in offset table
  offset_table_size=(uint64_t)(8 * output_res_y);
  chunk_start=(uint64_t)(header_size + offset_table_size); 
  offset=chunk_start;
  offset_buf[0]=0;
  for (output_y=0; output_y < output_res_y; output_y++) {
    storeU64LE(offset_buf, offset);
    if (bsr_config->cgi_mode == 1) {
      fwrite(offset_buf, 8, 1, stdout);
    } else {
      fwrite(offset_buf, 8, 1, output_file);
    }
    offset+=bytes_per_line;
  }

  return(0);
}

int outputEXRChunk(bsr_config_t *bsr_config, bsr_state_t *bsr_state, FILE *output_file) {
  int output_res_x;
  int output_res_y;
  int pixel_data_size=0;
  unsigned char pixel_data_size_buf[4];
  int output_y;
  unsigned char y_coordinate_buf[4];
  unsigned char *image_output_p;
  
  output_res_x=bsr_state->current_image_res_x;
  output_res_y=bsr_state->current_image_res_y;

  if (bsr_config->bits_per_color == 8) {
    pixel_data_size=(3 * output_res_x);
  } else if (bsr_config->bits_per_color == 16) {
    pixel_data_size=(6 * output_res_x);
  } else if (bsr_config->bits_per_color == 32) {
    pixel_data_size=(12 * output_res_x);
  }
  pixel_data_size_buf[0]=0;
  storeI32LE(pixel_data_size_buf, pixel_data_size);

  y_coordinate_buf[0]=0;
  image_output_p=bsr_state->image_output_buf;
  for (output_y=0; output_y < output_res_y; output_y++) {
    // y coordinate
    storeI32LE(y_coordinate_buf, output_y);
    if (bsr_config->cgi_mode == 1) {
      fwrite(y_coordinate_buf, 4, 1, stdout);
    } else {
      fwrite(y_coordinate_buf, 4, 1, output_file);
    }

    // pixel data size
    if (bsr_config->cgi_mode == 1) {
      fwrite(pixel_data_size_buf, 4, 1, stdout);
    } else {
      fwrite(pixel_data_size_buf, 4, 1, output_file);
    }

    // pixel data
    if (bsr_config->cgi_mode == 1) {
      fwrite(image_output_p, pixel_data_size, 1, stdout);
    } else {
      fwrite(image_output_p, pixel_data_size, 1, output_file);
    }

    image_output_p+=pixel_data_size;
  } // end for output_y

  return(0);
}

int outputEXR(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  FILE *output_file=NULL;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  int i;
  int header_size;

  //
  // main thread: display status update if not in CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->main_pid) && (bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Writing %s...", bsr_config->output_file_name);
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->main_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_IMAGE_OUTPUT_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_IMAGE_OUTPUT_BEGIN;
    }
  } // end if not main thread

//
// incomplete - work in progress
// multi-thread compression step goes here
//

  //
  // main thread: output EXR image to file or stdout 
  //
  if (bsr_state->perthread->my_pid == bsr_state->main_pid) {
    if (bsr_config->cgi_mode != 1) {
      output_file=fopen(bsr_config->output_file_name, "wb");
      if (output_file == NULL) {
        printf("Error: could not open %s for writing\n", bsr_config->output_file_name);
        fflush(stdout);
        exit(1);
      }
    } // end if not cgi_mode

    // output file components
    header_size=outputEXRHeader(bsr_config, bsr_state, output_file);
    outputEXROffsetTable(bsr_config, bsr_state, output_file, header_size);
    outputEXRChunk(bsr_config, bsr_state, output_file);

    if (bsr_config->cgi_mode != 1) {
      fclose(output_file);
    }
  } // end if main thread

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->main_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_IMAGE_OUTPUT_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_IMAGE_OUTPUT_CONTINUE);
  } else {
    waitForWorkerThreads(bsr_state, THREAD_STATUS_IMAGE_OUTPUT_COMPLETE);
    // ready to continue, set all worker thread status to continue
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_IMAGE_OUTPUT_CONTINUE;
    }
  } // end if not main thread

  //
  // main thread: display status message and close output file if not CGI mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->main_pid) && (bsr_config->cgi_mode != 1)) {
    if (bsr_config->print_status == 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.3fs)\n", elapsed_time);
      fflush(stdout);
    }
  }

  return(0);
}

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
#include "util.h"
#include "cgi.h"
#include "icc-profiles.h"

#ifdef BSR_USE_PNG
#define PNG_SETJMP_NOT_SUPPORTED
#include <png.h>
#endif

int outputPNG(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {

#ifdef BSR_USE_PNG

  FILE *output_file=NULL;
  png_structp png_ptr;
  png_infop info_ptr;
  unsigned char color_type=PNG_COLOR_TYPE_RGB;
  unsigned char bit_depth;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;

  //
  // display status update if not in CGI mode
  //
  if ((bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Writing %s...", bsr_config->output_file_name);
    fflush(stdout);
  }

  //
  // if not CGI mode, open output file
  //
  if (bsr_config->cgi_mode != 1) {
    output_file=fopen(bsr_config->output_file_name, "wb");
    if (output_file == NULL) {
      printf("Error: could not open %s for writing\n", bsr_config->output_file_name);
      fflush(stdout);
      exit(1);
    }
  }
  
  //
  // initialize PNG ptr and info_ptr
  //
  png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info_ptr=png_create_info_struct(png_ptr);
  if (bsr_config->cgi_mode == 1) {
    png_init_io(png_ptr, stdout);
  } else {
    png_init_io(png_ptr, output_file);
  } 
  if (bsr_config->bits_per_color == 16) {
    bit_depth=16;
  } else {
    bit_depth=8;
  }
  png_set_IHDR(png_ptr, info_ptr, bsr_state->current_image_res_x, bsr_state->current_image_res_y, bit_depth, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  
  //
  // set color profile
  //
  if (bsr_config->color_profile == 0) {
    // no ICC profile and linear gamma
    png_set_gAMA(png_ptr, info_ptr, 1.0);
  } else if (bsr_config->color_profile == 2) {
    // Display-P3
    png_set_iCCP(png_ptr, info_ptr, "Display-P3", 0, DisplayP3Compat_v4_icc, DisplayP3Compat_v4_icc_len);
  } else if (bsr_config->color_profile == 3) {
    // Rec. 2020
    png_set_iCCP(png_ptr, info_ptr, "Rec 2020", 0, Rec2020Compat_v4_icc, Rec2020Compat_v4_icc_len);
  } else if (bsr_config->color_profile == 4) {
    // Rec. 601 NTSC
    png_set_iCCP(png_ptr, info_ptr, "Rec 601 NTSC", 0, Rec601NTSC_v4_icc, Rec601NTSC_v4_icc_len);
  } else if (bsr_config->color_profile == 5) {
    // Rec. 601 PAL
    png_set_iCCP(png_ptr, info_ptr, "Re 601 PAL", 0, Rec601PAL_v4_icc, Rec601PAL_v4_icc_len);
  } else if (bsr_config->color_profile == 6) {
    // Rec. 709
    png_set_iCCP(png_ptr, info_ptr, "Rec 709", 0, Rec709_v4_icc, Rec709_v4_icc_len);
  } else if (bsr_config->color_profile == 7) {
    // no ICC profile and flat 2.0 gamma
    png_set_gAMA(png_ptr, info_ptr, 0.5);
  } else if (bsr_config->color_profile == 8) {
    // Rec. 2100 PQ
    png_set_iCCP(png_ptr, info_ptr, "Rec 2100 PQ", 0, Rec2100PQ_v4_icc, Rec2100PQ_v4_icc_len);
  } else {
    // default is sRGB
    png_set_iCCP(png_ptr, info_ptr, "sRGB", 0, sRGB_v4_icc, sRGB_v4_icc_len);
  }

  //
  // write PNG header and image data
  //
  png_write_info(png_ptr, info_ptr);
  png_write_image(png_ptr, bsr_state->row_pointers);
  png_write_end(png_ptr, NULL);

  //
  // display status message and close output file if not CGI mode
  //
  if (bsr_config->cgi_mode != 1) {
    if (bsr_config->print_status == 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.3fs)\n", elapsed_time);
      fflush(stdout);
    }

    // close output file
    fclose(output_file);
  }

#endif // BSR_USE_PNG

  return(0);
}

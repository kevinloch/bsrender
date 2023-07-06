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

#ifdef BSR_USE_JPEG
#include <jpeglib.h>
#include <jerror.h>
#endif

int outputJpeg(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {

#ifdef BSR_USE_JPEG

  FILE *output_file=NULL;
  struct jpeg_compress_struct jpeg_info;
  struct jpeg_error_mgr jpeg_err;
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
  // initialize jpeg_info
  //
  jpeg_info.err=jpeg_std_error(&jpeg_err);
  jpeg_create_compress(&jpeg_info);
  if (bsr_config->cgi_mode == 1) {
    jpeg_stdio_dest(&jpeg_info, stdout);
  } else {
    jpeg_stdio_dest(&jpeg_info, output_file);
  }
  jpeg_info.image_width=bsr_state->current_image_res_x;
  jpeg_info.image_height=bsr_state->current_image_res_y;
  jpeg_info.input_components=3;
  jpeg_info.in_color_space=JCS_RGB;
  jpeg_set_defaults(&jpeg_info);
  jpeg_set_quality(&jpeg_info, bsr_config->compression_quality, 1);
  jpeg_start_compress(&jpeg_info, 1);

  //
  // set color profile
  //
  if (bsr_config->color_profile == 2) {
    // Display-P3
    jpeg_write_icc_profile(&jpeg_info, DisplayP3Compat_v4_icc, DisplayP3Compat_v4_icc_len);
  } else if (bsr_config->color_profile == 3) {
    // Rec. 2020
    jpeg_write_icc_profile(&jpeg_info, Rec2020Compat_v4_icc, Rec2020Compat_v4_icc_len);
  } else if (bsr_config->color_profile == 4) {
    // Rec. 601 NTSC
    jpeg_write_icc_profile(&jpeg_info, Rec601NTSC_v4_icc, Rec601NTSC_v4_icc_len);
  } else if (bsr_config->color_profile == 5) {
    // Rec. 601 PAL
    jpeg_write_icc_profile(&jpeg_info, Rec601PAL_v4_icc, Rec601PAL_v4_icc_len);
  } else if (bsr_config->color_profile == 6) {
    // Rec. 709
    jpeg_write_icc_profile(&jpeg_info, Rec709_v4_icc, Rec709_v4_icc_len);
/*
  } else if (bsr_config->color_profile == 8) {
    // Rec. 2100 PQ
    jpeg_write_icc_profile(&jpeg_info, Rec2100PQ_v4_icc, Rec2100PQ_v4_icc_len);
*/
  } else {
    // default is sRGB
    jpeg_write_icc_profile(&jpeg_info, sRGB_v4_icc, sRGB_v4_icc_len);
  }

  //
  // compress and output jpeg
  //
  jpeg_write_scanlines(&jpeg_info, bsr_state->row_pointers, bsr_state->current_image_res_y);
  jpeg_finish_compress(&jpeg_info);

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

  // clean up libjpeg
  jpeg_destroy_compress(&jpeg_info);

#endif // BSR_USE_JPEG

  return(0);
}

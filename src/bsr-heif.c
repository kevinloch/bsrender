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
#include <time.h>
#include "util.h"
#include "cgi.h"
#include "icc-profiles.h"

#ifdef BSR_USE_HEIF
#include <libheif/heif.h>
#endif

int outputHeif(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {

#ifdef BSR_USE_HEIF

  FILE *output_file=NULL;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  struct heif_context *h_ctx;
  struct heif_encoder *h_encoder;
  struct heif_image *h_image;
  int h_chroma=heif_chroma_interleaved_RGB;
  struct heif_error h_error;
  uint8_t *h_image_plane_p;
  int bytes_per_pixel=0;
  size_t current_image_size;
  struct heif_color_profile_nclx h_nclx;
  struct heif_image_handle *h_out_image_handle;

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
  // initialize encoder
  //
  h_ctx=heif_context_alloc();
  heif_context_get_encoder_for_format(h_ctx, heif_compression_HEVC, &h_encoder);
  heif_encoder_set_lossy_quality(h_encoder, bsr_config->compression_quality);
  if (bsr_config->bits_per_color == 8) {
    h_chroma=heif_chroma_interleaved_RGB;
  } else if ((bsr_config->bits_per_color == 10) || (bsr_config->bits_per_color == 12)) {
    h_chroma=heif_chroma_interleaved_RRGGBB_LE;
  }
  h_error=heif_image_create(bsr_state->current_image_res_x, bsr_state->current_image_res_y, heif_colorspace_RGB, h_chroma, &h_image);
  if ((h_error.code) != heif_error_Ok) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: heif_image_create() failed with code %d.%d: %s\n", h_error.code, h_error.subcode, h_error.message);
    }
    exit(1);
  }
  h_error=heif_image_add_plane(h_image, heif_channel_interleaved, bsr_state->current_image_res_x, bsr_state->current_image_res_y, bsr_config->bits_per_color);
  if (h_error.code != heif_error_Ok) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: heif_add_plane() failed with code %d.%d: %s\n", h_error.code, h_error.subcode, h_error.message);
    }
    exit(1);
  }

  //
  // set color profile
  //
  if (bsr_config->color_profile == 0) {
    // Linear transfer function
    h_nclx.color_primaries=heif_color_primaries_ITU_R_BT_709_5;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_linear;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_ITU_R_BT_709_5;
  } else if (bsr_config->color_profile == 2) {
    // Display-P3
    h_nclx.color_primaries=heif_color_primaries_SMPTE_EG_432_1;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_IEC_61966_2_4;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_ITU_R_BT_709_5;
  } else if (bsr_config->color_profile == 3) {
    // Rec. 2020
    h_nclx.color_primaries=heif_color_primaries_ITU_R_BT_2020_2_and_2100_0;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_ITU_R_BT_2020_2_10bit;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_ITU_R_BT_2020_2_non_constant_luminance;
  } else if (bsr_config->color_profile == 4) {
    // Rec. 601 NTSC
    h_nclx.color_primaries=heif_color_primaries_ITU_R_BT_601_6;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_ITU_R_BT_601_6;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_ITU_R_BT_601_6;
  } else if (bsr_config->color_profile == 5) {
    // Rec. 601 PAL
    h_nclx.color_primaries=heif_color_primaries_ITU_R_BT_601_6;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_ITU_R_BT_601_6;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_ITU_R_BT_601_6;
  } else if (bsr_config->color_profile == 6) {
    // Rec. 709
    h_nclx.color_primaries=heif_color_primaries_ITU_R_BT_709_5;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_ITU_R_BT_709_5;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_ITU_R_BT_709_5;
  } else if (bsr_config->color_profile == 8) {
    // Rec. 2100 PQ
    h_nclx.color_primaries=heif_color_primaries_ITU_R_BT_2020_2_and_2100_0;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_ITU_R_BT_2100_0_PQ;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_SMPTE_ST_2085;
  } else {
    // default is sRGB
    h_nclx.color_primaries=heif_color_primaries_ITU_R_BT_709_5;
    h_nclx.transfer_characteristics=heif_transfer_characteristic_IEC_61966_2_4;
    h_nclx.matrix_coefficients=heif_matrix_coefficients_ITU_R_BT_709_5;
  }
  h_error=heif_image_set_nclx_color_profile(h_image, &h_nclx);
  if (h_error.code != heif_error_Ok) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: heif_image_set_nclx_color_profile() failed with code %d.%d: %s\n", h_error.code, h_error.subcode, h_error.message);
    }
    exit(1);
  }

  //
  // copy imgae data to h_image 
  //
  // Note to libheif maintainers: Why can't we just set a pointer to where the image data already is?
  // Having to needlessly allocate memory and copy a large (possibly > 1TB) image is sub-optimal
  //
  h_image_plane_p=heif_image_get_plane(h_image, heif_channel_interleaved, NULL);
  if (h_image_plane_p == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: heif_image_get_plane() returned NULL");
    }
    exit(1);
  }
  if (bsr_config->bits_per_color == 8) {
    bytes_per_pixel=3;
  } else if ((bsr_config->bits_per_color == 10) || (bsr_config->bits_per_color == 12)) {
    bytes_per_pixel=6;
  }
  current_image_size=(size_t)bytes_per_pixel * (size_t)bsr_state->current_image_res_x * (size_t)bsr_state->current_image_res_y;
  memcpy(h_image_plane_p, bsr_state->image_output_buf, current_image_size);

  //
  // compress image
  //
  h_error=heif_context_encode_image(h_ctx, h_image, h_encoder, NULL, &h_out_image_handle);
  if (h_error.code != heif_error_Ok) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: heif_context_encode_image() failed with code %d.%d: %s\n", h_error.code, h_error.subcode, h_error.message);
    }
    exit(1);
  }
  heif_encoder_release(h_encoder);

  //
  // output image data
  //
  if (bsr_config->cgi_mode != 1) {
    h_error=heif_context_write_to_file(h_ctx, bsr_config->output_file_name);
/*
  } else { // output to stdout is not supported by libheif?!
*/
  }
  if (h_error.code != heif_error_Ok) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: heif_context_write_to_file() failed with code %d.%d: %s\n", h_error.code, h_error.subcode, h_error.message);
    }
    exit(1);
  }

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
//    fclose(output_file);
  }

  // clean up
  heif_context_free(h_ctx);

#endif // BSR_USE_HEIF

  return(0);
}

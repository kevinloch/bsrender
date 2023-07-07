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

#ifdef BSR_USE_AVIF
#include <avif/avif.h>
#endif

int outputAvif(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {

#ifdef BSR_USE_AVIF

  FILE *output_file=NULL;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  avifEncoder *avif_encoder;
  avifRWData avif_output=AVIF_DATA_EMPTY;
  avifRGBImage avif_rgb;
  avifImage *avif_image;
  int bytes_per_pixel=3;
  avifResult avif_result;
  int avif_quantizer;

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
  // initialize avif_image
  //
  avif_image=avifImageCreate(bsr_state->current_image_res_x, bsr_state->current_image_res_y, bsr_config->bits_per_color, AVIF_PIXEL_FORMAT_YUV444);
  if (avif_image == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: avifImageCreate() failed (invalid arguments or memory allocation failed)\n");
    }
    exit(1);
  }

  //
  // set color profile
  //
  if (bsr_config->color_profile == 0) {
    // Linear transfer function
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_IEC61966_2_4;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_LINEAR;
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT709;
  } else if (bsr_config->color_profile == 2) {
    // Display-P3
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_SMPTE432;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_SRGB;
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT709;
  } else if (bsr_config->color_profile == 3) {
    // Rec. 2020
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_BT2020;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_BT2020_10BIT;
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT2020_NCL;
  } else if (bsr_config->color_profile == 4) {
    // Rec. 601 NTSC
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_BT601;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_BT601;
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT601;
  } else if (bsr_config->color_profile == 5) {
    // Rec. 601 PAL
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_BT601;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_BT601;
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT601;
  } else if (bsr_config->color_profile == 6) {
    // Rec. 709
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_BT709;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_BT709;
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT709;
  } else if (bsr_config->color_profile == 8) {
    // Rec. 2100 PQ
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_BT2020;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_SMPTE2084;
    //avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_SMPTE2085; // fails with "no content" error
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT709;
  } else {
    // default is sRGB
    avif_image->colorPrimaries=AVIF_COLOR_PRIMARIES_IEC61966_2_4;
    avif_image->matrixCoefficients=AVIF_MATRIX_COEFFICIENTS_BT709;
    avif_image->transferCharacteristics=AVIF_TRANSFER_CHARACTERISTICS_SRGB;
  }

  //
  // convert RGB to YUV
  //
  memset(&avif_rgb, 0, sizeof(avif_rgb));
  avifRGBImageSetDefaults(&avif_rgb, avif_image);
  avif_rgb.width=bsr_state->current_image_res_x;
  avif_rgb.height=bsr_state->current_image_res_y;
  avif_rgb.depth=bsr_config->bits_per_color;
  avif_rgb.format=AVIF_RGB_FORMAT_RGB;
  avif_rgb.pixels=bsr_state->image_output_buf;
  if (bsr_config->bits_per_color == 8) {
    bytes_per_pixel=3;
  } else if ((bsr_config->bits_per_color == 10) || (bsr_config->bits_per_color == 12) || (bsr_config->bits_per_color == 16)) {
    bytes_per_pixel=6;
    if (bsr_config->bits_per_color == 16) {
      avif_rgb.isFloat=AVIF_TRUE;
    }       
  }
  avif_rgb.rowBytes=bytes_per_pixel * bsr_state->current_image_res_x;
  avifImageRGBToYUV(avif_image, &avif_rgb);

  //
  // compress image
  //
  avif_encoder=avifEncoderCreate();
  if (avif_image == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: avifEncoderCreate() failed (memory allocation failed)\n");
    }
    exit(1);
  }
  //  avif_encoder->quality=bsr_config->compression_quality;
  // for backwards compatibility with older versions that do not take quality directly, set min/max quantizer level
  avif_quantizer=63 - (int)(((float)bsr_config->compression_quality * 0.63f) + 0.5f);
  if (avif_quantizer < 0) {
    avif_quantizer=0;
  } else if (avif_quantizer > 63) {
    avif_quantizer=63;
  }
  avif_encoder->minQuantizer=avif_quantizer;
  avif_encoder->maxQuantizer=avif_quantizer;
  avif_encoder->maxThreads=bsr_config->num_threads;
  avif_result=avifEncoderAddImage(avif_encoder, avif_image, 1, AVIF_ADD_IMAGE_FLAG_SINGLE);
  if (avif_result != AVIF_RESULT_OK) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: avifEncoderAddImage() failed: %s\n", avifResultToString(avif_result));
    }
    exit(1);
  }
  avif_result=avifEncoderFinish(avif_encoder, &avif_output);
  if (avif_result != AVIF_RESULT_OK) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: avifEncoderFinish() failed: %s\n", avifResultToString(avif_result));
    }
    exit(1);
  }


  //
  // output image data
  //
  if (bsr_config->cgi_mode != 1) {
    fwrite(avif_output.data, 1, avif_output.size, output_file);
  } else {
    fwrite(avif_output.data, 1, avif_output.size, stdout);
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
    fclose(output_file);
  }

  // clean up
  avifImageDestroy(avif_image);
  avifEncoderDestroy(avif_encoder);

#endif // BSR_USE_AVIF

  return(0);
}

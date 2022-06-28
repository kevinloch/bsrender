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

int drawCrossHairs(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  int i;
  pixel_composition_t *image_composition_p;
  double res_x;
  double res_y;
  double half_res_x;
  double half_res_y;
  pixel_composition_t *image_buf;

  //
  // get current image buffer and resolution
  //
  image_buf=bsr_state->current_image_buf;
  res_x=(double)bsr_state->current_image_res_x;
  res_y=(double)bsr_state->current_image_res_y;
  half_res_x=res_x / 2.0;
  half_res_y=res_y / 2.0;

  //
  // draw crosshairs
  //
  for (i=(half_res_x - (res_y * 0.02)); i < (half_res_x - (res_y * 0.005)); i++) {
    image_composition_p=image_buf + ((int)res_x * (int)half_res_y) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(half_res_x + (res_y * 0.005)); i < (half_res_x + (res_y * 0.02)); i++) {
    image_composition_p=image_buf + ((int)res_x * (int)half_res_y) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(half_res_y - (res_y * 0.02)); i < (half_res_y - (res_y * 0.005)); i++) {
    image_composition_p=image_buf + ((int)res_x * i) + (int)half_res_x;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(half_res_y + (res_y * 0.005)); i < (half_res_y + (res_y * 0.02)); i++) {
    image_composition_p=image_buf + ((int)res_x * i) + (int)half_res_x;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  return(0);
}

int drawGridLines(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  int i;
  pixel_composition_t *image_composition_p;
  double res_x;
  double res_y;
  double half_res_x;
  double half_res_y;
  pixel_composition_t *image_buf;

  //
  // get current image buffer and resolution
  //
  image_buf=bsr_state->current_image_buf;
  res_x=(double)bsr_state->current_image_res_x;
  res_y=(double)bsr_state->current_image_res_y;
  half_res_x=res_x / 2.0;
  half_res_y=res_y / 2.0;

  //
  // draw select raster lines
  //
  for (i=0; i < res_x; i++) {
    image_composition_p=image_buf + (int)res_x * (int)(res_y * 0.25) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_x; i++) {
    image_composition_p=image_buf + (int)res_x * (int)half_res_y + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_x; i++) {
    image_composition_p=image_buf + (int)res_x * (int)(res_y * 0.75) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_y; i++) {
    image_composition_p=image_buf + (int)res_x * i + (int)(res_x * 0.25);
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_y; i++) {
    image_composition_p=image_buf + (int)res_x * i + (int)(half_res_x);
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_y; i++) {
    image_composition_p=image_buf + (int)res_x * i + (int)(res_x * 0.75);
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  return(0);
}

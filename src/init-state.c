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
#include <sys/mman.h>
#include <math.h>

bsr_state_t *initState(bsr_config_t *bsr_config) {
  bsr_state_t *bsr_state;
  double pi_over_360=M_PI / 360.0;
  double pi_over_180=M_PI / 180.0;
  double camera_icrs_ra_rad;
  double camera_icrs_dec_rad;
  double target_icrs_ra_rad;
  double target_icrs_dec_rad;

  int mmap_protection;
  int mmap_visibility;
  size_t bsr_state_size=0;

  //
  // allocate shared memory for bsr_state
  //
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
  bsr_state_size=sizeof(bsr_state_t);
  bsr_state=(bsr_state_t *)mmap(NULL, bsr_state_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate shared memory for bsr_state\n");
    }
    return(NULL);
  }

  //
  // calculate number of rendering threads to be forked
  //
  bsr_state->num_worker_threads=bsr_config->num_threads-1;
  if (bsr_state->num_worker_threads < 1) {
    bsr_state->num_worker_threads=1;
  }

  //
  // process user-supplied arguments
  //
  bsr_state->camera_hfov=bsr_config->camera_fov * pi_over_360; // includes divide by 2
  bsr_state->camera_half_res_x=(double)bsr_config->camera_res_x / 2.0;
  bsr_state->camera_half_res_y=(double)bsr_config->camera_res_y / 2.0;
  bsr_state->pixels_per_radian=bsr_state->camera_half_res_x / bsr_state->camera_hfov;
  bsr_state->camera_3az_yz=bsr_config->camera_rotation * pi_over_180;
  bsr_state->camera_3az_xy=bsr_config->camera_pan * pi_over_180;
  bsr_state->camera_3az_xz=bsr_config->camera_tilt * -pi_over_180;

  //
  // select per thread buffer size
  //
  if (bsr_config->Airy_disk == 1) {
    bsr_state->per_thread_buffers = bsr_config->per_thread_buffer_Airy;
  } else {
    bsr_state->per_thread_buffers = bsr_config->per_thread_buffer;
  }

  //
  // optionally transform spherical icrs to euclidian icrs if x,y,z are zero
  //
  if (((bsr_config->camera_icrs_ra != 0.0) || (bsr_config->camera_icrs_dec != 0.0) || (bsr_config->camera_icrs_r != 0.0))\
    && (bsr_config->camera_icrs_x == 0.0) && (bsr_config->camera_icrs_y == 0.0) && (bsr_config->camera_icrs_z == 0.0)) {
    camera_icrs_ra_rad=bsr_config->camera_icrs_ra * pi_over_180;
    camera_icrs_dec_rad=bsr_config->camera_icrs_dec * pi_over_180;
    bsr_config->camera_icrs_x=bsr_config->camera_icrs_r * cos(camera_icrs_dec_rad) * cos(camera_icrs_ra_rad);
    bsr_config->camera_icrs_y=bsr_config->camera_icrs_r * cos(camera_icrs_dec_rad) * sin(camera_icrs_ra_rad);
    bsr_config->camera_icrs_z=bsr_config->camera_icrs_r * sin(camera_icrs_dec_rad);
  }
  if (((bsr_config->target_icrs_ra != 0.0) || (bsr_config->target_icrs_dec != 0.0) || (bsr_config->target_icrs_r != 0.0))\
    && (bsr_config->target_icrs_x == 0.0) && (bsr_config->target_icrs_y == 0.0) && (bsr_config->target_icrs_z == 0.0)) {
    target_icrs_ra_rad=bsr_config->target_icrs_ra * pi_over_180;
    target_icrs_dec_rad=bsr_config->target_icrs_dec * pi_over_180;
    bsr_config->target_icrs_x=bsr_config->target_icrs_r * cos(target_icrs_dec_rad) * cos(target_icrs_ra_rad);
    bsr_config->target_icrs_y=bsr_config->target_icrs_r * cos(target_icrs_dec_rad) * sin(target_icrs_ra_rad);
    bsr_config->target_icrs_z=bsr_config->target_icrs_r * sin(target_icrs_dec_rad);
  }

  //
  // translate original target x,y,z to new coordinates as seen by camera position
  //
  bsr_state->target_x=bsr_config->target_icrs_x - bsr_config->camera_icrs_x;
  bsr_state->target_y=bsr_config->target_icrs_y - bsr_config->camera_icrs_y;
  bsr_state->target_z=bsr_config->target_icrs_z - bsr_config->camera_icrs_z;

  //
  // initialize target xy angle used in star rotations
  //
  bsr_state->target_3az_xy=atan2(bsr_state->target_y, bsr_state->target_x);

  //
  // initialize target xz angle used in star rotations by setting target xy angle to 0
  //
  bsr_state->target_3az_xy_r=sqrt((bsr_state->target_x * bsr_state->target_x) + (bsr_state->target_y * bsr_state->target_y)); 
  bsr_state->target_x=bsr_state->target_3az_xy_r; // xy=0
  bsr_state->target_3az_xz=atan2(bsr_state->target_z, bsr_state->target_x);

  return(bsr_state);
}

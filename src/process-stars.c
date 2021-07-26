
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
#include <math.h>
#include "util.h"

//
// note: all functions in this file should remain in this file for compiler optimization
//

int sendPixelToMainThread(bsr_state_t *bsr_state, long long image_offset, double r, double g, double b) {
  //
  // This function inserts a single pixel (image_offset, r, g, b) into this thread's section of the main threaad buffer.
  // If the current slot of the main thread buffer is occupied it will wait until it is free, periodically checking if
  // the main thread is still alive.
  //

  volatile int success; // gcc optimization breaks algorithm without volatile keyword
  int idle_count;

  //
  // rewind if we have reached the end of our section of the main thread buffer
  //
  if (bsr_state->perthread->thread_buffer_index == bsr_state->per_thread_buffers) {
    // end of our section of thread_buf, rewind
    bsr_state->perthread->thread_buffer_index=0;
    bsr_state->perthread->thread_buf_p-=bsr_state->per_thread_buffers;
  }

  //
  // Attempt to insert into main thread buffer.  If current slot has not been cleared by main thread, wait until it is
  // clear.  Periodically check if main thread is still alive.
  //
  success=0;
  idle_count=0;
  while (success == 0) {
    // check if we've written a pixel to this location before and it has not been read/cleared by main thread yet
    if ((bsr_state->perthread->thread_buf_p->status_left == 0) && (bsr_state->perthread->thread_buf_p->status_right == 0)) {
      bsr_state->perthread->thread_buf_p->status_left=1;
      bsr_state->perthread->thread_buf_p->image_offset=image_offset;
      bsr_state->perthread->thread_buf_p->r=r;
      bsr_state->perthread->thread_buf_p->g=g;
      bsr_state->perthread->thread_buf_p->b=b;
      bsr_state->perthread->thread_buf_p->status_right=1;
      bsr_state->perthread->thread_buf_p++;
      bsr_state->perthread->thread_buffer_index++;
      success=1;
    } else {
      idle_count++;
      if (idle_count > 10000) {
        // check if main thread is still alive
        checkExceptions(bsr_state);
        idle_count=0;
      }
    } // end if buffer slot is available
  } // end while success=0

  return(0);
}

int sendDedupBufferToMainThread(bsr_state_t *bsr_state) {
  //
  // This function sends pixels from the dedup buffer to the main thread.
  // As pixels are sent it clears them from the dedup buffer and dedup buffer index.
  //

  dedup_buffer_t *dedup_buf_p;
  dedup_index_t *dedup_index_p;
  int dedup_buf_i;
  long long dedup_index_offset;

  dedup_buf_p=bsr_state->dedup_buf; // start at beginning of dedup buffer
  for (dedup_buf_i=0; dedup_buf_i < bsr_state->perthread->dedup_count; dedup_buf_i++) {
    //
    // send to main thread
    //
    sendPixelToMainThread(bsr_state, dedup_buf_p->image_offset, dedup_buf_p->r, dedup_buf_p->g, dedup_buf_p->b);

    //
    // set dedup index depending on dedup index mode
    //
    if (bsr_state->dedup_index_mode == 0) {
      dedup_index_offset=dedup_buf_p->image_offset;
    } else {
      dedup_index_offset=(dedup_buf_p->image_offset & 0xffffff); // truncate image_offset to 24 bits
    }
    dedup_index_p=bsr_state->dedup_index + dedup_index_offset;

    // 
    // clear dedup buffer index record
    //
    if (dedup_index_p->dedup_record_p == dedup_buf_p) {
      dedup_index_p->dedup_record_p=NULL;
    }

    //
    // clear dedup buffer record
    //
    dedup_buf_p->image_offset=-1;
    dedup_buf_p->r=0.0;
    dedup_buf_p->g=0.0;
    dedup_buf_p->b=0.0;
    dedup_buf_p++;
  } // end for dedup_buf_i

  //
  // set dedup_count to 0 indicating dedup buffer is empty
  //
  bsr_state->perthread->dedup_count=0;

  return(0);
}

int sendPixelToDedupBuffer(bsr_state_t *bsr_state, long long image_offset, double r, double g, double b) {
  //
  // This function attempts to insert a pixel into the dedup buffer.  If a record for this image_offset does not
  // exist yet, then the insert is made.  If a record for this image_offset already exists then the record is
  // updated with the new value added to the existing value.  If a dedup index collision is detected then the
  // pixel is sent directly to the main thread instead.  Finally, it checks if dedup buffer is full and if so
  // sends dedup buffer contents to main thread.
  //

  dedup_buffer_t *dedup_buf_p;
  dedup_index_t *dedup_index_p;
  long long dedup_index_offset;

  //
  // set dedup index depending on dedup index mode
  //
  if (bsr_state->dedup_index_mode == 0) {
    dedup_index_offset=(int)image_offset;
  } else {
    dedup_index_offset=(int)(image_offset & 0xffffff); // truncate image_offset to 24 bits
  }
  dedup_index_p=bsr_state->dedup_index + dedup_index_offset;

  //
  // try to insert pixel into dedup buffer
  //
  if (dedup_index_p->dedup_record_p == NULL) {
    // no dup yet, just store value in next dedup buffer record and update index
    dedup_buf_p=bsr_state->dedup_buf + bsr_state->perthread->dedup_count;
    bsr_state->perthread->dedup_count++;
    dedup_buf_p->image_offset=image_offset;
    dedup_buf_p->r=r;
    dedup_buf_p->g=g;
    dedup_buf_p->b=b;
    dedup_index_p->dedup_record_p=dedup_buf_p;
  } else if (dedup_index_p->dedup_record_p->image_offset == image_offset) {
    // duplicate pixel locaiton, add to existing dedup buffer record values
    dedup_buf_p=dedup_index_p->dedup_record_p;
    dedup_buf_p->r+=r;
    dedup_buf_p->g+=g;
    dedup_buf_p->b+=b;
  } else {
    // dedup index collision, send this pixel directly to main thread
    sendPixelToMainThread(bsr_state, image_offset, r, g, b);
  } // end if dedup_index_p->dedup_record_p

  //
  // check if dedup buffer is full and if it is send pixels to main thread
  //
  if (bsr_state->perthread->dedup_count == bsr_state->per_thread_buffers) {
    sendDedupBufferToMainThread(bsr_state);
  } // end if dedup buffer full

  return(0);
}

int processStars(bsr_config_t *bsr_config, bsr_state_t *bsr_state, FILE *input_file) {
  //
  // This function handles the most expensive operations in bsrender.  It performs the following:
  //
  // - reads stars from the supplied input file
  // - filters stars by distance from target or camera, and color temperature
  // - translates position relative to camera position
  // - rotates stars to center on target (and optional pan/tilt away from target)
  // - applies selected raster projection
  // - optionally maps Airy disk pixels
  // - deduplicates pixels to reduce load on memory bandwidth, which is typically the limiting performance factor on large servers with many cpus
  // - sends pixels to main thread for integration into the image composition buffer
  //

  int i;
  star_record_t star_record;
  int star_record_size=sizeof(star_record_t);
  double star_x;
  double star_y;
  double star_z;
  double star_r;
  double star_r2; // squared
  double star_3az_xy_r;
  double star_3az_xz_r;
  double star_3az_yz_r;
  double render_distance2; // distance from selected point to star (squared)
  double star_3az_xy;
  double star_3az_xz;
  double star_3az_yz;
  uint64_t color_temperature;
  float star_linear_1pc_intensity;
  double star_linear_intensity; // star linear intensity as viewed from camera
  double output_az;
  double output_az_by2;
  double output_el;
  int output_x=0;
  int output_y=0;
  double spherical_distance;
  double spherical_angle;
  double two_mollewide_angle;
  double mollewide_angle;
  const double pi_over_2=M_PI / 2.0;
  int Airymap_x;
  int Airymap_y;
  int Airymap_half_xy;
  int Airymap_xy;
  int Airymap_output_x;
  int Airymap_output_y;
  double *Airymap_red_p;
  double *Airymap_green_p;
  double *Airymap_blue_p;
  double star_rgb_red;
  double star_rgb_green;
  double star_rgb_blue;
  long long image_offset;
  double r;
  double g;
  double b;

  //
  // init shortcut variables
  //
  Airymap_half_xy=bsr_config->Airy_disk_max_extent;
  Airymap_xy=(Airymap_half_xy * 2) + 1;

  //
  // reset dedup buffer (probably not needed as it is cleared each time dedup buffer is emptied)
  //
  //bsr_state->perthread->dedup_count=0;

  //
  // read and process each line of input file
  // 
  fread(&star_record, star_record_size, 1, input_file);
  while ((feof(input_file) == 0) && (ferror(input_file) == 0)) {

    //
    // translate original star x,y,z to new coordinates as seen by camera position
    //
    star_x=star_record.icrs_x - bsr_config->camera_icrs_x;
    star_y=star_record.icrs_y - bsr_config->camera_icrs_y;
    star_z=star_record.icrs_z - bsr_config->camera_icrs_z;
    star_r2=(star_x * star_x) + (star_y * star_y) + (star_z * star_z); // leave squared for now for better performance

    //
    // extract color_temperature from combined field
    //
    color_temperature=star_record.intensity_and_temperature & 0x00000000fffffffful;

    // 
    // only continue if star distance is within min-max range from selected point (0=camera, 1=target),
    // greater than zero, and color temperature is within allowed range
    //
    if (bsr_config->render_distance_selector == 0) { // selected point is camera
      render_distance2=star_r2; // star distance from camera
    } else { // selected point is target
      // leave squared
      render_distance2=((star_record.icrs_x - bsr_config->target_icrs_x) * (star_record.icrs_x - bsr_config->target_icrs_x))\
                     + ((star_record.icrs_y - bsr_config->target_icrs_y) * (star_record.icrs_y - bsr_config->target_icrs_y))\
                     + ((star_record.icrs_z - bsr_config->target_icrs_z) * (star_record.icrs_z - bsr_config->target_icrs_z)); // important, use un-translated/rotated coordinates
    }
    if ((star_r2 > 0.0)\
     && (render_distance2 >= bsr_config->render_distance_min2) && (render_distance2 <= bsr_config->render_distance_max2)\
     && (color_temperature >= bsr_config->star_color_min) && (color_temperature <= bsr_config->star_color_max)) {

      //
      // extract star intensity from combined field and adjust for distance from camera
      //
      star_record.intensity_and_temperature >>=32;
      star_linear_1pc_intensity=*((float*)&star_record.intensity_and_temperature);
      star_linear_intensity=star_linear_1pc_intensity / star_r2;

      //
      // Note: 3D rotations for the camera, target, and stars are done with a 'triple-azimuth' (3az) method.  
      // Three full 360 degree orthogonal angles (x->y), (x->z), and (y->z) provide redundancy for constant precision and eliminates gimbal lock.
      // All angles range from -pi to +pi.  From the camera's perspective +x=forward, +y=left, and +z=up.  Thus the camera is
      // aimed at the target when the target is at 3az_xy=0, 3az_xz=0. Camera rotation (about the x axis) is controlled by the 3az_yz angle.
      // Optional after-aim pan and tilt is done by rotating the 3az_xy and 3az_xz angles respectively away from the target.
      //

      //
      // rotate star xy angle by target xy angle
      //
      // initialize 3az_xy and 3az_xy_r
      star_3az_xy=atan2(star_y, star_x);
      star_3az_xy_r=sqrt((star_x * star_x) + (star_y * star_y));
      // rotate
      star_3az_xy-=bsr_state->target_3az_xy;
      // update x and y
      star_x=star_3az_xy_r * cos(star_3az_xy);
      star_y=star_3az_xy_r * sin(star_3az_xy);

      //
      // rotate star xz angle by (rotated) target xz angle
      //
      // initialize 3az_xz and 3az_xz_r
      star_3az_xz=atan2(star_z, star_x);
      star_3az_xz_r=sqrt((star_x * star_x) + (star_z * star_z));
      // rotate
      star_3az_xz-=bsr_state->target_3az_xz;
      // update x and z
      star_x=star_3az_xz_r * cos(star_3az_xz);
      star_z=star_3az_xz_r * sin(star_3az_xz);

      //
      // rotate star yz angle by camera rotation angle
      //
      // initialize 3az_xy and 3az_xy_r
      star_3az_yz=atan2(star_z, star_y);
      star_3az_yz_r=sqrt((star_y * star_y) + (star_z * star_z));
      // rotate
      star_3az_yz+=bsr_state->camera_3az_yz;
      // update y and z
      star_y=star_3az_yz_r * cos(star_3az_yz);
      star_z=star_3az_yz_r * sin(star_3az_yz);

      //
      // optionally pan camera left/right (xy angle)
      //
      if (bsr_config->camera_pan != 0.0) {
        // initialize 3az_xy and 3az_xy_r
        star_3az_xy=atan2(star_y, star_x);
        star_3az_xy_r=sqrt((star_x * star_x) + (star_y * star_y));
        // rotate
        star_3az_xy+=bsr_state->camera_3az_xy;
        // update x and y
        star_x=star_3az_xy_r * cos(star_3az_xy);
        star_y=star_3az_xy_r * sin(star_3az_xy);
      }

      //
      // optionally tilt camera up/down (xz_angle)
      //
      if (bsr_config->camera_tilt != 0.0) {
        // initialize 3az_xz and 3az_xz_r
        star_3az_xz=atan2(star_z, star_x);
        star_3az_xz_r=sqrt((star_x * star_x) + (star_z * star_z));
        // rotate
        star_3az_xz+=bsr_state->camera_3az_xz;
        // update x and z
        star_x=star_3az_xz_r * cos(star_3az_xz);
        star_z=star_3az_xz_r * sin(star_3az_xz);
      }

      //
      // project star onto output raster x,y
      //
      if (bsr_config->camera_projection == 0) {
        // lat/lon 
        star_r=sqrt(star_r2);
        star_3az_xy=atan2(star_y, star_x);
        output_az=star_3az_xy;
        output_el=asin(star_z / star_r);
        output_x=(int)((-bsr_state->pixels_per_radian * output_az) + bsr_state->camera_half_res_x - 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * output_el) + bsr_state->camera_half_res_y - 0.5);
      } else if (bsr_config->camera_projection == 1) {
        // spherical
        star_r=sqrt(star_r2);
        star_3az_yz_r=sqrt((star_y * star_y) + (star_z * star_z));
        spherical_distance=asin(star_3az_yz_r / star_r);
        star_3az_yz=atan2(star_z, star_y);
        spherical_angle=star_3az_yz;
        output_az=spherical_distance * cos(spherical_angle);
        output_el=spherical_distance * sin(spherical_angle);
        if (bsr_config->spherical_orientation == 1) { // side by side orientation
          if (star_x > 0.0) { // star is in front of camera, draw on left side
            output_az+=pi_over_2;
          } else { // star is behind camera, draw on right
            output_az=-pi_over_2-output_az;
          }
        } else { // front=center orientation
          if (star_x < 0.0) { // star is behind camera we need to move to sides of front spherical frame
            if (star_y > 0.0) { // left
              output_az=M_PI-output_az;
            } else { // right
              output_az=-M_PI-output_az;
            } 
          }
        }
        output_x=(int)((-bsr_state->pixels_per_radian * output_az) + bsr_state->camera_half_res_x - 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * output_el) + bsr_state->camera_half_res_y - 0.5);
      } else if (bsr_config->camera_projection == 2) {
        // Hammer
        star_r=sqrt(star_r2);
        star_3az_xy=atan2(star_y, star_x);
        output_az_by2=star_3az_xy / 2.0;
        output_el=asin(star_z / star_r);
        output_x=(int)((-bsr_state->pixels_per_radian * M_PI * cos(output_el) * sin(output_az_by2) / (sqrt(1.0 + (cos(output_el) * cos(output_az_by2))))) + bsr_state->camera_half_res_x - 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * pi_over_2 * sin(output_el) / (sqrt(1.0 + (cos(output_el) * cos(output_az_by2))))) + bsr_state->camera_half_res_y - 0.5);
      } else if (bsr_config->camera_projection == 3) {
        // Mollewide
        star_r=sqrt(star_r2);
        star_3az_xy=atan2(star_y, star_x);
        output_az=star_3az_xy;
        output_el=asin(star_z / star_r);
        two_mollewide_angle=2.0 * asin(2.0 * output_el / M_PI);
        for (i=0; i < bsr_config->Mollewide_iterations; i++) {
          two_mollewide_angle-=(two_mollewide_angle + sin(two_mollewide_angle) - (M_PI * sin(output_el))) / (1.0 + cos(two_mollewide_angle));
        }
        mollewide_angle=two_mollewide_angle * 0.5;
        output_x=(int)((-bsr_state->pixels_per_radian * output_az * cos(mollewide_angle)) + bsr_state->camera_half_res_x - 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * pi_over_2 * sin(mollewide_angle)) + bsr_state->camera_half_res_y - 0.5);
      }

      //
      // if star is within raster bounds, send star (or Airy disk pixels) to dedup buffer
      //
      if ((output_x >= 0) && (output_x < bsr_config->camera_res_x) && (output_y >= 0) && (output_y < bsr_config->camera_res_y)) {
        if (bsr_config->Airy_disk == 1) {
          //
          // Airy disk mode, use Airy disk maps to find all pixel values for this star and send to dedup buffer
          //
          Airymap_red_p=bsr_state->Airymap_red;
          Airymap_green_p=bsr_state->Airymap_green;
          Airymap_blue_p=bsr_state->Airymap_blue;
          star_rgb_red=bsr_state->rgb_red[color_temperature];
          star_rgb_green=bsr_state->rgb_green[color_temperature];
          star_rgb_blue=bsr_state->rgb_blue[color_temperature];
          for (Airymap_y=0; Airymap_y < Airymap_xy; Airymap_y++) {
            for (Airymap_x=0; Airymap_x < Airymap_xy; Airymap_x++) {
              Airymap_output_x=output_x - Airymap_half_xy + Airymap_x;
              Airymap_output_y=output_y - Airymap_half_xy + Airymap_y;
              if ((Airymap_output_x >= 0) && (Airymap_output_x < bsr_config->camera_res_x) && (Airymap_output_y >= 0) && (Airymap_output_y < bsr_config->camera_res_y)
                && (*Airymap_red_p > 0.0) && (*Airymap_green_p > 0.0) && (*Airymap_blue_p > 0.0)) {
                //
                // Airymap pixel is within image raster, send to dedup buffer
                //
                image_offset=((long long)bsr_config->camera_res_x * (long long)Airymap_output_y) + (long long)Airymap_output_x;
                r=(star_linear_intensity * *Airymap_red_p * star_rgb_red);
                g=(star_linear_intensity * *Airymap_green_p * star_rgb_green);
                b=(star_linear_intensity * *Airymap_green_p * star_rgb_blue);
                sendPixelToDedupBuffer(bsr_state, image_offset, r, g, b);
              } // end if Airymap pixel is within image raster
              Airymap_red_p++;
              Airymap_green_p++;
              Airymap_blue_p++;
            } // end for Airymap_x
          } // end for Airymap_y
        } else {
          //
          // not Airy disk mode, send star pixel to dedup buffer
          //
          image_offset=((long long)bsr_config->camera_res_x * (long long)output_y) + (long long)output_x;
          r=(star_linear_intensity * bsr_state->rgb_red[color_temperature]);
          g=(star_linear_intensity * bsr_state->rgb_green[color_temperature]);
          b=(star_linear_intensity * bsr_state->rgb_blue[color_temperature]);
          sendPixelToDedupBuffer(bsr_state, image_offset, r, g, b);
        } // end if Airy disk mode
      } // end if star is within image raster
    } // end if within distance ranges

    //
    // read next star record from input file
    //
    fread(&star_record, star_record_size, 1, input_file);
  } // end input loop

  //
  // done with input file, check for any remaining pixels in dedup buffer and send to main thread
  //
  if (bsr_state->perthread->dedup_count > 0) {
    sendDedupBufferToMainThread(bsr_state);
  } // end if dedup buffer has remaining entries

  return(0);
}

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <math.h>
#include "util.h"

int processStars(bsr_config_t *bsr_config, bsr_state_t *bsr_state, FILE *input_file) {
  int i;
  volatile int success; // gcc optimization breaks code without volatile keyword
  star_record_t star_record;
  int star_record_size=sizeof(star_record_t);
  double star_x;
  double star_y;
  double star_z;
  double star_r;
  double star_xy_r;
  double star_xz_r;
  double star_yz_r;
  double render_distance;            // distance from selected point to star
  double star_3az_xy;
  double star_3az_xz;
  double star_3az_yz;
  const double two_pi=2.0 * M_PI;
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
  dedup_buffer_t *dedup_buf_p;
  dedup_index_t *dedup_index_p;
  long int dedup_index_offset;
  int dedup_count;
  int dedup_buf_i;
  int idle_count;
  long long input_count;

  Airymap_half_xy=bsr_config->Airy_disk_max_extent;
  Airymap_xy=(Airymap_half_xy * 2) + 1;

  //
  // read and process each line of input file and render stars to image
  //
  input_count=0;
  dedup_count=0;
  fread(&star_record, star_record_size, 1, input_file);
  while ((feof(input_file) == 0) && (ferror(input_file) == 0)) {
    input_count++;

    //
    // periodically check if main thread is still alive
    //
    if (input_count % 100000) {
      checkExceptions(bsr_state);
    } 

    //
    // convert original star x,y,z to new coordinates as seen by camera position
    //
    star_x=star_record.icrs_x - bsr_config->camera_icrs_x;
    star_y=star_record.icrs_y - bsr_config->camera_icrs_y;
    star_z=star_record.icrs_z - bsr_config->camera_icrs_z;
    star_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0) + pow(star_z, 2.0));

    //
    // extract color_temperature from combined field
    //
    color_temperature=star_record.intensity_and_temperature & 0x00000000fffffffful;

    // 
    // only continue if star distance is within min-max range from selected point (0=camera, 1=target),
    // greater than zero, and color temperature is within allowed range
    //
    if (bsr_config->render_distance_selector == 0) { // selected point is camera
      render_distance=star_r; // star distance from camera
    } else { // selected point is target
      render_distance=sqrt(pow((star_record.icrs_x - bsr_config->target_icrs_x), 2.0) + pow((star_record.icrs_y - bsr_config->target_icrs_y), 2.0) + pow((star_record.icrs_z - bsr_config->target_icrs_z), 2.0)); // important, use un-rotated coordinates
    }
    if ((star_r > 0.0)\
     && (render_distance >= bsr_config->render_distance_min) && (render_distance <= bsr_config->render_distance_max)\
     && (color_temperature >= bsr_config->star_color_min) && (color_temperature <= bsr_config->star_color_max)) {

      //
      // extract star intensity from combined field and adjust for distance from camera
      //
      star_record.intensity_and_temperature >>=32;
      star_linear_1pc_intensity=*((float*)&star_record.intensity_and_temperature);
      star_linear_intensity=star_linear_1pc_intensity / pow(star_r, 2.0);

      //
      // setup initial triple-azimuth (3az) variables for this star
      //
      star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
      //star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_y == 0.0)) {
        star_3az_xy=0.0;
      } else {
        star_3az_xy=atan2(star_y, star_x);
      }
      if ((star_x == 0.0) && (star_z == 0.0)) {
        star_3az_xz=0.0;
      } else {
        star_3az_xz=atan2(star_z, star_x);
      }
      if ((star_y == 0.0) && (star_z == 0.0)) {
        star_3az_yz=0.0;
      } else {
        star_3az_yz=atan2(star_z, star_y);
      }

      //
      // rotate star xy angle by target xy angle
      //
      star_3az_xy-=bsr_state->target_3az_xy;
      if (star_3az_xy > M_PI) {
        star_3az_xy = -two_pi + star_3az_xy;
      } else if (star_3az_xy < -M_PI) {
        star_3az_xy = two_pi + star_3az_xy;
      }
      // apply star xy rotation angle to star xz angle, only x has changed
      star_x=star_xy_r * cos(star_3az_xy);
      star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_z == 0.0)) {
        star_3az_xz=0.0;
      } else {
        star_3az_xz=atan2(star_z, star_x);
      }
      // apply star xy rotation angle to star yz angle, only y has changed
      star_y=star_xy_r * sin(star_3az_xy);
      //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_y == 0.0) && (star_z == 0.0)) {
        star_3az_yz=0.0;
      } else {
        star_3az_yz=atan2(star_z, star_y);
      }

      //
      // rotate star xz angle by (rotated) target xz angle
      //
      star_3az_xz-=bsr_state->target_3az_xz;
      if (star_3az_xz > M_PI) {
        star_3az_xz = -two_pi + star_3az_xz;
      } else if (star_3az_xz < -M_PI) {
        star_3az_xz = two_pi + star_3az_xz;
      }
      // apply star xz rotation to star xy angle, only x has changed
      star_x=star_xz_r * cos(star_3az_xz);
      //star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_y == 0.0)) {
        star_3az_xy=0.0;
      } else {
        star_3az_xy=atan2(star_y, star_x);
      }
      // apply star xz rotation to star yz angle, only z has changed
      star_z=star_xz_r * sin(star_3az_xz);
      star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_y == 0.0) && (star_z == 0.0)) {
        star_3az_yz=0.0;
      } else {
        star_3az_yz=atan2(star_z, star_y);
      }

      //
      // rotate star yz angle by camera rotation angle
      //
      star_3az_yz+=bsr_state->camera_3az_yz;
      if (star_3az_yz > M_PI) {
        star_3az_yz = -two_pi + star_3az_yz;
      } else if (star_3az_yz < -M_PI) {
        star_3az_yz = two_pi + star_3az_yz;
      }
      // apply star yz rotation to star xy angle, only y has changed
      star_y=star_yz_r * cos(star_3az_yz);
      star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_y == 0.0)) {
        star_3az_xy=0.0;
      } else {
        star_3az_xy=atan2(star_y, star_x);
      }
      // apply star yz rotation to star xz angle, only z has changed
      star_z=star_yz_r * sin(star_3az_yz);
      star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_z == 0.0)) {
        star_3az_xz=0.0;
      } else {
        star_3az_xz=atan2(star_z, star_x);
      }

      //
      // optionally pan camera left/right
      //
      if (bsr_config->camera_pan != 0) {
        // optionally rotate star xy angle by camera pan angle
        star_3az_xy+=bsr_state->camera_3az_xy;
        if (star_3az_xy > M_PI) {
          star_3az_xy = -two_pi + star_3az_xy;
        } else if (star_3az_xy < -M_PI) {
          star_3az_xy = two_pi + star_3az_xy;
        }
        // apply star xy rotation angle to star xz angle, only x has changed
        star_x=star_xy_r * cos(star_3az_xy);
        star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
        if ((star_x == 0.0) && (star_z == 0.0)) {
          star_3az_xz=0.0;
        } else {
          star_3az_xz=atan2(star_z, star_x);
        }
        // apply star xy rotation angle to star yz angle, only y has changed
        star_y=star_xy_r * sin(star_3az_xy);
        //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
        if ((star_y == 0.0) && (star_z == 0.0)) {
          star_3az_yz=0.0;
        } else {
          star_3az_yz=atan2(star_z, star_y);
        }
      }

      //
      // optionally pan camera up/down
      //
      if (bsr_config->camera_tilt != 0.0) {
        // optionally rotate star xz angle by camera tilt angle
        star_3az_xz+=bsr_state->camera_3az_xz;
        if (star_3az_xz > M_PI) {
          star_3az_xz = -two_pi + star_3az_xz;
        } else if (star_3az_xz < -M_PI) {
          star_3az_xz = two_pi + star_3az_xz;
        }
        // apply star xz rotation to star xy angle, only x has changed
        star_x=star_xz_r * cos(star_3az_xz);
        star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
        if ((star_x == 0.0) && (star_y == 0.0)) {
          star_3az_xy=0.0;
        } else {
          star_3az_xy=atan2(star_y, star_x);
        }
        // apply star xz rotation to star yz angle, only z has changed
        star_z=star_xz_r * sin(star_3az_xz);
        //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
        if ((star_y == 0.0) && (star_z == 0.0)) {
          star_3az_yz=0.0;
        } else {
          star_3az_yz=atan2(star_z, star_y);
        }
      }

      //
      // project star onto output raster x,y
      //
      if (bsr_config->camera_projection == 0) {
        // lat/lon 
        output_az=star_3az_xy;
        output_el=atan2(star_z, star_xy_r);
        output_x=(int)((-bsr_state->pixels_per_radian * output_az) + bsr_state->camera_half_res_x - 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * output_el) + bsr_state->camera_half_res_y - 0.5);
      } else if (bsr_config->camera_projection == 1) {
        // spherical
        star_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0) + pow(star_z, 2.0));
        spherical_distance=asin(sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)) / star_r);
        spherical_angle=star_3az_yz;
        output_az=spherical_distance * cos(spherical_angle);
        output_el=spherical_distance * sin(spherical_angle);
        if (bsr_config->spherical_orientation == 1) { // side by side orientation
          if (star_x > 0) { // star is in front of camera, draw on left side
            output_az+=pi_over_2;
          } else { // star is behind camera, draw on right
            output_az=-pi_over_2-output_az;
          }
        } else { // front=center orientation
          if (star_x < 0) { // star is behind camera we need to move to sides of front spherical frame
            if (star_y > 0) { // left
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
        output_az_by2=star_3az_xy / 2.0;
        output_el=atan2(star_z, star_xy_r);
        output_x=(int)((-bsr_state->pixels_per_radian * M_PI * cos(output_el) * sin(output_az_by2) / (sqrt(1.0 + (cos(output_el) * cos(output_az_by2))))) + bsr_state->camera_half_res_x - 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * pi_over_2 * sin(output_el) / (sqrt(1.0 + (cos(output_el) * cos(output_az_by2))))) + bsr_state->camera_half_res_y - 0.5);
      } else if (bsr_config->camera_projection == 3) {
        // Mollewide
        output_az=star_3az_xy;
        output_el=atan2(star_z, star_xy_r);
        two_mollewide_angle=2.0 * asin(2.0 * output_el / M_PI);
        for (i=0; i < bsr_config->Mollewide_iterations; i++) {
          two_mollewide_angle-=(two_mollewide_angle + sin(two_mollewide_angle) - (M_PI * sin(output_el))) / (1.0 + cos(two_mollewide_angle));
        }
        mollewide_angle=two_mollewide_angle * 0.5;
        output_x=(int)((-bsr_state->pixels_per_radian * output_az * cos(mollewide_angle)) + bsr_state->camera_half_res_x - 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * pi_over_2 * sin(mollewide_angle)) + bsr_state->camera_half_res_y - 0.5);
      }

      //
      // if star is within raster bounds, write to thread buffer to send to main thread for integration in final image
      //
      if ((output_x >= 0) && (output_x < bsr_config->camera_res_x) && (output_y >= 0) && (output_y < bsr_config->camera_res_y)) {
        if (bsr_config->Airy_disk == 1) {
          //
          // Airy disk mode, use Airy disk maps to find all pixel values for this star
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
                // Airymap pixel is within image raster, put in my thread buffer
                //
                dedup_index_offset=(bsr_config->camera_res_x * Airymap_output_y) + Airymap_output_x;
                dedup_index_p=bsr_state->dedup_index + dedup_index_offset;
                if (dedup_index_p->dedup_record == NULL) {
                  // no dup yet, just store value in next dedup buffer record and update index
                  dedup_count++;
                  dedup_buf_p=bsr_state->dedup_buf + (dedup_count - 1);
                  dedup_buf_p->image_offset=dedup_index_offset;
                  dedup_buf_p->r=(star_linear_intensity * *Airymap_red_p * star_rgb_red);
                  dedup_buf_p->g=(star_linear_intensity * *Airymap_green_p * star_rgb_green);
                  dedup_buf_p->b=(star_linear_intensity * *Airymap_green_p * star_rgb_blue);
                  dedup_index_p->dedup_record=dedup_buf_p;
                } else {
                  // duplicate pixel locaiton, add to existing dedup buffer record values
                  dedup_buf_p=dedup_index_p->dedup_record;
                  dedup_buf_p->r+=(star_linear_intensity * *Airymap_green_p * star_rgb_red);
                  dedup_buf_p->g+=(star_linear_intensity * *Airymap_green_p * star_rgb_green);
                  dedup_buf_p->b+=(star_linear_intensity * *Airymap_green_p * star_rgb_blue);
                }

                //
                // check if dedup buffer is full and if it is send pixels to main thread
                //
                if (dedup_count == bsr_state->per_thread_buffers) {
                  dedup_buf_p=bsr_state->dedup_buf;
                  for (dedup_buf_i=0; dedup_buf_i < dedup_count; dedup_buf_i++) {
                    dedup_index_p=bsr_state->dedup_index + dedup_buf_p->image_offset;
                    if (bsr_state->perthread->thread_buffer_index == bsr_state->per_thread_buffers) {
                      // end of our section of thread_buf, rewind
                      bsr_state->perthread->thread_buffer_index=0;
                      bsr_state->perthread->thread_buf_p-=bsr_state->per_thread_buffers;
                    }
                    success=0;
                    idle_count=0;
                    while (success == 0) {
                      // check if we've written a pixel to this location before and it has not been read/cleared by main thread yet
                      if ((bsr_state->perthread->thread_buf_p->status_left == 0) && (bsr_state->perthread->thread_buf_p->status_right == 0)) {
                        bsr_state->perthread->thread_buf_p->status_left=1;
                        bsr_state->perthread->thread_buf_p->image_offset=dedup_buf_p->image_offset;
                        bsr_state->perthread->thread_buf_p->r=dedup_buf_p->r;
                        bsr_state->perthread->thread_buf_p->g=dedup_buf_p->g;
                        bsr_state->perthread->thread_buf_p->b=dedup_buf_p->b;
                        bsr_state->perthread->thread_buf_p->status_right=1;
                        bsr_state->perthread->thread_buf_p++;
                        bsr_state->perthread->thread_buffer_index++;
                        // clear buffer location and reset index
                        dedup_buf_p->image_offset=-1;
                        dedup_buf_p->r=0.0;
                        dedup_buf_p->g=0.0;
                        dedup_buf_p->b=0.0;
                        dedup_index_p->dedup_record=NULL;
                        dedup_buf_p++;
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
                  } // end for dedup_buf_i
                  dedup_count=0;
                  dedup_buf_p=bsr_state->dedup_buf;
                } // end if dedup buffer full
              } // end if Airymap pixel is within image raster
              Airymap_red_p++;
              Airymap_green_p++;
              Airymap_blue_p++;
            } // end for Airymap_x
          } // end for Airymap_y
        } else {
          //
          // not Airy disk mode, check if this pixel in the dedup buffer
          //
          dedup_index_offset=(bsr_config->camera_res_x * output_y) + output_x;
          dedup_index_p=bsr_state->dedup_index + dedup_index_offset;
          if (dedup_index_p->dedup_record == NULL) {
            // no dup yet, just store value in next dedup buffer record and update index
            dedup_count++;
            dedup_buf_p=bsr_state->dedup_buf + (dedup_count - 1);
            dedup_buf_p->image_offset=dedup_index_offset;
            dedup_buf_p->r=(star_linear_intensity * bsr_state->rgb_red[color_temperature]);
            dedup_buf_p->g=(star_linear_intensity * bsr_state->rgb_green[color_temperature]);
            dedup_buf_p->b=(star_linear_intensity * bsr_state->rgb_blue[color_temperature]);
            dedup_index_p->dedup_record=dedup_buf_p;
          } else {
            // duplicate pixel locaiton, add to existing dedup buffer record values
            dedup_buf_p=dedup_index_p->dedup_record;
            dedup_buf_p->r+=(star_linear_intensity * bsr_state->rgb_red[color_temperature]);
            dedup_buf_p->g+=(star_linear_intensity * bsr_state->rgb_green[color_temperature]);
            dedup_buf_p->b+=(star_linear_intensity * bsr_state->rgb_blue[color_temperature]);
          }

          //
          // check if dedup buffer is full and if it is send pixels to main thread
          //
          if (dedup_count == bsr_state->per_thread_buffers) {
            dedup_buf_p=bsr_state->dedup_buf;
            for (dedup_buf_i=0; dedup_buf_i < dedup_count; dedup_buf_i++) { 
              dedup_index_p=bsr_state->dedup_index + dedup_buf_p->image_offset;
              if (bsr_state->perthread->thread_buffer_index == bsr_state->per_thread_buffers) {
                // end of our section of thread_buf, rewind
                bsr_state->perthread->thread_buffer_index=0;
                bsr_state->perthread->thread_buf_p-=bsr_state->per_thread_buffers;
              }
              success=0;
              idle_count=0;
              while (success == 0) {
                // check if we've written a pixel to this location before and it has not been read/cleared by main thread yet
                if ((bsr_state->perthread->thread_buf_p->status_left == 0) && (bsr_state->perthread->thread_buf_p->status_right == 0)) {
                  bsr_state->perthread->thread_buf_p->status_left=1;
                  bsr_state->perthread->thread_buf_p->image_offset=dedup_buf_p->image_offset;
                  bsr_state->perthread->thread_buf_p->r=dedup_buf_p->r;
                  bsr_state->perthread->thread_buf_p->g=dedup_buf_p->g;
                  bsr_state->perthread->thread_buf_p->b=dedup_buf_p->b;
                  bsr_state->perthread->thread_buf_p->status_right=1;
                  bsr_state->perthread->thread_buf_p++;
                  bsr_state->perthread->thread_buffer_index++;
                  // clear buffer location and reset index
                  dedup_buf_p->image_offset=-1;
                  dedup_buf_p->r=0.0;
                  dedup_buf_p->g=0.0;
                  dedup_buf_p->b=0.0;
                  dedup_index_p->dedup_record=NULL;
                  dedup_buf_p++;
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
            } // end for dedup_buf_i
            dedup_count=0;
            dedup_buf_p=bsr_state->dedup_buf;
          } // end if dedup buffer full
        } // end if Airy disk mode
      } // end if star is within image raster
    } // end if within distance ranges

    //
    // read next star record from input file
    //
    fread(&star_record, star_record_size, 1, input_file);
  } // end input loop

  //
  // check for remaining pixels in dedup buffer and send to main thread
  //
  if (dedup_count > 0) {
    dedup_buf_p=bsr_state->dedup_buf;
    for (dedup_buf_i=0; dedup_buf_i < dedup_count; dedup_buf_i++) {
      dedup_index_p=bsr_state->dedup_index + dedup_buf_p->image_offset;
      if (bsr_state->perthread->thread_buffer_index == bsr_state->per_thread_buffers) {
        // end of our section of thread_buf, rewind
        bsr_state->perthread->thread_buffer_index=0;
        bsr_state->perthread->thread_buf_p-=bsr_state->per_thread_buffers;
      }
      success=0;
      while (success == 0) {
        // check if we've written a pixel to this location before and it has not been read/cleared by main thread yet
        if ((bsr_state->perthread->thread_buf_p->status_left == 0) && (bsr_state->perthread->thread_buf_p->status_right == 0)) {
          bsr_state->perthread->thread_buf_p->status_left=1;
          bsr_state->perthread->thread_buf_p->image_offset=dedup_buf_p->image_offset;
          bsr_state->perthread->thread_buf_p->r=dedup_buf_p->r;
          bsr_state->perthread->thread_buf_p->g=dedup_buf_p->g;
          bsr_state->perthread->thread_buf_p->b=dedup_buf_p->b;
          bsr_state->perthread->thread_buf_p->status_right=1;
          bsr_state->perthread->thread_buf_p++;
          bsr_state->perthread->thread_buffer_index++;
          // clear buffer location and reset index
          dedup_buf_p->image_offset=-1;
          dedup_buf_p->r=0.0;
          dedup_buf_p->g=0.0;
          dedup_buf_p->b=0.0;
          dedup_index_p->dedup_record=NULL;
          dedup_buf_p++;
          success=1;
        } // end if buffer slot is available
      } // end while success=0
    } // end for dedup_buf_i
  } // end if dedup buffer has remaining entries

  return(0);
}

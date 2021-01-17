#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <math.h>

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
  double output_el;
  int output_x=0;
  int output_y=0;
  double spherical_distance;
  double spherical_angle;
  double two_mollewide_angle;
  double mollewide_angle;
  const double pi_over_2=M_PI / 2.0;

  //
  // read and process each line of input file and render stars to image
  //
  fread(&star_record, star_record_size, 1, input_file);
  while ((feof(input_file) == 0) && (ferror(input_file) == 0)) {

    //
    // convert original star x,y,z to new coordinates as seen by camera position
    //
    star_x=star_record.icrs_x - bsr_config->camera_icrs_x;
    star_y=star_record.icrs_y - bsr_config->camera_icrs_y;
    star_z=star_record.icrs_z - bsr_config->camera_icrs_z;
    star_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0) + pow(star_z, 2.0));

    // 
    // only continue if star distance is within min-max range from selected point (0=camera, 1=target)
    //
    if (bsr_config->render_distance_selector == 0) { // selected point is camera
      render_distance=star_r; // star distance from camera
    } else { // selected point is target
      render_distance=sqrt(pow((star_record.icrs_x - bsr_config->target_icrs_x), 2.0) + pow((star_record.icrs_y - bsr_config->target_icrs_y), 2.0) + pow((star_record.icrs_z - bsr_config->target_icrs_z), 2.0)); // important, use un-rotated coordinates
    }
    if ((render_distance >= bsr_config->render_distance_min) && (render_distance <= bsr_config->render_distance_max)) {
      //
      // extract color_temperature from combined field
      //
      color_temperature=star_record.intensity_and_temperature & 0x00000000fffffffful;

      //
      // extract star intensity from combined field and adjust for distance from camera
      //
      star_record.intensity_and_temperature >>=32;
      star_linear_1pc_intensity=*((float*)&star_record.intensity_and_temperature);
      star_linear_intensity=star_linear_1pc_intensity / pow(star_r, 2.0);

      //
      // setup initial triple-azimuth (3az) variables for this star
      // 3az is a hybrid of quaternions/versors and Euler angles.   It uses three full-circle orthoganal angles for redundancy to avoid problems around poles (gimbal lock) but can be
      // very quickly converted into lat/lon raster when done rotating.
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
        output_x=(int)((-bsr_state->pixels_per_radian * output_az) + bsr_state->camera_half_res_x + 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * output_el) + bsr_state->camera_half_res_y + 0.5);
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
        output_x=(int)((-bsr_state->pixels_per_radian * output_az) + bsr_state->camera_half_res_x + 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * output_el) + bsr_state->camera_half_res_y + 0.5);
      } else if (bsr_config->camera_projection == 2) {
        // Mollewide
        output_az=star_3az_xy;
        output_el=atan2(star_z, star_xy_r);
        two_mollewide_angle=2.0 * asin(2.0 * output_el / M_PI);
        for (i=0; i < bsr_config->Mollewide_iterations; i++) {
          two_mollewide_angle-=(two_mollewide_angle + sin(two_mollewide_angle) - (M_PI * sin(output_el))) / (1.0 + cos(two_mollewide_angle));
        }
        mollewide_angle=two_mollewide_angle * 0.5;
        output_x=(int)((-bsr_state->pixels_per_radian * output_az * cos(mollewide_angle)) + bsr_state->camera_half_res_x + 0.5);
        output_y=(int)((-bsr_state->pixels_per_radian * pi_over_2 * sin(mollewide_angle)) + bsr_state->camera_half_res_y + 0.5);
      }

      //
      // if star is within raster bounds, write to thread buffer to send to main thread for integration in final image
      //
      if ((output_x >= 0) && (output_x < bsr_config->camera_res_x) && (output_y >= 0) && (output_y < bsr_config->camera_res_y)) {
        // put pixel in my thread buffer
        if (bsr_state->thread_buffer_index == bsr_config->per_thread_buffer) {
          // end of our section of thread_buf, rewind
          bsr_state->thread_buffer_index=0;
          bsr_state->thread_buf_p-=bsr_config->per_thread_buffer;
        }
        success=0;
        while (success == 0) {
          // check if we've written a pixel to this location before and it has not been read/cleared by main thread yet
          if ((bsr_state->thread_buf_p->status_left == 0) && (bsr_state->thread_buf_p->status_right == 0)) {
            bsr_state->thread_buf_p->status_left=1;
            bsr_state->thread_buf_p->image_offset=(bsr_config->camera_res_x * output_y) + output_x;
            bsr_state->thread_buf_p->r=(star_linear_intensity * bsr_state->rgb_red[color_temperature]);
            bsr_state->thread_buf_p->g=(star_linear_intensity * bsr_state->rgb_green[color_temperature]);
            bsr_state->thread_buf_p->b=(star_linear_intensity * bsr_state->rgb_blue[color_temperature]);
            bsr_state->thread_buf_p->status_right=1;
            bsr_state->thread_buf_p++;
            bsr_state->thread_buffer_index++;
            success=1;
          } // end if buffer slot is available
        } // end while success=0
      } // end if within image raster
    } // end if within distance ranges

    //
    // read next star record from input file
    //
    fread(&star_record, star_record_size, 1, input_file);

  } // end input loop

  return(0);
}

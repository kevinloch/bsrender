#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <math.h>

int initState(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  double pi_over_360=M_PI / 360.0;
  double pi_over_180=M_PI / 180.0;
  double camera_icrs_ra_rad;
  double camera_icrs_dec_rad;
  double target_icrs_ra_rad;
  double target_icrs_dec_rad;

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
  // convert camera target x,y,z to triple-azimuth (3az) coordinates as seen from camera
  //
  bsr_state->target_x=bsr_config->target_icrs_x - bsr_config->camera_icrs_x;
  bsr_state->target_y=bsr_config->target_icrs_y - bsr_config->camera_icrs_y;
  bsr_state->target_z=bsr_config->target_icrs_z - bsr_config->camera_icrs_z;
  //bsr_state->target_r=sqrt(pow(bsr_state->target_x, 2.0) + pow(bsr_state->target_y, 2.0) + pow(bsr_state.target_z, 2.0));
  bsr_state->target_xy_r=sqrt(pow(bsr_state->target_x, 2.0) + pow(bsr_state->target_y, 2.0)); // may be used in future rotations or raster projections
  //bsr_state->target_xz_r=sqrt(pow(bsr_state->target_x, 2.0) + pow(bsr_state->target_z, 2.0)); // may be used in future rotations or raster projections
  //bsr_state->target_yz_r=sqrt(pow(bsr_state->target_y, 2.0) + pow(bsr_state->target_z, 2.0)); // may be used in future rotations or raster projections
  if ((bsr_state->target_x == 0.0) && (bsr_state->target_y == 0.0)) {
    bsr_state->target_3az_xy=0.0;
  } else {
    bsr_state->target_3az_xy=atan2(bsr_state->target_y, bsr_state->target_x);
  }
  if ((bsr_state->target_x == 0.0) && (bsr_state->target_z == 0.0)) {
    bsr_state->target_3az_xz=0.0;
  } else {
    bsr_state->target_3az_xz=atan2(bsr_state->target_z, bsr_state->target_x);
  }
  //if ((bsr_state->target_y == 0.0) && (bsr_state->target_z == 0.0)) {
  //  bsr_state->target_3az_yz=0.0;
  //} else {
  //  bsr_state->target_3az_yz=atan2(bsr_state->target_z, bsr_state->target_y);
  //}

  //
  // calculate target xz rotation after setting xy=0
  //
  // apply target xy rotation angle to target xz angle, only x has changed
  bsr_state->target_x=bsr_state->target_xy_r; // xy=0
  //bsr_state->target_xz_r=sqrt(pow(bsr_state->target_x, 2.0) + pow(bsr_state->target_z, 2.0)); // may be used in future rotations or raster projections
  if ((bsr_state->target_x == 0.0) && (bsr_state->target_z == 0.0)) {
    bsr_state->target_3az_xz=0.0;
  } else {
    bsr_state->target_3az_xz=atan2(bsr_state->target_z, bsr_state->target_x);
  }

  return(0);
}

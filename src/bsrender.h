#ifndef NLE_LEPTON_H
#define NLE_LEPTON_H

#define BSR_VERSION "0.9-dev-05"

#define _GNU_SOURCE // needed for strcasestr in string.h
#include <stdint.h> // needed for uint64_t

typedef struct {
  double icrs_x;
  double icrs_y;
  double icrs_z;
  uint64_t intensity_and_temperature;
} star_record_t;

typedef struct {
  double r;
  double g;
  double b;
} pixel_composition_t;

typedef struct {
  char config_file_name[256];
  int num_processes;  
  int min_parallax_quality;
  double render_distance_min;
  double render_distance_max;
  int render_distance_selector;
  int draw_cross_hairs;
  int draw_grid_lines;
  int cgi_mode;
  int cgi_max_res_x;
  int cgi_max_res_y;
  int cgi_min_parallax_quality;
  int camera_res_x;
  int camera_res_y;
  double camera_fov;
  double camera_wb_temp;
  double camera_pixel_limit_mag;
  double camera_pixel_limit;
  int camera_pixel_limit_mode;
  double camera_color_saturation;
  int camera_projection;
  int Mollewide_iterations;
  double camera_gamma;
  int sRGB_gamma;
  double camera_icrs_x;
  double camera_icrs_y;
  double camera_icrs_z;
  double camera_icrs_ra;
  double camera_icrs_dec;
  double camera_icrs_r;
  double target_icrs_x;
  double target_icrs_y;
  double target_icrs_z;
  double target_icrs_ra;
  double target_icrs_dec;
  double target_icrs_r;
  double camera_rotation;
  double camera_pan;
  double camera_tilt;
} bsr_config_t;

#endif // NLE_LEPTON_H


#ifndef BSRENDER_H
#define BSRENDER_H

#define BSR_VERSION "0.9-dev-26"

#define _GNU_SOURCE // needed for strcasestr in string.h
#include <stdint.h> // needed for uint64_t

typedef struct {
  double icrs_x;
  double icrs_y;
  double icrs_z;
  uint64_t intensity_and_temperature;
} star_record_t;

typedef struct {
  int status_left;
  long int image_offset;
  double r;
  double g;
  double b;
  int status_right; // left/right status fields help us deal with situation where truct is split between difference cache lines and written back at different times
} thread_buffer_t;

typedef struct {
  double r;
  double g;
  double b;
} pixel_composition_t;

typedef struct {
  thread_buffer_t *thread_buf; // globally mmaped
  thread_buffer_t *thread_buf_p; // used locally by each thread
  int thread_buffer_index; // used locally by each thread
  pixel_composition_t *image_composition_buf; // globally mmaped
  double *rgb_red;
  double *rgb_green;
  double *rgb_blue;
  double *Airymap_red;
  double *Airymap_green;
  double *Airymap_blue;
  double camera_hfov;
  double camera_half_res_x;
  double camera_half_res_y;
  double pixels_per_radian;
  double camera_3az_xy;
  double camera_3az_xz;
  double camera_3az_yz;
  double target_x;
  double target_y;
  double target_z;
  //double target_r;
  double target_xy_r;
  double target_xz_r;
  double target_yz_r;
  double target_3az_xy;
  double target_3az_xz;
  //double target_3az_yz;
  int my_thread_id; // used locally by each thread
} bsr_state_t;

typedef struct {
  char config_file_name[256];
  char data_file_directory[256];
  int num_threads;
  int per_thread_buffer;
  int cgi_mode;
  int cgi_max_res_x;
  int cgi_max_res_y;
  int cgi_min_parallax_quality;
  int min_parallax_quality;
  double render_distance_min;
  double render_distance_max;
  int render_distance_selector;
  double star_color_min;
  double star_color_max;
  int draw_crosshairs;
  int draw_grid_lines;
  int sRGB_gamma;
  int camera_res_x;
  int camera_res_y;
  double camera_fov;
  double camera_pixel_limit_mag;
  double camera_pixel_limit;
  int camera_pixel_limit_mode;
  int camera_wb_enable;
  double camera_wb_temp;
  double camera_color_saturation;
  double camera_gamma;
  int camera_projection;
  int spherical_orientation;
  int Mollewide_iterations;
  int Airy_disk;
  double Airy_disk_first_null;
  int Airy_disk_max_extent;
  double red_filter_long_limit;
  double red_filter_short_limit;
  double green_filter_long_limit;
  double green_filter_short_limit;
  double blue_filter_long_limit;
  double blue_filter_short_limit;
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

#endif // BSRENDER_H

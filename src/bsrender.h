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

#ifndef BSRENDER_H
#define BSRENDER_H

#define BSR_VERSION "1.0-dev-21i"

//
// user configurable compile-time ptions
//
#define BSR_USE_JPEG
#define BSR_USE_PNG
#define BSR_USE_EXR
#define BSR_USE_AVIF
#define BSR_USE_HEIF
#define BSR_32BIT_BUFFERS // use 32-bit floats in image composition, blur, and resize buffers. This reduces the size of these buffers by half which may
                          // be useful for extremely large image resolutions at the expense of summation precision within these buffers. This does not
                          // change the main thread, or dedup buffers which are relatively small and always double precision.


//
// Binary data file details
//
// As of v1.0, bsrender data files have a fixed-length 256 byte header. The header contains a single line of ascii text terminated with
// LF (0x0a) and NULL (0x0). The remaining unused bytes in the header must be set to 0x0. The first 11 bytes of the header must contain
// the appropriate magic number file identifier (BSRENDER_LE or BSRENDER_BE). The header contents can be viewed with 'tail -1 <filename>'
//
// The header is followed by a variable number of 33 byte star records.
//
// Each star record includes a 64-bit unsigned integer for Gaia DR3 'source_id', three 40-bit truncaed doubles for x,y,z,
// a 24-bit truncated float for linear_1pc_intensity, a 24-bit truncated float for linear_1pc_intensity_undimmed,
// a 16-bit unsigned int for color_temperature, and a 16-bit unsigned int for color_temperature_unreddened.
// these are closely packed into a 33 byte star record with each field encoded in the selected byte order.
//
// +---------------+---------+---------+---------+-----+-----+---+---+
// |   source_id   |    x    |    y    |    z    | li  |li-u | c |c-u|
// +---------------+---------+---------+---------+-----+-----+---+---+
// |      8        |    5    |    5    |    5    |  3  |  3  | 2 | 2 |
//                               bytes
//
// mkgalaxy and mkexternal have options to generate either little-endian or big-endian files but the bye-order of the file
// must match the architecture it is used on with bsrender. This is to avoid unnessary operations in the performance
// critical inner-loop of processStars() which iterates over potentially billions of star records. The filenames contain '-le'
// or '-be' to indicate byte order. Byte order is also indicated with the file identifier in the first 11 bytes of the header:
// BSRENDER_LE for little-endian and BSRENDER_BE for big-endian.
//

//
// "advanced" compile-time options
//
#define BSR_EXTERNAL_PREFIX "galaxy-external"
#define BSR_GDR3_PREFIX "galaxy-gdr3"
#define BSR_LE_SUFFIX "le" // filename suffix for little-endian files
#define BSR_BE_SUFFIX "be" // filename suffix for big-endian files
#define BSR_EXTENSION "bsr" // file extension
#define BSR_FILE_HEADER_SIZE 256 // bytes, ascii including magic number
#define BSR_MAGIC_NUMBER_LE "BSRENDER_LE" // file identifier for little-endian files, included in file header size
#define BSR_MAGIC_NUMBER_BE "BSRENDER_BE" // file identifier for big-endian files, included in file header size
#define BSR_STAR_RECORD_SIZE 33  // bytes
#define BSR_BLUR_RESCALE 16777216.0 // pixel values are divided by this number before Gaussian blur to help keep values between [0..1]
#define BSR_RESIZE_LOG_OFFSET 1.0E-6 // pixel values are converted to log(BSR_LOG_OFFSET + pixel value) before Lanczos scaline to minimize clipping artifacts

#define _GNU_SOURCE // needed for strcasestr in string.h
#include <stdint.h> // needed for uint64_t
#include <unistd.h>
#include <sys/stat.h>

//
// For most things we detect endianness runtime with littleEndianTest(). For certain expensive
// operations we use these compile time conditionals instead.
// A check in printVersion() will detect a mismatch from compiled and runtime modes
//
// from https://stackoverflow.com/questions/4239993/determining-endianness-at-compile-time
//
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__ || \
    defined(__BIG_ENDIAN__) || \
    defined(__ARMEB__) || \
    defined(__THUMBEB__) || \
    defined(__AARCH64EB__) || \
    defined(_MIBSEB) || defined(__MIBSEB) || defined(__MIBSEB__)
// It's a big-endian target architecture
#define BSR_BIG_ENDIAN_COMPILE
#elif defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ || \
    defined(__LITTLE_ENDIAN__) || \
    defined(__ARMEL__) || \
    defined(__THUMBEL__) || \
    defined(__AARCH64EL__) || \
    defined(_MIPSEL) || defined(__MIPSEL) || defined(__MIPSEL__)
// It's a little-endian target architecture
#define BSR_LITTLE_ENDIAN_COMPILE
#endif

//
// these checkpoints are used to monitor and control worker thread progress
// they need to be in correct numerical order corresponding to the order in which each
// step might be executed.
//
typedef enum {
  THREAD_STATUS_INVALID                           = -1,
  THREAD_STATUS_AIRY_MAP_BEGIN                    = 10,
  THREAD_STATUS_AIRY_MAP_COMPLETE                 = 11,
  THREAD_STATUS_AIRY_MAP_CONTINUE                 = 12,
  THREAD_STATUS_INIT_IMAGECOMP_BEGIN              = 20,
  THREAD_STATUS_INIT_IMAGECOMP_COMPLETE           = 21,
  THREAD_STATUS_INIT_IMAGECOMP_CONTINUE           = 22,
  THREAD_STATUS_PROCESS_STARS_BEGIN               = 30,
  THREAD_STATUS_PROCESS_STARS_COMPLETE            = 31,
  THREAD_STATUS_PROCESS_STARS_CONTINUE            = 32,
  THREAD_STATUS_POST_PROCESS_BEGIN                = 40,
  THREAD_STATUS_POST_PROCESS_COMPLETE             = 41,
  THREAD_STATUS_POST_PROCESS_CONTINUE             = 42,
  THREAD_STATUS_GAUSSIAN_BLUR_PREP_BEGIN          = 50,
  THREAD_STATUS_GAUSSIAN_BLUR_PREP_COMPLETE       = 51,
  THREAD_STATUS_GAUSSIAN_BLUR_HORIZONTAL_BEGIN    = 52,
  THREAD_STATUS_GAUSSIAN_BLUR_HORIZONTAL_COMPLETE = 53,
  THREAD_STATUS_GAUSSIAN_BLUR_VERTICAL_BEGIN      = 54,
  THREAD_STATUS_GAUSSIAN_BLUR_VERTICAL_COMPLETE   = 55,
  THREAD_STATUS_GAUSSIAN_BLUR_CONTINUE            = 56,
  THREAD_STATUS_LANCZOS_PREP_BEGIN                = 60,
  THREAD_STATUS_LANCZOS_PREP_COMPLETE             = 61,
  THREAD_STATUS_LANCZOS_RESAMPLE_BEGIN            = 62,
  THREAD_STATUS_LANCZOS_RESAMPLE_COMPLETE         = 63,
  THREAD_STATUS_LANCZOS_POINTERS_BEGIN            = 64,
  THREAD_STATUS_LANCZOS_POINTERS_COMPLETE         = 65,
  THREAD_STATUS_LANCZOS_CONTINUE                  = 66,
  THREAD_STATUS_SEQUENCE_PIXELS_BEGIN             = 70,
  THREAD_STATUS_SEQUENCE_PIXELS_COMPLETE          = 71,
  THREAD_STATUS_SEQUENCE_PIXELS_CONTINUE          = 72,
  THREAD_STATUS_IMAGE_COMPRESS_BEGIN              = 80,
  THREAD_STATUS_IMAGE_COMPRESS_COMPLETE           = 81,
  THREAD_STATUS_IMAGE_OUTPUT_BEGIN                = 82,
  THREAD_STATUS_IMAGE_OUTPUT_COMPLETE             = 83,
  THREAD_STATUS_IMAGE_OUTPUT_CONTINUE             = 84,
} bsr_thread_status_t;

typedef struct {
  float redX;
  float redY;
  float greenX;
  float greenY;
  float blueX;
  float blueY;
  float whiteX;
  float whiteY;
} chromaticities_t;

typedef struct {
  double r;
  double i;
  double j;
  double k;
} quaternion_t;

typedef struct {
  pid_t pid;
  int status;
} bsr_status_t;

typedef struct {
  int status_left;
  uint64_t image_offset;
  double r;
  double g;
  double b;
  int status_right; // left/right status fields help us deal with situation where struct is split between difference cache lines and written back at different times
} thread_buffer_t;

typedef struct {
  uint64_t image_offset;
  double r;
  double g;
  double b;
} dedup_buffer_t;

typedef struct {
  dedup_buffer_t *dedup_record_p;
} dedup_index_t;

typedef struct {
#ifdef BSR_32BIT_BUFFERS
  float r;
  float g;
  float b;
#else
  double r;
  double g;
  double b;
#endif
} pixel_composition_t;

typedef struct {
  //
  // these are not globally mmapped so they can be set differently by each thread after fork()
  //
  thread_buffer_t *thread_buf_p;
  int thread_buffer_index;
  int my_thread_id;
  pid_t my_pid;
  int dedup_count;
} bsr_thread_state_t;

typedef struct {
  int fd;
  struct stat sb;
  char *buf; // pointer to large input file, globally mmapped
  size_t buf_size;
} input_file_t;

typedef struct {
  //
  // bsr_state is globally mmapped so all of these variables will be the sync'ed between threads
  // the remaining comments in this struct refer to the objects the pointers point to, which may or may not be mmapped
  // depending on if they are initialized and/or updated by multiple threads
  //
  thread_buffer_t *thread_buf;                // updated by all threads, globally mmaped
  pixel_composition_t *image_composition_buf; // updated by all threads, globally mmaped
  unsigned char *image_output_buf;            // updated by all threads, globally mmaped
  unsigned char **row_pointers;               // updated by all threads, globally mmaped
  int *compressed_sizes;                      // updated by all threads, globally mmaped
  pixel_composition_t *image_blur_buf;        // updated by all threads, globally mmaped
  pixel_composition_t *image_resize_buf;      // updated by all threads, globally mmaped
  dedup_buffer_t *dedup_buf;        // thread-specific buffer, malloc'ed so each thread get's it's own local buffer when fork()'ed
  dedup_index_t *dedup_index;       // thread-specific buffer, malloc'ed so each thread get's it's own local buffer when fork()'ed
  unsigned char *compression_buf1;  // thread-specific buffer, malloc'ed so each thread get's it's own local buffer when fork()'ed
  unsigned char *compression_buf2;  // thread-specific buffer, malloc'ed so each thread get's it's own local buffer when fork()'ed
  input_file_t input_file_external;
  input_file_t input_file_pq100;
  input_file_t input_file_pq050;
  input_file_t input_file_pq030;
  input_file_t input_file_pq020;
  input_file_t input_file_pq010;
  input_file_t input_file_pq005;
  input_file_t input_file_pq003;
  input_file_t input_file_pq002;
  input_file_t input_file_pq001;
  input_file_t input_file_pq000;
  int dedup_index_mode;
  int resize_res_x;
  int resize_res_y;
  pixel_composition_t *current_image_buf; // just a pointer to one of the real image buffers which are all globally mmapped
  int current_image_res_x;
  int current_image_res_y;
  int num_worker_threads;
  pid_t main_pid;
  pid_t main_pgid;
  pid_t httpd_pid;
  bsr_thread_state_t *perthread; // thread-specific variables, not globally mmapped
  int per_thread_buffers;
  int thread_buffer_count;
  bsr_status_t *status_array;    // updated by all threads, globally mmaped
  double rgb_red[32768];
  double rgb_green[32768];
  double rgb_blue[32768];
  double *Airymap_red;           // multi-thread initialization, globally mmapped
  double *Airymap_green;         // multi-thread initialization, globally mmapped
  double *Airymap_blue;          // multi-thread initialization, globally mmapped
  double camera_hfov;
  double camera_half_res_x;
  double camera_half_res_y;
  double pixels_per_radian;
  double render_distance_min2;
  double render_distance_max2;
  double camera_pixel_limit;
  double linear_star_intensity_min;
  double linear_star_intensity_max;
  double anti_alias_per_pixel;
  quaternion_t target_rotation;
  int little_endian;
  size_t composition_buffer_size;
  size_t output_buffer_size;
  size_t row_pointers_size;
  size_t compressed_sizes_size;
  size_t blur_buffer_size;
  size_t resize_buffer_size;
  size_t thread_buffer_size;
  size_t status_array_size;
  size_t dedup_buffer_size;
  size_t dedup_index_size;
  size_t compression_buf_size;
  size_t Airymap_size;
  size_t bsr_state_size;
} bsr_state_t;

typedef struct {
  int use_bandpass_ratios;
  int use_gspphot_distance;
  int calibrate_parallax;
  int enable_maximum_distance;
  double maximum_distance;
  int output_little_endian;
} mkg_config_t;

typedef struct {
  char bsrender_cfg_version[256];
  char *QUERY_STRING_p;
  char config_file_name[256];
  char data_file_directory[256];
  char output_file_name[256];
  int print_status;
  int num_threads;
  int per_thread_buffer;
  int per_thread_buffer_Airy;
  int cgi_mode;
  int cgi_max_res_x;
  int cgi_max_res_y;
  int cgi_Gaia_min_parallax_quality;
  int cgi_allow_Airy_disk;
  double cgi_min_Airy_disk_first_null;
  int cgi_max_Airy_disk_max_extent;
  int cgi_max_Airy_disk_min_extent;
  int cgi_allow_anti_alias;
  int Gaia_db_enable;
  int Gaia_min_parallax_quality;
  int external_db_enable;
  double render_distance_min;
  double render_distance_max;
  int render_distance_selector;
  double star_intensity_min;
  double star_intensity_max;
  int star_intensity_selector;
  double star_color_min;
  double star_color_max;
  int extinction_dimming_undo;
  int extinction_reddening_undo;
  int camera_res_x;
  int camera_res_y;
  double camera_fov;
  double camera_pixel_limit_mag;
  int camera_pixel_limit_mode;
  int camera_wb_enable;
  double camera_wb_temp;
  double camera_color_saturation;
  double camera_gamma;
  int camera_projection;
  int spherical_orientation;
  int Mollewide_iterations;
  double red_filter_long_limit;
  double red_filter_short_limit;
  double green_filter_long_limit;
  double green_filter_short_limit;
  double blue_filter_long_limit;
  double blue_filter_short_limit;
  int Airy_disk_enable;
  double Airy_disk_first_null;
  int Airy_disk_max_extent;
  int Airy_disk_min_extent;
  double Airy_disk_obstruction;
  int anti_alias_enable;
  double anti_alias_radius;
  int skyglow_enable;
  double skyglow_temp;
  double skyglow_per_pixel_mag;
  int pre_limit_intensity;
  double Gaussian_blur_radius;
  double output_scaling_factor;
  int Lanczos_order;
  int draw_crosshairs;
  int draw_grid_lines;
  int output_format;
  int color_profile;
  int exr_compression;
  int compression_quality;
  int image_format;
  int hdr_neutral_white_ref;
  int bits_per_color;
  int image_number_format;
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

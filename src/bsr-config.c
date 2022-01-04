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
#include <stdlib.h>
#include <string.h>
#include <math.h>

void initConfig(bsr_config_t *bsr_config) {
  strcpy(bsr_config->config_file_name, "./bsrender.cfg");
  strcpy(bsr_config->data_file_directory, "./galaxydata");
  bsr_config->num_threads=16;
  bsr_config->per_thread_buffer=10000;
  bsr_config->per_thread_buffer_Airy=200000;
  bsr_config->cgi_mode=0;
  bsr_config->cgi_max_res_x=999999;
  bsr_config->cgi_max_res_y=999999;
  bsr_config->cgi_Gaia_min_parallax_quality=0;
  bsr_config->cgi_allow_Airy_disk=1;
  bsr_config->cgi_min_Airy_disk_first_null=0.3;
  bsr_config->cgi_max_Airy_disk_max_extent=1000;
  bsr_config->enable_Gaia=1;
  bsr_config->Gaia_min_parallax_quality=0;
  bsr_config->enable_external=1;
  bsr_config->render_distance_min=0.0;
  bsr_config->render_distance_min2=bsr_config->render_distance_min * bsr_config->render_distance_min;
  bsr_config->render_distance_max=1.0E99;
  bsr_config->render_distance_max2=bsr_config->render_distance_max * bsr_config->render_distance_max;
  bsr_config->render_distance_selector=0;
  bsr_config->star_color_min=0.0;
  bsr_config->star_color_max=1.0E99;
  bsr_config->camera_res_x=2000;
  bsr_config->camera_res_y=1000;
  bsr_config->camera_fov=360.0;
  bsr_config->camera_pixel_limit_mag=6.5;
  bsr_config->camera_pixel_limit=pow(100.0, (-bsr_config->camera_pixel_limit_mag / 5.0));
  bsr_config->camera_pixel_limit_mode=0;
  bsr_config->camera_wb_enable=1;
  bsr_config->camera_wb_temp=4300.0;
  bsr_config->camera_color_saturation=1.0;
  bsr_config->camera_gamma=1.0;
  bsr_config->camera_projection=0;
  bsr_config->spherical_orientation=0;
  bsr_config->Mollewide_iterations=5;
  bsr_config->red_filter_long_limit=705.0;
  bsr_config->red_filter_short_limit=550.0;
  bsr_config->green_filter_long_limit=600.0;
  bsr_config->green_filter_short_limit=445.0;
  bsr_config->blue_filter_long_limit=465.0;
  bsr_config->blue_filter_short_limit=395.0;
  bsr_config->Airy_disk=0;
  bsr_config->Airy_disk_first_null=0.75;
  bsr_config->Airy_disk_max_extent=100;
  bsr_config->Airy_disk_min_extent=1;
  bsr_config->Airy_disk_obstruction=0.0;
  bsr_config->skyglow_enable=0;
  bsr_config->skyglow_temp=4500.0;
  bsr_config->skyglow_per_pixel_mag=11.0;
  bsr_config->Gaussian_blur_radius=0.0;
  bsr_config->output_scaling_factor=1.0;
  bsr_config->draw_crosshairs=0;
  bsr_config->draw_grid_lines=0;
  bsr_config->sRGB_gamma=1;
  bsr_config->bits_per_color=8;
  bsr_config->camera_icrs_x=0.0;
  bsr_config->camera_icrs_y=0.0;
  bsr_config->camera_icrs_z=0.0;
  bsr_config->camera_icrs_ra=0.0;
  bsr_config->camera_icrs_dec=0.0;
  bsr_config->camera_icrs_r=0.0;
  bsr_config->target_icrs_x=0.0;
  bsr_config->target_icrs_y=0.0;
  bsr_config->target_icrs_z=0.0;
  bsr_config->target_icrs_ra=266.4168371;
  bsr_config->target_icrs_dec=-29.0078106;
  bsr_config->target_icrs_r=8178.0;
  bsr_config->camera_rotation=-58.6;
  bsr_config->camera_pan=0.0;
  bsr_config->camera_tilt=0.0;
}

void cleanupValueStr(char *value) {
  int i;
  int j;
  int start;
  int end;
  char tmpvalue[256];
  size_t value_length;

  //
  // safety bounds check
  //
  value_length=strlen(value);
  if (value_length > 254) {
    value_length=254;
  }
  
  //
  // find beginning of value, ignoring leading spaces or quotes
  //
  start=-1;
  for (i=0; ((start == -1) && (i < value_length)); i++) {
    if ((value[i] != 32) && (value[i] != 34) && (value[i] != 39)) { 
      start=i;
    }
  }

  //
  // find end of value, ignoring trailing spaces or quotes
  //
  end=-1;
  for (i=(value_length-1); ((end == -1) && (i >= 0)); i--) {
    if ((value[i] != 32) && (value[i] != 34) && (value[i] != 39)) {
      end=i;
    }
  }

  //
  // trim before and after value
  //
  j=0;
  for (i=start; i <= end; i++) {
    tmpvalue[j]=value[i];
    j++;
  }
  tmpvalue[j]=0;

  strcpy(value, tmpvalue);
}

void checkOptionBool(int *config_int, char *option, char *value, char *matchstr) {
  size_t matchstr_length;

  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) == option) && (option[matchstr_length] != '_')) {
    if (strcasestr(value, "yes") != NULL) {
      *config_int=1;
    } else {
      *config_int=0;
    }
  }
}

void checkOptionInt(int *config_int, char *option, char *value, char *matchstr) {
  size_t matchstr_length;
  
  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) == option) && (option[matchstr_length] != '_')) {
    *config_int=strtol(value, NULL, 10);
  }
}

void checkOptionDouble(double *config_double, char *option, char *value, char *matchstr) {
  size_t matchstr_length;
  
  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) == option) && (option[matchstr_length] != '_')) {
    *config_double=strtod(value, NULL);
  }
}

void checkOptionStr(char *config_str,  char *option, char *value, char *matchstr) {
  size_t matchstr_length;
  
  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) == option) && (option[matchstr_length] != '_')) {
    strcpy(config_str, value);
  }
}

void setOptionValue(bsr_config_t *bsr_config, char *option, char *value, int from_cgi) {
  //
  // privileged values that can be set from config file or command line only (not cgi query_string)
  //
  if (from_cgi == 0) {
    checkOptionStr(bsr_config->data_file_directory, option, value, "data_file_directory");
    checkOptionInt(&bsr_config->num_threads, option, value, "num_threads");
    checkOptionInt(&bsr_config->per_thread_buffer, option, value, "per_thread_buffer");
    checkOptionInt(&bsr_config->per_thread_buffer_Airy, option, value, "per_thread_buffer_Airy");
    checkOptionBool(&bsr_config->cgi_mode, option, value, "cgi_mode");
    checkOptionInt(&bsr_config->cgi_max_res_x, option, value, "cgi_max_res_x");
    checkOptionInt(&bsr_config->cgi_max_res_y, option, value, "cgi_max_res_y");
    checkOptionInt(&bsr_config->cgi_Gaia_min_parallax_quality, option, value, "cgi_Gaia_min_parallax_quality");
    checkOptionBool(&bsr_config->cgi_allow_Airy_disk, option, value, "cgi_allow_Airy_disk");
    checkOptionDouble(&bsr_config->cgi_min_Airy_disk_first_null, option, value, "cgi_min_Airy_disk_first_null");
    checkOptionInt(&bsr_config->cgi_max_Airy_disk_max_extent, option, value, "cgi_max_Airy_disk_max_extent");
  }

  //
  // values that can be set from config file, command line, or cgi query_string
  //
  checkOptionBool(&bsr_config->enable_Gaia, option, value, "enable_Gaia");
  checkOptionInt(&bsr_config->Gaia_min_parallax_quality, option, value, "Gaia_min_parallax_quality");
  checkOptionBool(&bsr_config->enable_external, option, value, "enable_external");
  checkOptionDouble(&bsr_config->render_distance_min, option, value, "render_distance_min");
  bsr_config->render_distance_min2=bsr_config->render_distance_min * bsr_config->render_distance_min;
  checkOptionDouble(&bsr_config->render_distance_max, option, value, "render_distance_max");
  bsr_config->render_distance_max2=bsr_config->render_distance_max * bsr_config->render_distance_max;
  checkOptionInt(&bsr_config->render_distance_selector, option, value, "render_distance_selector");
  checkOptionDouble(&bsr_config->star_color_min, option, value, "star_color_min");
  checkOptionDouble(&bsr_config->star_color_max, option, value, "star_color_max");
  checkOptionInt(&bsr_config->camera_res_x, option, value, "camera_res_x");
  checkOptionInt(&bsr_config->camera_res_y, option, value, "camera_res_y");
  checkOptionDouble(&bsr_config->camera_fov, option, value, "camera_fov");
  checkOptionDouble(&bsr_config->camera_pixel_limit_mag, option, value, "camera_pixel_limit_mag");
  bsr_config->camera_pixel_limit=pow(100.0, (-bsr_config->camera_pixel_limit_mag / 5.0));
  checkOptionInt(&bsr_config->camera_pixel_limit_mode, option, value, "camera_pixel_limit_mode");
  checkOptionBool(&bsr_config->camera_wb_enable, option, value, "camera_wb_enable");
  checkOptionDouble(&bsr_config->camera_wb_temp, option, value, "camera_wb_temp");
  checkOptionDouble(&bsr_config->camera_color_saturation, option, value, "camera_color_saturation");
  checkOptionDouble(&bsr_config->camera_gamma, option, value, "camera_gamma");
  checkOptionInt(&bsr_config->camera_projection, option, value, "camera_projection");
  checkOptionInt(&bsr_config->spherical_orientation, option, value, "spherical_orientation");
  checkOptionInt(&bsr_config->Mollewide_iterations, option, value, "Mollewide_iterations");
  checkOptionDouble(&bsr_config->red_filter_long_limit, option, value, "red_filter_long_limit");
  checkOptionDouble(&bsr_config->red_filter_short_limit, option, value, "red_filter_short_limit");
  checkOptionDouble(&bsr_config->green_filter_long_limit, option, value, "green_filter_long_limit");
  checkOptionDouble(&bsr_config->green_filter_short_limit, option, value, "green_filter_short_limit");
  checkOptionDouble(&bsr_config->blue_filter_long_limit, option, value, "blue_filter_long_limit");
  checkOptionDouble(&bsr_config->blue_filter_short_limit, option, value, "blue_filter_short_limit");
  checkOptionBool(&bsr_config->Airy_disk, option, value, "Airy_disk");
  checkOptionDouble(&bsr_config->Airy_disk_first_null, option, value, "Airy_disk_first_null");
  checkOptionInt(&bsr_config->Airy_disk_max_extent, option, value, "Airy_disk_max_extent");
  checkOptionInt(&bsr_config->Airy_disk_min_extent, option, value, "Airy_disk_min_extent");
  checkOptionDouble(&bsr_config->Airy_disk_obstruction, option, value, "Airy_disk_obstruction");
  checkOptionBool(&bsr_config->skyglow_enable, option, value, "skyglow_enable");
  checkOptionDouble(&bsr_config->skyglow_temp, option, value, "skyglow_temp");
  checkOptionDouble(&bsr_config->skyglow_per_pixel_mag, option, value, "skyglow_per_pixel_mag");
  checkOptionDouble(&bsr_config->Gaussian_blur_radius, option, value, "Gaussian_blur_radius");
  checkOptionDouble(&bsr_config->output_scaling_factor, option, value, "output_scaling_factor");
  checkOptionBool(&bsr_config->draw_crosshairs, option, value, "draw_crosshairs");
  checkOptionBool(&bsr_config->draw_grid_lines, option, value, "draw_grid_lines");
  checkOptionBool(&bsr_config->sRGB_gamma, option, value, "sRGB_gamma");
  checkOptionInt(&bsr_config->bits_per_color, option, value, "bits_per_color");
  checkOptionDouble(&bsr_config->camera_icrs_x, option, value, "camera_icrs_x");
  checkOptionDouble(&bsr_config->camera_icrs_y, option, value, "camera_icrs_y");
  checkOptionDouble(&bsr_config->camera_icrs_z, option, value, "camera_icrs_z");
  checkOptionDouble(&bsr_config->camera_icrs_ra, option, value, "camera_icrs_ra");
  checkOptionDouble(&bsr_config->camera_icrs_dec, option, value, "camera_icrs_dec");
  checkOptionDouble(&bsr_config->camera_icrs_r, option, value, "camera_icrs_r");
  checkOptionDouble(&bsr_config->target_icrs_x, option, value, "target_icrs_x");
  checkOptionDouble(&bsr_config->target_icrs_y, option, value, "target_icrs_y");
  checkOptionDouble(&bsr_config->target_icrs_z, option, value, "target_icrs_z");
  checkOptionDouble(&bsr_config->target_icrs_ra, option, value, "target_icrs_ra");
  checkOptionDouble(&bsr_config->target_icrs_dec, option, value, "target_icrs_dec");
  checkOptionDouble(&bsr_config->target_icrs_r, option, value, "target_icrs_r");
  checkOptionDouble(&bsr_config->camera_rotation, option, value, "camera_rotation");
  checkOptionDouble(&bsr_config->camera_pan, option, value, "camera_pan");
  checkOptionDouble(&bsr_config->camera_tilt, option, value, "camera_tilt");
}

int loadConfigFromFile(bsr_config_t *bsr_config) {
  FILE *config_file;
  char *input_line_p;
  char *symbol_p;
  size_t input_line_length;
  char input_line[256];
  char input_line_trimmed[256];
  size_t input_line_trimmed_length;
  char option[256];
  size_t option_length;
  char value[256];
  size_t value_length;
  
  // attempt to open config file
/*
  // need to add cgi auto-detect here
  printf("Loading configuration from %s\n", bsr_config->config_file_name);
  fflush(stdout);
*/
  config_file=fopen(bsr_config->config_file_name, "r");
  if (config_file == NULL) {
/*
    printf("Warning: could not open %s\n", bsr_config->config_file_name);
    fflush(stdout);
*/
    return(0);
  }

  //
  // read and process each line of config file
  //
  input_line_p=fgets(input_line, 256, config_file);
  while(input_line_p != NULL) {
    input_line_length=strlen(input_line);

    //
    // search for comment symbol and remove any comments, otherwise just remove newline
    //
    symbol_p=strchr(input_line, '#');
    if (symbol_p != NULL) {
      input_line_trimmed_length=(symbol_p - input_line);
    } else {
      input_line_trimmed_length=input_line_length-1;
    } 
    strncpy(input_line_trimmed, input_line, input_line_trimmed_length);
    input_line_trimmed[input_line_trimmed_length]=0;

    //
    // search for option/value delimiter and split option and value strings
    //
    symbol_p=strchr(input_line_trimmed, '=');
    if (symbol_p != NULL) {
      option_length=(symbol_p - input_line_trimmed);
      value_length=(input_line_trimmed_length - option_length);
      // enforce range restrictions on option and value
      if ((option_length < 254) && (value_length < 254)) {
        strncpy(option, input_line_trimmed, option_length);
        option[option_length]=0;
        strncpy(value, (symbol_p+=1), value_length);
        value[value_length]=0;

        // remove non-alphanumeric characters before and after value
        cleanupValueStr(value);

        // send to option value processing fucntion
        setOptionValue(bsr_config, option, value, 0); // 0 == not from cgi

      } // end option_length and value_length checks
    } // end symbol_p check
    input_line_p=fgets(input_line, 256, config_file);
  } // end while input_line_raw

  return(0);
}

int loadConfigFromQueryString(bsr_config_t *bsr_config, char *query_string) {
  int done;
  char *query_p;
  char segment[2048];
  int segment_length;
  char *symbol_p;
  char option[256];
  size_t option_length;
  char value[256];
  size_t value_length;

  //
  // load first segment from query_string
  //
  done=0;
  query_p=query_string;
  if (query_p[0] == 0) {
    done=1;
  } else {
    symbol_p=strchr(query_p, '&');
    if (symbol_p != NULL) {
      segment_length=symbol_p - query_p;
    } else {
      segment_length=strlen(query_p);
    }
  }
  while (done == 0) {
    strncpy(segment, query_p, segment_length);
    segment[segment_length]=0;

    //
    // search for option/value delimiter and split option and value strings
    //
    symbol_p=strchr(segment, '=');
    if ((segment_length > 0) && (symbol_p != NULL)) {
      option_length=(symbol_p - segment);
      value_length=(segment_length - option_length);

      //
      // enforce range restrictions on option and value
      //
      if ((option_length < 254) && (value_length < 254)) {
        strncpy(option, segment, option_length);
        option[option_length]=0;
        strncpy(value, (symbol_p+=1), value_length);
        value[value_length]=0;

        //
        // remove non-alphanumeric characters before and after value
        //
        cleanupValueStr(value);

        //
        // send to option value processing fucntion
        //
        setOptionValue(bsr_config, option, value, 1); // 1 == from cgi

      } // end option_length and value_length checks
    } // end symbol_p check

    //
    // load next segment
    //
    query_p+=segment_length;
    if (query_p[0] == 0) {
      done=1;
    } else {
      query_p++;
      symbol_p=strchr(query_p, '&');
      if (symbol_p != NULL) {
        segment_length=symbol_p - query_p;
      } else {
        segment_length=strlen(query_p);
      }
    }
  } // end while not done

  return(0);
}

int validateCGIOptions(bsr_config_t *bsr_config) {
  if (bsr_config->camera_res_x < 1) {
    bsr_config->camera_res_x=1;
  }
  if (bsr_config->camera_res_x > bsr_config->cgi_max_res_x) {
    bsr_config->camera_res_x=bsr_config->cgi_max_res_x;
  }
  if (bsr_config->camera_res_y < 1) {
    bsr_config->camera_res_y=1;
  }
  if (bsr_config->camera_res_y > bsr_config->cgi_max_res_y) {
    bsr_config->camera_res_y=bsr_config->cgi_max_res_y;
  }
  if (bsr_config->Gaia_min_parallax_quality < bsr_config->cgi_Gaia_min_parallax_quality) {
    bsr_config->Gaia_min_parallax_quality=bsr_config->cgi_Gaia_min_parallax_quality;
  }
  if (bsr_config->cgi_allow_Airy_disk == 0) {
    bsr_config->Airy_disk=0;
  }
  if (bsr_config->Airy_disk_first_null < bsr_config->cgi_min_Airy_disk_first_null) {
    bsr_config->Airy_disk_first_null=bsr_config->cgi_min_Airy_disk_first_null;
  }
  if (bsr_config->Airy_disk_max_extent > bsr_config->cgi_max_Airy_disk_max_extent) {
    bsr_config->Airy_disk_max_extent=bsr_config->cgi_max_Airy_disk_max_extent;
  }
  if (bsr_config->Airy_disk_min_extent > bsr_config->cgi_max_Airy_disk_max_extent) {
    bsr_config->Airy_disk_min_extent=bsr_config->cgi_max_Airy_disk_max_extent;
  }

  return(0);
}

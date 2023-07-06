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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "usage.h"

void initConfig(bsr_config_t *bsr_config) {
  bsr_config->bsrender_cfg_version[0]=0;
  bsr_config->QUERY_STRING_p=NULL;
  strncpy(bsr_config->config_file_name, "bsrender.cfg", 255);
  bsr_config->config_file_name[255]=0;
  strncpy(bsr_config->data_file_directory, "galaxydata", 255);
  bsr_config->data_file_directory[255]=0;
  strncpy(bsr_config->output_file_name, "galaxy.png", 255);
  bsr_config->output_file_name[255]=0;
  bsr_config->print_status=1;
  bsr_config->num_threads=16;
  bsr_config->per_thread_buffer=1000;
  bsr_config->per_thread_buffer_Airy=100000;
  bsr_config->cgi_mode=0;
  bsr_config->cgi_max_res_x=999999;
  bsr_config->cgi_max_res_y=999999;
  bsr_config->cgi_Gaia_min_parallax_quality=0;
  bsr_config->cgi_allow_Airy_disk=1;
  bsr_config->cgi_min_Airy_disk_first_null=0.3;
  bsr_config->cgi_max_Airy_disk_max_extent=1000;
  bsr_config->cgi_max_Airy_disk_min_extent=3;
  bsr_config->cgi_allow_anti_alias=1;
  bsr_config->Gaia_db_enable=1;
  bsr_config->Gaia_min_parallax_quality=0;
  bsr_config->external_db_enable=1;
  bsr_config->render_distance_min=0.0;
  bsr_config->render_distance_max=1.0E99;
  bsr_config->render_distance_selector=0;
  bsr_config->star_intensity_min=1.0E99;
  bsr_config->star_intensity_max=-1.0E99;
  bsr_config->star_intensity_selector=0;
  bsr_config->star_color_min=0.0;
  bsr_config->star_color_max=1.0E99;
  bsr_config->extinction_dimming_undo=0;
  bsr_config->extinction_reddening_undo=0;
  bsr_config->camera_res_x=4000;
  bsr_config->camera_res_y=2000;
  bsr_config->camera_fov=360.0;
  bsr_config->camera_pixel_limit_mag=8.0;
  bsr_config->camera_pixel_limit_mode=-1;
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
  bsr_config->Airy_disk_enable=0;
  bsr_config->Airy_disk_first_null=0.75;
  bsr_config->Airy_disk_max_extent=100;
  bsr_config->Airy_disk_min_extent=1;
  bsr_config->Airy_disk_obstruction=0.0;
  bsr_config->anti_alias_enable=0;
  bsr_config->anti_alias_radius=1.0;
  bsr_config->skyglow_enable=0;
  bsr_config->skyglow_temp=4500.0;
  bsr_config->skyglow_per_pixel_mag=14.0;
  bsr_config->pre_limit_intensity=1;
  bsr_config->Gaussian_blur_radius=0.0;
  bsr_config->output_scaling_factor=1.0;
  bsr_config->Lanczos_order=3;
  bsr_config->draw_crosshairs=0;
  bsr_config->draw_grid_lines=0;
  bsr_config->output_format=0;
  bsr_config->color_profile=-1;
  bsr_config->exr_compression=3;
  bsr_config->compression_quality=80;
  bsr_config->image_format=0;
  bsr_config->hdr_neutral_white_ref=100;
  bsr_config->bits_per_color=8;
  bsr_config->image_number_format=0;
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
  strncpy(value, tmpvalue, 255);
  value[255]=0;
}

int checkOptionBool(int *config_int, char *option, char *value, char *matchstr) {
  int match=0;

  if ((strcasestr(option, matchstr) == option) && (strlen(option) == strlen(matchstr))) {
    match=1;
    if (strcasestr(value, "yes") != NULL) {
      *config_int=1;
    } else {
      *config_int=0;
    }
  }

  return(match);
}

int checkOptionInt(int *config_int, char *option, char *value, char *matchstr) {
  int match=0;

  if ((strcasestr(option, matchstr) == option) && (strlen(option) == strlen(matchstr))) {
    match=1;
    *config_int=strtol(value, NULL, 10);
  }

  return(match);
}

int checkOptionDouble(double *config_double, char *option, char *value, char *matchstr) {
  int match=0;

  if ((strcasestr(option, matchstr) == option) && (strlen(option) == strlen(matchstr))) {
    match=1;
    *config_double=strtod(value, NULL);
  }

  return(match);
}

int checkOptionStr(char *config_str,  char *option, char *value, char *matchstr) {
  int match=0;

  if ((strcasestr(option, matchstr) == option) && (strlen(option) == strlen(matchstr))) {
    match=1;
    strncpy(config_str, value, 255);
    config_str[255]=0;
  }

  return(match);
}

void setOptionValue(bsr_config_t *bsr_config, char *option, char *value, int from_cgi) {
  int match_count=0;

  //
  // privileged values that can be set from config file or command line only (not CGI QUERY_STRING)
  // note: config file name is additionally privileged in that it can only be set from command line
  // so it is not even processed in this function.
  //
  if (from_cgi == 0) {
    match_count+=checkOptionStr(bsr_config->bsrender_cfg_version, option, value, "bsrender_cfg_version");
    match_count+=checkOptionStr(bsr_config->data_file_directory, option, value, "data_file_directory");
    match_count+=checkOptionStr(bsr_config->output_file_name, option, value, "output_file_name");
    match_count+=checkOptionBool(&bsr_config->print_status, option, value, "print_status");
    match_count+=checkOptionInt(&bsr_config->num_threads, option, value, "num_threads");
    match_count+=checkOptionInt(&bsr_config->per_thread_buffer, option, value, "per_thread_buffer");
    match_count+=checkOptionInt(&bsr_config->per_thread_buffer_Airy, option, value, "per_thread_buffer_Airy");
    match_count+=checkOptionBool(&bsr_config->cgi_mode, option, value, "cgi_mode");
    match_count+=checkOptionInt(&bsr_config->cgi_max_res_x, option, value, "cgi_max_res_x");
    match_count+=checkOptionInt(&bsr_config->cgi_max_res_y, option, value, "cgi_max_res_y");
    match_count+=checkOptionInt(&bsr_config->cgi_Gaia_min_parallax_quality, option, value, "cgi_Gaia_min_parallax_quality");
    match_count+=checkOptionBool(&bsr_config->cgi_allow_Airy_disk, option, value, "cgi_allow_Airy_disk");
    match_count+=checkOptionDouble(&bsr_config->cgi_min_Airy_disk_first_null, option, value, "cgi_min_Airy_disk_first_null");
    match_count+=checkOptionInt(&bsr_config->cgi_max_Airy_disk_max_extent, option, value, "cgi_max_Airy_disk_max_extent");
    match_count+=checkOptionInt(&bsr_config->cgi_max_Airy_disk_min_extent, option, value, "cgi_max_Airy_disk_min_extent");
    match_count+=checkOptionBool(&bsr_config->cgi_allow_anti_alias, option, value, "cgi_allow_anti_alias");
  }

  //
  // values that can be set from config file, command line, or CGI QUERY_STRING
  //
  match_count+=checkOptionBool(&bsr_config->Gaia_db_enable, option, value, "Gaia_db_enable");
  match_count+=checkOptionInt(&bsr_config->Gaia_min_parallax_quality, option, value, "Gaia_min_parallax_quality");
  match_count+=checkOptionBool(&bsr_config->external_db_enable, option, value, "external_db_enable");
  match_count+=checkOptionDouble(&bsr_config->render_distance_min, option, value, "render_distance_min");
  match_count+=checkOptionDouble(&bsr_config->render_distance_max, option, value, "render_distance_max");
  match_count+=checkOptionInt(&bsr_config->render_distance_selector, option, value, "render_distance_selector");
  match_count+=checkOptionDouble(&bsr_config->star_intensity_min, option, value, "star_intensity_min");
  match_count+=checkOptionDouble(&bsr_config->star_intensity_max, option, value, "star_intensity_max");
  match_count+=checkOptionInt(&bsr_config->star_intensity_selector, option, value, "star_intensity_selector");
  match_count+=checkOptionDouble(&bsr_config->star_color_min, option, value, "star_color_min");
  match_count+=checkOptionDouble(&bsr_config->star_color_max, option, value, "star_color_max");
  match_count+=checkOptionBool(&bsr_config->extinction_dimming_undo, option, value, "extinction_dimming_undo");
  match_count+=checkOptionBool(&bsr_config->extinction_reddening_undo, option, value, "extinction_reddening_undo");
  match_count+=checkOptionInt(&bsr_config->camera_res_x, option, value, "camera_res_x");
  match_count+=checkOptionInt(&bsr_config->camera_res_y, option, value, "camera_res_y");
  match_count+=checkOptionDouble(&bsr_config->camera_fov, option, value, "camera_fov");
  match_count+=checkOptionDouble(&bsr_config->camera_pixel_limit_mag, option, value, "camera_pixel_limit_mag");
  match_count+=checkOptionInt(&bsr_config->camera_pixel_limit_mode, option, value, "camera_pixel_limit_mode");
  match_count+=checkOptionBool(&bsr_config->camera_wb_enable, option, value, "camera_wb_enable");
  match_count+=checkOptionDouble(&bsr_config->camera_wb_temp, option, value, "camera_wb_temp");
  match_count+=checkOptionDouble(&bsr_config->camera_color_saturation, option, value, "camera_color_saturation");
  match_count+=checkOptionDouble(&bsr_config->camera_gamma, option, value, "camera_gamma");
  match_count+=checkOptionInt(&bsr_config->camera_projection, option, value, "camera_projection");
  match_count+=checkOptionInt(&bsr_config->spherical_orientation, option, value, "spherical_orientation");
  match_count+=checkOptionInt(&bsr_config->Mollewide_iterations, option, value, "Mollewide_iterations");
  match_count+=checkOptionDouble(&bsr_config->red_filter_long_limit, option, value, "red_filter_long_limit");
  match_count+=checkOptionDouble(&bsr_config->red_filter_short_limit, option, value, "red_filter_short_limit");
  match_count+=checkOptionDouble(&bsr_config->green_filter_long_limit, option, value, "green_filter_long_limit");
  match_count+=checkOptionDouble(&bsr_config->green_filter_short_limit, option, value, "green_filter_short_limit");
  match_count+=checkOptionDouble(&bsr_config->blue_filter_long_limit, option, value, "blue_filter_long_limit");
  match_count+=checkOptionDouble(&bsr_config->blue_filter_short_limit, option, value, "blue_filter_short_limit");
  match_count+=checkOptionBool(&bsr_config->Airy_disk_enable, option, value, "Airy_disk_enable");
  match_count+=checkOptionDouble(&bsr_config->Airy_disk_first_null, option, value, "Airy_disk_first_null");
  match_count+=checkOptionInt(&bsr_config->Airy_disk_max_extent, option, value, "Airy_disk_max_extent");
  match_count+=checkOptionInt(&bsr_config->Airy_disk_min_extent, option, value, "Airy_disk_min_extent");
  match_count+=checkOptionDouble(&bsr_config->Airy_disk_obstruction, option, value, "Airy_disk_obstruction");
  match_count+=checkOptionBool(&bsr_config->anti_alias_enable, option, value, "anti_alias_enable");
  match_count+=checkOptionDouble(&bsr_config->anti_alias_radius, option, value, "anti_alias_radius");
  match_count+=checkOptionBool(&bsr_config->skyglow_enable, option, value, "skyglow_enable");
  match_count+=checkOptionDouble(&bsr_config->skyglow_temp, option, value, "skyglow_temp");
  match_count+=checkOptionDouble(&bsr_config->skyglow_per_pixel_mag, option, value, "skyglow_per_pixel_mag");
  match_count+=checkOptionBool(&bsr_config->pre_limit_intensity, option, value, "pre_limit_intensity");
  match_count+=checkOptionDouble(&bsr_config->Gaussian_blur_radius, option, value, "Gaussian_blur_radius");
  match_count+=checkOptionDouble(&bsr_config->output_scaling_factor, option, value, "output_scaling_factor");
  match_count+=checkOptionInt(&bsr_config->Lanczos_order, option, value, "Lanczos_order");
  match_count+=checkOptionBool(&bsr_config->draw_crosshairs, option, value, "draw_crosshairs");
  match_count+=checkOptionBool(&bsr_config->draw_grid_lines, option, value, "draw_grid_lines");
  match_count+=checkOptionInt(&bsr_config->output_format, option, value, "output_format");
  match_count+=checkOptionInt(&bsr_config->color_profile, option, value, "color_profile");
  match_count+=checkOptionInt(&bsr_config->exr_compression, option, value, "exr_compression");
  match_count+=checkOptionInt(&bsr_config->compression_quality, option, value, "compression_quality");
  match_count+=checkOptionInt(&bsr_config->hdr_neutral_white_ref, option, value, "hdr_neutral_white_ref");
  match_count+=checkOptionDouble(&bsr_config->camera_icrs_x, option, value, "camera_icrs_x");
  match_count+=checkOptionDouble(&bsr_config->camera_icrs_y, option, value, "camera_icrs_y");
  match_count+=checkOptionDouble(&bsr_config->camera_icrs_z, option, value, "camera_icrs_z");
  match_count+=checkOptionDouble(&bsr_config->camera_icrs_ra, option, value, "camera_icrs_ra");
  match_count+=checkOptionDouble(&bsr_config->camera_icrs_dec, option, value, "camera_icrs_dec");
  match_count+=checkOptionDouble(&bsr_config->camera_icrs_r, option, value, "camera_icrs_r");
  match_count+=checkOptionDouble(&bsr_config->target_icrs_x, option, value, "target_icrs_x");
  match_count+=checkOptionDouble(&bsr_config->target_icrs_y, option, value, "target_icrs_y");
  match_count+=checkOptionDouble(&bsr_config->target_icrs_z, option, value, "target_icrs_z");
  match_count+=checkOptionDouble(&bsr_config->target_icrs_ra, option, value, "target_icrs_ra");
  match_count+=checkOptionDouble(&bsr_config->target_icrs_dec, option, value, "target_icrs_dec");
  match_count+=checkOptionDouble(&bsr_config->target_icrs_r, option, value, "target_icrs_r");
  match_count+=checkOptionDouble(&bsr_config->camera_rotation, option, value, "camera_rotation");
  match_count+=checkOptionDouble(&bsr_config->camera_pan, option, value, "camera_pan");
  match_count+=checkOptionDouble(&bsr_config->camera_tilt, option, value, "camera_tilt");

  if ((bsr_config->QUERY_STRING_p == NULL) && (bsr_config->print_status == 1)) {
    if (match_count == 0) {
      printf("Unknown configuration option: %s\n", option);
      fflush(stdout);
    } else if (match_count > 1) {
      printf("Warning, ambiguous configuration option: %s\n", option); // should never happen
      fflush(stdout);
    }
  }
}

int processConfigSegment(bsr_config_t *bsr_config, char *segment, int from_cgi) {
  int segment_length;
  char *symbol_p;
  char option[256];
  size_t option_length;
  char value[256];
  size_t value_length;

  //
  // special handling for --help
  //
  if ((strcasestr(segment, "help") == segment) && (strlen(segment) == strlen("help"))) {
    printUsage();
    exit(0);
  }

  //
  // search for option/value delimiter and split option and value strings
  //
  segment_length=strlen(segment);
  symbol_p=strchr(segment, '=');
  if ((segment_length > 0) && (symbol_p != NULL)) {
    option_length=(symbol_p - segment);
    value_length=(segment_length - option_length);

    //
    // enforce range restrictions on option and value
    //
    if (option_length > 255) {
      option_length=255;
    }
    strncpy(option, segment, option_length);
    option[option_length]=0;
    if (value_length > 255) {
      value_length=255;
    }
    strncpy(value, (symbol_p+=1), value_length);
    value[value_length]=0;

    //
    // remove non-alphanumeric characters before and after value
    //
    cleanupValueStr(value);

    //
    // send to option value processing fucntion
    //
    setOptionValue(bsr_config, option, value, from_cgi);
  } else if ((segment_length > 0) && (segment[0] != 32) && (bsr_config->QUERY_STRING_p == NULL) && (bsr_config->print_status == 1)) {
    printf("Unknown configuration option: %s\n", segment);
    fflush(stdout);
  } // end symbol_p check

  return(0);
}

int loadConfigFromFile(bsr_config_t *bsr_config) {
  int from_cgi;
  FILE *config_file;
  char *input_line_p;
  char *symbol_p;
  size_t input_line_length;
  char input_line[256];
  char segment[256];
  size_t segment_length;
  
  //
  // attempt to open config file
  //
  config_file=fopen(bsr_config->config_file_name, "r");
  if ((config_file == NULL) && (bsr_config->QUERY_STRING_p == NULL)) {
    printf("Warning: could not open %s\n", bsr_config->config_file_name);
    fflush(stdout);
    return(0);
  }
  if ((bsr_config->QUERY_STRING_p == NULL) && (bsr_config->print_status == 1)) {
    printf("Loading configuration file %s\n", bsr_config->config_file_name);
    fflush(stdout);
  }

  //
  // attempt to open config file
  //
  config_file=fopen(bsr_config->config_file_name, "r");
  if ((config_file == NULL) && (bsr_config->QUERY_STRING_p == NULL)) {
    printf("Warning: could not open %s\n", bsr_config->config_file_name);
    fflush(stdout);
    return(0);
  }

  //
  // read and process each line of config file
  //
  input_line_p=fgets(input_line, 256, config_file);
  while(input_line_p != NULL) {
    //
    // search for comment symbol and remove any comments, otherwise just remove newline
    //
    input_line_length=strlen(input_line);
    symbol_p=strchr(input_line, '#');
    if (symbol_p != NULL) {
      segment_length=(symbol_p - input_line);
    } else {
      segment_length=input_line_length-1;
    } 
    if (segment_length > 255) {
      segment_length=255;
    }
    strncpy(segment, input_line, segment_length);
    segment[segment_length]=0;

    //
    // process config segment
    //
    from_cgi=0;
    processConfigSegment(bsr_config, segment, from_cgi);

    //
    // load next line from config file
    //
    input_line_p=fgets(input_line, 256, config_file);
  } // end while input_line_raw

  return(0);
}

int loadConfigFromQueryString(bsr_config_t *bsr_config, char *query_string_2048) {
  int done;
  int from_cgi;
  char *query_p;
  char segment[2048];
  int segment_length;
  char *symbol_p;

  //
  // load first segment from query_string_2048
  //
  done=0;
  query_p=query_string_2048;
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

  //
  // loop for processing segments from query string
  //
  while (done == 0) {
    //
    // process config segment
    //
    if (segment_length > 2047) {
      segment_length=2047;
    }
    strncpy(segment, query_p, segment_length);
    segment[segment_length]=0;
    from_cgi=1;
    processConfigSegment(bsr_config, segment, from_cgi);

    //
    // load next segment from query string
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

int processCmdArgs(bsr_config_t *bsr_config, int argc, char **argv) {
  int i;
  int from_cgi;
  char *option_start;

  if (argc == 1) {
    return(0);
  } else {
    for (i=1; i <= (argc - 1); i++) {
      if (argv[i][1] == 'c') {
        // configuration file name
        if (argv[i][2] != 0) {
          // option concatenated onto switch
          option_start=argv[i];
          strncpy(bsr_config->config_file_name, (option_start + (size_t)2), 255);
          bsr_config->config_file_name[255]=0;
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          strncpy(bsr_config->config_file_name, option_start, 255);
          bsr_config->config_file_name[255]=0;
        } // end if no space
      } else if (argv[i][1] == 'd') {
        // data files directory
        if (argv[i][2] != 0) {
          // option concatenated onto switch
          option_start=argv[i];
          strncpy(bsr_config->data_file_directory, (option_start + (size_t)2), 255);
          bsr_config->data_file_directory[255]=0;
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          strncpy(bsr_config->data_file_directory, option_start, 255);
          bsr_config->data_file_directory[255]=0;
        } // end if no space
      } else if (argv[i][1] == 'h') {
        // print help
        printUsage();
        exit(0);
      } else if (argv[i][1] == 'o') {
        // output filename
        if (argv[i][2] != 0) {
          // option concatenated onto switch
          option_start=argv[i];
          strncpy(bsr_config->output_file_name, (option_start + (size_t)2), 255);
          bsr_config->output_file_name[255]=0;
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          strncpy(bsr_config->output_file_name, option_start, 255);
          bsr_config->output_file_name[255]=0;
        } // end if no space
      } else if (argv[i][1] == 'q') {
        // quite mode - suppress non-error status messages
        bsr_config->print_status=0;
      } else if (argv[i][1] == '-') {
        // long config option
        from_cgi=0;
        processConfigSegment(bsr_config, (argv[i]+2), from_cgi);
      } else if ((bsr_config->QUERY_STRING_p == NULL) && (bsr_config->print_status == 1)) {
        printf("Unknown configuration option: %s\n", argv[i]);
        fflush(stdout);
      } // end which option
    } // end for argc
  } // end if any options
  return(0);
}

int validateConfig(bsr_config_t *bsr_config) {
  //
  // minimum required threads is 2 (one main thread and one worker thread)
  //
  if (bsr_config->num_threads < 2) {
    bsr_config->num_threads=2;
  }

  //
  // translate output_format to internal config variables
  // 0 = PNG 8-bit unsigned integer per color
  // 1 = PNG 16-bit unsigned integer per color
  // 2 = EXR 16-bit floating-point per color
  // 3 = EXR 32-bit floating-point per color
  // 4 = EXR 32-bit unsigned integer per color
  // 5 = JPG 8-bit unsigned integer per color
  // 6 = AVIF 8-bit unsigned integer per color
  // 7 = AVIF 10-bit unsigned integer per color
  // 8 = AVIF 12-bit unsigned integer per color
  // 9 = AVIF 16-bit floating-point per color
  // 10 = HEIF 8-bit unsigned integer per color
  // 11 = HEIF 10-bit unsigned integer per color
  // 12 = HEIF 12-bit unsigned integer per color
  //
  if (bsr_config->output_format == 0) {
#ifndef BSR_USE_PNG
    printf("Error: not compiled with PNG support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=0;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=8;
  } else if (bsr_config->output_format == 1) {
#ifndef BSR_USE_PNG 
    printf("Error: not compiled with PNG support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=0;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=16;
  } else if (bsr_config->output_format == 2) {
#ifndef BSR_USE_EXR
    printf("Error: not compiled with EXR support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=1;
    bsr_config->image_number_format=1;
    bsr_config->bits_per_color=16;
  } else if (bsr_config->output_format == 3) {
#ifndef BSR_USE_EXR
    printf("Error: not compiled with EXR support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=1;
    bsr_config->image_number_format=1;
    bsr_config->bits_per_color=32;
  } else if (bsr_config->output_format == 4) {
#ifndef BSR_USE_EXR
    printf("Error: not compiled with EXR support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=1;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=32;
  } else if (bsr_config->output_format == 5) {
#ifndef BSR_USE_JPEG
    printf("Error: not compiled with JPEG support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=2;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=8;
  } else if (bsr_config->output_format == 6) {
#ifndef BSR_USE_AVIF
    printf("Error: not compiled with AVIF support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=3;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=8;
  } else if (bsr_config->output_format == 7) {
#ifndef BSR_USE_AVIF
    printf("Error: not compiled with AVIF support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=3;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=10;
  } else if (bsr_config->output_format == 8) {
#ifndef BSR_USE_AVIF
    printf("Error: not compiled with AVIF support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=3;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=12;
/*
  } else if (bsr_config->output_format == 9) {
#ifndef BSR_USE_AVIF
    printf("Error: not compiled with AVIF support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=3;
    bsr_config->image_number_format=1;
    bsr_config->bits_per_color=16;
*/
  } else if (bsr_config->output_format == 10) {
#ifndef BSR_USE_HEIF
    printf("Error: not compiled with HEIF support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=4;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=8;
  } else if (bsr_config->output_format == 11) {
#ifndef BSR_USE_HEIF
    printf("Error: not compiled with HEIF support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=4;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=10;
  } else if (bsr_config->output_format == 12) {
#ifndef BSR_USE_HEIF
    printf("Error: not compiled with HEIF support\n");
    fflush(stdout);
    exit(1);
#endif
    bsr_config->image_format=4;
    bsr_config->image_number_format=0;
    bsr_config->bits_per_color=12;
  } else {
    if ((bsr_config->QUERY_STRING_p == NULL) && (bsr_config->print_status == 1)) {
      printf("Error: invalid output_format (%d). See --help (Output section) for output format codes.\n", bsr_config->output_format);
      fflush(stdout);
    }
    exit(1);
  }

  //
  // libheif does not appear to support writing to stdout
  //
  if ((bsr_config->cgi_mode == 1) && ((bsr_config->output_format == 10) || (bsr_config->output_format == 11) || (bsr_config->output_format == 12))) {
    printf("Error: HEIF output format is not supported in CGI mode\n");
    fflush(stdout);
    exit(1);
  }

  //
  // change output filename if still default "galaxy.png" and non-PNG output format is configured
  //
  if ((bsr_config->image_format == 1)\
       && ((strstr(bsr_config->output_file_name, "galaxy.png") == bsr_config->output_file_name))\
       && (strlen(bsr_config->output_file_name) == strlen("galaxy.png"))) {
    strncpy(bsr_config->output_file_name, "galaxy.exr", 255);
    bsr_config->output_file_name[255]=0;
  } else if ((bsr_config->image_format == 2)\
       && ((strstr(bsr_config->output_file_name, "galaxy.png") == bsr_config->output_file_name))\
       && (strlen(bsr_config->output_file_name) == strlen("galaxy.png"))) {
    strncpy(bsr_config->output_file_name, "galaxy.jpg", 255);
    bsr_config->output_file_name[255]=0;
  } else if ((bsr_config->image_format == 3)\
       && ((strstr(bsr_config->output_file_name, "galaxy.png") == bsr_config->output_file_name))\
       && (strlen(bsr_config->output_file_name) == strlen("galaxy.png"))) {
    strncpy(bsr_config->output_file_name, "galaxy.avif", 255);
    bsr_config->output_file_name[255]=0;
  } else if ((bsr_config->image_format == 4)\
       && ((strstr(bsr_config->output_file_name, "galaxy.png") == bsr_config->output_file_name))\
       && (strlen(bsr_config->output_file_name) == strlen("galaxy.png"))) {
    strncpy(bsr_config->output_file_name, "galaxy.heif", 255);
    bsr_config->output_file_name[255]=0;
  }

  //
  // check camera_pixel_limit_mode
  //
  if (bsr_config->camera_pixel_limit_mode == -1) {
    // process default setting
    if (bsr_config->image_number_format == 0) { // integer format
      bsr_config->camera_pixel_limit_mode=0;
    } else if (bsr_config->image_number_format == 1) { // floating-point format
      bsr_config->camera_pixel_limit_mode=2;
    } 
  } else if ((bsr_config->camera_pixel_limit_mode == 2) && (bsr_config->image_number_format == 0)) {
    // integer formats must be limited to range [0..1]
    if ((bsr_config->QUERY_STRING_p == NULL) && (bsr_config->print_status == 1)) {
      printf("Warning: integer image formats require clipping pixel values above 1.0. Setting camera_pixel_limit_mode=0\n");
      fflush(stdout);
    }
    bsr_config->camera_pixel_limit_mode=0;
  }

  //
  // check color_profile
  //
  if (bsr_config->color_profile == -1) {
    // process default setting
    if ((bsr_config->image_format == 0) || (bsr_config->image_format == 2) || (bsr_config->image_format == 3) || (bsr_config->image_format == 4)) {
      bsr_config->color_profile=1;  // PNG,JPG,AVIF default is sRGB
    } else if (bsr_config->image_format == 1) {
      bsr_config->color_profile=0;  // EXR default is none
    }
  } else if (bsr_config->color_profile == 8) {
    if ((bsr_config->QUERY_STRING_p == NULL) && (bsr_config->print_status == 1)) {
      printf("HDR color profile selected: disabling pre_limit_intensity\n");
      fflush(stdout);
    }
    // disable pre_limit_intensity if using HDR color profile
    bsr_config->pre_limit_intensity=0;
  }

  return(0);
}

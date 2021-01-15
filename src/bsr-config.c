#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void initConfig(bsr_config_t *bsr_config) {
  strcpy(bsr_config->config_file_name, "./bsrender.cfg");
  strcpy(bsr_config->data_file_directory, "./galaxydata");
  bsr_config->num_threads=16;
  bsr_config->per_thread_buffer=20000;
  bsr_config->min_parallax_quality=10;
  bsr_config->render_distance_min=0.0;
  bsr_config->render_distance_max=1.0E99;
  bsr_config->render_distance_selector=0;
  bsr_config->draw_cross_hairs=0;
  bsr_config->draw_grid_lines=0;
  bsr_config->cgi_mode=0;
  bsr_config->cgi_max_res_x=32768;
  bsr_config->cgi_max_res_y=16384;
  bsr_config->cgi_min_parallax_quality=0;
  bsr_config->camera_res_x=1920;
  bsr_config->camera_res_y=1080;
  bsr_config->camera_fov=360.0;
  bsr_config->camera_wb_temp=4200.0;
  bsr_config->camera_pixel_limit_mag=7.0;
  bsr_config->camera_pixel_limit=pow(100.0, (-bsr_config->camera_pixel_limit_mag / 5.0));
  bsr_config->camera_pixel_limit_mode=0;
  bsr_config->camera_color_saturation=4;
  bsr_config->camera_projection=0;
  bsr_config->Mollewide_iterations=5;
  bsr_config->camera_gamma=1.0;
  bsr_config->sRGB_gamma=1;
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
  bsr_config->camera_rotation=-60.2;
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

  value_length=strlen(value);
  if (value_length > 254) {
    value_length=254;
  }
  
  // find beginning of real value
  start=-1;
  for (i=0; ((start == -1) && (i < value_length)); i++) {
    if ((value[i] != 32) && (value[i] != 34) && (value[i] != 39)) { 
      start=i;
    }
  }

  // find end of real value
  end=-1;
  for (i=(value_length-1); ((end == -1) && (i >= 0)); i--) {
    if ((value[i] != 32) && (value[i] != 34) && (value[i] != 39)) {
      end=i;
    }
  }

  // trim before and after real value
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

void setOptionValue(bsr_config_t *bsr_config, char *option, char *value) {
  checkOptionStr(bsr_config->data_file_directory, option, value, "data_file_directory");
  checkOptionInt(&bsr_config->num_threads, option, value, "num_threads");
  checkOptionInt(&bsr_config->per_thread_buffer, option, value, "per_thread_buffer");
  checkOptionInt(&bsr_config->min_parallax_quality, option, value, "min_parallax_quality");
  checkOptionDouble(&bsr_config->render_distance_min, option, value, "render_distance_min");
  checkOptionDouble(&bsr_config->render_distance_max, option, value, "render_distance_max");
  checkOptionInt(&bsr_config->render_distance_selector, option, value, "render_distance_selector");
  checkOptionBool(&bsr_config->draw_cross_hairs, option, value, "draw_cross_hairs");
  checkOptionBool(&bsr_config->draw_grid_lines, option, value, "draw_grid_lines");
  checkOptionBool(&bsr_config->cgi_mode, option, value, "cgi_mode");
  checkOptionInt(&bsr_config->cgi_max_res_x, option, value, "cgi_max_res_x");
  checkOptionInt(&bsr_config->cgi_max_res_y, option, value, "cgi_max_res_y");
  checkOptionInt(&bsr_config->cgi_min_parallax_quality, option, value, "cgi_min_parallax_quality");
  checkOptionInt(&bsr_config->camera_res_x, option, value, "camera_res_x");
  checkOptionInt(&bsr_config->camera_res_y, option, value, "camera_res_y");
  checkOptionDouble(&bsr_config->camera_fov, option, value, "camera_fov");
  checkOptionDouble(&bsr_config->camera_wb_temp, option, value, "camera_wb_temp");
  checkOptionDouble(&bsr_config->camera_pixel_limit_mag, option, value, "camera_pixel_limit_mag");
  bsr_config->camera_pixel_limit=pow(100.0, (-bsr_config->camera_pixel_limit_mag / 5.0));
  checkOptionInt(&bsr_config->camera_pixel_limit_mode, option, value, "camera_pixel_limit_mode");
  checkOptionDouble(&bsr_config->camera_color_saturation, option, value, "camera_color_saturation");
  checkOptionInt(&bsr_config->camera_projection, option, value, "camera_projection");
  checkOptionInt(&bsr_config->Mollewide_iterations, option, value, "Mollewide_iterations");
  checkOptionDouble(&bsr_config->camera_gamma, option, value, "camera_gamma");
  checkOptionBool(&bsr_config->sRGB_gamma, option, value, "sRGB_gamma");
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

int loadConfig(bsr_config_t *bsr_config) {
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
  printf("Loading configuration from %s\n", bsr_config->config_file_name);
  fflush(stdout);
  config_file=fopen(bsr_config->config_file_name, "r");
  if (config_file == NULL) {
    printf("Warning: could not open %s\n", bsr_config->config_file_name);
    fflush(stdout);
    return(0);
  }

  // read and process each line of config file
  input_line_p=fgets(input_line, 256, config_file);
  while(input_line_p != NULL) {
    input_line_length=strlen(input_line);

    // search for comment symbol and remove any comments, otherwise just remove newline
    symbol_p=strchr(input_line, '#');
    if (symbol_p != NULL) {
      input_line_trimmed_length=(symbol_p - input_line);
    } else {
      input_line_trimmed_length=input_line_length-1;
    } 
    strncpy(input_line_trimmed, input_line, input_line_trimmed_length);
    input_line_trimmed[input_line_trimmed_length]=0;

    // search for option/value delimiter and split option and value strings
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
        setOptionValue(bsr_config, option, value);

      } // end option_length and value_length checks
    } // end symbol_p check
    input_line_p=fgets(input_line, 256, config_file);
  } // end while input_line_raw

  return(0);
}
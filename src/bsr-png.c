#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdlib.h>
#include <math.h>
#define PNG_SETJMP_NOT_SUPPORTED
#include <png.h>
#include "util.h"
#include "cgi.h"

int writePNGFile(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  pixel_composition_t *current_image_p;
  png_byte *image_output_buf;
  png_byte *image_output_p;
  double pixel_r;
  double pixel_g;
  double pixel_b;
  int i;
  int output_x;
  int output_y;
  FILE *output_file=NULL;
  png_structp png_ptr;
  png_infop info_ptr;
  png_byte color_type=PNG_COLOR_TYPE_RGB;
  png_byte bit_depth=8;
  png_bytep *row_pointers;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  double one_over_2dot4;
  int output_res_x;
  int output_res_y;

  //
  // display status message if not in cgi mode
  //
  if (bsr_config->cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Converting to 8-bits per color...");
    fflush(stdout);
  }

  //
  // get current image pointer and resolution
  //
  current_image_p=bsr_state->current_image_buf;
  output_res_x=bsr_state->current_image_res_x;
  output_res_y=bsr_state->current_image_res_y;


  //
  // allocate memory for image output buffer (8 bits per color rgb)
  //
  image_output_buf=(png_byte *)malloc(output_res_x * output_res_y * 3 * sizeof(png_byte));
  if (image_output_buf == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate memory for image output buffer\n");
      fflush(stdout);
    }
    return(1);
  }
  image_output_p=image_output_buf;
/*
  for (i=0; i < (output_res_x * output_res_y); i++) {
    *image_output_p=0;
    image_output_p++;
    *image_output_p=0;
    image_output_p++;
    *image_output_p=0;
    image_output_p++;
  }
*/

  //
  // allocate memory for row_pointers
  //
  row_pointers=(png_bytep *)malloc(output_res_y * sizeof(png_bytep));
  if (row_pointers == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate memory for libpng row_pointers\n");
      fflush(stdout);
    }
    return(1);
  }

  //
  // convert double precision image buffer to 8-bits per color
  //
  image_output_p=image_output_buf;
  output_x=0;
  output_y=0;
  row_pointers[0]=image_output_p;
  one_over_2dot4=1.0 / 2.4;
  for (i=0; i < (output_res_x * output_res_y); i++) {

    //
    // set new row pointer if we have reached end of row
    //
    if (output_x == output_res_x) {
      output_x=0;
      output_y++;
      row_pointers[output_y]=image_output_p;
    }

    //
    // copy pixel data from current_image_buf
    //
    pixel_r=current_image_p->r;
    pixel_g=current_image_p->g;
    pixel_b=current_image_p->b;

    //
    // limit pixel intensity to range [0..1]
    //
    if (bsr_config->camera_pixel_limit_mode == 0) {
      limitIntensity(&pixel_r, &pixel_g, &pixel_b);
    } else if (bsr_config->camera_pixel_limit_mode == 1) {
      limitIntensityPreserveColor(&pixel_r, &pixel_g, &pixel_b);
    }

    //
    // optionally apply sRGB gamma
    //
    if (bsr_config->sRGB_gamma == 1) {
      // apply sRGB gamma
      if (pixel_r <= 0.0031308) {
        pixel_r=pixel_r * 12.92;
      } else {
        pixel_r=(1.055 * pow(pixel_r, one_over_2dot4) - 0.055);
      }
      if (pixel_g <= 0.0031308) {
        pixel_g=pixel_g * 12.92;
      } else {
        pixel_g=(1.055 * pow(pixel_g, one_over_2dot4) - 0.055);
      }
      if (pixel_b <= 0.0031308) {
        pixel_b=pixel_b * 12.92;
      } else {
        pixel_b=(1.055 * pow(pixel_b, one_over_2dot4) - 0.055);
      }
    }    

    //
    // convert r,g,b to 8 bit values and copy to 8-bit output buffer
    //
    *image_output_p=(unsigned char)((pixel_r * 255.0) + 0.5);
    image_output_p++;
    *image_output_p=(unsigned char)((pixel_g * 255.0) + 0.5);
    image_output_p++;
    *image_output_p=(unsigned char)((pixel_b * 255.0) + 0.5);
    image_output_p++;

    output_x++;
    current_image_p++;
  }

  //
  // display status update if not in cgi mode
  //
  if (bsr_config->cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.4fs)\n", elapsed_time);
    fflush(stdout);

    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Writing galaxy.png...");
    fflush(stdout);
  }

  //
  // if not cgi mode, open output file
  //
  if (bsr_config->cgi_mode != 1) {
    output_file=fopen("galaxy.png", "wb");
    if (output_file == NULL) {
      printf("Error: could not open galaxy.png for writing\n");
      fflush(stdout);
      return(1);
    }
  }
  
  //
  // initialize PNG library
  //
  png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info_ptr=png_create_info_struct(png_ptr);
  if (bsr_config->cgi_mode == 1) {
    png_init_io(png_ptr, stdout);
  } else {
    png_init_io(png_ptr, output_file);
  } 
  
  //
  // override default sRGB gamma header info if sRGB is disabled in config
  //
  if (bsr_config->sRGB_gamma == 0) {
    png_set_gAMA(png_ptr, info_ptr, 1.0);
  }

  //
  // write PNG header
  //
  png_set_IHDR(png_ptr, info_ptr, output_res_x, output_res_y, bit_depth, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  png_write_info(png_ptr, info_ptr);
  // wringe PNG image data
  png_write_image(png_ptr, row_pointers);
  // end write
  png_write_end(png_ptr, NULL);

  //
  // if not cgi mode display status message and close output file
  //
  if (bsr_config->cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.4fs)\n", elapsed_time);
    fflush(stdout);

    // clean up
    fclose(output_file);
  }

  //
  // clean up
  //
  free(image_output_buf);
  free(row_pointers);

  return(0);
}

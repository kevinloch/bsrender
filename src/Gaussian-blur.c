#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int GaussianBlur(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  double radius;
  int sample_size;
  int half_sample_size;
  int blur_res_x;
  int blur_res_y;
  int blur_i;
  int blur_x;
  int blur_y;
  int source_x;
  int source_y;
  int kernel_i;
  int current_image_res_x;
  int current_image_res_y;
  int current_image_offset;
  double *G_kernel_array;
  double *G_kernel_p;
  double G_kernel_sum;
  double G;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  pixel_composition_t *current_image_p;
  pixel_composition_t *image_blur_p;
  size_t blur_buffer_size;
  double G_r;
  double G_g;
  double G_b;

  //
  // determine Gaussian kernel size
  //
  radius=bsr_config->Gaussian_blur_radius;
  sample_size=((int)ceil(radius) * 6) + 1;
  half_sample_size=((int)ceil(radius) * 3) + 1;

  //
  // display status message if not in cgi mode
  //
  if (bsr_config->cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Applying Gaussian blur with radius %.3e...", radius);
    fflush(stdout);
  }

  //
  // allocate memory for 1D Gaussian kernel
  //
  G_kernel_array=(double *)malloc(sample_size * sizeof(double));
  if (G_kernel_array == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate memory for Gaussian blur kernel\n");
      fflush(stdout);
    }
    return(1);
  }

  //
  // generate 1D Gaussian kernel
  //
  G_kernel_sum=0.0;
  G_kernel_p=G_kernel_array;
  for (kernel_i=-half_sample_size+1; kernel_i < half_sample_size; kernel_i++) {
    G=exp(-(double)kernel_i * (double)kernel_i / (2.0 * radius * radius)) / sqrt(2.0 * M_PI * radius * radius);
    G_kernel_sum+=G;
    *G_kernel_p=G;
    G_kernel_p++;
  } // end for kernel

  //
  // normalize Gaussian kernel
  //
  G_kernel_p=G_kernel_array;
  for (kernel_i=-half_sample_size+1; kernel_i < half_sample_size; kernel_i++) {
    *G_kernel_p/=G_kernel_sum;
    G_kernel_p++;
  } // end for kernel

  //
  // get current image pointer and resolution
  //
  current_image_p=bsr_state->current_image_buf;
  current_image_res_x=bsr_state->current_image_res_x;
  current_image_res_y=bsr_state->current_image_res_y;
  blur_res_x=current_image_res_x;
  blur_res_y=current_image_res_y;

  //
  // allocate memory for image blur buffer
  //
  blur_buffer_size=blur_res_x * blur_res_y * sizeof(pixel_composition_t);
  bsr_state->image_blur_buf=(pixel_composition_t *)malloc(blur_buffer_size);
  if (bsr_state->image_blur_buf == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate memory for image blur buffer\n");
      fflush(stdout);
    }
    return(1);
  }

  //
  // apply Gaussian 1D kernel to each pixel horizontally and put output in blur buffer
  //
  blur_x=0;
  blur_y=0;
  image_blur_p=bsr_state->image_blur_buf;
  for (blur_i=0; blur_i < (blur_res_x * blur_res_y); blur_i++) {
    if (blur_x == blur_res_x) {
      blur_x=0;
      blur_y++;
    }

    //
    // apply Gaussian kernel to this pixel horizontally
    //
    G_r=0.0;
    G_g=0.0;
    G_b=0.0;
    G_kernel_p=G_kernel_array;
    for (kernel_i=-half_sample_size+1; kernel_i < half_sample_size; kernel_i++) {
      source_x=blur_x + kernel_i;
      current_image_offset=(blur_y * blur_res_x) + source_x;
      current_image_p=bsr_state->current_image_buf + current_image_offset;
      if ((source_x >= 0) && (source_x < current_image_res_x)) {
        G_r+=(current_image_p->r * *G_kernel_p);
        G_g+=(current_image_p->g * *G_kernel_p);
        G_b+=(current_image_p->b * *G_kernel_p);
        G_kernel_p++;
      } // end if within current image bounds
    } // end for kernel

    //
    // copy blurred pixel to blur buffer
    //
    image_blur_p->r=G_r;
    image_blur_p->g=G_g;
    image_blur_p->b=G_b;

    blur_x++;
    image_blur_p++;
  } // end for blur_i

  //
  // apply Gaussian 1D kernel to each pixel vertically and put output back in 'current_image_buffer'
  //
  blur_x=0;
  blur_y=0;
  // note swapping buffers to go back to 'current_image_buffer' so some variable names are backwards
  current_image_p=bsr_state->image_blur_buf;
  image_blur_p=bsr_state->current_image_buf;
  for (blur_i=0; blur_i < (blur_res_x * blur_res_y); blur_i++) {
    if (blur_x == blur_res_x) {
      blur_x=0;
      blur_y++;
    }

    //
    // apply Gaussian kernel to this pixel vertically
    //
    G_r=0.0;
    G_g=0.0;
    G_b=0.0;
    G_kernel_p=G_kernel_array;
    for (kernel_i=-half_sample_size+1; kernel_i < half_sample_size; kernel_i++) {
      source_y=blur_y + kernel_i;
      current_image_offset=(source_y * blur_res_x) + blur_x;
      current_image_p=bsr_state->image_blur_buf + current_image_offset;
      if ((source_y >= 0) && (source_y < current_image_res_y)) {
        G_r+=(current_image_p->r * *G_kernel_p);
        G_g+=(current_image_p->g * *G_kernel_p);
        G_b+=(current_image_p->b * *G_kernel_p);
        G_kernel_p++;
      } // end if within current image bounds
    } // end for kernel

    //
    // copy blurred pixel to blur buffer
    //
    image_blur_p->r=G_r;
    image_blur_p->g=G_g;
    image_blur_p->b=G_b;

    blur_x++;
    image_blur_p++;
  } // end for blur_i

  //
  // free image blur buffer
  //
  free(bsr_state->image_blur_buf);

  //
  // display execution time if not in cgi mode
  //
  if (bsr_config->cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.4fs)\n", elapsed_time);
    fflush(stdout);
  }

  return(0);
}

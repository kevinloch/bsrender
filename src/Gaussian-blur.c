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
  int matrix_x;
  int matrix_y;
  int current_image_res_x;
  int current_image_res_y;
  int current_image_offset;
  double *G_matrix_array;
  double *G_matrix_p;
  double G_matrix_sum;
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
  // determine Gaussian matrix size
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
  // allocate memory for Gaussian matrix
  //
  G_matrix_array=(double *)malloc(sample_size * sample_size * sizeof(double));
  if (G_matrix_array == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate memory for image blur buffer\n");
      fflush(stdout);
    }
    return(1);
  }

  //
  // generate Gaussian matrix
  //
  G_matrix_sum=0.0;
  G_matrix_p=G_matrix_array;
  for (matrix_y=-half_sample_size+1; matrix_y < half_sample_size; matrix_y++) {
    for (matrix_x=-half_sample_size+1; matrix_x < half_sample_size; matrix_x++) {
      G=exp(-(((double)matrix_x * (double)matrix_x) + ((double)matrix_y * (double)matrix_y)) / (2.0 * radius * radius) / (2.0 * M_PI * radius * radius)  );
      G_matrix_sum+=G;
      *G_matrix_p=G;
      G_matrix_p++;
    } // end for matrix_x
  } // end for matrix_y

  //
  // normalize Gaussian matrix
  //
  G_matrix_p=G_matrix_array;
  for (matrix_y=-half_sample_size+1; matrix_y < half_sample_size; matrix_y++) {
    for (matrix_x=-half_sample_size+1; matrix_x < half_sample_size; matrix_x++) {
      *G_matrix_p/=G_matrix_sum;
      G_matrix_p++;
    } // end for matrix_x
  } // end for matrix_y

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
  // apply Gaussian matrix to each pixel and put output in blur buffer
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
    // apply Gaussian matrix to this pixel
    //
    G_r=0.0;
    G_g=0.0;
    G_b=0.0;
    G_matrix_p=G_matrix_array;
    for (matrix_y=-half_sample_size+1; matrix_y < half_sample_size; matrix_y++) {
      source_y=blur_y + matrix_y;
      for (matrix_x=-half_sample_size+1; matrix_x < half_sample_size; matrix_x++) {
        source_x=blur_x + matrix_x;
        current_image_offset=(source_y * blur_res_x) + source_x;
        current_image_p=bsr_state->current_image_buf + current_image_offset;
        if ((source_x >= 0) && (source_x < current_image_res_x) && (source_y >= 0) && (source_y < current_image_res_y)) {
          G_r+=(current_image_p->r * *G_matrix_p);
          G_g+=(current_image_p->g * *G_matrix_p);
          G_b+=(current_image_p->b * *G_matrix_p);
          G_matrix_p++;
        } // end if within current image bounds
      } // end for matrix_x
    } // end for matrix_y

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
  // free old image buffer and update current_image_buf pointer
  //
  free(bsr_state->current_image_buf);
  bsr_state->current_image_buf=bsr_state->image_blur_buf;
  bsr_state->current_image_res_x=blur_res_x;
  bsr_state->current_image_res_y=blur_res_y;

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

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int resizeLanczos(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  pixel_composition_t *current_image_p;
  pixel_composition_t *image_resize_p;
  int current_image_offset;
  int resize_res_x;
  int resize_res_y;
  double source_w;
  double half_source_w;
  double source_x_center;
  double source_y_center;
  int source_x;
  int source_y;
  int resize_x;
  int resize_y;
  double L_kernel;
  double L_x_r;
  double L_x_g;
  double L_x_b;
  double L_y_r;
  double L_y_g;
  double L_y_b;
  double L_distance_x;
  double L_distance_y;
  int resize_i;
  int i_x;
  int i_y;
  int Lanczos_order;
  int current_image_res_x;
  int current_image_res_y;
  int lines_per_thread;
  volatile int cont;
  int all_workers_done;
  int i;

  //
  // all threads: get current image resolution and calculate resize resolution
  //
  current_image_res_x=bsr_state->current_image_res_x;
  current_image_res_y=bsr_state->current_image_res_y;
  resize_res_x=bsr_state->resize_res_x;
  resize_res_y=bsr_state->resize_res_y;
  source_w=1.0 / bsr_config->output_scaling_factor;
  half_source_w=source_w / 2.0;
  lines_per_thread=(int)ceil(((double)resize_res_y / (double)(bsr_state->num_worker_threads + 1)));

  //
  // main thread: display status message if not in cgi mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config->cgi_mode != 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Resizing image from %dx%d to %dx%d...", current_image_res_x, current_image_res_y, resize_res_x, resize_res_y);
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    cont=0;
    while (cont == 0) {
      if (bsr_state->status_array[bsr_state->perthread->my_thread_id] >= 20) {
        cont=1;
      }
    }
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i]=20;
    }
  } // end if not main thread

  //
  // all threads: copy rendered image to resize buffer using Lanczos interpolation
  //
  Lanczos_order=2;
  resize_x=0;
  resize_y=(bsr_state->perthread->my_thread_id * lines_per_thread);
  image_resize_p=bsr_state->image_resize_buf + (bsr_state->perthread->my_thread_id * lines_per_thread * resize_res_x);
  for (resize_i=0; ((resize_i < (resize_res_x * lines_per_thread)) && (resize_y < (resize_res_y - 1))); resize_i++) {
    if (resize_x == resize_res_x) {
      resize_x=0;
      resize_y++;
    }
    source_x_center=((double)resize_x * source_w) + half_source_w + 0.5;
    source_y_center=((double)resize_y * source_w) + half_source_w + 0.5;
    L_y_r=0.0;
    L_y_g=0.0;
    L_y_b=0.0;
    for (i_y=((int)source_y_center - Lanczos_order + 1); i_y <= ((int)source_y_center + Lanczos_order); i_y++) {
      source_y=i_y;

      //
      // L_x
      //
      L_x_r=0.0;
      L_x_g=0.0;
      L_x_b=0.0;
      for (i_x=((int)source_x_center - Lanczos_order + 1); i_x <= ((int)source_x_center + Lanczos_order); i_x++) {
        source_x=i_x;
        current_image_offset=(source_y * current_image_res_x) + source_x;
        current_image_p=bsr_state->current_image_buf + current_image_offset;
        if ((source_x >= 0) && (source_x < current_image_res_x) && (source_y >= 0) && (source_y < current_image_res_y)) {
          L_distance_x=source_x_center - (double)i_x;
          if (L_distance_x == 0.0) {
            L_x_r+=current_image_p->r;
            L_x_g+=current_image_p->g;
            L_x_b+=current_image_p->b;
          } else if ((L_distance_x >= -Lanczos_order) && (L_distance_x < Lanczos_order)) {
            L_kernel=Lanczos_order * sin(M_PI * L_distance_x) * sin(M_PI * L_distance_x / Lanczos_order) / (M_PI * M_PI * L_distance_x * L_distance_x);
            L_x_r+=(current_image_p->r * L_kernel);
            L_x_g+=(current_image_p->g * L_kernel);
            L_x_b+=(current_image_p->b * L_kernel);
          } else {
            // L_x == 0
          } // end if L_distance_x
        } // end if within composition buffer bounds
      } // end for i_x

      //
      // L_y
      //
      L_distance_y=source_y_center - (double)i_y;
      if (L_distance_y == 0.0) {
        L_y_r+=L_x_r;
        L_y_g+=L_x_g;
        L_y_b+=L_x_b;
      } else if ((L_distance_y >= -Lanczos_order) && (L_distance_y < Lanczos_order)) {
        L_kernel=Lanczos_order * sin(M_PI * L_distance_y) * sin(M_PI * L_distance_y / Lanczos_order) / (M_PI * M_PI * L_distance_y * L_distance_y);
        L_y_r+=(L_x_r * L_kernel);
        L_y_g+=(L_x_g * L_kernel);
        L_y_b+=(L_x_b * L_kernel);
      } else {
        // L_y = 0
      } // end if L_distance_y
    } // end for i_y

    //
    // handle clipping
    //
    if (L_y_r < 0.0) {
      L_y_r=0.0;
    }
    if (L_y_g < 0.0) {
      L_y_g=0.0;
    }
    if (L_y_b < 0.0) {
      L_y_b=0.0;
    }

    //
    // copy to output buffer
    //
    image_resize_p->r=L_y_r;
    image_resize_p->g=L_y_g;
    image_resize_p->b=L_y_b;

    resize_x++;
    image_resize_p++;
  }

  //
  // worker threads: signal this thread is done and wait until main thread says we can continue to next step.
  // main thread: wait until all other threads are done and then signal that they can continue to next step.
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    bsr_state->status_array[bsr_state->perthread->my_thread_id]=21;
    cont=0;
    while (cont == 0) {
      if (bsr_state->status_array[bsr_state->perthread->my_thread_id] >= 22) {
        cont=1;
      }
    }
  } else {
    all_workers_done=0;
    while (all_workers_done == 0) {
      all_workers_done=1;
      for (i=1; i <= bsr_state->num_worker_threads; i++) {
        if (bsr_state->status_array[i] < 21) {
          all_workers_done=0;
        }
      }
    }
    //
    // ready to continue, set all worker thread status to 22
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i]=22;
    }
  } // end if not main thread


  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    //
    // main thread: update current_image_buf pointer
    //
    bsr_state->current_image_buf=bsr_state->image_resize_buf;
    bsr_state->current_image_res_x=resize_res_x;
    bsr_state->current_image_res_y=resize_res_y;

    //
    // main thread: output execution time if not in cgi mode
    //
    if (bsr_config->cgi_mode != 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.4fs)\n", elapsed_time);
      fflush(stdout);
    }
  } // end if main thread

  return(0);
}

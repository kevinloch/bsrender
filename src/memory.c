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
#include <sys/mman.h>
#include <time.h>

int freeMemory(bsr_state_t *bsr_state) {
  if (bsr_state->image_composition_buf != NULL) {
    munmap(bsr_state->image_composition_buf, bsr_state->composition_buffer_size);
  }
  if (bsr_state->image_blur_buf != NULL) {
    munmap(bsr_state->image_blur_buf, bsr_state->blur_buffer_size);
  }
  if (bsr_state->image_resize_buf != NULL) {
    munmap(bsr_state->image_resize_buf, bsr_state->resize_buffer_size);
  }
  if (bsr_state->thread_buf != NULL) {
    munmap(bsr_state->thread_buf, bsr_state->thread_buffer_size);
  }
  if (bsr_state->status_array != NULL) {
    munmap(bsr_state->status_array, bsr_state->status_array_size);
  }
  if (bsr_state->Airymap_red != NULL) {
    munmap(bsr_state->Airymap_red, bsr_state->Airymap_size);
  }
  if (bsr_state->Airymap_green != NULL) {
    munmap(bsr_state->Airymap_green, bsr_state->Airymap_size);
  }
  if (bsr_state->Airymap_blue != NULL) {
    munmap(bsr_state->Airymap_blue, bsr_state->Airymap_size);
  }
  if (bsr_state->dedup_buf != NULL) {
    free(bsr_state->dedup_buf);
  }
  if (bsr_state->dedup_index != NULL) {
    free(bsr_state->dedup_index);
  }
  if (bsr_state->image_output_buf != NULL) {
    munmap(bsr_state->image_output_buf, bsr_state->output_buffer_size);
  }
  if (bsr_state->row_pointers != NULL) {
    munmap(bsr_state->row_pointers, bsr_state->row_pointers_size);
  }
  if (bsr_state->compressed_sizes != NULL) {
    munmap(bsr_state->compressed_sizes, bsr_state->compressed_sizes_size);
  }
  if (bsr_state->compression_buf1 != NULL) {
    free(bsr_state->compression_buf1);
  }
  if (bsr_state->compression_buf2 != NULL) {
    free(bsr_state->compression_buf2);
  }
  // must be freed last
  if (bsr_state != NULL) {
    munmap(bsr_state, bsr_state->bsr_state_size);
  }
  return(0);
}

int allocateMemory(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  int mmap_protection;
  int mmap_visibility;
  int Airymap_width;
  int dedup_index_count;
  dedup_buffer_t *dedup_buf_p;
  dedup_index_t *dedup_index_p;
  int i;
  thread_buffer_t *main_thread_buf_p;
  int output_res_x;
  int output_res_y;
  int lines_per_block=0;
  int pixel_data_size=0;

  //
  // allocate shared memory for Airy disk maps if Airy disk mode enabled
  //
  if (bsr_config->Airy_disk_enable == 1) {
    mmap_protection=PROT_READ | PROT_WRITE;
    mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
    Airymap_width=bsr_config->Airy_disk_max_extent + 1;
    bsr_state->Airymap_size=(size_t)Airymap_width * (size_t)Airymap_width * sizeof(double);
    bsr_state->Airymap_red=(double *)mmap(NULL, bsr_state->Airymap_size, mmap_protection, mmap_visibility, -1, 0);
    bsr_state->Airymap_green=(double *)mmap(NULL, bsr_state->Airymap_size, mmap_protection, mmap_visibility, -1, 0);
    bsr_state->Airymap_blue=(double *)mmap(NULL, bsr_state->Airymap_size, mmap_protection, mmap_visibility, -1, 0);
  }

  //
  // allocate shared memory for image composition buffer (floating-point rgb)
  //
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
  bsr_state->composition_buffer_size=(size_t)bsr_config->camera_res_x * (size_t)bsr_config->camera_res_y * sizeof(pixel_composition_t);
  bsr_state->image_composition_buf=(pixel_composition_t *)mmap(NULL, bsr_state->composition_buffer_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->image_composition_buf == MAP_FAILED) {
    if ((bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
      printf("Error: could not allocate shared memory for image composition buffer\n");
    }
    exit(1);
  }
  bsr_state->current_image_buf=bsr_state->image_composition_buf;
  bsr_state->current_image_res_x=bsr_config->camera_res_x;
  bsr_state->current_image_res_y=bsr_config->camera_res_y;

  //
  // allocate shared memory for image blur buffer if needed
  //
  if (bsr_config->Gaussian_blur_radius > 0.0) {
    mmap_protection=PROT_READ | PROT_WRITE;
    mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
    bsr_state->blur_buffer_size=(size_t)bsr_config->camera_res_x * (size_t)bsr_config->camera_res_y * sizeof(pixel_composition_t);
    bsr_state->image_blur_buf=(pixel_composition_t *)mmap(NULL, bsr_state->blur_buffer_size, mmap_protection, mmap_visibility, -1, 0);
    if (bsr_state->image_blur_buf == MAP_FAILED) {
      if (bsr_config->cgi_mode != 1) {
        printf("Error: could not allocate shared memory for image blur buffer\n");
        fflush(stdout);
      }
      exit(1);
    }
  }

  //
  // allocate shared memory for image resize buffer if needed
  //
  if (bsr_config->output_scaling_factor != 1.0) {
    bsr_state->resize_res_x=(int)(((double)bsr_config->camera_res_x * bsr_config->output_scaling_factor) + 0.5);
    bsr_state->resize_res_y=(int)(((double)bsr_config->camera_res_y * bsr_config->output_scaling_factor) + 0.5);
    mmap_protection=PROT_READ | PROT_WRITE;
    mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
    bsr_state->resize_buffer_size=(size_t)bsr_state->resize_res_x * (size_t)bsr_state->resize_res_y * sizeof(pixel_composition_t);
    bsr_state->image_resize_buf=(pixel_composition_t *)mmap(NULL, bsr_state->resize_buffer_size, mmap_protection, mmap_visibility, -1, 0);
    if (bsr_state->image_resize_buf == MAP_FAILED) {
      if (bsr_config->cgi_mode != 1) {
        printf("Error: could not allocate shared memory for image resize buffer\n");
        fflush(stdout);
      }
      exit(1);
    }
  }

  //
  // allocate non-shared memory for dedup buffer and initialize
  //
  if ((bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing dedup buffers and indexes...");
    fflush(stdout);
  }
  bsr_state->dedup_buffer_size=(size_t)bsr_state->per_thread_buffers * sizeof(dedup_buffer_t);
  bsr_state->dedup_buf=(dedup_buffer_t *)malloc(bsr_state->dedup_buffer_size);
  if (bsr_state->dedup_buf == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate memory for dedup buffer\n");
    }
    exit(1);
  }
  // initialize dedup buffer
  dedup_buf_p=bsr_state->dedup_buf;
  for (i=0; i < (bsr_state->per_thread_buffers); i++) {
    dedup_buf_p->image_offset=-1;
    dedup_buf_p->r=0.0;
    dedup_buf_p->g=0.0;
    dedup_buf_p->b=0.0;
    dedup_buf_p++;
  }
  bsr_state->perthread->dedup_count=0;

  //
  // allocate non-shared memory for dedup index and initialize
  //
  if ((uint64_t)bsr_config->camera_res_x * (uint64_t)bsr_config->camera_res_y <= 16777216) {
    bsr_state->dedup_index_mode=0; // use image_offset for dedup index
    dedup_index_count=bsr_config->camera_res_x * bsr_config->camera_res_y;
    bsr_state->dedup_index_size=(size_t)dedup_index_count * sizeof(dedup_index_t);
  } else {
    bsr_state->dedup_index_mode=1; // use lowest 24 bits of image_offset for dedup index
    dedup_index_count=0xffffff;
    bsr_state->dedup_index_size=(size_t)dedup_index_count * sizeof(dedup_index_t);
  }
  bsr_state->dedup_index=(dedup_index_t *)malloc(bsr_state->dedup_index_size);
  if (bsr_state->dedup_index == NULL) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate shared memory for dedup index\n");
    }
    exit(1);
  }
  // initialize dedup index
  dedup_index_p=bsr_state->dedup_index;
  for (i=0; i < dedup_index_count; i++) {
    dedup_index_p->dedup_record_p=NULL;
    dedup_index_p++;
  }
  if ((bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  //
  // allocate shared memory for main thread buffer and status array
  //
  if ((bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing main thread buffer...");
    fflush(stdout);
  }
  if (bsr_state->per_thread_buffers < 1) {
    bsr_state->per_thread_buffers=1;
  }
  bsr_state->thread_buffer_count=bsr_state->num_worker_threads * bsr_state->per_thread_buffers;
  bsr_state->thread_buffer_size=(size_t)bsr_state->thread_buffer_count * sizeof(thread_buffer_t);
  bsr_state->thread_buf=(thread_buffer_t *)mmap(NULL, bsr_state->thread_buffer_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->thread_buf == MAP_FAILED) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate shared memory for main thread buffer\n");
    }
    exit(1);
  }
  // initialize thread buffer
  main_thread_buf_p=bsr_state->thread_buf;
  for (i=0; i < (bsr_state->num_worker_threads * bsr_state->per_thread_buffers); i++) {
    main_thread_buf_p->status_left=0;
    main_thread_buf_p->status_right=0;
    main_thread_buf_p++;
  }
  // allocate shared memory for thread status array
  bsr_state->status_array_size=(size_t)(bsr_state->num_worker_threads + 1) * sizeof(bsr_status_t);
  bsr_state->status_array=(bsr_status_t *)mmap(NULL, bsr_state->status_array_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->status_array == MAP_FAILED) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate shared memory for thread status array\n");
    }
    exit(1);
  }
  if ((bsr_config->cgi_mode != 1) && (bsr_config->print_status == 1)) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  //
  // allocate shared memory for image output buffer
  //
  if (bsr_config->output_scaling_factor != 1.0) {
    output_res_x=bsr_state->resize_res_x;
    output_res_y=bsr_state->resize_res_y;
  } else {
    output_res_x=bsr_config->camera_res_x;
    output_res_y=bsr_config->camera_res_y;
  }
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
  if (bsr_config->bits_per_color == 32) {
    bsr_state->output_buffer_size=(size_t)output_res_x * (size_t)output_res_y * (size_t)12 * sizeof(unsigned char);
  } else if ((bsr_config->bits_per_color == 10) || (bsr_config->bits_per_color == 12) || (bsr_config->bits_per_color == 16)) {
    bsr_state->output_buffer_size=(size_t)output_res_x * (size_t)output_res_y * (size_t)6 * sizeof(unsigned char);
  } else { // default 8 bits per color
    bsr_state->output_buffer_size=(size_t)output_res_x * (size_t)output_res_y * (size_t)3 * sizeof(unsigned char);
  }
  bsr_state->image_output_buf=(unsigned char *)mmap(NULL, bsr_state->output_buffer_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->image_output_buf == MAP_FAILED) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate shared memory for image output buffer\n");
      fflush(stdout);
    }
    exit(1);
  }

  //
  // allocate shared memory for row_pointers used in various ways by different image output encoders
  //
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
  bsr_state->row_pointers_size=(size_t)output_res_y * sizeof(unsigned char *);
  bsr_state->row_pointers=(unsigned char **)mmap(NULL, bsr_state->row_pointers_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->row_pointers == MAP_FAILED) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not allocate shared memory for row_pointers\n");
      fflush(stdout);
    }
    exit(1);
  }

  //
  // allocate memory for image compression buffers if required
  //
  if ((bsr_config->image_format == 1) && ((bsr_config->exr_compression == 2) || (bsr_config->exr_compression == 3))) {
    // allocate shared memory for compressed_sizes table
    mmap_protection=PROT_READ | PROT_WRITE;
    mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
    bsr_state->compressed_sizes_size=(size_t)output_res_y * sizeof(int);
    bsr_state->compressed_sizes=(int *)mmap(NULL, bsr_state->compressed_sizes_size, mmap_protection, mmap_visibility, -1, 0);
    if (bsr_state->row_pointers == MAP_FAILED) {
      if (bsr_config->cgi_mode != 1) {
        printf("Error: could not allocate shared memory for compressed_sizes array\n");
        fflush(stdout);
      }
      exit(1);
    }

    // allocate non-shared memory for compression_buf1 and 2
    if (bsr_config->exr_compression == 2) {
      // deflate, 1 line per block
      lines_per_block=1;
    } else if (bsr_config->exr_compression == 3) {
      // deflate, 16 line per block
      lines_per_block=16;
    }
    if (bsr_config->bits_per_color == 16) {
      pixel_data_size=6 * output_res_x * lines_per_block;
    } else if (bsr_config->bits_per_color == 32) {
      pixel_data_size=12 * output_res_x * lines_per_block;
    }
    // allocate non-shared memory for compression_buf1
    bsr_state->compression_buf_size=(size_t)pixel_data_size * sizeof(unsigned char);
    bsr_state->compression_buf1=(unsigned char *)malloc(bsr_state->compression_buf_size);
    if (bsr_state->compression_buf1 == NULL) {
      if (bsr_config->cgi_mode != 1) {
        printf("Error: could not allocate memory for compression buffer 1\n");
      }
      exit(1);
    }
    // allocate non-shared memory for compression_buf2
    bsr_state->compression_buf2=(unsigned char *)malloc(bsr_state->compression_buf_size);
    if (bsr_state->compression_buf2== NULL) {
      if (bsr_config->cgi_mode != 1) {
        printf("Error: could not allocate memory for compression buffer 2\n");
      }
      exit(1);
    }
  } // end if exr_compression

  return(0);
}

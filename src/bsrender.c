//
// Billion Star 3D Rendering Engine
// Kevin M. Loch
//
// 3D rendering engine for the ESA Gaia EDR3 star dataset
//

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <math.h>
#include <time.h>
#include "bsr-config.h"
#include "usage.h"
#include "util.h"
#include "rgb.h"
#include "bsr-png.h"
#include "process-stars.h"
#include "init-state.h"
#include "cgi.h"
#include "diffraction.h"
#include "post-process.h"

int processCmdArgs(bsr_config_t *bsr_config, int argc, char **argv) {
  int i;
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
          strcpy(bsr_config->config_file_name, (option_start + (size_t)2)); 
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          strcpy(bsr_config->config_file_name, option_start);
        } // end if no space

      } else if (argv[i][1] == 'd') {
        // data files directory
        if (argv[i][2] != 0) {
          // option concatenated onto switch
          option_start=argv[i];
          strcpy(bsr_config->data_file_directory, (option_start + (size_t)2));
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          strcpy(bsr_config->data_file_directory, option_start);
        } // end if no space

      } else if (argv[i][1] == 'h') {
        // print help
        printUsage();
        exit(0);
      } // end which option
    } // end for argc
  } // end if any options
  return(0);
}

int main(int argc, char **argv) {
  bsr_config_t bsr_config;
  bsr_state_t *bsr_state;
  bsr_thread_state_t perthread;
  char file_path[1024];
  FILE *input_file_external=NULL;
  FILE *input_file_pq100=NULL;
  FILE *input_file_pq050=NULL;
  FILE *input_file_pq030=NULL;
  FILE *input_file_pq020=NULL;
  FILE *input_file_pq010=NULL;
  FILE *input_file_pq005=NULL;
  FILE *input_file_pq003=NULL;
  FILE *input_file_pq002=NULL;
  FILE *input_file_pq001=NULL;
  FILE *input_file_pq000=NULL;
  int mmap_protection;
  int mmap_visibility;
  struct timespec overall_starttime;
  struct timespec overall_endtime;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  size_t bsr_state_size=0;
  size_t composition_buffer_size=0;
  size_t blur_buffer_size=0;
  size_t resize_buffer_size=0;
  size_t status_array_size=0;
  size_t dedup_buffer_size=0;
  size_t dedup_index_size=0;
  size_t Airymap_size=0;
  size_t thread_buffer_size=0;
  int thread_buffer_count;
  int all_workers_done;
  int i;
  double rgb_red[32768];
  double rgb_green[32768];
  double rgb_blue[32768];
  pixel_composition_t *image_composition_p;
  dedup_buffer_t *dedup_buf_p;
  dedup_index_t *dedup_index_p;
  thread_buffer_t *main_thread_buf_p;
  int main_thread_buffer_index;
  int buffer_is_empty;
  int empty_passes;
  int Airymap_xy;

  //
  // initialize bsr_config to default values
  //
  initConfig(&bsr_config);

  //
  // process command line arguments
  //
  processCmdArgs(&bsr_config, argc, argv);

  //
  // load config file
  //
  loadConfigFromFile(&bsr_config);

  //
  // if cgi mode, print CGI output header and read query string, otherwise print version
  //
  if (bsr_config.cgi_mode == 1) {
    printCGIheader();
    getCGIOptions(&bsr_config);
    validateCGIOptions(&bsr_config);
  } else {
    printf("bsrender version %s\n", BSR_VERSION);
  }

  //
  // allocate shared memory for bsr_state
  //
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
  bsr_state_size=sizeof(bsr_state_t);
  bsr_state=(bsr_state_t *)mmap(NULL, bsr_state_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state == NULL) {
    if (bsr_config.cgi_mode != 1) {
      printf("Error: could not allocate shared memory for bsr_state\n");
    }
    return(1);
  }
  bsr_state->perthread=&perthread;

  //
  // initialize bsr_state (store processed settings and setup camera target)
  //
  initState(&bsr_config, bsr_state);

  //
  // calculate number of rendering threads to be forked
  //
  bsr_state->num_worker_threads=bsr_config.num_threads-1;
  if (bsr_state->num_worker_threads < 1) {
    bsr_state->num_worker_threads=1;
  }

  //
  // initialize total execution timer and display major performance affecting options
  //
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &overall_starttime);
    if (bsr_config.enable_Gaia == 1) {
      printf("Minimum Gaia parallax quality: %d, total threads: %d, buffers per rendering thread: %d pixels\n", bsr_config.Gaia_min_parallax_quality, (bsr_state->num_worker_threads + 1), bsr_state->per_thread_buffers);
    } else {
      printf("Total threads: %d, buffers per rendering thread: %d pixels\n", (bsr_state->num_worker_threads + 1), bsr_state->per_thread_buffers);
    }
    fflush(stdout);
  }

  //
  // initialize RGB color lookup tables
  //
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing rgb color tables...");
    fflush(stdout);
  }
  bsr_state->rgb_red=rgb_red;
  bsr_state->rgb_green=rgb_green;
  bsr_state->rgb_blue=rgb_blue;
  initRGBTables(&bsr_config, bsr_state);
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  //
  // allocate memory for Airy disk maps and initialize
  //
  if (bsr_config.Airy_disk == 1) {
    if (bsr_config.cgi_mode != 1) {
      clock_gettime(CLOCK_REALTIME, &starttime);
      printf("Initializing Airy disk maps...");
      fflush(stdout);
    }
    Airymap_xy=(bsr_config.Airy_disk_max_extent * 2) + 1;
    Airymap_size=Airymap_xy * Airymap_xy * sizeof(double);
    bsr_state->Airymap_red=(double *)malloc(Airymap_size);
    bsr_state->Airymap_green=(double *)malloc(Airymap_size);
    bsr_state->Airymap_blue=(double *)malloc(Airymap_size);
    initAiryMaps(&bsr_config, bsr_state);
    if (bsr_config.cgi_mode != 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.3fs)\n", elapsed_time);
      fflush(stdout);
    }
  }

  //
  // allocate shared memory for image composition buffer (floating point rgb)
  //
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing image composition buffer...");
    fflush(stdout);
  }
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
  composition_buffer_size=bsr_config.camera_res_x * bsr_config.camera_res_y * sizeof(pixel_composition_t);
  bsr_state->image_composition_buf=(pixel_composition_t *)mmap(NULL, composition_buffer_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->image_composition_buf == NULL) {
    if (bsr_config.cgi_mode != 1) {
      printf("Error: could not allocate shared memory for image composition buffer\n");
    }
    return(1);
  }
  // initialize (clear) image composition buffer
  image_composition_p=bsr_state->image_composition_buf;
  for (i=0; i < (bsr_config.camera_res_x * bsr_config.camera_res_y); i++) {
    image_composition_p->r=0.0;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
    image_composition_p++;
  }
  bsr_state->current_image_buf=bsr_state->image_composition_buf;
  bsr_state->current_image_res_x=bsr_config.camera_res_x;
  bsr_state->current_image_res_y=bsr_config.camera_res_y;
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  //
  // allocate shared memory for image blur buffer if needed
  //
  if (bsr_config.Gaussian_blur_radius > 0.0) {
    mmap_protection=PROT_READ | PROT_WRITE;
    mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
    blur_buffer_size=bsr_config.camera_res_x * bsr_config.camera_res_y * sizeof(pixel_composition_t);
    bsr_state->image_blur_buf=(pixel_composition_t *)mmap(NULL, blur_buffer_size, mmap_protection, mmap_visibility, -1, 0);
    if (bsr_state->image_blur_buf == NULL) {
      if (bsr_config.cgi_mode != 1) {
        printf("Error: could not allocate shared memory for image blur buffer\n");
        fflush(stdout);
      }
      return(1);
    }
  }

  //
  // allocate shared memory for image resize buffer if needed
  //
  if (bsr_config.output_scaling_factor != 1.0) {
    bsr_state->resize_res_x=(int)(((double)bsr_config.camera_res_x * bsr_config.output_scaling_factor) + 0.5);
    bsr_state->resize_res_y=(int)(((double)bsr_config.camera_res_y * bsr_config.output_scaling_factor) + 0.5);
    mmap_protection=PROT_READ | PROT_WRITE;
    mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
    resize_buffer_size=bsr_state->resize_res_x * bsr_state->resize_res_y * sizeof(pixel_composition_t);
    bsr_state->image_resize_buf=(pixel_composition_t *)mmap(NULL, resize_buffer_size, mmap_protection, mmap_visibility, -1, 0);
    if (bsr_state->image_resize_buf == NULL) {
      if (bsr_config.cgi_mode != 1) {
        printf("Error: could not allocate shared memory for image resize buffer\n");
        fflush(stdout);
      }
      return(1);
    }
  }

  //
  // allocate non-shared memory for dedup buffer and initialize
  //
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing dedup buffer and index...");
    fflush(stdout);
  }
  dedup_buffer_size=bsr_state->per_thread_buffers * sizeof(dedup_buffer_t);
  bsr_state->dedup_buf=(dedup_buffer_t *)malloc(dedup_buffer_size);
  if (bsr_state->dedup_buf == NULL) {
    if (bsr_config.cgi_mode != 1) {
      printf("Error: could not allocate shared memory for dedup buffer\n");
    }
    return(1);
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

  //
  // allocate non-shared memory for dedup index and initialize
  //
  dedup_index_size=bsr_config.camera_res_x * bsr_config.camera_res_y * sizeof(dedup_index_t);
  bsr_state->dedup_index=(dedup_index_t *)malloc(dedup_index_size);
  if (bsr_state->dedup_index == NULL) {
    if (bsr_config.cgi_mode != 1) {
      printf("Error: could not allocate shared memory for dedup index\n");
    }
    return(1);
  }
  // initialize dedup index
  dedup_index_p=bsr_state->dedup_index;
  for (i=0; i < (bsr_config.camera_res_x * bsr_config.camera_res_y); i++) {
    dedup_index_p->dedup_record=NULL;
    dedup_index_p++;
  }
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  //
  // allocate shared memory for thread buffer and status array
  //
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Initializing rendering thread buffers...");
    fflush(stdout);
  }
  if (bsr_state->per_thread_buffers < 1) {
    bsr_state->per_thread_buffers=1;
  }
  thread_buffer_count=bsr_state->num_worker_threads * bsr_state->per_thread_buffers;
  thread_buffer_size=thread_buffer_count * sizeof(thread_buffer_t);
  bsr_state->thread_buf=(thread_buffer_t *)mmap(NULL, thread_buffer_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->thread_buf == NULL) {
    if (bsr_config.cgi_mode != 1) {
      printf("Error: could not allocate shared memory for thread buffer\n");
    }
    return(1);
  }
  // initialize thread buffer
  main_thread_buf_p=bsr_state->thread_buf;
  for (i=0; i < (bsr_state->num_worker_threads * bsr_state->per_thread_buffers); i++) {
    main_thread_buf_p->status_left=0;
    main_thread_buf_p->status_right=0;
    main_thread_buf_p++;
  }
  // allocate shared memory for thread status array
  status_array_size=(bsr_state->num_worker_threads + 1) * sizeof(int);
  bsr_state->status_array=(int *)mmap(NULL, status_array_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state->status_array == NULL) {
    if (bsr_config.cgi_mode != 1) {
      printf("Error: could not allocate shared memory for thread status array\n");
    }
    return(1);
  }
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);
  }

  //
  // attempt to open 'external' input file
  //
  if (bsr_config.enable_external == 1) {
    sprintf(file_path, "%s/galaxy-external.dat", bsr_config.data_file_directory);
    input_file_external=fopen(file_path, "rb");
    if (input_file_external == NULL) {
      if (bsr_config.cgi_mode != 1) {
        printf("Error: could not open %s\n", file_path);
        fflush(stdout);
      }
      return(1);
    }
  } // end if enable external

  //
  // attempt to open Gaia input file(s)
  //
  if (bsr_config.enable_Gaia == 1) {
    sprintf(file_path, "%s/galaxy-pq100.dat", bsr_config.data_file_directory);
    input_file_pq100=fopen(file_path, "rb");
    if (input_file_pq100 == NULL) {
      if (bsr_config.cgi_mode != 1) {
        printf("Error: could not open %s\n", file_path);
        fflush(stdout);
      }
      return(1);
    }
    if (bsr_config.Gaia_min_parallax_quality < 100) {
      sprintf(file_path, "%s/galaxy-pq050.dat", bsr_config.data_file_directory);
      input_file_pq050=fopen(file_path, "rb");
      if (input_file_pq050 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 50) {
      sprintf(file_path, "%s/galaxy-pq030.dat", bsr_config.data_file_directory);
      input_file_pq030=fopen(file_path, "rb");
      if (input_file_pq030 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 30) {
      sprintf(file_path, "%s/galaxy-pq020.dat", bsr_config.data_file_directory);
      input_file_pq020=fopen(file_path, "rb");
      if (input_file_pq020 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 20) {
      sprintf(file_path, "%s/galaxy-pq010.dat", bsr_config.data_file_directory);
      input_file_pq010=fopen(file_path, "rb");
      if (input_file_pq010 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 10) {
      sprintf(file_path, "%s/galaxy-pq005.dat", bsr_config.data_file_directory);
      input_file_pq005=fopen(file_path, "rb");
      if (input_file_pq005 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 5) {
      sprintf(file_path, "%s/galaxy-pq003.dat", bsr_config.data_file_directory);
      input_file_pq003=fopen(file_path, "rb");
      if (input_file_pq003 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 3) {
      sprintf(file_path, "%s/galaxy-pq002.dat", bsr_config.data_file_directory);
      input_file_pq002=fopen(file_path, "rb");
      if (input_file_pq002 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 2) {
      sprintf(file_path, "%s/galaxy-pq001.dat", bsr_config.data_file_directory);
      input_file_pq001=fopen(file_path, "rb");
      if (input_file_pq001 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
    if (bsr_config.Gaia_min_parallax_quality < 1) {
      sprintf(file_path, "%s/galaxy-pq000.dat", bsr_config.data_file_directory);
      input_file_pq000=fopen(file_path, "rb");
      if (input_file_pq000 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        return(1);
      }
    }
  } // end if enable Gaia

  //
  // display begin rendering status if not in cgi mode
  //
  if (bsr_config.cgi_mode != 1) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Rendering stars to image composition buffer...");
    fflush(stdout);
  }

  //
  // fork rendering threads
  //
  bsr_state->master_pid=getpid();
  if (bsr_state->num_worker_threads > 0) {
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->perthread->my_pid=getpid();
      if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
        bsr_state->perthread->my_thread_id=i; // this gets inherited by forked process
        bsr_state->status_array[i]=0;
        fork();
      }
    }
  }
  bsr_state->perthread->my_pid=getpid();
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    bsr_state->perthread->my_thread_id=0;
  }

  // 
  // begin thread specific processing.  Rendering threads read data files and send pixles to main thread.  Main integrates these pixels into
  // the image composition buffer
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {

    //
    // rendering thread: set our thread buffer postion to the beginning of this threads block
    //
    bsr_state->perthread->thread_buf_p=bsr_state->thread_buf + ((bsr_state->perthread->my_thread_id - 1) * bsr_state->per_thread_buffers);
    bsr_state->perthread->thread_buffer_index=0; // index within this threads block

    //
    // worker thread: set input file and send to rendering function
    //
    if (bsr_config.enable_external == 1) {
      processStars(&bsr_config, bsr_state, input_file_external);
    } // end if enable external
    if (bsr_config.enable_Gaia == 1) {
      processStars(&bsr_config, bsr_state, input_file_pq100);
      if (bsr_config.Gaia_min_parallax_quality < 100) {
        processStars(&bsr_config, bsr_state, input_file_pq050);
      }
      if (bsr_config.Gaia_min_parallax_quality < 50) {
        processStars(&bsr_config, bsr_state, input_file_pq030);
      }
      if (bsr_config.Gaia_min_parallax_quality < 30) {
        processStars(&bsr_config, bsr_state, input_file_pq020);
      }
      if (bsr_config.Gaia_min_parallax_quality < 20) {
        processStars(&bsr_config, bsr_state, input_file_pq010);
      }
      if (bsr_config.Gaia_min_parallax_quality < 10) {
        processStars(&bsr_config, bsr_state, input_file_pq005);
      }
      if (bsr_config.Gaia_min_parallax_quality < 05) {
        processStars(&bsr_config, bsr_state, input_file_pq003);
      }
      if (bsr_config.Gaia_min_parallax_quality < 03) {
        processStars(&bsr_config, bsr_state, input_file_pq002);
      }
      if (bsr_config.Gaia_min_parallax_quality < 02) {
        processStars(&bsr_config, bsr_state, input_file_pq001);
      }
      if (bsr_config.Gaia_min_parallax_quality < 01) {
        processStars(&bsr_config, bsr_state, input_file_pq000);
      }
    } // end if enable Gaia

    //
    // let main thread know we are done
    //
    bsr_state->status_array[bsr_state->perthread->my_thread_id]=1;
  } else {

    //
    // main thread: scan rendering thread buffers for pixels to integrate into image until all rendering threads are done
    //
    empty_passes=0;
    while (empty_passes < 2) { // do second pass once empty
      main_thread_buf_p=bsr_state->thread_buf;
      buffer_is_empty=1;
      // scan buffer for new pixel data
      for (main_thread_buffer_index=0; main_thread_buffer_index < thread_buffer_count; main_thread_buffer_index++) {
        if ((main_thread_buf_p->status_left == 1) && (main_thread_buf_p->status_right == 1)) {
          // buffer location has new pixel data, add to image composition buffer
          if (buffer_is_empty == 1) {
            buffer_is_empty=0; 
          }
          image_composition_p=bsr_state->image_composition_buf + main_thread_buf_p->image_offset;
          image_composition_p->r+=main_thread_buf_p->r;
          image_composition_p->g+=main_thread_buf_p->g;
          image_composition_p->b+=main_thread_buf_p->b;
          // set this buffer location to free
          main_thread_buf_p->status_left=0;
	  main_thread_buf_p->status_right=0;
        }
        main_thread_buf_p++;
      } // end for thread_buffer_index
      // if buffer is completely empty, check if all threads are done
      if (buffer_is_empty == 1) {
        all_workers_done=1;
        for (i=1; i <= bsr_state->num_worker_threads; i++) {
          if (bsr_state->status_array[i] < 1) {
            all_workers_done=0;
          }
        }
        if (all_workers_done == 1) {
          // if buffer is empty and all worker threads are done, increment empty_passes
          empty_passes++;
        }
      } 
    } // end while not done

    //
    // main thread: report rendering time if not in cgi mode
    //
    if (bsr_config.cgi_mode != 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf(" (%.3fs)\n", elapsed_time);
      fflush(stdout);
    }
  } // end if main thread

  //
  // all threads: post processing
  //
  postProcess(&bsr_config, bsr_state);

  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    //
    // main thread: output png file
    //
    writePNGFile(&bsr_config, bsr_state);

    //
    // main thread: clean up
    //
    if (input_file_external != NULL) {
      fclose(input_file_external);
    }
    if (input_file_pq100 != NULL) {
      fclose(input_file_pq100);
    }
    if (input_file_pq050 != NULL) {
      fclose(input_file_pq050);
    }
    if (input_file_pq030 != NULL) {
      fclose(input_file_pq030);
    }
    if (input_file_pq020 != NULL) {
      fclose(input_file_pq020);
    }
    if (input_file_pq010 != NULL) {
      fclose(input_file_pq010);
    }
    if (input_file_pq005 != NULL) {
      fclose(input_file_pq005);
    }
    if (input_file_pq003 != NULL) {
      fclose(input_file_pq003);
    }
    if (input_file_pq002 != NULL) {
      fclose(input_file_pq002);
    }
    if (input_file_pq001 != NULL) {
      fclose(input_file_pq001);
    }
    if (input_file_pq000 != NULL) {
      fclose(input_file_pq000);
    }
    if (bsr_state->image_composition_buf != NULL) {
      munmap(bsr_state->image_composition_buf, composition_buffer_size);
    }
    if (bsr_state->image_blur_buf != NULL) {
      munmap(bsr_state->image_blur_buf, blur_buffer_size);
    }
    if (bsr_state->image_resize_buf != NULL) {
      munmap(bsr_state->image_resize_buf, resize_buffer_size);
    }
    munmap(bsr_state->thread_buf, thread_buffer_size);
    munmap(bsr_state->status_array, status_array_size);
    free(bsr_state->Airymap_red);
    free(bsr_state->Airymap_green);
    free(bsr_state->Airymap_blue);

    //
    // main thread: output total runtime
    //
    if (bsr_config.cgi_mode != 1) {
      clock_gettime(CLOCK_REALTIME, &overall_endtime);
      elapsed_time=((double)(overall_endtime.tv_sec - 1500000000) + ((double)overall_endtime.tv_nsec / 1.0E9)) - ((double)(overall_starttime.tv_sec - 1500000000) + ((double)overall_starttime.tv_nsec) / 1.0E9);
      printf("Total execution time: %.3fs\n", elapsed_time);
      fflush(stdout);
    }

    //
    // clean up bsr_state
    //
    munmap(bsr_state, bsr_state_size);
  } // end if main thread

  return(0);
}

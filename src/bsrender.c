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
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
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
#include "post-process.h"
#include "image-composition.h"
#include "memory.h"
#include "quantize.h"

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
  struct timespec overall_starttime;
  struct timespec overall_endtime;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  int all_workers_done;
  int i;
  double rgb_red[32768];
  double rgb_green[32768];
  double rgb_blue[32768];
  pixel_composition_t *image_composition_p;
  int main_thread_buffer_index;
  int buffer_is_empty;
  int empty_passes;
  thread_buffer_t *main_thread_buf_p;

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
  // allocate memory for and initialize bsr_state
  //
  bsr_state=initState(&bsr_config);
  if (bsr_state == NULL) {
    if (bsr_config.cgi_mode != 1) {
      printf("Error: could initialize bsr_state\n");
      fflush(stdout);
    }
    exit(1);
  }
  bsr_state->perthread=&perthread;

  //
  // calculate number of worker threads to be forked
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
  // allocate memory and initialize various buffers that get attached to bsr_state
  //
  allocateMemory(&bsr_config, bsr_state);

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
      exit(1);
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
      exit(1);
    }
    if (bsr_config.Gaia_min_parallax_quality < 100) {
      sprintf(file_path, "%s/galaxy-pq050.dat", bsr_config.data_file_directory);
      input_file_pq050=fopen(file_path, "rb");
      if (input_file_pq050 == NULL) {
        if (bsr_config.cgi_mode != 1) {
          printf("Error: could not open %s\n", file_path);
          fflush(stdout);
        }
        exit(1);
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
        exit(1);
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
        exit(1);
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
        exit(1);
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
        exit(1);
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
        exit(1);
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
        exit(1);
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
        exit(1);
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
        exit(1);
      }
    }
  } // end if enable Gaia

  //
  // fork worker threads
  //
  bsr_state->master_pid=getpid();
  bsr_state->master_pgid=getpgrp();
  bsr_state->httpd_pid=getppid();
  if (bsr_state->num_worker_threads > 0) {
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->perthread->my_pid=getpid();
      if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
        bsr_state->perthread->my_thread_id=i; // this gets inherited by forked process
        bsr_state->status_array[i].status=THREAD_STATUS_INVALID;
        fork();
      }
    }
  }
  bsr_state->perthread->my_pid=getpid();
  bsr_state->status_array[bsr_state->perthread->my_thread_id].pid=bsr_state->perthread->my_pid;
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    bsr_state->perthread->my_thread_id=0;
  }

// 
// begin thread specific processing.
//

  //
  // all threads: initialize (clear) image composition buffer
  //
  initImageCompositionBuffer(&bsr_config, bsr_state);

  //
  // main thread: display begin rendering status if not in cgi mode
  //
  if ((bsr_state->perthread->my_pid == bsr_state->master_pid) && (bsr_config.cgi_mode != 1)) {
    clock_gettime(CLOCK_REALTIME, &starttime);
    printf("Rendering stars to image composition buffer...");
    fflush(stdout);
  }

  //
  // worker threads:  wait for main thread to say go
  // main thread: tell worker threads to go
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    waitForMainThread(bsr_state, THREAD_STATUS_PROCESS_STARS_BEGIN);
  } else {
    // main thread
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_PROCESS_STARS_BEGIN;
    }
  } // end if not main thread

  //
  // worker threads: process stars from binary data files
  //
  if (bsr_state->perthread->my_pid != bsr_state->master_pid) {
    //
    // worker threads: set our thread buffer postion to the beginning of this threads block
    //
    bsr_state->perthread->thread_buf_p=bsr_state->thread_buf + ((bsr_state->perthread->my_thread_id - 1) * bsr_state->per_thread_buffers);
    bsr_state->perthread->thread_buffer_index=0; // index within this threads block

    //
    // worker threads: set input file and send to rendering function
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
    // let main thread know we are done, then wait until main thread says ok to continue
    //
    bsr_state->status_array[bsr_state->perthread->my_thread_id].status=THREAD_STATUS_PROCESS_STARS_COMPLETE;
    waitForMainThread(bsr_state, THREAD_STATUS_PROCESS_STARS_CONTINUE);
  } else {
    //
    // main thread: scan worker thread buffers for pixels to integrate into image until all worker threads are done
    //
    empty_passes=0;
    while (empty_passes < 2) { // do second pass once empty
      // check if any worker threads have died
      checkExceptions(bsr_state);

      // scan buffer for new pixel data
      main_thread_buf_p=bsr_state->thread_buf;
      buffer_is_empty=1;
      for (main_thread_buffer_index=0; main_thread_buffer_index < bsr_state->thread_buffer_count; main_thread_buffer_index++) {
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
          if (bsr_state->status_array[i].status < THREAD_STATUS_PROCESS_STARS_COMPLETE) {
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
    // main thread: tell worker threads it's ok to continue
    //
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      bsr_state->status_array[i].status=THREAD_STATUS_PROCESS_STARS_CONTINUE;
    }

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

  //
  // all threads: convert image to 8 or 16 bits per pixel
  //
  quantize(&bsr_config, bsr_state);  

  //
  // main thread: output png image and cleanup
  //
  if (bsr_state->perthread->my_pid == bsr_state->master_pid) {
    //
    // main thread: output png file
    //
    outputPNG(&bsr_config, bsr_state);

    //
    // main thread: clean up file pointers
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

    //
    // main thread: clean up memory allocations
    //
    freeMemory(bsr_state);

    //
    // main thread: output total runtime
    //
    if (bsr_config.cgi_mode != 1) {
      clock_gettime(CLOCK_REALTIME, &overall_endtime);
      elapsed_time=((double)(overall_endtime.tv_sec - 1500000000) + ((double)overall_endtime.tv_nsec / 1.0E9)) - ((double)(overall_starttime.tv_sec - 1500000000) + ((double)overall_starttime.tv_nsec) / 1.0E9);
      printf("Total execution time: %.3fs\n", elapsed_time);
      fflush(stdout);
    }
  } // end if main thread

  return(0);
}

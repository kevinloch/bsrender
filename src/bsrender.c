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
#include "overlay.h"
#include "process-stars.h"
#include "init-state.h"

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
  bsr_state_t bsr_state;
  char file_path[1024];
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
  pid_t master_pid;
  pid_t mypid;
  int mmap_protection;
  int mmap_visibility;
  struct timespec overall_starttime;
  struct timespec overall_endtime;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  size_t composition_buffer_size;
  size_t status_array_size;
  int *status_array;
  int i;
  double rgb_red[32768];
  double rgb_green[32768];
  double rgb_blue[32768];
  int num_worker_threads; 
  pixel_composition_t *image_composition_p;
  int thread_buffer_size;
  int thread_buffer_count;
  int all_workers_done;
  thread_buffer_t *main_thread_buf_p;
  int main_thread_buffer_index;
  int buffer_is_empty;
  int empty_passes;
  long star_count;

  // initialize total execution timer
  clock_gettime(CLOCK_REALTIME, &overall_starttime);

  // initialize bsr_config to default values
  initConfig(&bsr_config);

  // process command line arguments
  processCmdArgs(&bsr_config, argc, argv);

  // load config file
  loadConfig(&bsr_config);

  num_worker_threads=bsr_config.num_threads-1;
  if (num_worker_threads < 1) {
    num_worker_threads=1;
  }

  // display major performance affecting options
  printf("Minimum parallax quality: %d, total threads: %d, buffers per rendering thread: %d stars\n", bsr_config.min_parallax_quality, (num_worker_threads + 1), bsr_config.per_thread_buffer);
  fflush(stdout);

  // initialize bsr_state (setup camera target)
  initState(&bsr_config, &bsr_state);

  // initialize RGB color lookup tables
  bsr_state.rgb_red=rgb_red;
  bsr_state.rgb_green=rgb_green;
  bsr_state.rgb_blue=rgb_blue;
  initRGBTables(&bsr_config, &bsr_state);

  // allocate memory for image composition buffer (floating point rgb)
  clock_gettime(CLOCK_REALTIME, &starttime);
  printf("Initializing image composition buffer...");
  fflush(stdout);
  composition_buffer_size=bsr_config.camera_res_x * bsr_config.camera_res_y * sizeof(pixel_composition_t);
  bsr_state.image_composition_buf=(pixel_composition_t *)malloc(composition_buffer_size);
  if (bsr_state.image_composition_buf == NULL) {
    printf("Error: could not allocate shared memory for image composition buffer\n");
    return(1);
  }
  // initialize (clear) image composition buffer
  image_composition_p=bsr_state.image_composition_buf;
  for (i=0; i < (bsr_config.camera_res_x * bsr_config.camera_res_y); i++) {
    image_composition_p->r=0.0;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
    image_composition_p++;
  }
  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
  printf(" (%.3fs)\n", elapsed_time);
  fflush(stdout);

  // allocate shared memory for thread buffer
  clock_gettime(CLOCK_REALTIME, &starttime);
  printf("Initializing rendering thread buffers...");
  fflush(stdout);
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;
  if (bsr_config.per_thread_buffer < 1) {
    bsr_config.per_thread_buffer=1;
  }
  thread_buffer_count=num_worker_threads * bsr_config.per_thread_buffer;
  thread_buffer_size=thread_buffer_count * sizeof(thread_buffer_t);
  bsr_state.thread_buf=(thread_buffer_t *)mmap(NULL, thread_buffer_size, mmap_protection, mmap_visibility, -1, 0);
  if (bsr_state.thread_buf == NULL) {
    printf("Error: could not allocate shared memory for thread buffer\n");
  }
  // initialize thread buffer
  main_thread_buf_p=bsr_state.thread_buf;
  for (i=0; i < (num_worker_threads * bsr_config.per_thread_buffer); i++) {
    main_thread_buf_p->status_left=0;
    main_thread_buf_p->status_right=0;
    main_thread_buf_p++;
  }

  // allocate shared memory for thread status array
  status_array_size=num_worker_threads * sizeof(int);
  status_array=(int *)mmap(NULL, status_array_size, mmap_protection, mmap_visibility, -1, 0);
  if (status_array == NULL) {
    printf("Error: could not allocate shared memory for thread status array\n");
  }

  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
  printf(" (%.3fs)\n", elapsed_time);
  fflush(stdout);

  // attempt to open input file(s)
  sprintf(file_path, "%s/galaxy-pq100.dat", bsr_config.data_file_directory);
  input_file_pq100=fopen(file_path, "rb");
  if (input_file_pq100 == NULL) {
    printf("Error: could not open %s\n", file_path);
    fflush(stdout);
    return(1);
  }
  if (bsr_config.min_parallax_quality < 100) {
    sprintf(file_path, "%s/galaxy-pq050.dat", bsr_config.data_file_directory);
    input_file_pq050=fopen(file_path, "rb");
    if (input_file_pq050 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 50) {
    sprintf(file_path, "%s/galaxy-pq030.dat", bsr_config.data_file_directory);
    input_file_pq030=fopen(file_path, "rb");
    if (input_file_pq030 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 30) {
    sprintf(file_path, "%s/galaxy-pq020.dat", bsr_config.data_file_directory);
    input_file_pq020=fopen(file_path, "rb");
    if (input_file_pq020 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 20) {
    sprintf(file_path, "%s/galaxy-pq010.dat", bsr_config.data_file_directory);
    input_file_pq010=fopen(file_path, "rb");
    if (input_file_pq010 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 10) {
    sprintf(file_path, "%s/galaxy-pq005.dat", bsr_config.data_file_directory);
    input_file_pq005=fopen(file_path, "rb");
    if (input_file_pq005 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 5) {
    sprintf(file_path, "%s/galaxy-pq003.dat", bsr_config.data_file_directory);
    input_file_pq003=fopen(file_path, "rb");
    if (input_file_pq003 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 3) {
    sprintf(file_path, "%s/galaxy-pq002.dat", bsr_config.data_file_directory);
    input_file_pq002=fopen(file_path, "rb");
    if (input_file_pq002 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 2) {
    sprintf(file_path, "%s/galaxy-pq001.dat", bsr_config.data_file_directory);
    input_file_pq001=fopen(file_path, "rb");
    if (input_file_pq001 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }
  if (bsr_config.min_parallax_quality < 1) {
    sprintf(file_path, "%s/galaxy-pq000.dat", bsr_config.data_file_directory);
    input_file_pq000=fopen(file_path, "rb");
    if (input_file_pq000 == NULL) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
      return(1);
    }
  }

  clock_gettime(CLOCK_REALTIME, &starttime);
  printf("Rendering stars to image composition buffer...");
  fflush(stdout);

  // fork parallel rendering processes
  master_pid=getpid();
  if (num_worker_threads > 0) {
    for (i=0; i < num_worker_threads; i++) {
      mypid=getpid();
      if (mypid == master_pid) {
        bsr_state.my_thread_id=i; // this gets inherited by forked process
        status_array[i]=0;
        fork();
      }
    }
  }
  mypid=getpid();

  if (mypid != master_pid) {
    // we are a worker thread
    // set our thread buffer postion to the beginning of this threads block
    bsr_state.thread_buf_p=bsr_state.thread_buf + (bsr_state.my_thread_id * bsr_config.per_thread_buffer);
    bsr_state.thread_buffer_index=0; // index within this threads block

    // main rendering function reads stars from input file and renders to image composition buffer
    processStars(&bsr_config, &bsr_state, input_file_pq100);
    if (bsr_config.min_parallax_quality < 100) {
      processStars(&bsr_config, &bsr_state, input_file_pq050);
    }
    if (bsr_config.min_parallax_quality < 50) {
      processStars(&bsr_config, &bsr_state, input_file_pq030);
    }
    if (bsr_config.min_parallax_quality < 30) {
      processStars(&bsr_config, &bsr_state, input_file_pq020);
    }
    if (bsr_config.min_parallax_quality < 20) {
      processStars(&bsr_config, &bsr_state, input_file_pq010);
    }
    if (bsr_config.min_parallax_quality < 10) {
      processStars(&bsr_config, &bsr_state, input_file_pq005);
    }
    if (bsr_config.min_parallax_quality < 05) {
      processStars(&bsr_config, &bsr_state, input_file_pq003);
    }
    if (bsr_config.min_parallax_quality < 03) {
      processStars(&bsr_config, &bsr_state, input_file_pq002);
    }
    if (bsr_config.min_parallax_quality < 02) {
      processStars(&bsr_config, &bsr_state, input_file_pq001);
    }
    if (bsr_config.min_parallax_quality < 01) {
      processStars(&bsr_config, &bsr_state, input_file_pq000);
    }
    // let main thread know we are done
    status_array[bsr_state.my_thread_id]=1;
  } else {
    // we are the master thread, process thread buffer entries until all worker threads have completed and buffer is empty
    star_count=0;
    empty_passes=0;
    while (empty_passes < 2) { // do second pass once empty
      main_thread_buf_p=bsr_state.thread_buf;
      buffer_is_empty=1;
      // scan buffer for new pixel data
      for (main_thread_buffer_index=0; main_thread_buffer_index < thread_buffer_count; main_thread_buffer_index++) {
        if ((main_thread_buf_p->status_left == 1) && (main_thread_buf_p->status_right == 1)) {
          // buffer location has new pixel data, add to image composition buffer
          star_count++;
          if (buffer_is_empty == 1) {
            buffer_is_empty=0; 
          }
          image_composition_p=bsr_state.image_composition_buf + main_thread_buf_p->image_offset;
          image_composition_p->r+=main_thread_buf_p->r;
          image_composition_p->g+=main_thread_buf_p->g;
          image_composition_p->b+=main_thread_buf_p->b;
          // set this buffer location to free
          main_thread_buf_p->status_left=0;
	  main_thread_buf_p->status_right=0;
        }
        main_thread_buf_p++;
      } // end for thread_buffer_index

      if (buffer_is_empty == 1) {
        // if buffer is completely empty, check if all threads are done
        all_workers_done=1;
        for (i=0; i < num_worker_threads; i++) {
          if (status_array[i] == 0) {
            all_workers_done=0;
          }
        }
        if (all_workers_done == 1) {
          // if buffer is empty and all worker threads are done, increment empty_passes
          empty_passes++;
        }
      } 
    } // end while 

    // optionally draw overlays
    if (bsr_config.draw_cross_hairs == 1) {
      drawCrossHairs(&bsr_config, &bsr_state);
    }
    if (bsr_config.draw_grid_lines == 1) {
      drawGridLines(&bsr_config, &bsr_state);
    }

    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);

    // output png file
    writePNGFile(&bsr_config, &bsr_state);

    // clean up
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
    free(bsr_state.image_composition_buf);
    munmap(bsr_state.thread_buf, thread_buffer_size);
    munmap(status_array, status_array_size);

    // output total runtime
    clock_gettime(CLOCK_REALTIME, &overall_endtime);
    elapsed_time=((double)(overall_endtime.tv_sec - 1500000000) + ((double)overall_endtime.tv_nsec / 1.0E9)) - ((double)(overall_starttime.tv_sec - 1500000000) + ((double)overall_starttime.tv_nsec) / 1.0E9);
    printf("Rendered %ld stars to %.2f megapixels, total runtime: %.3fs\n", star_count, ((double)(bsr_config.camera_res_x * bsr_config.camera_res_y) / 1.0E6), elapsed_time);
    fflush(stdout);

  } // end if master process
}

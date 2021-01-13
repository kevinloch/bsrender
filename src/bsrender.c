//
// Billion Star 3D Rendering Engine Proof of Concept
// Kevin M. Loch
//
// 3D rendering engine for the ESA Gaia EDR3 star dataset
//
// This rendering engine reads star data from binary file 'galaxy.dat' and outputs an ascii PPM file as output.
// See comments in mkgalaxy.c for the binary data file format.
// In this proof of concept all options must be changed here in the source code and recompiled.
// The output ppm file can be converted to png with the program pnmtopng from the netpbm collection
// Support for direct png file output, config file and command line options will be added as development progresses
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
  FILE *input_file;
  double two_pi=2.0 * M_PI;
  double pi_over_2=M_PI / 2.0;
  double pi_over_180=M_PI / 180.0;
  double pi_over_360=M_PI / 360.0;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  pid_t master_pid;
  pid_t child_pid;
  pid_t mypid;
  int mmap_protection;
  int mmap_visibility;
  size_t composition_buffer_size;
  size_t status_array_size;
  int *status_array;
  int my_thread_id=0;
  int done;
  int i;
  double star_x;
  double star_y;
  double star_z;
  double star_r;
  double star_xy_r;
  double star_xz_r;
  double star_yz_r;
  double target_x;
  double target_y;
  double target_z;
  //double target_r;
  double target_xy_r;
  //double target_xz_r;
  //double target_yz_r;
  double output_az;
  double output_el;
  int output_x=0;
  int output_y=0;
  double star_linear_intensity; // star linear intensity as viewed from camera
  int image_offset;
  float star_linear_1pc_intensity;
  uint64_t color_temperature;
  double rgb_red[32768];
  double rgb_green[32768];
  double rgb_blue[32768];
  double mollewide_angle;
  double two_mollewide_angle;
  double camera_icrs_ra_rad;
  double camera_icrs_dec_rad;
  double target_icrs_ra_rad;
  double target_icrs_dec_rad;
  double camera_half_res_x;   // half of x-axis resolution
  double camera_half_res_y;   // half of y-axis resolution
  double camera_hfov;         // half of fov angle, to speed calculations (rad)
  double pixels_per_radian;   // pixels per radian for x and y axis
  double camera_3az_yz;       // rotation of camera (radians)
  double camera_3az_xy;       // pan of camera (radians)
  double camera_3az_xz;       // tilt of camera (radians)
  double render_distance;            // distance from selected point to star
  int num_worker_threads;            // number of worker threads to fork from master

  // input file record
  star_record_t star_record;
  int star_record_size=sizeof(star_record_t);

  // image buffers
  pixel_composition_t *image_composition_buf;
  pixel_composition_t *image_composition_p;

  // variables for triple-azimuth (3az) coordinates are derived from x,y,z positions
  // triple-azimuth coordinates allow for fast and accurate rotations.  They are a hybrid of Euler angles and quaternions/versors.
  // Three 360 degree angles provide redundancy to avoid gimbal lock or degredation around poles (like quaternions/versors).
  // Using orthoganal angles (like Eeuler rotations) allows for simple and fast rendering once rotations are complete.

  // 3az coordinates of the camera target direction as seen from camera.   Distance is not used and not calculated
  // these values are also the rotation factors applied to stars
  double target_3az_xy;
  double target_3az_xz;
  //double target_3az_yz;

  // 3az coordinates of star as seen from camera (before coordinate rotation)
  double star_3az_xy;
  double star_3az_xz;
  double star_3az_yz;

  // initialize bsr_config to default values
  initConfig(&bsr_config);

  // process command line arguments
  processCmdArgs(&bsr_config, argc, argv);

  // load config file
  if (loadConfig(&bsr_config) == 1) {
    exit(1);
  }

  // get master process id for when we fork
  master_pid=getpid();
  mmap_protection=PROT_READ | PROT_WRITE;
  mmap_visibility=MAP_SHARED | MAP_ANONYMOUS;

  // allocate shared memory for image composition buffer (floating point rgb)
  clock_gettime(CLOCK_REALTIME, &starttime);
  printf("Initializing image composition buffer...");
  fflush(stdout);
  composition_buffer_size=bsr_config.camera_res_x * bsr_config.camera_res_y * sizeof(pixel_composition_t);
  image_composition_buf=(pixel_composition_t *)mmap(NULL, composition_buffer_size, mmap_protection, mmap_visibility, -1, 0);
  if (image_composition_buf == NULL) {
    printf("Error: could not allocate shared memory for image composition buffer\n");
    return(1);
  }
  // initialize image composition buffer
  image_composition_p=image_composition_buf;
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

  // allocate shared memory for thread status array
  num_worker_threads=bsr_config.num_processes - 1;
  if (num_worker_threads > 0) {
    status_array_size=num_worker_threads * sizeof(int);
    status_array=(int *)mmap(NULL, status_array_size, mmap_protection, mmap_visibility, -1, 0);
    if (status_array == NULL) {
      printf("Error: could not allocate shared memory for thread status array\n");
    }
  }

  // initialize RGB color lookup tables
  initRGBTables(&bsr_config, rgb_red, rgb_green, rgb_blue);

/*
  printf("Opening galaxy.dat\n");
  fflush(stdout);
*/

  // attempt to open input file(s)
  input_file=fopen("galaxy.dat", "rb");
  if (input_file == NULL) {
    printf("Error: could not open galaxy.dat\n");
    fflush(stdout);
    return(1);
  }

  clock_gettime(CLOCK_REALTIME, &starttime);
  printf("Rendering stars to image composition buffer...");
  fflush(stdout);

  // process user-supplied arguments
  camera_hfov=bsr_config.camera_fov * pi_over_360; // includes divide by 2
  camera_half_res_x=(double)bsr_config.camera_res_x / 2.0;
  camera_half_res_y=(double)bsr_config.camera_res_y / 2.0;
  pixels_per_radian=camera_half_res_x / camera_hfov;
  camera_3az_yz=bsr_config.camera_rotation * pi_over_180;
  camera_3az_xy=bsr_config.camera_pan * pi_over_180;
  camera_3az_xz=bsr_config.camera_tilt * -pi_over_180;

  // optionally transform spherical icrs to euclidian icrs
  if ((bsr_config.camera_icrs_dec != 0) || (bsr_config.camera_icrs_ra != 0)) {
    camera_icrs_ra_rad=bsr_config.camera_icrs_ra * pi_over_180;
    camera_icrs_dec_rad=bsr_config.camera_icrs_dec * pi_over_180;
    bsr_config.camera_icrs_x=bsr_config.camera_icrs_r * cos(camera_icrs_dec_rad) * cos(camera_icrs_ra_rad);
    bsr_config.camera_icrs_y=bsr_config.camera_icrs_r * cos(camera_icrs_dec_rad) * sin(camera_icrs_ra_rad);
    bsr_config.camera_icrs_z=bsr_config.camera_icrs_r * sin(camera_icrs_dec_rad);
  }
  if ((bsr_config.target_icrs_dec != 0) || (bsr_config.target_icrs_ra != 0)) {
    target_icrs_ra_rad=bsr_config.target_icrs_ra * pi_over_180;
    target_icrs_dec_rad=bsr_config.target_icrs_dec * pi_over_180;
    bsr_config.target_icrs_x=bsr_config.target_icrs_r * cos(target_icrs_dec_rad) * cos(target_icrs_ra_rad);
    bsr_config.target_icrs_y=bsr_config.target_icrs_r * cos(target_icrs_dec_rad) * sin(target_icrs_ra_rad);
    bsr_config.target_icrs_z=bsr_config.target_icrs_r * sin(target_icrs_dec_rad);
  }

/*
  printf("debug, camera_res_x: %d, camera_res_y: %d, camera_fov: %.4e, camera hfov: %.4e, pixels_per_radian: %.4e\n", bsr_config.camera_res_x, bsr_config.camera_res_y, bsr_config.camera_fov, camera_hfov, pixels_per_radian);
  fflush(stdout);
*/

  // convert camera target x,y,z to 3az coordinates as seen from camera
  target_x=bsr_config.target_icrs_x - bsr_config.camera_icrs_x;
  target_y=bsr_config.target_icrs_y - bsr_config.camera_icrs_y;
  target_z=bsr_config.target_icrs_z - bsr_config.camera_icrs_z;
  //target_r=sqrt(pow(target_x, 2.0) + pow(target_y, 2.0) + pow(target_z, 2.0));
  target_xy_r=sqrt(pow(target_x, 2.0) + pow(target_y, 2.0)); // may be used in future rotations or raster projections
  //target_xz_r=sqrt(pow(target_x, 2.0) + pow(target_z, 2.0)); // may be used in future rotations or raster projections
  //target_yz_r=sqrt(pow(target_y, 2.0) + pow(target_z, 2.0)); // may be used in future rotations or raster projections

  if ((target_x == 0.0) && (target_y == 0.0)) {
    target_3az_xy=0.0;
  } else {
    target_3az_xy=atan2(target_y, target_x);
  }
  if ((target_x == 0.0) && (target_z == 0.0)) {
    target_3az_xz=0.0;
  } else {
    target_3az_xz=atan2(target_z, target_x);
  }
  //if ((target_y == 0.0) && (target_z == 0.0)) {
  //  target_3az_yz=0.0;
  //} else {
  //  target_3az_yz=atan2(target_z, target_y);
  //}

  // calculate target xz rotation after setting xy=0
  // apply target xy rotation angle to target xz angle, only x has changed
  target_x=target_xy_r; // xy=0
  //target_xz_r=sqrt(pow(target_x, 2.0) + pow(target_z, 2.0)); // may be used in future rotations or raster projections
  if ((target_x == 0.0) && (target_z == 0.0)) {
    target_3az_xz=0.0;
  } else {
    target_3az_xz=atan2(target_z, target_x);
  }

  // fork parallel rendering processes
  if (num_worker_threads > 0) {
    for (i=0; i < num_worker_threads; i++) {
      mypid=getpid();
      if (mypid == master_pid) {
        my_thread_id=i; // this gets inherited by forked process
        status_array[i]=0;
        child_pid=fork();
      }
    }
  }

  // read star record from input file
  fread(&star_record, star_record_size, 1, input_file);

  // read and process each line of input file and render stars to image
  while ((feof(input_file) == 0) && (ferror(input_file) == 0)) {

    // convert star x,y,z to 3az coordinates as seen by camera position
    star_x=star_record.icrs_x - bsr_config.camera_icrs_x;
    star_y=star_record.icrs_y - bsr_config.camera_icrs_y;
    star_z=star_record.icrs_z - bsr_config.camera_icrs_z;
    star_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0) + pow(star_z, 2.0));

/*
    printf("debug, star_icrs_x: %.9f, star_icrs_y: %.9f, star_icrs_z: %.9f, star_x: %.9f, star_y: %.9f, star_z: %.9f, star_r: %.9f\n", star_record.icrs_x, star_record.icrs_y, star_record.icrs_z, star_x, star_y, star_z, star_r);
    fflush(stdout);
*/

    // only continue if star distance is within min-max range from selected point (0=camera, 1=target)
    if (bsr_config.render_distance_selector == 0) { // selected point is camera
      render_distance=star_r; // star distance from camera
    } else { // selected point is target
      render_distance=sqrt(pow((star_record.icrs_x - bsr_config.target_icrs_x), 2.0) + pow((star_record.icrs_y - bsr_config.target_icrs_y), 2.0) + pow((star_record.icrs_z - bsr_config.target_icrs_z), 2.0)); // important, use un-rotated coordinates
    }
    if ((render_distance >= bsr_config.render_distance_min) && (render_distance <= bsr_config.render_distance_max)) {
      star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
      //star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections

/*
      printf("debug, star_icrs_x: %.9f, star_icrs_y: %.9f, star_icrs_z: %.9f, star_x: %.9f, star_y: %.9f, star_z: %.9f, star_r: %.9f, star_xy_r: %.9f, star_xz_r: %.9f, star_yz_r: %.9f\n", star_record.icrs_x, star_record.icrs_y, star_record.icrs_z, star_x, star_y, star_z, star_r, star_xy_r, star_xz_r, star_yz_r);
      fflush(stdout);
*/

      if ((star_x == 0.0) && (star_y == 0.0)) {
        star_3az_xy=0.0;
      } else {
        star_3az_xy=atan2(star_y, star_x);
      }
      if ((star_x == 0.0) && (star_z == 0.0)) {
        star_3az_xz=0.0;
      } else {
        star_3az_xz=atan2(star_z, star_x);
      }
      if ((star_y == 0.0) && (star_z == 0.0)) {
        star_3az_yz=0.0;
      } else {
        star_3az_yz=atan2(star_z, star_y);
      }

      // rotate star by camera target in xy and xz planes so camera target = (xy=0,xz=0) and yz is camera rotation angle
      // This will allow for direct mapping of xy,xz angles as polar coordinates to our camera raster x,y axis and we can also easily tell if a star is within fov

      // rotate star xy angle by target xy angle
      star_3az_xy-=target_3az_xy;
      if (star_3az_xy > M_PI) {
        star_3az_xy = -two_pi + star_3az_xy;
      } else if (star_3az_xy < -M_PI) {
        star_3az_xy = two_pi + star_3az_xy;
      }
      // apply star xy rotation angle to star xz angle, only x has changed
      star_x=star_xy_r * cos(star_3az_xy);
      star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_z == 0.0)) {
        star_3az_xz=0.0;
      } else {
        star_3az_xz=atan2(star_z, star_x);
      }
      // apply star xy rotation angle to star yz angle, only y has changed
      star_y=star_xy_r * sin(star_3az_xy);
      //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_y == 0.0) && (star_z == 0.0)) {
        star_3az_yz=0.0;
      } else {
        star_3az_yz=atan2(star_z, star_y);
      }

      // rotate star xz angle by (rotated) target xz angle
      star_3az_xz-=target_3az_xz;
      if (star_3az_xz > M_PI) {
        star_3az_xz = -two_pi + star_3az_xz;
      } else if (star_3az_xz < -M_PI) {
        star_3az_xz = two_pi + star_3az_xz;
      }
      // apply star xz rotation to star xy angle, only x has changed
      star_x=star_xz_r * cos(star_3az_xz);
      //star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_y == 0.0)) {
        star_3az_xy=0.0;
      } else {
        star_3az_xy=atan2(star_y, star_x);
      }
      // apply star xz rotation to star yz angle, only z has changed
      star_z=star_xz_r * sin(star_3az_xz);
      star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_y == 0.0) && (star_z == 0.0)) {
        star_3az_yz=0.0;
      } else {
        star_3az_yz=atan2(star_z, star_y);
      }

      // rotate star yz angle by camera rotation angle
      star_3az_yz+=camera_3az_yz;
      if (star_3az_yz > M_PI) {
        star_3az_yz = -two_pi + star_3az_yz;
      } else if (star_3az_yz < -M_PI) {
        star_3az_yz = two_pi + star_3az_yz;
      }
      // apply star yz rotation to star xy angle, only y has changed
      star_y=star_yz_r * cos(star_3az_yz);
      star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_y == 0.0)) {
        star_3az_xy=0.0;
      } else {
        star_3az_xy=atan2(star_y, star_x);
      }
      // apply star yz rotation to star xz angle, only z has changed
      star_z=star_yz_r * sin(star_3az_yz);
      //star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
      if ((star_x == 0.0) && (star_z == 0.0)) {
        star_3az_xz=0.0;
      } else {
        star_3az_xz=atan2(star_z, star_x);
      }

      if (bsr_config.camera_pan != 0) {
        // optionally rotate star xy angle by camera pan angle
        star_3az_xy+=camera_3az_xy;
        if (star_3az_xy > M_PI) {
          star_3az_xy = -two_pi + star_3az_xy;
        } else if (star_3az_xy < -M_PI) {
          star_3az_xy = two_pi + star_3az_xy;
        }
        // apply star xy rotation angle to star xz angle, only x has changed
        star_x=star_xy_r * cos(star_3az_xy);
        star_xz_r=sqrt(pow(star_x, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
        if ((star_x == 0.0) && (star_z == 0.0)) {
          star_3az_xz=0.0;
        } else {
          star_3az_xz=atan2(star_z, star_x);
        }
        // apply star xy rotation angle to star yz angle, only y has changed
        star_y=star_xy_r * sin(star_3az_xy);
        //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
        if ((star_y == 0.0) && (star_z == 0.0)) {
          star_3az_yz=0.0;
        } else {
          star_3az_yz=atan2(star_z, star_y);
        }
      }

      if (bsr_config.camera_tilt != 0.0) {
        // optionally rotate star xz angle by camera tilt angle
        star_3az_xz+=camera_3az_xz;
        if (star_3az_xz > M_PI) {
          star_3az_xz = -two_pi + star_3az_xz;
        } else if (star_3az_xz < -M_PI) {
          star_3az_xz = two_pi + star_3az_xz;
        }
        // apply star xz rotation to star xy angle, only x has changed
        star_x=star_xz_r * cos(star_3az_xz);
        star_xy_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0)); // may be used in future rotations or raster projections
        if ((star_x == 0.0) && (star_y == 0.0)) {
          star_3az_xy=0.0;
        } else {
          star_3az_xy=atan2(star_y, star_x);
        }
        // apply star xz rotation to star yz angle, only z has changed
        star_z=star_xz_r * sin(star_3az_xz);
        //star_yz_r=sqrt(pow(star_y, 2.0) + pow(star_z, 2.0)); // may be used in future rotations or raster projections
        if ((star_y == 0.0) && (star_z == 0.0)) {
          star_3az_yz=0.0;
        } else {
          star_3az_yz=atan2(star_z, star_y);
        }

      }

      // at this point we can filter by stars that are within our fov
      if ((fabs(star_3az_xy) < camera_hfov) && (fabs(star_3az_xz) < camera_hfov)) {

/*
        printf("intensity_and_temperature: %16lx\n", star_record.intensity_and_temperature);
        fflush(stdout);
*/

        // extract color_temperature from combined field
        color_temperature=star_record.intensity_and_temperature & 0x00000000fffffffful;

        // extract star intensity from combined field and adjust for distance from camera
        star_record.intensity_and_temperature >>=32;
        star_linear_1pc_intensity=*((float*)&star_record.intensity_and_temperature);
        star_linear_intensity=star_linear_1pc_intensity / pow(star_r, 2.0);

/*
        printf("star_linear_1pc_intensity: %.4e, distance: %.4e, star_linear_intensity: %.4e, color: %ld\n", star_linear_1pc_intensity, star_r, star_linear_intensity, color_temperature);
        fflush(stdout);
*/

        //
        // render star to floating point image composition buffer
        //
        output_el=atan2(star_z, star_xy_r);
        output_az=star_3az_xy;

        if (bsr_config.camera_projection == 0) {
          // lat/lon 
          output_x=(int)((-pixels_per_radian * output_az) + camera_half_res_x + 0.5);
          output_y=(int)((-pixels_per_radian * output_el) + camera_half_res_y + 0.5);
        } else if (bsr_config.camera_projection == 1) {
          // spherical TBD
          output_x=0;
          output_y=0;
        } else if (bsr_config.camera_projection == 2) {
          // Mollewide
          two_mollewide_angle=2.0 * asin(2.0 * output_el / M_PI);
          for (i=0; i < bsr_config.Mollewide_iterations; i++) {
            two_mollewide_angle-=(two_mollewide_angle + sin(two_mollewide_angle) - (M_PI * sin(output_el))) / (1.0 + cos(two_mollewide_angle));
          }
          mollewide_angle=two_mollewide_angle * 0.5;
          output_x=(int)((-pixels_per_radian * output_az * cos(mollewide_angle)) + camera_half_res_x + 0.5);
          output_y=(int)((-pixels_per_radian * pi_over_2 * sin(mollewide_angle)) + camera_half_res_y + 0.5);
        }

        if ((output_x >= 0) && (output_x < bsr_config.camera_res_x) && (output_y >= 0) && (output_y < bsr_config.camera_res_y)) {
          image_offset=(bsr_config.camera_res_x * output_y) + output_x;
          image_composition_p=image_composition_buf + image_offset;
          image_composition_p->r+=(star_linear_intensity * rgb_red[color_temperature]);
          image_composition_p->g+=(star_linear_intensity * rgb_green[color_temperature]);
          image_composition_p->b+=(star_linear_intensity * rgb_blue[color_temperature]);

/*
          printf("debug, x: %d, y: %d, image_offset: %d, flux: %.4e, temp: %d, r: %.4e, g: %.4e, b: %.4e, rgb_red: %d, rgb_green: %d, rgb_blue: %d\n", output_x, output_y, image_offset, star_linear_intensity, star_record.color_temperature, image_composition_p->r, image_composition_p->g, image_composition_p->b, rgb_red[star_record.color_temperature], rgb_green[star_record.color_temperature], rgb_blue[star_record.color_temperature]);
          fflush(stdout);
*/

        } // end if within image raster
      } // end if within fov
    } // end if within distance ranges

    // read star record from input file
    fread(&star_record, star_record_size, 1, input_file);
  } // end input loop

  mypid=getpid();
  if (mypid != master_pid) {
    // if i'm a worker thread and I made it this far set my thread status to complete (1)
    status_array[my_thread_id]=1;
  } else {
    // wait until all worker threads have completed
    if (num_worker_threads > 0) {
      done=waitForWorkerThreads(status_array, num_worker_threads);
    } // end if num_worker_threads

    if (bsr_config.draw_cross_hairs == 1) {
      drawCrossHairs(&bsr_config, image_composition_buf, camera_half_res_x, camera_half_res_y);
    }

    if (bsr_config.draw_grid_lines == 1) {
      drawGridLines(&bsr_config, image_composition_buf, camera_half_res_x, camera_half_res_y);
    }

    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    printf(" (%.3fs)\n", elapsed_time);
    fflush(stdout);

    writePNGFile(&bsr_config, image_composition_buf);

    // clean up
    fclose(input_file);
    munmap(image_composition_buf, composition_buffer_size);
    if (num_worker_threads > 0) {
      munmap(status_array, status_array_size);
    }
  } // end if master process
}

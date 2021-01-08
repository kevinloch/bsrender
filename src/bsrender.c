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

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define OUTPUT_PPM

// generate RGB color values for a given blackbody temperature
// set temp to desired white balance temperature

int initRGBTables(double camera_wb_temp, double camera_color_saturation, double rgb_red[], double rgb_green[], double rgb_blue[]) {
  int i;

  // CIE wavelengths
  double red_center=700.0E-9;
  double green_center=546.1E-9;
  double blue_center=435.8E-9;

  double red_freq;
  double green_freq;
  double blue_freq;
  double red_intensity;
  double green_intensity;
  double blue_intensity;
  double kb=1.380649E-23;
  double h=6.62607015E-34;
  double c=299792458.0;
  double temp;
  double red_wb_factor;
  double green_wb_factor;
  double blue_wb_factor;
  double normalization_factor;
  double color_max;
  double color_min;
  double color_mid;

  red_freq=c / red_center;
  green_freq=c / green_center;
  blue_freq=c / blue_center;

  // calculate white balance factors
  red_intensity=  (2.0 * h * pow(red_freq,   3.0) / pow(c, 2.0)) / (exp(h * red_freq   / (kb * camera_wb_temp)) - 1);
  green_intensity=(2.0 * h * pow(green_freq, 3.0) / pow(c, 2.0)) / (exp(h * green_freq / (kb * camera_wb_temp)) - 1);
  blue_intensity= (2.0 * h * pow(blue_freq,  3.0) / pow(c, 2.0)) / (exp(h * blue_freq  / (kb * camera_wb_temp)) - 1);
  green_wb_factor=1.0 / green_intensity;
  red_wb_factor=1.0 / red_intensity;
  blue_wb_factor=1.0 / blue_intensity;
  
  // calculate rgb values for each integer Kelving temp from 0 - 32767K
  for (i=0; i < 32768; i++) {
    temp=(double)i;
    red_intensity=red_wb_factor     * (2.0 * h * pow(red_freq,   3.0) / pow(c, 2.0)) / (exp(h * red_freq   / (kb * temp)) - 1);
    green_intensity=green_wb_factor * (2.0 * h * pow(green_freq, 3.0) / pow(c, 2.0)) / (exp(h * green_freq / (kb * temp)) - 1);
    blue_intensity= blue_wb_factor  * (2.0 * h * pow(blue_freq,  3.0) / pow(c, 2.0)) / (exp(h * blue_freq  / (kb * temp)) - 1);

    // calculate normalization factor so r+g+b=1
    if (red_intensity == 0) {
      red_intensity=1.0;
    }
    normalization_factor=1.0 / (red_intensity + green_intensity + blue_intensity);
    red_intensity=red_intensity * normalization_factor;
    green_intensity=green_intensity * normalization_factor;
    blue_intensity=blue_intensity * normalization_factor;

    // apply camera color saturation adjustment
    color_max=-1.0;
    if (red_intensity > color_max) {
      color_max=red_intensity;
    } 
    if (green_intensity > color_max) {
      color_max=green_intensity;
    } 
    if (blue_intensity > color_max) {
      color_max=blue_intensity;
    } 
    color_min=2.0;
    if (red_intensity < color_min) {
      color_min=red_intensity;
    }
    if (green_intensity < color_min) {
      color_min=green_intensity;
    }
    if (blue_intensity < color_min) {
      color_min=blue_intensity;
    } 
    color_mid=(color_max + color_min) / 2.0;
    red_intensity=color_mid + (camera_color_saturation * (red_intensity - color_mid));
    if (red_intensity < 0) {
      red_intensity=0;
    }
    green_intensity=color_mid + (camera_color_saturation * (green_intensity - color_mid));
    if (green_intensity < 0) {
      green_intensity=0;
    }
    blue_intensity=color_mid + (camera_color_saturation * (blue_intensity - color_mid));
    if (blue_intensity < 0) {
      blue_intensity=0;
    }

    // re-normalize and store in rgb arrays
    normalization_factor=1.0 / (red_intensity + green_intensity + blue_intensity);
    rgb_red[i]=red_intensity * normalization_factor;
    rgb_green[i]=green_intensity * normalization_factor;
    rgb_blue[i]=blue_intensity * normalization_factor;
  } // end for i
  return(0);
}

int limitIntensityPreserveColor(double *pixel_r, double *pixel_g, double *pixel_b) {
  double pixel_max;

  // limit pixel to range 0.0-1.0 while maintaining color (max channel=1.0)
  if (*pixel_r < 0) {
    *pixel_r=0;
  }
  if (*pixel_g < 0) {
    *pixel_g=0;
  }
  if (*pixel_b < 0) {
    *pixel_b=0;
  }
  if ((*pixel_r > 1.0) || (*pixel_g > 1.0) || (*pixel_b > 1.0)) {
    pixel_max=0;
    if (*pixel_r > pixel_max) {
      pixel_max=*pixel_r;
    }
    if (*pixel_g > pixel_max) {
      pixel_max=*pixel_g;
    }
    if (*pixel_b > pixel_max) {
      pixel_max=*pixel_b;
    }
    *pixel_r=*pixel_r / pixel_max;
    *pixel_g=*pixel_g / pixel_max;
    *pixel_b=*pixel_b / pixel_max;
  }
  return(0);
}

int main(int argc, char **argv) {
  FILE *input_file;
  double two_pi=2.0 * M_PI;
  //double pi_over_2=M_PI / 2.0;
  double pi_over_180=M_PI / 180.0;
  double pi_over_360=M_PI / 360.0;
  
  // temp working variables
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
  int output_x;
  int output_y;
  double star_linear_intensity; // star linear intensity as viewed from camera
  int image_offset;
  double pixel_r;
  double pixel_g;
  double pixel_b;
  int raster_r;
  int raster_g;
  int raster_b;
  float star_linear_1pc_intensity;
  uint64_t color_temperature;
  double rgb_red[32768];
  double rgb_green[32768];
  double rgb_blue[32768];

  // input file record
  typedef struct {
    double icrs_x;
    double icrs_y;
    double icrs_z;
    uint64_t intensity_and_temperature;
  } star_record_t;
  star_record_t star_record;
  int star_record_size=sizeof(star_record_t);

  // user-supplied/derived camera parameters
  int min_parallax_quality;
  int camera_res_x;           // camera image resolution x-axis
  double camera_half_res_x;   // half of x-axis resolution
  int camera_res_y;           // camera image resolution y-axix
  double camera_half_res_y;   // half of y-axis resolution
  double camera_fov;          // user-supplied field of view angle for x-axis (deg)
  double camera_hfov;         // half of fov angle, to speed calculations (rad)
  double pixels_per_radian;   // pixels per radian for x and y axis
  double camera_saturation_mag;     // Mangnitude at white individual camera pixels saturate 
  double pixel_intensity_limit;         // maximum allowed flux value for a single pixel r,g, or b channel
  double camera_color_saturation;  // camera chroma saturation level
  double camera_wb_temp;      // camera white balance temperature
  double camera_gamma;        // gamma adjustment for camera output
  double camera_icrs_x;        // user-supplied camera position (parsec)
  double camera_icrs_y;        // user-supplied camera position (parsec)
  double camera_icrs_z;        // user-supplied camera position (parsec)
  double target_icrs_x;     // user-supplied camera aiming position (parsec)
  double target_icrs_y;     // user-supplied camera aiming position (parsec)
  double target_icrs_z;     // user-supplied camera aiming position (parsec)
  double target_icrs_ra;
  double target_icrs_ra_rad;
  double target_icrs_dec;
  double target_icrs_dec_rad;
  double target_icrs_r;
  double camera_rotation;     // rotation of camera relative to standardized "up" direction.   This standard up direction varies with the 3az direction of the camera target
  double camera_pan;          // pan left-right of camera after rotating to target
  double camera_tilt;
  double camera_3az_yz;       // rotation of camera (radians)
  double camera_3az_xy;       // pan of camera (radians)
  double camera_3az_xz;       // tilt of camera (radians)
  double render_distance;            // distance from selected point to star
  double render_distance_min;        // only render stars at least this far away from selected point
  double render_distance_max;        // only render stars within this distance from selected point
  int render_distance_selector;         //  min/max render distance from: 0=camera, 1=target

  // image buffers
  typedef struct {
    double r;
    double g;
    double b;
  } pixel_composition_t;
  typedef struct {
    char r;
    char g;
    char b;
  } pixel_output_t;
  pixel_composition_t *image_composition_buf;
  pixel_composition_t *image_composition_p;
  pixel_output_t *image_output_buf;
  pixel_output_t *image_output_p;

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

  // rendering options
  min_parallax_quality=10;

  // camera options
  camera_fov=360.0;
  camera_saturation_mag=8.2;
  //pixel_intensity_limit=5.0E-4;
  camera_wb_temp=4100.0;
  camera_color_saturation=4.0;
  camera_gamma=1.00;
  camera_pan=0.0;
  camera_tilt=0.0;

  // camera resolution
  // 4K 16:9
  //camera_res_x=3840;
  //camera_res_y=2160;
  // 4K 2:1
  camera_res_x=4096;
  camera_res_y=2048;
  // 8K 16:9
  //camera_res_x=7680;
  //camera_res_y=4320;
  // 8K 2:1
  //camera_res_x=8192;
  //camera_res_y=4096;
  // 16K
  //camera_res_x=15360;
  //camera_res_y=8640;
  // 16K 2:1
  //camera_res_x=16384;
  //camera_res_y=8192;
  // 20K 2:1
  //camera_res_x=20480;
  //camera_res_y=10240;
  // 20K 16:9
  //camera_res_x=19200;
  //camera_res_y=10800;
  // 20K 1:1
  //camera_res_x=20480;
  //camera_res_y=20480;
  // 32K 2:1
  //camera_res_x=32768;
  //camera_res_y=16384;
  //
  //camera_res_x=4096;
  //camera_res_y=4096;

  // camera position in icrs x,y,z
  //camera_icrs_x=-0.974923292E6;
  //camera_icrs_y=-0.222541177E6;
  //camera_icrs_z=+0.45598507E6;

  // camera position in icrs x,y,z
  camera_icrs_x=+0.0E0;
  camera_icrs_y=+0.0E0;
  camera_icrs_z=+0.0E0;

  // camera target in icrs x,y,z
  //target_icrs_x=+0.0E0;
  //target_icrs_y=+0.0E0;
  //target_icrs_z=+0.0E0;
  //camera_rotation=-33.0;

  // camera target icrs ra,dec,r
  //target_icrs_dec=0.0E0;
  //target_icrs_ra=0.0E0;
  //target_icrs_r=1.0E0;
  //camera_rotation=0.0;

  // camera target galactic south pole in icrs x,y,z
  //target_icrs_x=+0.974923292E0;
  //target_icrs_y=+0.222541177E0;
  //target_icrs_z=-0.45598507E0;
  //camera_rotation=-33.0;

  // camera target galactic center
  target_icrs_dec=-29.0078106;
  target_icrs_ra=266.4168371;
  target_icrs_r=8178.0;
  camera_rotation=-60.2;
  //camera_pan=21.5;
  //camera_tilt=0.0;

  // min/max render distance
  render_distance_min=0.0;
  render_distance_max=1.0E99;
  //
  // 1.099055223
  //render_distance_min=1.099050E6;
  //render_distance_max=1.099060E6; 
  //
  //render_distance_min=2.0E3;
  //render_distance_max=1.0E99;

  render_distance_selector=0; // 0 = camera, 1 = target

/*
  printf("init, initializing image buffers\n");
  fflush(stdout);
*/

  // allocate memory for image composition buffer (floating point rgb) and initialize
  image_composition_buf = (pixel_composition_t *)malloc(camera_res_x * camera_res_y * sizeof(pixel_composition_t));
  if (image_composition_buf == NULL) {
    printf("Error: could not allocate memory for image composition buffer\n");
    return(1);
  }
  image_composition_p=image_composition_buf;
  for (i=0; i < (camera_res_x * camera_res_y); i++) {
    image_composition_p->r=0.0;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
    image_composition_p++;
  }

  // allocate memory for image output buffer (8 bits per color rgb) and initialize
  image_output_buf = (pixel_output_t *)malloc(camera_res_x * camera_res_y * sizeof(pixel_output_t));
  if (image_output_buf == NULL) {
    printf("Error: could not allocate memory for image output buffer\n");
    return(1);
  }
  image_output_p=image_output_buf;
  for (i=0; i < (camera_res_x * camera_res_y); i++) {
    image_output_p->r=(char)0;
    image_output_p->g=(char)0;
    image_output_p->b=(char)0;
    image_output_p++;
  }

  // initialize RGB color lookup tables
  initRGBTables(camera_wb_temp, camera_color_saturation, rgb_red, rgb_green, rgb_blue);

/*
  printf("init, opening input file\n");
  fflush(stdout);
*/

  // attempt to open input file(s)
  input_file=fopen("galaxy.dat", "r");
  if (input_file == NULL) {
    printf("init, Error: could not open galaxy.dat\n");
    fflush(stdout);
    return(1);
  }

  // process user-supplied arguments
  camera_hfov=camera_fov * pi_over_360; // includes divide by 2
  camera_half_res_x=(double)camera_res_x / 2.0;
  camera_half_res_y=(double)camera_res_y / 2.0;
  pixels_per_radian=camera_half_res_x / camera_hfov;
  camera_3az_yz=camera_rotation * pi_over_180;
  camera_3az_xy=camera_pan * pi_over_180;
  camera_3az_xz=camera_tilt * pi_over_180;
  pixel_intensity_limit=pow(100.0, (-camera_saturation_mag / 5.0));

  // transform spherical icrs to euclidian icrs
  target_icrs_ra_rad=target_icrs_ra * pi_over_180;
  target_icrs_dec_rad=target_icrs_dec * pi_over_180;
  target_icrs_x=target_icrs_r * cos(target_icrs_dec_rad) * cos(target_icrs_ra_rad);
  target_icrs_y=target_icrs_r * cos(target_icrs_dec_rad) * sin(target_icrs_ra_rad);
  target_icrs_z=target_icrs_r * sin(target_icrs_dec_rad);

/*
  printf("debug, camera_res_x: %d, camera_res_y: %d, camera_fov: %.4e, camera hfov: %.4e, pixels_per_radian: %.4e\n", camera_res_x, camera_res_y, camera_fov, camera_hfov, pixels_per_radian);
  fflush(stdout);
*/

  // convert camera target x,y,z to 3az coordinates as seen from camera
  target_x=target_icrs_x - camera_icrs_x;
  target_y=target_icrs_y - camera_icrs_y;
  target_z=target_icrs_z - camera_icrs_z;
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

  // read star record from input file
  fread(&star_record, star_record_size, 1, input_file);

  // read and process each line of input file and render stars to image
  while ((feof(input_file) == 0) && (ferror(input_file) == 0)) {

    // convert star x,y,z to 3az coordinates as seen by camera position
    star_x=star_record.icrs_x - camera_icrs_x;
    star_y=star_record.icrs_y - camera_icrs_y;
    star_z=star_record.icrs_z - camera_icrs_z;
    star_r=sqrt(pow(star_x, 2.0) + pow(star_y, 2.0) + pow(star_z, 2.0));

/*
    printf("debug, star_icrs_x: %.9f, star_icrs_y: %.9f, star_icrs_z: %.9f, star_x: %.9f, star_y: %.9f, star_z: %.9f, star_r: %.9f\n", star_record.icrs_x, star_record.icrs_y, star_record.icrs_z, star_x, star_y, star_z, star_r);
    fflush(stdout);
*/

    // only continue if star distance is within min-max range from selected point (0=camera, 1=target)
    if (render_distance_selector == 0) { // selected point is camera
      render_distance=star_r; // star distance from camera
    } else { // selected point is target
      render_distance=sqrt(pow((star_record.icrs_x - target_icrs_x), 2.0) + pow((star_record.icrs_y - target_icrs_y), 2.0) + pow((star_record.icrs_z - target_icrs_z), 2.0)); // important, use un-rotated coordinates
    }
    if ((render_distance >= render_distance_min) && (render_distance <= render_distance_max)) {
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

/*
      printf("debug, initial, star_3az_xy: %.9f, star_3az_xz: %.9f, star_3az_yz: %.9f\n", star_3az_xy, star_3az_xz, star_3az_yz);
      fflush(stdout);
*/
    
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

/*
      printf("debug, rotated once, star_3az_xy: %.9f, star_3az_xz: %.9f, star_3az_yz: %.9f\n", star_3az_xy, star_3az_xz, star_3az_yz);
      fflush(stdout);
*/

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

/*
      printf("debug, rotated twice, star_3az_xy: %.9f, star_3az_xz: %.9f, star_3az_yz: %.9f\n", star_3az_xy, star_3az_xz, star_3az_yz);
      fflush(stdout);
*/

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

      if (camera_pan != 0) {
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

      if (camera_tilt != 0.0) {
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
/*
      printf("debug, camera rotated on yz axis, star_3az_xy: %.9f, star_3az_xz: %.9f, star_3az_yz: %.9f\n", star_3az_xy, star_3az_xz, star_3az_yz);
      fflush(stdout);
*/

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

        // cylindrical projection
        output_x=(int)((-output_az * pixels_per_radian) + camera_half_res_x + 0.5);
        output_y=(int)((-output_el * pixels_per_radian) + camera_half_res_y + 0.5);

        if ((output_x >= 0) && (output_x < camera_res_x) && (output_y >= 0) && (output_y < camera_res_y)) {
          image_offset=(camera_res_x * output_y) + output_x;
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


/*
  // optionally draw cross hairs
  for (i=(camera_half_res_x - (camera_res_y * 0.02)); i < (camera_half_res_x - (camera_res_y * 0.005)); i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * camera_half_res_y) + i;
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(camera_half_res_x + (camera_res_y * 0.005)); i < (camera_half_res_x + (camera_res_y * 0.02)); i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * camera_half_res_y) + i;
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(camera_half_res_y - (camera_res_y * 0.02)); i < (camera_half_res_y - (camera_res_y * 0.005)); i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * i) + (int)camera_half_res_x;
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(camera_half_res_y + (camera_res_y * 0.005)); i < (camera_half_res_y + (camera_res_y * 0.02)); i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * i) + (int)camera_half_res_x;
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
*/

/*
  // optionally select raster lines
  for (i=0; i < camera_res_x; i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * (camera_res_y * 0.25)) + i;
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < camera_res_x; i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * camera_half_res_y) + i;
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < camera_res_x; i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * (camera_res_y * 0.75)) + i;
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < camera_res_y; i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * i) + (int)(camera_res_x * 0.25);
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < camera_res_y; i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * i) + (int)(camera_half_res_x);
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < camera_res_y; i++) {
    image_composition_p=image_composition_buf + (int)(camera_res_x * i) + (int)(camera_res_x * 0.75);
    image_composition_p->r=(pixel_intensity_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
*/

#ifdef OUTPUT_PPM
  // output PPM header
  printf("P3\n%d %d\n255\n", camera_res_x, camera_res_y);
  image_composition_p=image_composition_buf;
  for (i=0; i < (camera_res_x * camera_res_y); i++) {
    // convert flux values to output range ~0-1.0 with camera sensitivity reference level = 1.0
    pixel_r=image_composition_p->r / pixel_intensity_limit;
    pixel_g=image_composition_p->g / pixel_intensity_limit;
    pixel_b=image_composition_p->b / pixel_intensity_limit;

    // apply camera gamma setting
    pixel_r=pow(pixel_r, camera_gamma);
    pixel_g=pow(pixel_g, camera_gamma);
    pixel_b=pow(pixel_b, camera_gamma);

    limitIntensityPreserveColor(&pixel_r, &pixel_g, &pixel_b);

    // apply sRGB gamma
    if (pixel_r <= 0.0031308) {
      pixel_r=pixel_r * 12.92;
    } else {
      pixel_r=(1.055 * pow(pixel_r, (1.0 / 2.4)) - 0.055);
    }
    if (pixel_g <= 0.0031308) {
      pixel_g=pixel_g * 12.92;
    } else {
      pixel_g=(1.055 * pow(pixel_g, (1.0 / 2.4)) - 0.055);
    }
    if (pixel_b <= 0.0031308) {
      pixel_b=pixel_b * 12.92;
    } else {
      pixel_b=(1.055 * pow(pixel_b, (1.0 / 2.4)) - 0.055);
    }

    // convert r,g,b to 8 bit values with safety clippings
    raster_r=(int)(pixel_r * 255.0);
    if (raster_r > 255) {
      raster_r=255;
    } else if (raster_r < 0) {
      raster_r=0;
    }
    raster_g=(int)(pixel_g * 255.0);
    if (raster_g > 255) {
      raster_g=255;
    } else if (raster_g < 0) {
      raster_g=0;
    }
    raster_b=(int)(pixel_b * 255.0);
    if (raster_b > 255) {
      raster_b=255;
    } else if (raster_b < 0) {
      raster_b=0;
    }
    printf("%d %d %d\n", raster_r, raster_g, raster_b);
    image_composition_p++;
  }
#endif

  // clean up
  fclose(input_file);
  free(image_composition_buf);
  free(image_output_buf);

  return(0);
}

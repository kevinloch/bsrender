//
// Billion Star 3D Rendering Engine
// Kevin M. Loch
//
// pre-processor for external/manual star databases
// This program creates binary data files for use by the rendering engine
//
// The format of the binary data files is given by the star_record data type
//
//   typedef struct {
//    double icrs_x;
//    double icrs_y;
//    double icrs_z;
//    uint64_t intensity_and_temperature;
//  } star_record_t;
//
// where intensity_and_temperature has a single precision value for the the normalized (to vega at 1 parsec) flux in the 32 MSB,
// and the color temperature of the star as an unsigned integer between 0..32767 in the 32 LSB.
// For speed in the rendering engine no byte-order checking/setting is done so the data files must be generated on a compatible platform for the rendering engine.
//

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
  FILE *input_file;
  FILE *output_file;
  char *input_line_p;
  char input_line[256];
  char *field_start;
  char *field_end;
  size_t field_length;
  char tmpstr[32];
  long long input_count;
  long long output_count;
  float linear_1pc_intensity;
  uint64_t color_temperature;
  star_record_t star_record;
  int star_record_size=sizeof(star_record_t);

  // fields imported from external data csv
  double ra;
  double dec;
  double distance;
  double temperature;
  double apparent_magnitude;

  // temp working vars
  double ra_rad;
  double dec_rad;
  double linear_intensity; // intensity relative to vega

  //
  // print version
  //
  printf("mkexternal version %s\n", BSR_VERSION);

  //
  // attempt to open input file
  //
  printf("init, Opening input file external.csv\n");
  input_file=fopen("external.csv", "rb");
  if (input_file == NULL) {
    printf("init, Error: could not open external.csv\n");
    fflush(stdout);
    return(1);
  }

  //
  // attempt to open ouptut file
  //
  printf("init, Opening output file galaxy-external.dat\n");
  output_file=fopen("galaxy-external.dat", "wb");
  if (output_file == NULL) {
    printf("init, Error: could not open galaxy-external.dat for writing\n");
    fflush(stdout);
    return(1);
  }

  //
  // read and process each line of input file
  //
  input_count=0;
  output_count=0;
  input_line_p=fgets(input_line, 256, input_file);
  while (input_line_p != NULL) {

/*
    printf("input_line: %s", input_line_p);
    fflush(stdout);
*/
    //
    // ra,dec,distance,temperature,apparent_magnitude,common_name,notes
    //
    if (input_line[0] != 'r') { // skip csv header line
      input_count++;
      field_start=input_line;

      //
      // ra
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      ra=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // dec
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      dec=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // distance
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      distance=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // temperature
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      temperature=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // apparent_magnitude
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      apparent_magnitude=strtod(tmpstr, NULL);

      //
      // transform spherical icrs to euclidian icrs
      //
      ra_rad=ra * M_PI / 180.0;
      dec_rad=dec * M_PI / 180.0;
      star_record.icrs_x=distance * cos(dec_rad) * cos(ra_rad);
      star_record.icrs_y=distance * cos(dec_rad) * sin(ra_rad);
      star_record.icrs_z=distance * sin(dec_rad);

      //
      // convert apparent_magnitude to intensity at 1pc
      //
      linear_intensity=pow(100.0, (-apparent_magnitude / 5.0));
      if (distance == 0.0) {
        // special handling for the Sun
        linear_1pc_intensity=(float)(linear_intensity * 2.3504E-11);
      } else {
        linear_1pc_intensity=(float)(linear_intensity * pow(distance, 2.0));
      }

      //
      // combine intensity and temperature into one 64 bit value
      //
      color_temperature=(int)(temperature + 0.5);
      if (color_temperature > 32767) {
        color_temperature=32767;
      }
      star_record.intensity_and_temperature=0;
      star_record.intensity_and_temperature=(uint64_t)*((uint32_t*)&linear_1pc_intensity);
      star_record.intensity_and_temperature <<= 32;
      star_record.intensity_and_temperature |= color_temperature;

printf("ra_rad: %.6e, dec_rad: %.6e, distance: %.6e, linear_intensity: %.6e, linear_1pc_intensity: %.6e, temperature: %ld, intensity_and_temperature: %16lx\n", ra_rad, dec_rad, distance, linear_intensity, linear_1pc_intensity, color_temperature, star_record.intensity_and_temperature);
fflush(stdout);

      //
      // output transformed fields to output dat file
      //
      fwrite(&star_record, star_record_size, 1, output_file);
      output_count++;
    } // end ignore csv header lines

    //
    // periodic status
    //
    if ((input_count > 0) && ((input_count % 1000000) == 0)) {
      printf("status, input records: %9lld, galaxy-external.dat: %8lld\n", input_count, output_count);
    }

    input_line_p=fgets(input_line, 256, input_file);
  }  // end while read file
      
  //
  // print final status 
  //
  printf("status, input records: %9lld, galaxy-external.dat: %8lld\n", input_count, output_count);

  //
  // clean up
  //
  fclose(input_file);
  fclose(output_file);
  return(0);
}

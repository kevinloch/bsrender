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

//
// pre-processor for external/manual star databases
// This program creates binary data files for use by the rendering engine
//

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

void printUsage() {
  printf("mkexternal version %s\n", BSR_VERSION);
  printf("\n\
NAME\n\
     mkexternal -- create binary data file for use with bsrender\n\
\n\
SYNOPSIS\n\
     mkexternal [-l] [-g] [-h]\n\
 \n\
OPTIONS:\n\
\n\
     -l\n\
          Force output to little-endian format (default is to match this platform)\n\
\n\
     -g\n\
          Force output big-endian format (default is to match this platform)\n\
\n\
     -h\n\
          Show help\n\
\n\
DESCRIPTON\n\
 mkgalaxy processes extracted fields from ESA's Gaia DR3 dataset for use with bsrender. Uses the output from 'gaia-dr3-extract.sh' in the bsrender package\n\
 \n");
}

int setDefaults(mkg_config_t *mkg_config) {
  int little_endian;

  little_endian=littleEndianTest();
  if (little_endian == 1) {
    mkg_config->output_little_endian=1;
  } else {
    mkg_config->output_little_endian=0;
  }

  return(0);
}

int processCmdArgs(mkg_config_t *mkg_config, int argc, char **argv) {
  int i;

  if (argc == 1) {
    return(0);
  } else {
    for (i=1; i <= (argc - 1); i++) {
      if (argv[i][1] == 'l') {
        // force output to littl-endian
        mkg_config->output_little_endian=1;
      } else if (argv[i][1] == 'g') {
        // force output to big-endian
        mkg_config->output_little_endian=0;
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
  int i;
  FILE *input_file;
  char file_name[256];
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
  float linear_1pc_intensity_undimmed;
  uint64_t color_temperature;
  uint64_t color_temperature_unreddened;
  mkg_config_t mkg_config;
  double icrs_x;
  double icrs_y;
  double icrs_z;
  char star_record[BSR_STAR_RECORD_SIZE];
  size_t star_record_size=(size_t)BSR_STAR_RECORD_SIZE;
  char file_header[BSR_FILE_HEADER_SIZE];
  size_t file_header_size;
  int little_endian;
  int same_endian; // 1 == output endianness is same as arch

  // fields imported from external data csv
  double ra;
  double dec;
  double distance;
  double apparent_magnitude;
  double undimmed_magnitude;
  double apparent_temperature;
  double unreddened_temperature;

  // temp working vars
  double ra_rad;
  double dec_rad;
  double linear_intensity; // intensity relative to vega
  double linear_intensity_undimmed;
  uint64_t tmp64;
  uint32_t tmp32;
  uint16_t tmp16;

  //
  // set default options
  //
  setDefaults(&mkg_config);

  //
  // proces command line options
  processCmdArgs(&mkg_config, argc, argv);

  //
  // print version and options
  //
  printf("mkgalaxy version %s\n", BSR_VERSION);
  if (mkg_config.output_little_endian == 1) {
    printf("Output data files will be in little-endian format\n");
  } else {
    printf("Output data files will be in big-endian format\n");
  }

  //
  // attempt to open input file
  //
  printf("Opening input file external.csv\n");
  input_file=fopen("external.csv", "rb");
  if (input_file == NULL) {
    printf("Error: could not open external.csv\n");
    fflush(stdout);
    return(1);
  }

  //
  // attempt to open ouptut file
  //
  if (mkg_config.output_little_endian == 1) {
    sprintf(file_name, "%s-%s.%s", BSR_EXTERNAL_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
  } else {
    sprintf(file_name, "%s-%s.%s", BSR_EXTERNAL_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
  }
  printf("Opening output file %s\n", file_name);
  output_file=fopen(file_name, "wb");
  if (output_file == NULL) {
    printf("Error: could not open %s for writing\n", file_name);
    fflush(stdout);
    return(1);
  }

  //
  // check endianness
  //
  little_endian=littleEndianTest();
  if ((mkg_config.output_little_endian ^ little_endian) == 0) {
    same_endian=1;
  } else {
    same_endian=0;
  }

  //
  // write file header
  //
  if (mkg_config.output_little_endian == 1) {
    snprintf(file_header, BSR_FILE_HEADER_SIZE, "%s, mkexternal version: %s\n", BSR_MAGIC_NUMBER_LE, BSR_VERSION);
  } else {
    snprintf(file_header, BSR_FILE_HEADER_SIZE, "%s, mkexternal version: %s\n", BSR_MAGIC_NUMBER_BE, BSR_VERSION);
  } 
  // pad the rest of file_header with zeros
  file_header_size=strnlen(file_header, (BSR_FILE_HEADER_SIZE - 1));
  for (i=(int)file_header_size; i < (BSR_FILE_HEADER_SIZE - 1); i++) {
    file_header[i]=0;
  }
  printf("Writing file headers\n");
  fflush(stdout);
  fwrite(file_header, BSR_FILE_HEADER_SIZE, 1, output_file);

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
    // ra,dec,distance,apparent_magnitude,undimmed_magnitude,apparent_temperature,unreddened_temperature,common_name,notes
    //
    if (input_line[0] != 'r') { // skip csv header line
      input_count++;
      field_start=input_line;

      //
      // ra
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      ra=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // dec
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      dec=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // distance
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      distance=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // apparent_magnitude
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      apparent_magnitude=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // undimmed_magnitude
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      undimmed_magnitude=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // apparent_temperature
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      apparent_temperature=strtod(tmpstr, NULL);
      field_start=(field_end+1);

      //
      // unreddened_temperature
      //
      field_end=strchr(field_start, ',');
      field_length=(field_end - field_start);
      if (field_length > 31) {
        field_length=31;
      }
      strncpy(tmpstr, field_start, field_length);
      tmpstr[field_length]=0;
      unreddened_temperature=strtod(tmpstr, NULL);

      //
      // transform spherical icrs to euclidian icrs
      //
      ra_rad=ra * M_PI / 180.0;
      dec_rad=dec * M_PI / 180.0;
      icrs_x=distance * cos(dec_rad) * cos(ra_rad);
      icrs_y=distance * cos(dec_rad) * sin(ra_rad);
      icrs_z=distance * sin(dec_rad);

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
      // convert undimmed_magnitude to intensity at 1pc
      //
      linear_intensity_undimmed=pow(100.0, (-undimmed_magnitude / 5.0));
      if (distance == 0.0) {
        // special handling for the Sun
        linear_1pc_intensity_undimmed=(float)(linear_intensity_undimmed * 2.3504E-11);
      } else {
        linear_1pc_intensity_undimmed=(float)(linear_intensity_undimmed * pow(distance, 2.0));
      }

      //
      // convert apparent_temperature to int and verify range
      //
      color_temperature=(uint64_t)(apparent_temperature + 0.5);
      if (color_temperature < 500ul) {
        color_temperature=500ul;
      } else if (color_temperature > 32767ul) {
        color_temperature=32767ul;
      }

      //
      // convert unreddened_temperature to int and verify range
      //
      color_temperature_unreddened=(uint64_t)(unreddened_temperature + 0.5);
      if (color_temperature_unreddened < 500ul) {
        color_temperature_unreddened=500ul;
      } else if (color_temperature_unreddened > 32767ul) {
        color_temperature_unreddened=32767ul;
      }

      // As of v1.0, bsrender data files have a fixed-length 256-bit ascii header (including the file identifier in the first 11 bytes),
      // followed by a variable number of 33-byte star records. The ascii header can be viewed with 'head -1 <filename>'.
      //
      // Each star record includes a 64-bit unsigned integer for Gaia DR3 'source_id', three 40-bit truncaed doubles for x,y,z,
      // a 24-bit truncated float for linear_1pc_intensity, a 24-bit truncated float for linear_1pc_intensity_undimmed,
      // a 16-bit unsigned int for color_temperature, and a 16-bit unsigned int for color_temperature_unreddened.
      // these are packed into a 33 byte star record with each field encoded in the selected byte order.
      //
      // +---------------+---------+---------+---------+-----+-----+---+---+
      // |   source_id   |    x    |    y    |    z    | li  |li-u | c |c-u|
      // +---------------+---------+---------+---------+-----+-----+---+---+
      // |      8        |    5    |    5    |    5    |  3  |  3  | 2 | 2 |
      //                               bytes
      //
      // mkgalaxy and mkexternal have options to generate either little-endian or big-endian files but the bye-order of the file
      // must match the architecture it is used on with bsrender. This is to avoid unnessary operations in the performance
      // critical inner-loop of processStars() which iterates over potentially billions of star records. The filenames contain '-le'
      // or '-be' to indicate byte order. Byte order is also indicated with the file identifier in the first 11 bytes of the header:
      // BSRENDER_LE for little-endian and BSRENDER_BE for big-endian.
      //

      //
      // pack star record fields into 33-byte star_record
      //

      // set source_id to 0
      star_record[0]=0;
      star_record[1]=0;
      star_record[2]=0;
      star_record[3]=0;
      star_record[4]=0;
      star_record[5]=0;
      star_record[6]=0;
      star_record[7]=0;

      if (same_endian == 1) {
        //
        // output endianness is same as this arch
        //

        // icrs_x
        tmp64=0;
        tmp64=*(uint64_t *)&icrs_x;
        if (little_endian == 1) {
          tmp64 >>= 24; // skip 24 lsb if source is little-endian
        }
        star_record[8]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[9]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[10]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[11]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[12]=*(char *)&tmp64;

        // icrs_y
        tmp64=0;
        tmp64=*(uint64_t *)&icrs_y;
        if (little_endian == 1) {
          tmp64 >>= 24; // skip 24 lsb if source is little-endian
        }
        star_record[13]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[14]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[15]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[16]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[17]=*(char *)&tmp64;

        // icrs_z
        tmp64=0;
        tmp64=*(uint64_t *)&icrs_z;
        if (little_endian == 1) {
          tmp64 >>= 24; // skip 24 lsb if source is little-endian
        }
        star_record[18]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[19]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[20]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[21]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[22]=*(char *)&tmp64;

        // linear_1pc_intensity
        tmp32=0;
        tmp32=*(uint32_t *)&linear_1pc_intensity;
        if (little_endian == 1) {
          tmp32 >>= 8; // skip 8 lsb if source is little-endian
        }
        star_record[23]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[24]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[25]=*(char *)&tmp32;

        // linear_1pc_intensity_undimmed
        tmp32=0;
        tmp32=*(uint32_t *)&linear_1pc_intensity_undimmed;
        if (little_endian == 1) {
          tmp32 >>= 8; // skip 8 lsb if source is little-endian
        }
        star_record[26]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[27]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[28]=*(char *)&tmp32;

        // color_temperature
        tmp16=0;
        tmp16=*(uint16_t *)&color_temperature;
        star_record[29]=*(char *)&tmp16;
        tmp16 >>= 8;
        star_record[30]=*(char *)&tmp16;

        // color_temperature_unreddened
        tmp16=0;
        tmp16=*(uint16_t *)&color_temperature_unreddened;
        star_record[31]=*(char *)&tmp16;
        tmp16 >>= 8;
        star_record[32]=*(char *)&tmp16;
      } else {
        // icrs_x
        tmp64=0;
        tmp64=*(uint64_t *)&icrs_x;
        if (little_endian == 1) {
          tmp64 >>= 24; // skip 24 lsb if source is little-endian
        }
        star_record[12]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[11]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[10]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[9]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[8]=*(char *)&tmp64;

        // icrs_y
        tmp64=0;
        tmp64=*(uint64_t *)&icrs_y;
        if (little_endian == 1) {
          tmp64 >>= 24; // skip 24 lsb if source is little-endian
        }
        star_record[17]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[16]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[15]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[14]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[13]=*(char *)&tmp64;

        // icrs_z
        tmp64=0;
        tmp64=*(uint64_t *)&icrs_z;
        if (little_endian == 1) {
          tmp64 >>= 24; // skip 24 lsb if source is little-endian
        }
        star_record[22]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[21]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[20]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[19]=*(char *)&tmp64;
        tmp64 >>= 8;
        star_record[18]=*(char *)&tmp64;

        // linear_1pc_intensity
        tmp32=0;
        tmp32=*(uint32_t *)&linear_1pc_intensity;
        if (little_endian == 1) {
          tmp32 >>= 8; // skip 8 lsb if source is little-endian
        }
        star_record[25]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[24]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[23]=*(char *)&tmp32;

        // linear_1pc_intensity_undimmed
        tmp32=0;
        tmp32=*(uint32_t *)&linear_1pc_intensity_undimmed;
        if (little_endian == 1) {
          tmp32 >>= 8; // skip 8 lsb if source is little-endian
        }
        star_record[28]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[27]=*(char *)&tmp32;
        tmp32 >>= 8;
        star_record[26]=*(char *)&tmp32;

        // color_temperature
        tmp16=0;
        tmp16=*(uint16_t *)&color_temperature;
        star_record[30]=*(char *)&tmp16;
        tmp16 >>= 8;
        star_record[29]=*(char *)&tmp16;

        // color_temperature_unreddened
        tmp16=0;
        tmp16=*(uint16_t *)&color_temperature_unreddened;
        star_record[32]=*(char *)&tmp16;
        tmp16 >>= 8;
        star_record[31]=*(char *)&tmp16;
      }

      //
      // output star_record to output dat file
      //
      fwrite(&star_record, star_record_size, 1, output_file);
      output_count++;

    } // end ignore csv header lines

    //
    // periodic status
    //
    if ((input_count > 0) && ((input_count % 1000000) == 0)) {
      printf("Input records: %9lld, %s: %8lld\n", input_count, file_name, output_count);
    }

    input_line_p=fgets(input_line, 256, input_file);
  }  // end while read file
      
  //
  // print final status 
  //
  printf("Input records: %9lld, %s: %8lld\n", input_count, file_name, output_count);

  //
  // clean up
  //
  fclose(input_file);
  fclose(output_file);
  return(0);
}

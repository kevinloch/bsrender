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
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>

int openInputFile(bsr_config_t *bsr_config, char *file_path, input_file_t *input_file, int little_endian) {
  int mmap_protection;
  int mmap_visibility;

  input_file->fd=open(file_path, O_RDONLY);
  if (input_file->fd < 0) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not open %s\n", file_path);
      fflush(stdout);
    }
    exit(1);
  }

  fstat(input_file->fd, &input_file->sb);
  if (input_file->sb.st_size == 0) {
    // mmap will not map zero length files but we don't want that to abort the entire program
    // processStars() will not try to read anything from this file so input_file->buf is irrelevant
    input_file->buf=NULL;
    return(0);
  }

  mmap_protection=PROT_READ;
  mmap_visibility=MAP_SHARED;
  input_file->buf=mmap(NULL, input_file->sb.st_size, mmap_protection, mmap_visibility, input_file->fd, 0);
  if (input_file->buf == MAP_FAILED) {
    if (bsr_config->cgi_mode != 1) {
      printf("Error: could not mmap file %s, errno: %d\n", file_path, errno);
      fflush(stdout);
    }
    exit(1);
  }

  //
  // verify file has correct endianness signature for this platform
  //
  if (little_endian == 1) {
    if (strncmp(input_file->buf, BSR_MAGIC_NUMBER_LE, 11) != 0) {
      if (bsr_config->cgi_mode != 1) {
        printf("Error: input file %s is not a bsrender data file or is not in little endian format as this platform requires\n", file_path);
      }
      exit(1);
    }
  } else {
    if (strncmp(input_file->buf, BSR_MAGIC_NUMBER_BE, 11) != 0) {
      if (bsr_config->cgi_mode != 1) {
        printf("Error: input file %s is not a bsrender data files or is not in big endian format as this platform requires\n", file_path);
      }
      exit(1);
    }
  }

  return(0);
}

int openInputFiles(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  char file_path[1024];
  input_file_t input_file;
  int little_endian;

  //
  // check endianness
  //
  little_endian=bsr_state->little_endian;

  //
  // open each input file
  //
  if (bsr_config->external_db_enable == 1) {
    if (little_endian == 1) {
      sprintf(file_path, "%s/%s-%s.%s", bsr_config->data_file_directory, BSR_EXTERNAL_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
    } else {
      sprintf(file_path, "%s/%s-%s.%s", bsr_config->data_file_directory, BSR_EXTERNAL_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
    }
    openInputFile(bsr_config, file_path, &input_file, little_endian);
    bsr_state->input_file_external.buf=input_file.buf;
    bsr_state->input_file_external.fd=input_file.fd;
    bsr_state->input_file_external.buf_size=input_file.sb.st_size;
  } // end if enable_external_db
  if (bsr_config->Gaia_db_enable == 1) {
    if (little_endian == 1) {
      sprintf(file_path, "%s/%s-pq100-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
    } else {
      sprintf(file_path, "%s/%s-pq100-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
    }
    openInputFile(bsr_config, file_path, &input_file, little_endian);
    bsr_state->input_file_pq100.buf=input_file.buf;
    bsr_state->input_file_pq100.fd=input_file.fd;
    bsr_state->input_file_pq100.buf_size=input_file.sb.st_size;
    if (bsr_config->Gaia_min_parallax_quality < 100) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq050-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq050-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq050.buf=input_file.buf;
      bsr_state->input_file_pq050.fd=input_file.fd;
      bsr_state->input_file_pq050.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 50) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq030-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq030-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq030.buf=input_file.buf;
      bsr_state->input_file_pq030.fd=input_file.fd;
      bsr_state->input_file_pq030.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 30) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq020-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq020-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq020.buf=input_file.buf;
      bsr_state->input_file_pq020.fd=input_file.fd;
      bsr_state->input_file_pq020.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 20) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq010-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq010-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq010.buf=input_file.buf;
      bsr_state->input_file_pq010.fd=input_file.fd;
      bsr_state->input_file_pq010.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 10) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq005-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq005-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq005.buf=input_file.buf;
      bsr_state->input_file_pq005.fd=input_file.fd;
      bsr_state->input_file_pq005.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 05) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq003-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq003-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq003.buf=input_file.buf;
      bsr_state->input_file_pq003.fd=input_file.fd;
      bsr_state->input_file_pq003.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 03) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq002-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq002-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq002.buf=input_file.buf;
      bsr_state->input_file_pq002.fd=input_file.fd;
      bsr_state->input_file_pq002.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 02) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq001-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq001-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq001.buf=input_file.buf;
      bsr_state->input_file_pq001.fd=input_file.fd;
      bsr_state->input_file_pq001.buf_size=input_file.sb.st_size;
    }
    if (bsr_config->Gaia_min_parallax_quality < 01) {
      if (little_endian == 1) {
        sprintf(file_path, "%s/%s-pq000-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_LE_SUFFIX, BSR_EXTENSION);
      } else {
        sprintf(file_path, "%s/%s-pq000-%s.%s", bsr_config->data_file_directory, BSR_GDR3_PREFIX, BSR_BE_SUFFIX, BSR_EXTENSION);
      }
      openInputFile(bsr_config, file_path, &input_file, little_endian);
      bsr_state->input_file_pq000.buf=input_file.buf;
      bsr_state->input_file_pq000.fd=input_file.fd;
      bsr_state->input_file_pq000.buf_size=input_file.sb.st_size;
    }
  } // end if enable Gaia_db

  return(0);
}

int closeInputFiles(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {

  if (bsr_config->external_db_enable == 1) {
    munmap(bsr_state->input_file_external.buf, bsr_state->input_file_external.buf_size);
    close(bsr_state->input_file_external.fd);
  } // end if external_db_enable
  if (bsr_config->Gaia_db_enable == 1) {
    munmap(bsr_state->input_file_pq100.buf, bsr_state->input_file_pq100.buf_size);
    close(bsr_state->input_file_pq100.fd);
    if (bsr_config->Gaia_min_parallax_quality < 100) {
      munmap(bsr_state->input_file_pq050.buf, bsr_state->input_file_pq050.buf_size);
      close(bsr_state->input_file_pq050.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 50) {
      munmap(bsr_state->input_file_pq030.buf, bsr_state->input_file_pq030.buf_size);
      close(bsr_state->input_file_pq030.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 30) {
      munmap(bsr_state->input_file_pq020.buf, bsr_state->input_file_pq020.buf_size);
      close(bsr_state->input_file_pq020.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 20) {
      munmap(bsr_state->input_file_pq010.buf, bsr_state->input_file_pq010.buf_size);
      close(bsr_state->input_file_pq010.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 10) {
      munmap(bsr_state->input_file_pq005.buf, bsr_state->input_file_pq005.buf_size);
      close(bsr_state->input_file_pq005.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 05) {
      munmap(bsr_state->input_file_pq003.buf, bsr_state->input_file_pq003.buf_size);
      close(bsr_state->input_file_pq003.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 03) {
      munmap(bsr_state->input_file_pq002.buf, bsr_state->input_file_pq002.buf_size);
      close(bsr_state->input_file_pq002.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 02) {
      munmap(bsr_state->input_file_pq001.buf, bsr_state->input_file_pq001.buf_size);
      close(bsr_state->input_file_pq001.fd);
    }
    if (bsr_config->Gaia_min_parallax_quality < 01) {
      munmap(bsr_state->input_file_pq000.buf, bsr_state->input_file_pq000.buf_size);
      close(bsr_state->input_file_pq000.fd);
    } 
  } // end if enable Gaia_db

  return(0);
}


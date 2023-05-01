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
#include "bsr-config.h"

int printCGIHeader(bsr_config_t *bsr_config) {
  if (bsr_config->image_format == 0) {
    printf("Content-type: image/png\n");
    printf("Content-Disposition: attachment; filename=\"galaxy.png\"\n");
  } else if (bsr_config->image_format == 1) {
    printf("Content-type: image/x-exr\n");
    printf("Content-Disposition: attachment; filename=\"galaxy.exr\"\n");
  }
//  printf("Expires: 0\n");
//  printf("Cache-control: no-cache\n");
//  printf("Pragma: no-cache\n");
  printf("\n");
  fflush(stdout);
  return(0);
}

int sanitizeQueryString(char *query_string_2048) {
  int i;
  int tmpstr_len;
  char tmpstr[2048];
  char tmpchar;
  int sanitized_len;
  char hextmp[3];

  //
  // make a working copy of query_string_2048
  //
  strncpy(tmpstr, query_string_2048, 2047);
  tmpstr[2047]=0;
  tmpstr_len=strlen(tmpstr);

  //
  // validate each character of query_string_2048 (tmpstr)
  //
  sanitized_len=0;
  for (i=0; i < tmpstr_len; i++) {
    //
    // check for and convert hex encoded characters (common in CGI QUERY_STRING)
    //
    if ((tmpstr[i] == 37) && (i <= (tmpstr_len - 4))) {
      hextmp[0]=tmpstr[i+1];
      hextmp[1]=tmpstr[i+2];
      hextmp[2]=0;
      tmpchar=strtol(hextmp, 0, 16);
      i+=2;
    } else {
      tmpchar=tmpstr[i];
    }

    //
    // check if tmpchar is a valid symbol and copy back to query_string_2048
    //
    if (strchr("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.-+&=_", tmpchar) == NULL) {
      query_string_2048[sanitized_len]=32; // convert invalid character to space
    } else {
      query_string_2048[sanitized_len]=tmpchar;
    }
    sanitized_len++;
  }
  query_string_2048[sanitized_len]=0;

  return(0);
}

int getCGIOptions(bsr_config_t *bsr_config) {
  char query_string_2048[2048];

  if (bsr_config->QUERY_STRING_p == NULL) {
    return(1);
  } else {
    strncpy(query_string_2048, bsr_config->QUERY_STRING_p, 2047);
    query_string_2048[2047]=0;
    sanitizeQueryString(query_string_2048);
    loadConfigFromQueryString(bsr_config, query_string_2048);
  }

  return(0);
}

int enforceCGILimits(bsr_config_t *bsr_config) {
  if (bsr_config->camera_res_x < 1) {
    bsr_config->camera_res_x=1;
  }
  if (bsr_config->camera_res_x > bsr_config->cgi_max_res_x) {
    bsr_config->camera_res_x=bsr_config->cgi_max_res_x;
  }
  if (bsr_config->camera_res_y < 1) {
    bsr_config->camera_res_y=1;
  }
  if (bsr_config->camera_res_y > bsr_config->cgi_max_res_y) {
    bsr_config->camera_res_y=bsr_config->cgi_max_res_y;
  }
  if (bsr_config->Gaia_min_parallax_quality < bsr_config->cgi_Gaia_min_parallax_quality) {
    bsr_config->Gaia_min_parallax_quality=bsr_config->cgi_Gaia_min_parallax_quality;
  }
  if (bsr_config->cgi_allow_Airy_disk == 0) {
    bsr_config->Airy_disk_enable=0;
  }
  if (bsr_config->Airy_disk_first_null < bsr_config->cgi_min_Airy_disk_first_null) {
    bsr_config->Airy_disk_first_null=bsr_config->cgi_min_Airy_disk_first_null;
  }
  if (bsr_config->Airy_disk_max_extent > bsr_config->cgi_max_Airy_disk_max_extent) {
    bsr_config->Airy_disk_max_extent=bsr_config->cgi_max_Airy_disk_max_extent;
  }
  if (bsr_config->Airy_disk_min_extent > bsr_config->cgi_max_Airy_disk_min_extent) {
    bsr_config->Airy_disk_min_extent=bsr_config->cgi_max_Airy_disk_min_extent;
  }
  if (bsr_config->cgi_allow_anti_alias == 0) {
    bsr_config->anti_alias_enable=0;
  }

  return(0);
}

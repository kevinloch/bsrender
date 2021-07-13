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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bsr-config.h"

int printCGIheader(){

  printf("Content-type: image/png\n");
//  printf("Expires: 0\n");
//  printf("Cache-control: no-cache\n");
//  printf("Pragma: no-cache\n");
  printf("\n");
  fflush(stdout);
  return(0);
}

int sanitizeQueryString(char *query_string) {
  int i;
  int tmpstr_len;
  char tmpstr[2048];
  char tmpchar;
  char *output_str;
  int output_len;
  char hextmp[3];

  strcpy(tmpstr, query_string);
  tmpstr[2047]=0;
  tmpstr_len=strlen(tmpstr);

  output_str=query_string;
  output_len=0;
  for (i=0; i < tmpstr_len; i++) {
    // check for and convert hex encoded characters
    if ((tmpstr[i] == 37) && (i <= (tmpstr_len - 4))) {
      hextmp[0]=tmpstr[i+1];
      hextmp[1]=tmpstr[i+2];
      hextmp[2]=0;
      tmpchar=strtol(hextmp, 0, 16);
      i+=2;
    } else {
      tmpchar=tmpstr[i];
    }

    // check if tmpchar is valid and copy to output_str
    if (strchr("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.-+&=_", tmpchar) == NULL) {
      output_str[output_len]=32; // convert invalid character to space
    } else {
      output_str[output_len]=tmpchar;
    }
    output_len++;
  }
  output_str[output_len]=0;

  return(0);
}

int getCGIOptions(bsr_config_t *bsr_config) {
  char *tmpstr;
  char query_string[2048];

  tmpstr=NULL;
  tmpstr=getenv("QUERY_STRING");
  if (tmpstr == NULL) {
    return(1);
  } else {
    strcpy(query_string, tmpstr);
    query_string[2047]=0;
    sanitizeQueryString(query_string);
    loadConfigFromQueryString(bsr_config, query_string);
  }

  return(0);
}

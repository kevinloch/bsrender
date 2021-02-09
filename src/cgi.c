#include "bsrender.h"
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

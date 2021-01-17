#include <stdio.h>

int printCGIheader(){
  printf("Content-type: image/png\n");
  printf("Expires: 0\n");
  printf("Cache-control: no-cache\n");
  printf("Pragma: no-cache\n");
  printf("\n");
  return(0);
}

int printCGIfooter() {
  printf("\n");
  return(0);
}

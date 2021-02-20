#include <stdio.h>
#include <math.h>

int main(int argc, char **cargv) {
  //
  // This program generates a table of the Bessel function of the first kind, order one J1(x), with 10 samples per 'x'.   This is used to make Airy disk pattern maps.
  // After the first null each Airy ring is approximately 3.2 'x' units wide for about 32 samples per ring. 
  //
  double bessel_J1;
  double x;
  double x_increment=0.1;
  double x_max=255.0;
  int i;
  int j;
  double theta;
  double theta_increment = 0.00000001; // smaller increments enable more accuracy for smaller values of bessel_J1 but take longer to compute
 
  j=0;
  for (x=0.0; x <= x_max; x+=x_increment) {
    //
    // numerical integration
    //
    bessel_J1=0.0;
    if (x != 0.0) {
      for (theta=0.0; theta <= M_PI; theta+=theta_increment) {
       bessel_J1+=(cos(theta - (x * sin(theta))) * theta_increment / M_PI);
      }
    }

    // output sample
    if (j == 16) {
      j=0;
      printf("\n");
    }
    printf(" %9.2e", bessel_J1);
    j++;
    if (x <= x_max) {
      printf(",");
    }
  } // end for x
  return(0);
}

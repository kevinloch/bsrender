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

#include <stdio.h>
#include <math.h>

int main(int argc, char **cargv) {
  //
  // This program generates a table of the Bessel function of the first kind, order one J1(x), with 10 samples per 'x'.   This is used to make Airy disk pattern maps.
  // After the first null each Airy ring is approximately 3.2 'x' units wide for about 32 samples per ring. 
  // This is included for documentation purposes, the output is already included in Bessel.h
  //
  double Bessel_J1;
  double x;
  double x_increment=0.1;
  double x_max=1023.0;
  int j;
  double theta;
  double theta_increment = 0.00000001; // smaller increments enable more accuracy for smaller values of Bessel_J1 but take longer to compute
 
  j=0;
  for (x=0.0; x <= x_max; x+=x_increment) {
    //
    // numerical integration
    //
    Bessel_J1=0.0;
    if (x != 0.0) {
      for (theta=0.0; theta <= M_PI; theta+=theta_increment) {
       Bessel_J1+=(cos(theta - (x * sin(theta))) * theta_increment / M_PI);
      }
    }

    // output sample
    if (j == 16) {
      j=0;
      printf("\n");
    }
    printf(" %9.2e", Bessel_J1);
    j++;
    if (x <= x_max) {
      printf(",");
    }
  } // end for x
  return(0);
}

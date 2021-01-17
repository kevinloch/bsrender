int limitIntensity(double *pixel_r, double *pixel_g, double *pixel_b) {

  //
  // limit pixel to range 0.0-1.0 without regard to color
  //
/*
  // pixel should never be negative
  if (*pixel_r < 0.0) {
    *pixel_r=0.0;
  }
  if (*pixel_g < 0.0) {
    *pixel_g=0.0;
  }
  if (*pixel_b < 0.0) {
    *pixel_b=0.0;
  }
*/
  if (*pixel_r > 1.0) {
    *pixel_r=1.0;
  }
  if (*pixel_g > 1.0) {
    *pixel_g=1.0;
  }
  if (*pixel_b > 1.0) {
    *pixel_b=1.0;
  }

  return(0);
}

int limitIntensityPreserveColor(double *pixel_r, double *pixel_g, double *pixel_b) {
  double pixel_max;

  //
  // limit pixel to range 0.0-1.0 while maintaining color (max channel=1.0)
  //
/*
  // pixel should never be negative 
  if (*pixel_r < 0) {
    *pixel_r=0;
  }
  if (*pixel_g < 0) {
    *pixel_g=0;
  }
  if (*pixel_b < 0) {
    *pixel_b=0;
  }
*/
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

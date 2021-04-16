#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr

int drawCrossHairs(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  int i;
  pixel_composition_t *image_composition_p;
  double res_x;
  double res_y;
  double half_res_x;
  double half_res_y;
  pixel_composition_t *image_buf;

  //
  // get current image buffer and resolution
  //
  image_buf=bsr_state->current_image_buf;
  res_x=(double)bsr_state->current_image_res_x;
  res_y=(double)bsr_state->current_image_res_y;
  half_res_x=res_x / 2.0;
  half_res_y=res_y / 2.0;

  //
  // draw crosshairs
  //
  for (i=(half_res_x - (res_y * 0.02)); i < (half_res_x - (res_y * 0.005)); i++) {
    image_composition_p=image_buf + (int)(res_x * half_res_y) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(half_res_x + (res_y * 0.005)); i < (half_res_x + (res_y * 0.02)); i++) {
    image_composition_p=image_buf + (int)(res_x * half_res_y) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(half_res_y - (res_y * 0.02)); i < (half_res_y - (res_y * 0.005)); i++) {
    image_composition_p=image_buf + (int)(res_x * (double)i) + (int)half_res_x;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(half_res_y + (res_y * 0.005)); i < (half_res_y + (res_y * 0.02)); i++) {
    image_composition_p=image_buf + (int)(res_x * (double)i) + (int)half_res_x;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  return(0);
}

int drawGridLines(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  int i;
  pixel_composition_t *image_composition_p;
  double res_x;
  double res_y;
  double half_res_x;
  double half_res_y;
  pixel_composition_t *image_buf;

  //
  // get current image buffer and resolution
  //
  image_buf=bsr_state->current_image_buf;
  res_x=(double)bsr_state->current_image_res_x;
  res_y=(double)bsr_state->current_image_res_y;
  half_res_x=res_x / 2.0;
  half_res_y=res_y / 2.0;

  //
  // draw select raster lines
  //
  for (i=0; i < res_x; i++) {
    image_composition_p=image_buf + (int)(res_x * (res_y * 0.25)) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_x; i++) {
    image_composition_p=image_buf + (int)(res_x * half_res_y) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_x; i++) {
    image_composition_p=image_buf + (int)(res_x * (res_y * 0.75)) + i;
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_y; i++) {
    image_composition_p=image_buf + (int)(res_x * (double)i) + (int)(res_x * 0.25);
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_y; i++) {
    image_composition_p=image_buf + (int)(res_x * (double)i) + (int)(half_res_x);
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < res_y; i++) {
    image_composition_p=image_buf + (int)(res_x * (double)i) + (int)(res_x * 0.75);
    image_composition_p->r=0.9;
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  return(0);
}

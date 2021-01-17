#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr

int drawCrossHairs(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  int i;
  pixel_composition_t *image_composition_p;

  //
  // optionally draw cross hairs
  //
  for (i=(bsr_state->camera_half_res_x - (bsr_config->camera_res_y * 0.02)); i < (bsr_state->camera_half_res_x - (bsr_config->camera_res_y * 0.005)); i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * bsr_state->camera_half_res_y) + i;
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(bsr_state->camera_half_res_x + (bsr_config->camera_res_y * 0.005)); i < (bsr_state->camera_half_res_x + (bsr_config->camera_res_y * 0.02)); i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * bsr_state->camera_half_res_y) + i;
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(bsr_state->camera_half_res_y - (bsr_config->camera_res_y * 0.02)); i < (bsr_state->camera_half_res_y - (bsr_config->camera_res_y * 0.005)); i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * i) + (int)bsr_state->camera_half_res_x;
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=(bsr_state->camera_half_res_y + (bsr_config->camera_res_y * 0.005)); i < (bsr_state->camera_half_res_y + (bsr_config->camera_res_y * 0.02)); i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * i) + (int)bsr_state->camera_half_res_x;
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  return(0);
}

int drawGridLines(bsr_config_t *bsr_config, bsr_state_t *bsr_state) {
  int i;
  pixel_composition_t *image_composition_p;

  //
  // optionally select raster lines
  //
  for (i=0; i < bsr_config->camera_res_x; i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * (bsr_config->camera_res_y * 0.25)) + i;
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < bsr_config->camera_res_x; i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * bsr_state->camera_half_res_y) + i;
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < bsr_config->camera_res_x; i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * (bsr_config->camera_res_y * 0.75)) + i;
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < bsr_config->camera_res_y; i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * i) + (int)(bsr_config->camera_res_x * 0.25);
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < bsr_config->camera_res_y; i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * i) + (int)(bsr_state->camera_half_res_x);
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  for (i=0; i < bsr_config->camera_res_y; i++) {
    image_composition_p=bsr_state->image_composition_buf + (int)(bsr_config->camera_res_x * i) + (int)(bsr_config->camera_res_x * 0.75);
    image_composition_p->r=(bsr_config->camera_pixel_limit * 0.9);
    image_composition_p->g=0.0;
    image_composition_p->b=0.0;
  }
  return(0);
}

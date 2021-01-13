#ifndef BSR_OVERLAY_H
#define BSR_OVERLAY_H

int drawCrossHairs(bsr_config_t *bsr_config, pixel_composition_t *image_composition_buf, int camera_half_res_x, int camera_half_res_y);
int drawGridLines(bsr_config_t *bsr_config, pixel_composition_t *image_composition_buf, int camera_half_res_x, int camera_half_res_y);

#endif // BSR_OVERLAY_H

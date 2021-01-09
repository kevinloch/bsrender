#ifndef NLE_LEPTON_H
#define NLE_LEPTON_H

#define BSRENDER_VERSION "0.9-dev-01"

typedef struct {
  double icrs_x;
  double icrs_y;
  double icrs_z;
  uint64_t intensity_and_temperature;
} star_record_t;

typedef struct {
  double r;
  double g;
  double b;
} pixel_composition_t;

#endif // NLE_LEPTON_H


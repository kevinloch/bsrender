#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include "Gaia-EDR3-transmissivity.h"

double getGaiaTransmissivityG(int wavelength) {
  double transmissivity;

  // 320-1100nm is the range of the Gaia EDR3 passband data file https://www.cosmos.esa.int/web/gaia/edr3-passbands
  if (wavelength < 320) {
    return(0.0);
  } else if (wavelength > 1100) {
    return(0.0);
  }

  transmissivity=Gaia_EDR3_transmissivity_G[wavelength - 320];

  return(transmissivity);
}

double getGaiaTransmissivityBp(int wavelength) {
  double transmissivity;

  // 320-1100nm is the range of the Gaia EDR3 passband data file https://www.cosmos.esa.int/web/gaia/edr3-passbands
  if (wavelength < 320) {
    return(0.0);
  } else if (wavelength > 1100) {
    return(0.0);
  }

  transmissivity=Gaia_EDR3_transmissivity_bp[wavelength - 320];

  return(transmissivity);
}

double getGaiaTransmissivityRp(int wavelength) {
  double transmissivity;

  // 320-1100nm is the range of the Gaia EDR3 passband data file https://www.cosmos.esa.int/web/gaia/edr3-passbands
  if (wavelength < 320) {
    return(0.0);
  } else if (wavelength > 1100) {
    return(0.0);
  }

  transmissivity=Gaia_EDR3_transmissivity_rp[wavelength - 320];

  return(transmissivity);
}

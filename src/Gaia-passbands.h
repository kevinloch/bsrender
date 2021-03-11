#ifndef BSR_GAIA_PASSBANDS_H
#define BSR_GAIA_PASSBANDS_H

const double Gaia_Gband_long_limit=1050.0;
const double Gaia_Gband_short_limit=320.0;
const double Gaia_Gband_scalar=6.77892686e-01; // G-band transmissivity at 550nm

double getGaiaTransmissivityG(int wavelength); 
double getGaiaTransmissivityBp(int wavelength); 
double getGaiaTransmissivityRp(int wavelength); 

#endif // BSR_GAIA_PASSBANDS_H

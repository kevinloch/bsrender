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

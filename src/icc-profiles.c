//
// Billion Star 3D Rendering Engine
// Kevin M. Loch
//
// 3D rendering engine for the ESA Gaia DR3 star dataset

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

//
// Basic chomaticities for EXR format
//
const chromaticities_t sRGB_c = { // https://en.wikipedia.org/wiki/SRGB
  .redX  = 0.6400f,
  .redY  = 0.3300f,
  .greenX= 0.3000f,
  .greenY= 0.6000f,
  .blueX = 0.1500f,
  .blueY = 0.0600f,
  .whiteX= 0.3127f,
  .whiteY= 0.3290f
};

const chromaticities_t DisplayP3_c = { // https://en.wikipedia.org/wiki/DCI-P3
  .redX  = 0.6800f,
  .redY  = 0.3200f,
  .greenX= 0.2650f,
  .greenY= 0.6900f,
  .blueX = 0.1500f,
  .blueY = 0.0600f,
  .whiteX= 0.3127f,
  .whiteY= 0.3290f
};

const chromaticities_t Rec2020_c = { // https://en.wikipedia.org/wiki/Rec._2020
  .redX  = 0.7080f,
  .redY  = 0.3290f,
  .greenX= 0.1700f,
  .greenY= 0.7970f,
  .blueX = 0.1310f,
  .blueY = 0.0460f,
  .whiteX= 0.3127f,
  .whiteY= 0.3290f
};

const chromaticities_t Rec601NTSC_c = { // https://en.wikipedia.org/wiki/Rec._601
  .redX  = 0.6300f,
  .redY  = 0.3400f,
  .greenX= 0.3100f,
  .greenY= 0.5950f,
  .blueX = 0.1550f,
  .blueY = 0.0700f,
  .whiteX= 0.3127f,
  .whiteY= 0.3290f
};

const chromaticities_t Rec601PAL_c = { // https://en.wikipedia.org/wiki/Rec._601
  .redX  = 0.6400f,
  .redY  = 0.3300f,
  .greenX= 0.2900f,
  .greenY= 0.6000f,
  .blueX = 0.1500f,
  .blueY = 0.0600f,
  .whiteX= 0.3127f,
  .whiteY= 0.3290f
};

const chromaticities_t Rec709_c = { // https://en.wikipedia.org/wiki/Rec._709
  .redX  = 0.6400f,
  .redY  = 0.3300f,
  .greenX= 0.3000f,
  .greenY= 0.6000f,
  .blueX = 0.1500f,
  .blueY = 0.0600f,
  .whiteX= 0.3127f,
  .whiteY= 0.3290f
};

//
// ICC Profiles from https://github.com/saucecontrol/Compact-ICC-Profiles
//
const unsigned char sRGB_v4_icc[] = {
  0x00, 0x00, 0x01, 0xe0, 0x6c, 0x63, 0x6d, 0x73, 0x04, 0x20, 0x00, 0x00,
  0x6d, 0x6e, 0x74, 0x72, 0x52, 0x47, 0x42, 0x20, 0x58, 0x59, 0x5a, 0x20,
  0x07, 0xe2, 0x00, 0x03, 0x00, 0x14, 0x00, 0x09, 0x00, 0x0e, 0x00, 0x1d,
  0x61, 0x63, 0x73, 0x70, 0x4d, 0x53, 0x46, 0x54, 0x00, 0x00, 0x00, 0x00,
  0x73, 0x61, 0x77, 0x73, 0x63, 0x74, 0x72, 0x6c, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x68, 0x61, 0x6e, 0x64,
  0xa3, 0xb2, 0xab, 0xdf, 0x5c, 0xa7, 0x03, 0x12, 0xa8, 0x55, 0xa4, 0xec,
  0x35, 0x7a, 0xd1, 0xf3, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x64, 0x65, 0x73, 0x63, 0x00, 0x00, 0x00, 0xfc, 0x00, 0x00, 0x00, 0x24,
  0x63, 0x70, 0x72, 0x74, 0x00, 0x00, 0x01, 0x20, 0x00, 0x00, 0x00, 0x22,
  0x77, 0x74, 0x70, 0x74, 0x00, 0x00, 0x01, 0x44, 0x00, 0x00, 0x00, 0x14,
  0x63, 0x68, 0x61, 0x64, 0x00, 0x00, 0x01, 0x58, 0x00, 0x00, 0x00, 0x2c,
  0x72, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x84, 0x00, 0x00, 0x00, 0x14,
  0x67, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x98, 0x00, 0x00, 0x00, 0x14,
  0x62, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0xac, 0x00, 0x00, 0x00, 0x14,
  0x72, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x67, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x62, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x08,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x73, 0x00, 0x52, 0x00, 0x47, 0x00, 0x42,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x43, 0x00, 0x43, 0x00, 0x30, 0x00, 0x21,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x73, 0x66, 0x33, 0x32,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x0c, 0x3f, 0x00, 0x00, 0x05, 0xdd,
  0xff, 0xff, 0xf3, 0x26, 0x00, 0x00, 0x07, 0x90, 0x00, 0x00, 0xfd, 0x92,
  0xff, 0xff, 0xfb, 0xa1, 0xff, 0xff, 0xfd, 0xa2, 0x00, 0x00, 0x03, 0xdc,
  0x00, 0x00, 0xc0, 0x71, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x6f, 0xa0, 0x00, 0x00, 0x38, 0xf2, 0x00, 0x00, 0x03, 0x8f,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x62, 0x96,
  0x00, 0x00, 0xb7, 0x89, 0x00, 0x00, 0x18, 0xda, 0x58, 0x59, 0x5a, 0x20,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0xa0, 0x00, 0x00, 0x0f, 0x85,
  0x00, 0x00, 0xb6, 0xc4, 0x70, 0x61, 0x72, 0x61, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x03, 0x00, 0x00, 0x00, 0x02, 0x66, 0x69, 0x00, 0x00, 0xf2, 0xa7,
  0x00, 0x00, 0x0d, 0x59, 0x00, 0x00, 0x13, 0xd0, 0x00, 0x00, 0x0a, 0x5b
};
const unsigned int sRGB_v4_icc_len = 480;

const unsigned char DisplayP3Compat_v4_icc[] = {
  0x00, 0x00, 0x01, 0xe0, 0x6c, 0x63, 0x6d, 0x73, 0x04, 0x20, 0x00, 0x00,
  0x6d, 0x6e, 0x74, 0x72, 0x52, 0x47, 0x42, 0x20, 0x58, 0x59, 0x5a, 0x20,
  0x07, 0xe2, 0x00, 0x03, 0x00, 0x14, 0x00, 0x09, 0x00, 0x0e, 0x00, 0x1d,
  0x61, 0x63, 0x73, 0x70, 0x4d, 0x53, 0x46, 0x54, 0x00, 0x00, 0x00, 0x00,
  0x73, 0x61, 0x77, 0x73, 0x63, 0x74, 0x72, 0x6c, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x68, 0x61, 0x6e, 0x64,
  0x1e, 0x87, 0x9c, 0x47, 0x82, 0xcf, 0x9b, 0xc1, 0x16, 0x6f, 0x39, 0x60,
  0xae, 0x45, 0x75, 0x1e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x64, 0x65, 0x73, 0x63, 0x00, 0x00, 0x00, 0xfc, 0x00, 0x00, 0x00, 0x24,
  0x63, 0x70, 0x72, 0x74, 0x00, 0x00, 0x01, 0x20, 0x00, 0x00, 0x00, 0x22,
  0x77, 0x74, 0x70, 0x74, 0x00, 0x00, 0x01, 0x44, 0x00, 0x00, 0x00, 0x14,
  0x63, 0x68, 0x61, 0x64, 0x00, 0x00, 0x01, 0x58, 0x00, 0x00, 0x00, 0x2c,
  0x72, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x84, 0x00, 0x00, 0x00, 0x14,
  0x67, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x98, 0x00, 0x00, 0x00, 0x14,
  0x62, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0xac, 0x00, 0x00, 0x00, 0x14,
  0x72, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x67, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x62, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x08,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x73, 0x00, 0x50, 0x00, 0x33, 0x00, 0x43,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x43, 0x00, 0x43, 0x00, 0x30, 0x00, 0x21,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x73, 0x66, 0x33, 0x32,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x0c, 0x42, 0x00, 0x00, 0x05, 0xde,
  0xff, 0xff, 0xf3, 0x25, 0x00, 0x00, 0x07, 0x93, 0x00, 0x00, 0xfd, 0x90,
  0xff, 0xff, 0xfb, 0xa1, 0xff, 0xff, 0xfe, 0x4e, 0x00, 0x00, 0x03, 0x9a,
  0x00, 0x00, 0xc0, 0x14, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x83, 0xdf, 0x00, 0x00, 0x3d, 0xbf, 0x00, 0x00, 0x00, 0x00,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4a, 0xbf,
  0x00, 0x00, 0xb1, 0x37, 0x00, 0x00, 0x0a, 0xb5, 0x58, 0x59, 0x5a, 0x20,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x28, 0x38, 0x00, 0x00, 0x11, 0x0a,
  0x00, 0x00, 0xc8, 0x78, 0x70, 0x61, 0x72, 0x61, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x03, 0x00, 0x00, 0x00, 0x02, 0x66, 0x69, 0x00, 0x00, 0xf2, 0xa7,
  0x00, 0x00, 0x0d, 0x59, 0x00, 0x00, 0x13, 0xd0, 0x00, 0x00, 0x0a, 0x5b
};
const unsigned int DisplayP3Compat_v4_icc_len = 480;

const unsigned char Rec2020Compat_v4_icc[] = {
  0x00, 0x00, 0x01, 0xe0, 0x6c, 0x63, 0x6d, 0x73, 0x04, 0x20, 0x00, 0x00,
  0x6d, 0x6e, 0x74, 0x72, 0x52, 0x47, 0x42, 0x20, 0x58, 0x59, 0x5a, 0x20,
  0x07, 0xe2, 0x00, 0x03, 0x00, 0x14, 0x00, 0x09, 0x00, 0x0e, 0x00, 0x1d,
  0x61, 0x63, 0x73, 0x70, 0x4d, 0x53, 0x46, 0x54, 0x00, 0x00, 0x00, 0x00,
  0x73, 0x61, 0x77, 0x73, 0x63, 0x74, 0x72, 0x6c, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x68, 0x61, 0x6e, 0x64,
  0xe9, 0x37, 0xbc, 0x8d, 0x23, 0xdd, 0x2a, 0x15, 0x5a, 0xb3, 0x46, 0xb4,
  0xc4, 0xdc, 0x56, 0xd9, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x64, 0x65, 0x73, 0x63, 0x00, 0x00, 0x00, 0xfc, 0x00, 0x00, 0x00, 0x24,
  0x63, 0x70, 0x72, 0x74, 0x00, 0x00, 0x01, 0x20, 0x00, 0x00, 0x00, 0x22,
  0x77, 0x74, 0x70, 0x74, 0x00, 0x00, 0x01, 0x44, 0x00, 0x00, 0x00, 0x14,
  0x63, 0x68, 0x61, 0x64, 0x00, 0x00, 0x01, 0x58, 0x00, 0x00, 0x00, 0x2c,
  0x72, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x84, 0x00, 0x00, 0x00, 0x14,
  0x67, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x98, 0x00, 0x00, 0x00, 0x14,
  0x62, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0xac, 0x00, 0x00, 0x00, 0x14,
  0x72, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x67, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x62, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x08,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x32, 0x00, 0x30, 0x00, 0x32, 0x00, 0x43,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x43, 0x00, 0x43, 0x00, 0x30, 0x00, 0x21,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x73, 0x66, 0x33, 0x32,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x0c, 0x42, 0x00, 0x00, 0x05, 0xde,
  0xff, 0xff, 0xf3, 0x25, 0x00, 0x00, 0x07, 0x93, 0x00, 0x00, 0xfd, 0x90,
  0xff, 0xff, 0xfb, 0xa1, 0xff, 0xff, 0xfe, 0x7c, 0x00, 0x00, 0x03, 0xac,
  0x00, 0x00, 0xbf, 0xdb, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0xac, 0x69, 0x00, 0x00, 0x47, 0x70, 0x00, 0x00, 0x00, 0x00,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x2a, 0x6a,
  0x00, 0x00, 0xac, 0xe3, 0x00, 0x00, 0x07, 0xa8, 0x58, 0x59, 0x5a, 0x20,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x20, 0x03, 0x00, 0x00, 0x0b, 0xad,
  0x00, 0x00, 0xcb, 0x85, 0x70, 0x61, 0x72, 0x61, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x03, 0x00, 0x00, 0x00, 0x02, 0x38, 0xe5, 0x00, 0x00, 0xe8, 0xe0,
  0x00, 0x00, 0x17, 0x20, 0x00, 0x00, 0x38, 0xe4, 0x00, 0x00, 0x14, 0xbd
};
const unsigned int Rec2020Compat_v4_icc_len = 480;

const unsigned char Rec601NTSC_v4_icc[] = {
  0x00, 0x00, 0x01, 0xe0, 0x6c, 0x63, 0x6d, 0x73, 0x04, 0x20, 0x00, 0x00,
  0x6d, 0x6e, 0x74, 0x72, 0x52, 0x47, 0x42, 0x20, 0x58, 0x59, 0x5a, 0x20,
  0x07, 0xe5, 0x00, 0x04, 0x00, 0x1b, 0x00, 0x0a, 0x00, 0x1b, 0x00, 0x00,
  0x61, 0x63, 0x73, 0x70, 0x4d, 0x53, 0x46, 0x54, 0x00, 0x00, 0x00, 0x00,
  0x73, 0x61, 0x77, 0x73, 0x63, 0x74, 0x72, 0x6c, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x68, 0x61, 0x6e, 0x64,
  0x54, 0x36, 0xe4, 0x73, 0xe8, 0xae, 0x23, 0x3e, 0x5a, 0x41, 0xfc, 0x2a,
  0x61, 0x3e, 0x84, 0x85, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x64, 0x65, 0x73, 0x63, 0x00, 0x00, 0x00, 0xfc, 0x00, 0x00, 0x00, 0x24,
  0x63, 0x70, 0x72, 0x74, 0x00, 0x00, 0x01, 0x20, 0x00, 0x00, 0x00, 0x22,
  0x77, 0x74, 0x70, 0x74, 0x00, 0x00, 0x01, 0x44, 0x00, 0x00, 0x00, 0x14,
  0x63, 0x68, 0x61, 0x64, 0x00, 0x00, 0x01, 0x58, 0x00, 0x00, 0x00, 0x2c,
  0x72, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x84, 0x00, 0x00, 0x00, 0x14,
  0x67, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x98, 0x00, 0x00, 0x00, 0x14,
  0x62, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0xac, 0x00, 0x00, 0x00, 0x14,
  0x72, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x67, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x62, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x08,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x52, 0x00, 0x36, 0x00, 0x30, 0x00, 0x31,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x43, 0x00, 0x43, 0x00, 0x30, 0x00, 0x21,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x73, 0x66, 0x33, 0x32,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x0c, 0x42, 0x00, 0x00, 0x05, 0xde,
  0xff, 0xff, 0xf3, 0x25, 0x00, 0x00, 0x07, 0x93, 0x00, 0x00, 0xfd, 0x90,
  0xff, 0xff, 0xfb, 0xa1, 0xff, 0xff, 0xfd, 0xa2, 0x00, 0x00, 0x03, 0xdc,
  0x00, 0x00, 0xc0, 0x6e, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x6a, 0x92, 0x00, 0x00, 0x38, 0xc0, 0x00, 0x00, 0x03, 0x7e,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x64, 0xa9,
  0x00, 0x00, 0xb4, 0x0a, 0x00, 0x00, 0x17, 0x61, 0x58, 0x59, 0x5a, 0x20,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x27, 0x9b, 0x00, 0x00, 0x13, 0x36,
  0x00, 0x00, 0xb8, 0x4e, 0x70, 0x61, 0x72, 0x61, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x03, 0x00, 0x00, 0x00, 0x02, 0x38, 0xe5, 0x00, 0x00, 0xe8, 0xf0,
  0x00, 0x00, 0x17, 0x10, 0x00, 0x00, 0x38, 0xe4, 0x00, 0x00, 0x14, 0xbd
};
const unsigned int Rec601NTSC_v4_icc_len = 480;

const unsigned char Rec601PAL_v4_icc[] = {
  0x00, 0x00, 0x01, 0xe0, 0x6c, 0x63, 0x6d, 0x73, 0x04, 0x20, 0x00, 0x00,
  0x6d, 0x6e, 0x74, 0x72, 0x52, 0x47, 0x42, 0x20, 0x58, 0x59, 0x5a, 0x20,
  0x07, 0xe5, 0x00, 0x04, 0x00, 0x1b, 0x00, 0x0a, 0x00, 0x1b, 0x00, 0x00,
  0x61, 0x63, 0x73, 0x70, 0x4d, 0x53, 0x46, 0x54, 0x00, 0x00, 0x00, 0x00,
  0x73, 0x61, 0x77, 0x73, 0x63, 0x74, 0x72, 0x6c, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x68, 0x61, 0x6e, 0x64,
  0x30, 0xef, 0x69, 0x97, 0xe0, 0x26, 0xdb, 0x92, 0x81, 0x95, 0xd0, 0x62,
  0x08, 0xbc, 0x3d, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x64, 0x65, 0x73, 0x63, 0x00, 0x00, 0x00, 0xfc, 0x00, 0x00, 0x00, 0x24,
  0x63, 0x70, 0x72, 0x74, 0x00, 0x00, 0x01, 0x20, 0x00, 0x00, 0x00, 0x22,
  0x77, 0x74, 0x70, 0x74, 0x00, 0x00, 0x01, 0x44, 0x00, 0x00, 0x00, 0x14,
  0x63, 0x68, 0x61, 0x64, 0x00, 0x00, 0x01, 0x58, 0x00, 0x00, 0x00, 0x2c,
  0x72, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x84, 0x00, 0x00, 0x00, 0x14,
  0x67, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x98, 0x00, 0x00, 0x00, 0x14,
  0x62, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0xac, 0x00, 0x00, 0x00, 0x14,
  0x72, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x67, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x62, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x08,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x36, 0x00, 0x30, 0x00, 0x31, 0x00, 0x50,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x43, 0x00, 0x43, 0x00, 0x30, 0x00, 0x21,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x73, 0x66, 0x33, 0x32,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x0c, 0x42, 0x00, 0x00, 0x05, 0xde,
  0xff, 0xff, 0xf3, 0x25, 0x00, 0x00, 0x07, 0x93, 0x00, 0x00, 0xfd, 0x90,
  0xff, 0xff, 0xfb, 0xa1, 0xff, 0xff, 0xfd, 0xa2, 0x00, 0x00, 0x03, 0xdc,
  0x00, 0x00, 0xc0, 0x6e, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x74, 0x8b, 0x00, 0x00, 0x3b, 0x77, 0x00, 0x00, 0x03, 0xb8,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x5e, 0x1b,
  0x00, 0x00, 0xb5, 0x34, 0x00, 0x00, 0x1a, 0xd9, 0x58, 0x59, 0x5a, 0x20,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x30, 0x00, 0x00, 0x0f, 0x55,
  0x00, 0x00, 0xb4, 0x9c, 0x70, 0x61, 0x72, 0x61, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x03, 0x00, 0x00, 0x00, 0x02, 0x38, 0xe5, 0x00, 0x00, 0xe8, 0xf0,
  0x00, 0x00, 0x17, 0x10, 0x00, 0x00, 0x38, 0xe4, 0x00, 0x00, 0x14, 0xbd
};
const unsigned int Rec601PAL_v4_icc_len = 480;

const unsigned char Rec709_v4_icc[] = {
  0x00, 0x00, 0x01, 0xe0, 0x6c, 0x63, 0x6d, 0x73, 0x04, 0x20, 0x00, 0x00,
  0x6d, 0x6e, 0x74, 0x72, 0x52, 0x47, 0x42, 0x20, 0x58, 0x59, 0x5a, 0x20,
  0x07, 0xe2, 0x00, 0x03, 0x00, 0x14, 0x00, 0x09, 0x00, 0x0e, 0x00, 0x1d,
  0x61, 0x63, 0x73, 0x70, 0x4d, 0x53, 0x46, 0x54, 0x00, 0x00, 0x00, 0x00,
  0x73, 0x61, 0x77, 0x73, 0x63, 0x74, 0x72, 0x6c, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x68, 0x61, 0x6e, 0x64,
  0x9e, 0x30, 0xbd, 0xa0, 0xd7, 0xde, 0xd7, 0xee, 0xa3, 0x0c, 0x38, 0x81,
  0xe9, 0x9d, 0xc7, 0xc9, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0a,
  0x64, 0x65, 0x73, 0x63, 0x00, 0x00, 0x00, 0xfc, 0x00, 0x00, 0x00, 0x24,
  0x63, 0x70, 0x72, 0x74, 0x00, 0x00, 0x01, 0x20, 0x00, 0x00, 0x00, 0x22,
  0x77, 0x74, 0x70, 0x74, 0x00, 0x00, 0x01, 0x44, 0x00, 0x00, 0x00, 0x14,
  0x63, 0x68, 0x61, 0x64, 0x00, 0x00, 0x01, 0x58, 0x00, 0x00, 0x00, 0x2c,
  0x72, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x84, 0x00, 0x00, 0x00, 0x14,
  0x67, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x98, 0x00, 0x00, 0x00, 0x14,
  0x62, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0xac, 0x00, 0x00, 0x00, 0x14,
  0x72, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x67, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x62, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0xc0, 0x00, 0x00, 0x00, 0x20,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x08,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x52, 0x00, 0x37, 0x00, 0x30, 0x00, 0x39,
  0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x06,
  0x00, 0x00, 0x00, 0x1c, 0x00, 0x43, 0x00, 0x43, 0x00, 0x30, 0x00, 0x21,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x73, 0x66, 0x33, 0x32,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x0c, 0x42, 0x00, 0x00, 0x05, 0xde,
  0xff, 0xff, 0xf3, 0x25, 0x00, 0x00, 0x07, 0x93, 0x00, 0x00, 0xfd, 0x90,
  0xff, 0xff, 0xfb, 0xa1, 0xff, 0xff, 0xfd, 0xa2, 0x00, 0x00, 0x03, 0xdc,
  0x00, 0x00, 0xc0, 0x6e, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x6f, 0xa0, 0x00, 0x00, 0x38, 0xf5, 0x00, 0x00, 0x03, 0x90,
  0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x62, 0x97,
  0x00, 0x00, 0xb7, 0x87, 0x00, 0x00, 0x18, 0xda, 0x58, 0x59, 0x5a, 0x20,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x9f, 0x00, 0x00, 0x0f, 0x84,
  0x00, 0x00, 0xb6, 0xc3, 0x70, 0x61, 0x72, 0x61, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x03, 0x00, 0x00, 0x00, 0x02, 0x38, 0xe5, 0x00, 0x00, 0xe8, 0xf0,
  0x00, 0x00, 0x17, 0x10, 0x00, 0x00, 0x38, 0xe4, 0x00, 0x00, 0x14, 0xbd
};
const unsigned int Rec709_v4_icc_len = 480;

// open source test profile for Rec. 2100 PQ. Derived partially from Rec2020Compat_v4_icc
// but with cicp and lumi tags added. Images will appear very low contrast/washed out if
// viewer does not use the cicp and lumi tags
const unsigned char Rec2100PQ_v4_icc[] = {
0x00, 0x00, 0x02, 0x18, 0x61, 0x70, 0x70, 0x6c, 0x04, 0x40, 0x00, 0x00, 0x6d, 0x6e, 0x74, 0x72,
0x52, 0x47, 0x42, 0x20, 0x58, 0x59, 0x5a, 0x20, 0x07, 0xe7, 0x00, 0x07, 0x00, 0x06, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x61, 0x63, 0x73, 0x70, 0x41, 0x50, 0x50, 0x4c, 0x00, 0x00, 0x00, 0x00,
0x62, 0x73, 0x72, 0x65, 0x6e, 0x64, 0x65, 0x72, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d,
0x00, 0x00, 0x00, 0x00, 0x93, 0xea, 0xf9, 0xc8, 0x3e, 0xb5, 0x57, 0x08, 0x0c, 0x5b, 0xf5, 0xb2,
0xfb, 0xe0, 0xc7, 0xb6, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x0c, 0x64, 0x65, 0x73, 0x63, 0x00, 0x00, 0x01, 0x14, 0x00, 0x00, 0x00, 0x24,
0x72, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x38, 0x00, 0x00, 0x00, 0x14, 0x67, 0x58, 0x59, 0x5a,
0x00, 0x00, 0x01, 0x4c, 0x00, 0x00, 0x00, 0x14, 0x62, 0x58, 0x59, 0x5a, 0x00, 0x00, 0x01, 0x60,
0x00, 0x00, 0x00, 0x14, 0x72, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0x74, 0x00, 0x00, 0x00, 0x20,
0x67, 0x54, 0x52, 0x43, 0x00, 0x00, 0x01, 0x74, 0x00, 0x00, 0x00, 0x20, 0x62, 0x54, 0x52, 0x43,
0x00, 0x00, 0x01, 0x74, 0x00, 0x00, 0x00, 0x20, 0x77, 0x74, 0x70, 0x74, 0x00, 0x00, 0x01, 0x94,
0x00, 0x00, 0x00, 0x14, 0x63, 0x70, 0x72, 0x74, 0x00, 0x00, 0x01, 0xa8, 0x00, 0x00, 0x00, 0x24,
0x63, 0x68, 0x61, 0x64, 0x00, 0x00, 0x01, 0xcc, 0x00, 0x00, 0x00, 0x2c, 0x63, 0x69, 0x63, 0x70,
0x00, 0x00, 0x01, 0xf8, 0x00, 0x00, 0x00, 0x0c, 0x6c, 0x75, 0x6d, 0x69, 0x00, 0x00, 0x02, 0x04,
0x00, 0x00, 0x00, 0x14, 0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, 0x1c,
0x00, 0x32, 0x00, 0x31, 0x00, 0x50, 0x00, 0x51, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0xac, 0x69, 0x00, 0x00, 0x47, 0x70, 0x00, 0x00, 0x00, 0x00, 0x58, 0x59, 0x5a, 0x20,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x2a, 0x6a, 0x00, 0x00, 0xac, 0xe3, 0x00, 0x00, 0x07, 0xa8,
0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x20, 0x03, 0x00, 0x00, 0x0b, 0xad,
0x00, 0x00, 0xcb, 0x85, 0x70, 0x61, 0x72, 0x61, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00,
0x00, 0x02, 0x38, 0xe5, 0x00, 0x00, 0xe8, 0xe0, 0x00, 0x00, 0x17, 0x20, 0x00, 0x00, 0x38, 0xe4,
0x00, 0x00, 0x14, 0xbd, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf6, 0xd6,
0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0xd3, 0x2d, 0x6d, 0x6c, 0x75, 0x63, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0c, 0x65, 0x6e, 0x55, 0x53, 0x00, 0x00, 0x00, 0x06,
0x00, 0x00, 0x00, 0x1c, 0x00, 0x43, 0x00, 0x43, 0x00, 0x30, 0x00, 0x00, 0x73, 0x66, 0x33, 0x32,
0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x0c, 0x42, 0x00, 0x00, 0x05, 0xde, 0xff, 0xff, 0xf3, 0x25,
0x00, 0x00, 0x07, 0x93, 0x00, 0x00, 0xfd, 0x90, 0xff, 0xff, 0xfb, 0xa1, 0xff, 0xff, 0xfe, 0x7c,
0x00, 0x00, 0x03, 0xac, 0x00, 0x00, 0xbf, 0xdb, 0x63, 0x69, 0x63, 0x70, 0x00, 0x00, 0x00, 0x00,
0x09, 0x10, 0x00, 0x01, 0x58, 0x59, 0x5a, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x64, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};
const unsigned int Rec2100PQ_v4_icc_len = 536;

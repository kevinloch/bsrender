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

#ifndef BSR_ICC_PROFILES_H
#define BSR_ICC_PROFILES_H

const chromaticities_t sRGB_c;
const chromaticities_t DisplayP3_c;
const chromaticities_t Rec2020_c;
const chromaticities_t Rec601NTSC_c;
const chromaticities_t Rec601PAL_c;
const chromaticities_t Rec709_c;
extern const unsigned char sRGB_v4_icc[];
extern const unsigned int sRGB_v4_icc_len;
extern const unsigned char DisplayP3Compat_v4_icc[];
extern const unsigned int DisplayP3Compat_v4_icc_len;
extern const unsigned char Rec2020Compat_v4_icc[];
extern const unsigned int Rec2020Compat_v4_icc_len;
extern const unsigned char Rec601NTSC_v4_icc[];
extern const unsigned int Rec601NTSC_v4_icc_len;
extern const unsigned char Rec601PAL_v4_icc[];
extern const unsigned int Rec601PAL_v4_icc_len;
extern const unsigned char Rec709_v4_icc[];
extern const unsigned int Rec709_v4_icc_len;

#endif // BSR_ICC_PROFILES_H

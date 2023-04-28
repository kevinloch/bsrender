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

#ifndef BSR_UTIL_H
#define BSR_UTIL_H

int littleEndianTest();
int storeStr32(unsigned char *dest, char *src);
int storeU8(unsigned char *dest, unsigned char src);
int storeU16LE(unsigned char *dest, uint16_t src);
int storeU16BE(unsigned char *dest, uint16_t src);
int storeI32LE(unsigned char *dest, int32_t src);
int storeU32LE(unsigned char *dest, uint32_t src);
int storeU64LE(unsigned char *dest, uint64_t src);
int storeHalfLE(unsigned char *dest, float src);
int storeFloatLE(unsigned char *dest, float src);
int getQueryString(bsr_config_t *bsr_config);
int printVersion(bsr_config_t *bsr_config);
int waitForWorkerThreads(bsr_state_t *bsr_state, int min_status);
int waitForMainThread(bsr_state_t *bsr_state, int min_status);
int checkExceptions(bsr_state_t *bsr_state);
int limitIntensity(bsr_config_t *bsr_config, double *pixel_r, double *pixel_g, double *pixel_b);
int limitIntensityPreserveColor(bsr_config_t *bsr_config, double *pixel_r, double *pixel_g, double *pixel_b);

#endif // BSR_UTIL_H

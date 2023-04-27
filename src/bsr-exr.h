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

#ifndef BSR_EXR_H
#define BSR_EXR_H

//
// EXR attribute definitions copied from openexr-3.1.7/src/lib/OpenEXRCore/openexr_attr.h
//
/*
** SPDX-License-Identifier: BSD-3-Clause
** Copyright Contributors to the OpenEXR Project.
*/

/** Enum declaring allowed values for \c uint8_t value stored in built-in compression type. */
typedef enum
{
    EXR_COMPRESSION_NONE  = 0,
    EXR_COMPRESSION_RLE   = 1,
    EXR_COMPRESSION_ZIPS  = 2,
    EXR_COMPRESSION_ZIP   = 3,
    EXR_COMPRESSION_PIZ   = 4,
    EXR_COMPRESSION_PXR24 = 5,
    EXR_COMPRESSION_B44   = 6,
    EXR_COMPRESSION_B44A  = 7,
    EXR_COMPRESSION_DWAA  = 8,
    EXR_COMPRESSION_DWAB  = 9,
    EXR_COMPRESSION_LAST_TYPE /**< Invalid value, provided for range checking. */
} exr_compression_t;

/** Enum declaring allowed values for \c uint8_t value stored in \c lineOrder type. */
typedef enum
{
    EXR_LINEORDER_INCREASING_Y = 0,
    EXR_LINEORDER_DECREASING_Y = 1,
    EXR_LINEORDER_RANDOM_Y     = 2,
    EXR_LINEORDER_LAST_TYPE /**< Invalid value, provided for range checking. */
} exr_lineorder_t;

/** Enum declaring allowed values for part type. */
typedef enum
{
    EXR_STORAGE_SCANLINE = 0,  /**< Corresponds to type of \c scanlineimage. */
    EXR_STORAGE_TILED,         /**< Corresponds to type of \c tiledimage. */
    EXR_STORAGE_DEEP_SCANLINE, /**< Corresponds to type of \c deepscanline. */
    EXR_STORAGE_DEEP_TILED,    /**< Corresponds to type of \c deeptile. */
    EXR_STORAGE_LAST_TYPE      /**< Invalid value, provided for range checking. */
} exr_storage_t;

/** @brief Enum capturing the underlying data type on a channel. */
typedef enum
{
    EXR_PIXEL_UINT  = 0,
    EXR_PIXEL_HALF  = 1,
    EXR_PIXEL_FLOAT = 2,
    EXR_PIXEL_LAST_TYPE
} exr_pixel_type_t;

/** Hint for lossy compression methods about how to treat values
 * (logarithmic or linear), meaning a human sees values like R, G, B,
 * luminance difference between 0.1 and 0.2 as about the same as 1.0
 * to 2.0 (logarithmic), where chroma coordinates are closer to linear
 * (0.1 and 0.2 is about the same difference as 1.0 and 1.1).
 */
typedef enum
{
    EXR_PERCEPTUALLY_LOGARITHMIC  = 0,
    EXR_PERCEPTUALLY_LINEAR   = 1
} exr_perceptual_treatment_t;

int storeI32LE(unsigned char *dest, int32_t src);
int storeU8(unsigned char *dest, unsigned char src);
int storeStr32(unsigned char *dest, char *src);
int outputEXR(bsr_config_t *bsr_config, bsr_state_t *bsr_state);

#endif // BSR_EXR_H

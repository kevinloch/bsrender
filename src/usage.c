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
#include <stdio.h>

void printUsage() {
  printf("bsrender version %s, render 2D image from 3D star database\n", BSR_VERSION);
  printf("\
\n\
Usage:\n\
     bsrender [OPTION]...\n\
\n\
     See sample conifguration file for default settings. Individual configuration options are applied in this order:\n\
     1. built-in defaults, 2. configuration file, 3. command line flags, 4. environment QUERY_STRING (if CGI mode)\n\
\n\
Command line only options:\n\
     -c FILE                              Set configuration file name (default: bsrender.cfg)\n\
     --help, -h                           Show usage\n\
\n\
Privileged options - these cannot be changed by remote users in CGI mode:\n\
     --data_file_directory=DIR, -d        Path to galaxy-* data files, limit 255 characters\n\
     --output_file_name=FILE, -o          Output filename, may include path, limit 255 characters\n\
     --print_status=BOOL, -q              yes = sppress non-error status messages (also -q)\n\
                                          no = will allow informational status messages\n\
                                          All messages are always suppressed in CGI mode\n\
     --num_threads=NUM                    Total number of threads including main thread and worker\n\
                                          threads (minimum 2)\n\
                                          For best performance set to number of vcpus\n\
     --per_thread_buffer=NUM              Number of stars to buffer between each worker thread and main thread\n\
                                          Also sets size of dedup buffer for each thread\n\
     --per_thread_buffer_Airy=NUM         Number of stars to buffer between each worker thread and main thread\n\
                                          when Airy disks are enabled\n\
                                          Also sets size of dedup buffer for each thread\n\
     --cgi_mode=BOOL                      yes = enable CGI mode (html headers and png data written to stdout)\n\
     --cgi_max_res_x=NUM                  Maximum allowed horizontal resolution for CGI users\n\
     --cgi_max_res_y=NUM                  Maximum allowed vertical resolution for CGI users\n\
     --cgi_Gaia_min_parallax_quality=NUM  Minimum allowed parallax quality of Gaia stars for CGI users\n\
     --cgi_allow_Airy_disk=BOOL           yes = Airy disk mode is allowed for CGI users\n\
     --cgi_allow_anti_alias=BOOL          yes = anti-aliasing mode is allowed for CGI users\n\
     --cgi_min_Airy_disk_first_null=FLOAT Minimum allowed first null distance for CGI users\n\
     --cgi_max_Airy_disk_min_extent=NUM   Maximum allowed Airy disk minimum extent for CGI users\n\
     --cgi_max_Airy_disk_max_extent=NUM   Maximum allowed Airy disk extent for CGI users\n\
\n\
Star filters:\n\
     --enable_Gaia=BOOL                   yes = Enable galaxy-pq*.dat with Gaia stars\n\
     --Gaia_min_parallax_quality=NUM      Minimum parallax quality of Gaia stars (GDR3 'parallax_over_error')\n\
                                          Valid values: 0,1,2,3,5,10,20,30,50,100\n\
     --enable_external=BOOL               yes = Enable galaxy-external.dat with non-Gaia stars\n\
     --render_distance_min=FLOAT          Minimum star distance\n\
     --render_distance_max=FLOAT          Maximum star distance\n\
     --render_distance_selector=NUM       min/max star distance is measured from 0=camera, 1=target\n\
     --star_intensity_min=FLOAT           Minimum star intensity (Vega scale magnitude). Note: this applies after\n\
                                          extinction undimming if enabled\n\
     --star_intensity_max=FLOAT           Maximum star intensity (Vega scale magnitude). Note: this applies after\n\
                                          extinction undimming if enabled\n\
     --star_intensity_selector=FLOAT      Min/max star intensity is measured from 0=camera, 1=Earth, 2=10 parsecs\n\
     --star_color_min=FLOAT               Minimum star apparent color temperature in Kelvin\n\
     --star_color_max=FLOAT               Maximum star apparent color temperature in Kelvin\n\
\n\
Extinction:\n\
     --extinction_dimming_undo=BOOL       yes = undo extinction dimming (based on Gaia DR3 AG_GSPPHOT)\n\
     --extinction_reddening_undo=BOOL     yes = undo extinction reddening (based on Gaia DR3 TEFF_GSPPHOT)\n\
\n\
Camera:\n\
     --camera_res_x=NUM                   Horizontal resolution\n\
     --camera_res_y=NUM                   Vertical resolution\n\
     --camera_fov=FLOAT                   Field of vew in decimal degrees\n\
     --camera_pixel_limit_mag=FLOAT       Pixel exposure limit in Vega scale magnitude\n\
     --camera_pixel_limit_mode=NUM        How to handle overexposed pixels: 0=saturate to white, 1=preserve color\n\
     --camera_wb_enable=BOOL              yes = Enable white balance correction\n\
     --camera_wb_temp=FLOAT               White balance color temperature in Kelvin\n\
     --camera_color_saturation=FLOAT      Chroma saturation level (4.0 = 4x crhoma)\n\
     --camera_gamma=FLOAT                 Image gamma adjustment. This option never changes PNG header gamma as it\n\
                                          is intended to modify the way the image looks\n\
     --camera_projection=NUM              Raster projection: 0 = lat/lon, 1 = spherical, 2 = Hammer, 3 = Mollewide\n\
     --spherical_orientation=NUM          Spherical projection orientation: 0 = forward centered, 1 = forward on\n\
                                          left, rear on right\n\
     --Mollewide_iterations=NUM           Number of iterations for Mollewide projection algorithm\n\
\n\
Camera bandpass filters:\n\
     --red_filter_long_limit=FLOAT        Red channel passpand long wavelength limit in nm\n\
     --red_filter_short_limit=FLOAT       Red channel passband short wavelength limit in nm\n\
     --green_filter_long_limit=FLOAT      Green channel passband long wavelength limit in nm\n\
     --green_filter_short_limit=FLOAT     Green channel passband short wavelength limit in nm\n\
     --blue_filter_long_limit=FLOAT       Blue channel passband long wavelength limit in nm\n\
     --blue_filter_short_limit=FLOAT      Blue channel passband short wavelength limit in nm\n\
\n\
Diffraction:\n\
     --Airy_disk=BOOL                     yes = spread star flux with Airy disk pattern\n\
                                          no = star flux is mapped to exactly one output pixel\n\
     --Airy_disk_first_null=FLOAT         Radius to the first Airy disk null (green channel) in pixels\n\
                                          This sets the scale for the Airy disk patterns\n\
     --Airy_disk_max_extent=NUM           Maximum extent of Airy disk pattern in pixels\n\
                                          To minimize rendering time actual extent is autoscaled between min-max\n\
                                          for each star depending on how bright it is\n\
     --Airy_disk_min_extent=NUM           Minimum extent of Airy disk pattern in pixels. Minimum extent for all stars\n\
                                          Larger values can dramatically increase rendering time\n\
     --Airy_disk_obstruction=FLOAT        Aperture obstruction ratio (secondary mirror for example). Set to 0.0\n\
                                          for unobstructed aperture. Hubble = 0.127\n\
\n\
Anti-aliasing:\n\
     --anti_alias_enable=BOOL             yes = spread pixel intensity to neighboring pixels\n\
                                          no = pixel intensity is mapped to nearest pixel\n\
                                          This also applies to each Airy disk pixel\n\
     --anti_alias_radius=FLOAT            Radius of anti-aliasing spread in pixels. Valid range 0.5 - 2.0\n\
\n\
Skyglow:\n\
     --skyglow_enable=BOOL                yes = Enable skyglow effect\n\
     --skyglow_temp=FLOAT                 Effective temperature of skyglow in Kelvin\n\
     --skyglow_per_pixel_mag=FLOAT        Intensity of skyglow per output pixel in Vega scale magnitude\n\
\n\
Post-processing:\n\
     --Gaussian_blur_radius=NUM           Optional Gaussian blur with this radius in pixels\n\
     --output_scaling_factor=FLOAT        Optional output scaling using Lanczos2 interpolation\n\
\n\
Overlays:\n\
     --draw_crosshairs=BOOL               yes = Draw small crosshairs in center of image. Note: This will not be\n\
                                          centered on target if pan, tilt, or side-by-side spherical mode is used\n\
     --draw_grid_lines=BOOL               yes = Draw horizontal and vertical lines at 25%%, 50%%, and 75%% of\n\
                                          width and height\n\
\n\
Output:\n\
     --output_format=NUM                  0 = PNG 8-bit unsigned integer per color\n\
                                          1 = PNG 16-bit unsigned integer per color\n\
                                          2 = EXR 16-bit floating-point per color\n\
                                          3 = EXR 32-bit floating-point per color\n\
                                          4 = EXR 32-bit unsigned integer per color\n\
     --icc_profile=NUM                    -1 = default: PNG=sRGB, EXR=None (EXR clients assume Rec. 709 colorspace)\n\
                                          0 = None - linear gamma, 1 = sRGB, 2 = Display-P3 (compatible Z),\n\
                                          3 = Rec. 2020 (compatible Z), 4 = Rec. 601 NTSC, 5 = Rec. 601 PAL,\n\
                                          6 = Rec. 709, 7 = None - flat 2.0 gamma\n\
                                          ICC profiles are v4 from https://github.com/saucecontrol/Compact-ICC-Profiles \n\
                                          The EXR format does not support ICC profiles or encoding gamma. Instead, the\n\
                                          chromaticities associated with the selected profile is included in the header\n\
\n\
Camera position in Euclidian ICRS coordinates:\n\
     --camera_icrs_x=FLOAT                x coordinate in parsecs\n\
     --camera_icrs_y=FLOAT                y coordinate in parsecs\n\
     --camera_icrs_z=FLOAT                z coordinate in parsecs\n\
\n\
Camera position in spherical ICRS coordinates. These override Euclidian if not zero:\n\
     --camera_icrs_ra=FLOAT               Right ascension in decimal degrees\n\
     --camera_icrs_dec=FLOAT              Declination in decimal degrees\n\
     --camera_icrs_r=FLOAT                Distance in parsecs\n\
\n\
Camera target in Euclidian ICRS coordinates:\n\
     --target_icrs_x=FLOAT                x coordinate in parsecs\n\
     --target_icrs_y=FLOAT                y coordinate in parsecs\n\
     --target_icrs_z=FLOAT                z coordinate in parsecs\n\
\n\
Camera target in spherical ICRS coordinates. These override Euclidian if not zero:\n\
     --target_icrs_ra=FLOAT              Right ascension in decimal degrees\n\
     --target_icrs_dec=FLOAT             Declination in decimal degrees\n\
     --target_icrs_r=FLOAT               Distance in parsecs\n\
\n\
Optional camera rotation/pan/tilt after aiming at target:\n\
     --camera_rotation=FLOAT             Camera rotation once aimed at target in decimal degrees\n\
     --camera_pan=FLOAT                  Camera left-right pan once aimed at target and rotated in decimal degrees\n\
     --camera_tilt=FLOAT                 Camera up/down tilt once aimed at arget, rotated, and panned in decimal\n\
                                         degrees\n\
\n");
}

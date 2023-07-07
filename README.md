# Billion Star 3D Rendering Engine

bsrender is a 3D rendering engine for large star databases such as the ESA's mission Gaia DR3 data set with over a billion stars. It generates images from any position inside or outside the galaxy and can be run from the command line or from a web server in CGI mode. Direct rendering is used for every star without the use of pre-rendered frames or low-res previews.

  Resource |  Link    
---|---
Project website|[bsrender.io](https://bsrender.io/)
Live demo|[bsrender.io/demo/](https://bsrender.io/demo/)
Source code|[github.com/kevinloch/bsrender](https://github.com/kevinloch/bsrender)
Sample binary data files|[bsrender.io/sample_data/1.0-dev/](https://bsrender.io/sample_data/1.0-dev/)

## Key features

  - Generate images in JPG, PNG, AVIF, HEIF, or OpenEXR formats with multiple bit depth depending on format. HEIF is supported in cli mode only (not CGI/web interface)
  - HDR is supported on all formats using an experimental oepn source Rec. 2100 PQ ICC profile (PNG, JPG), or in-header signaling (EXR, AVIF, HEIF). Known to work with Chrome browser on M1/M2 Macbooks.
  - 3D translations/rotations of camera position/aiming using ICRS equitorial or Euclidian coordinates. Camera can be placed anywhere in the Universe
  - Customizable camera resolution, field of view, sensitivity, white balance, color saturation, and gamma
  - Several raster projection modes are supported: lat/lon (equirectangular), Spherical (forward hemisphere centered), Spherical (front/rear hemispheres), Hammer, Mollewide
  - Stars can be filtered by parallax quality (parallax over error), distance from camera or target, and apparent temperature
  - Camera color is modeled with a Planck spectrum for the apparent temperature of each star and customizable bandpass filters for each color channel
  - Apparent star temperature (color) is derived from Gaia bp/G and/or rp/G flux ratios for most stars
  - Multithreading support with customizable number of threads
  
## Optional features

- After aiming at target, separate controls are provided to pan and tilt away from target for maximum flexibility with aiming. Camera can also be rotated about it's view axis for desired orientation
- Support for user-supplied stars. A sample external database is provided with the Sun and all stars brighter than magnitude 3 that are not included in the Gaia dataset
- Airy disks provide photorealistic renderings of individual stars and clusters when enabled
- Gaussian blur and/or Lanczos output scaling. This allows high resolution renderings to be smoothed and downsampled on a server before downloading
- Anti-aliasing, helpful when combining multiple frames into videos or simulating DSLR/MILC images
- Skyglow for simulating views through Earth's atmosphere
- A sample html/javascript interface includes presets for a few camera targets and several common Hubble bandpass filter settings (along with typical LRGB). Also allows copy/paste settings URL for sharing links to your rendering settings
  
## Memory requirements

Using the full Gaia dataset with 1.4B stars requires at least 64GB of ram to run as fast as possible. This is for the operating system to cache the 46GB dataset in memory in addition to ram used by bsrender. Larger resolutions and/or use of blur or output scaling will increase memory requirements. Full 64-bit support allows for extremely large resolutions, limtied only by available ram and CPU time. 128000x64000 downsampled to 3200x16000 has been rendered with 512GB ram.

## Installation

This program is written in C and requires gcc, GNU make, libpng, libjpeg, libavif, libheif, and zlib to compile. You can disable compiling in specific output formats by commenting out '#define BSR_USE_<format>' in bsrender.h and removing the associated -l<library> flag from BSR_LIBS in Makefile.

On Linux or Mac w/Xcode, go to the 'src' directory and type:

    make

On FreeBSD and other systmes with non-GNU make as default:

    gmake

There is no 'make install' feature yet so you will have to manually copy the executables to /usr/local/bin or wherever you want them to live. For example (on some systems this my need to be done as root):

    cp bsrender /usr/local/bin; chmod 755 /usr/local/bin/bsrender
    cp mkgalaxy /usr/local/bin; chmod 755 /usr/local/bin/mkgalaxy
    cp mkexternal /usr/local/bin; chmod 755 /usr/local/bin/mkexternal
    cp ../scripts/getgalaxydata.sh /usr/local/bin; chmod 755 /usr/local/bin/getgalaxydata.sh
    cp ../scripts/gaia-dr3-extract.sh /usr/local/bin; chmod 755 /usr/local/bin/gaia-dr3-extract.sh

To show the version and all command line options run:

    bsrender --help

## Configuration file (optional)

Most options can also be set with a configuration file. By default bsrender looks for 'bsrender.cfg' in the same directory it is run from. A sample configuration file is provided in the source distribution. The -c command line option can be used to specify an alternate configuration filename/location. The compiled-in defaults are the same as those in the sample bsrender.cfg.

## Data files (required)

The Gaia archive .csv files are not suitable for direct 3D rendering. They must be processed into a binary data format that includes the 3D position, apparent temperature, and normalized intensity (relative to Vega at one parsec) of each star. The fastest and easiest way to use bsrender is to download pre-generated data files from [bsrender.io/sample_data/1.0-dev/](https://bsrender.io/sample_data/1.0-dev/) (46GB).

By default bsrender expects these files to be located in the subdirectory 'galaxydata' relative to where bsrender is run from. The -d command line option or data_file_directory config file option can be used to specify an alternate filename/location. To create the subdirectory and download the sample data files to it:

    mkdir galaxydata
    cd galaxydata
    getalaxydata.sh

Alternatively, you can create a symlink to the data file directory. For example if you put the binary data files in /data and you are currently in the directory you want to run bsrender from:

    ln -s /data galaxydata

### Generating a custom external data file (optional)

An "external" star database of user-supplied stars is supported. By default bsrender expects galaxy-external.dat to be in the data files directory but this can be disabled with the use_external_db=no configuration option. You can use a custom generated galaxy-external.dat with the sample Gaia data files downloaded from bsrender.io, it is not necessary to recreate the Gaia data files just to use a custom galaxy-external.dat. It is also possible to use only the external data file by setting use_external_db=yes and use_Gaia_db=no. In this mode bsrender can be used for any arbitrary star database instead of the Gaia dataset and there is no need to download or create the Gaia data files.

A sample external.csv source and binary galaxy-external.dat are provided in the source distribution and in the sample data files downloaded by getgalaxydata.sh. The sample file includes the Sun and all other stars that are too bright to be included in the Gaia dataset, plus a few fainter stars obscured by bright stars. To customize, edit external.csv to add/delete/modify any stars you want and then run mkexternal to generate galaxy-external.dat:

    mkexternal

Be sure to copy the new galaxy-external.dat to your data files directory.

### Generating the Gaia data files manually (optional)

Generating your own Gaia data files requires 805GB of disk space and approximately 24 hours of CPU time. You might want to do this if you want to use non-default settings of mkgalaxy, or you are using bsrender on a system architecture incompatible with the sample data files which are in little-endian (x86/arm) format.

Download the entire 'gaia_source' directory (644GB) from the Gaia archive website [http://cdn.gea.esac.esa.int/Gaia/gdr3/](http://cdn.gea.esac.esa.int/Gaia/gdr3/). Note that this version of bsrender will only work with gdr3. Do not uncompress the downloaded csv.gz files. 

The script 'gaia-edr3-extract.sh' will extract the columns bsrender uses from the compressed csv source files and create a single uncompressed 'gaia-edr3-extracted.csv' (115GB). From the directory above gaia_source containing the compressed csv files run the script:

    gaia-edr3-extract.sh

The utility 'mkgalaxy' is used to process gaia-edr3-extracted.csv into the binary data files used by bsrender. There are a few command line options for mkgalaxy that can be shown with:

    mkgalaxy --help

To generate the Gaia binary data files run mkgalaxy in the same directory as gaia-edr3-extracted.csv:

    mkgalaxy

This may take up to 24 hours to complete, depending on system and disk speed. Be sure to copy the new files to your data files directory. Note that binary data files created on similar but different systems may not be identical due to different non-significant bits of floating point values. This has no effect on the precision or operation of bsrender.

## Operation

  By default it will render a 360 degree lat/lon (equirectangular) projected panorama of the entire sky from the sun with the default settings in the sample bsrender.cfg. Descriptions of configuration options and their function are provided in the sample configuration file. Options are set in the following sequence:

  1. Compiled-in defaults
  2. Command line options -c and -h
  3. Configuration file options
  4. Other command line options
  5. CGI options from environment variable QUERY_STRING, if in CGI mode

  Some options are privileged and cannot be set by CGI users.

### Helpful hints

  - Reducing the 'camera\_fov' (zooming in) will generally require increasing 'camera\_pixel\_limit\_mag' (pixel intensity limit in the web interface) which makes camera more sensitive to maintain the same subjective iamge brightness. This is because dense star fields aggregate to brighter individual pixels with wider field of view. For very narrow fields of view with individual stars, enabling Airy disks is Highly recommended. Otherwise the individual star pixels can be very hard to see and increasing 'camera\_pixel\_limit\_mag' may just saturate those pixels without increasing subjective brightness.
  - Similarly, increasing the camera resolution will generally require increasing 'camera\_pixel\_limit\_mag' to maintain the same subjective image brightness. Use caution with increasing 'camera\_pixel\_limit\_mag' too high with very high resolutions and/or narrow fields of view. Colors will desaturate as pixel intensity is saturated unless 'camera\_pixel\_limit\_mode' is set to 1 (preserve color) and even then unnatural colors will result. The key is to remain aware of when stars start to map to individual pixels and the approximate magnitude of those stars. Enabling Airy disks provides significant freedom to "overexpose" pixels as overexposed stars will appear larger and still preserve some of their color in the outer parts of the Airy disk.
  - Rendering time depends on many factors. It is essential that there is enough ram for the operating system to cache the entire binary dataset. Enabling airy disks has minimal impact on rendering time unless there are a large number of highly overexposed stars or with a large setting for 'Airy\_disk\_min\_extent'. Wider fields of view contain more stars and take longer to render. Very large image resolutions take longer, mainly due to the time spent initializing and processing the image buffers, but also in image generation. Optional Gaussian blur and Lanczos2 resizing add minimal time but are also slower at larger resolutions.
  - When resizing with Lanczos2 resampling, best results are obtained by also using Gaussing blur at 1/4 the downscaling factor. If reducing by 2x, set blur radius to 0.5. if reducing by 8x set blur radius to 2.0 etc.
  - Star 'temperature' is apparent temperature not actual star temperature, except for supplemental stars in he external.csv dataset. This apparent temperature corresponds to a Planck blackbody spectrum that is the closest fit to the Gaia rp, bp and G flux data. Despite ignoring the distortion of stellar spectra by extinction this produces amazingly accurate star colors, often indistinguishable from Hubble photographs when Airy disks are enabled and the correct simulated Hubble passband filters are selected.
  - Due to uncertainty in the parallax data of approximately 20 microarcseconds, things start to look weird as the camera is positioned more than a short distance away from the sun. This is a limitation of the source data and not any bug or problem with the rendering engine. If override parallax is enabled in mkgalaxy (by setting -p > 0), there will be a spherical shell of residual stars at 1000 / minimum\_parallax parsecs from the Sun. This is of course artificial but is better than having some stars (like LMC and SMC) much farther away from the galaxy than they really are. The sample data files were generated with a 20 microarcsecond minimum parallax enforced and a 50 kpc artifical shell of distance-limited stars.
  - Color profiles tell an image viewer information about how the image was encoded (color space, gamma, etc.). If a viewer ignores the color profile it will most likely assume it was encoded with the sRGB color space and gamma. For this reason the sRGB profile is the safest and most compatible profile to use. Note that while bsrender applies the encoding gamma specified in the selected standard, it does not otherwise change the colors saved to the output image. This is because the configurable camera bandpass filters do not necessarily repersent human vision so color calibration beyond white balance is purely subjective. On a color managed viewer a wide-gamut profile like Rec. 2020 will render more highly saturated colors for the same RGB values than a narrow-gamut profile like sRGB. Some of the Hubble and the LRGB camera bandpass presets in sample-frontend.html will give natural looking colors with the sRGB profile. Presets based on the IEC 1931 standard observer RGB color matching functions (representing human vision) give natural looking colors with the Rec. 2020 profile. Of course false or oversaturated colors are sometimes desirable and overall color saturation can be adjusted with any profile.
 - When generating images for use with ffmpeg to make videos, a flat 2.0 encoding gamma should be used due to the way ffmpeg handles image import.

### CGI mode

when cgi\_mode=yes is set in the config file html headers and png data will be output to stdout, with all other output suppressed (unless run with -h).
  CGI requests should be made with http GET requests using the same key/value pairs as in the config file. Some options (data\_file\_directory, num\_threads, per\_thread\_buffer, cgi\_) cannot be overridden via CGI and some are limited by the cgi\_ options in config file.

## Methodology

Gaia source data is pre-processed with 'gaia-edr3-extract.sh' and then 'mkgalaxy' to tranform the relevant source data feilds into the most efficient form for direct renderng. Spherical ICRS coordinates ('ra', 'dec', and r derived from 'parallax') are transformed into double precision Euclidian x,y,z coordinates. The star's linear intensity (relative to Vega at 1pc) is derived from 'phot\_G\_mean\_flux'. The apparent star color temperature is derived by finding the best match (to the closest integer Kelvin) for bp/G and/or rp/G flux ratios to a Planck spectrum integrated within the Gaia rp, bp, and G passbands. If reliable bp and rp flux are not available the color wavenumber ('nu\_eff\_used\_in\_astrometry' or 'pseudocolor') is treated as the peak wavelength of a Planck spectrum and converted into an apparent color temperature. These five derived fields (x, y, z, color\_temperature, linear\_1pc\_intensity) are encoded into binary data files for use by 'bsrender'. The Gaia source 'parralax\_over\_error' value is used to split stars into 10 data files by "parallax quality".

The rendering engine 'bsrender' uses these data files to generate an image file (or stream in CGI mode). Extensive options for configuring rendering are provided through a configuration file 'bsrender.cfg' and most of those can also be set via CGI GET request.

Before rendering begins an rgb table is initialized using the 'camera\_wb\_temp', 'camera\_color\_saturation', and camera color channel passband options. A Planck spectrum is simulated for the white balance temperature to generate white balance factors for each color channel. Then r,g,b values (normalized to the integrated wide band flux) are calculated for each temperature between 0-32767K. Later during rendering the r,g,b values for a star's apparent temperature is multiplied by the linear star intensity (adjusted for distance) to generate the star's contribution to a pixel or Airy disk map of pixels. If Airy disks are enabled then a pre-computed map of Airy disk pixel factors is generated. The first null of the Airy disk in the green channel is scaled to match 'Airy\_disk\_first\_null'. The included Bessel function table supports an Airy\_disk\_max\_extent of 1000 pixels which is over 3000 orders of diffraction with Airy\_disk\_first\_null set to the minimum of 0.3.

Stars can be filtered by parallax quality, distance from either camera or target, and apparent color temperature.

The coordinate system used internally by bsrender is Euclidian x,y,z with equitorial orientation. From the camera's perspective +x=forward, +y=left, and +z=up. Quaternion algebra is used for 3D rotations of stars which provides maximum speed, consistent precision, and avoids gimbal lock. Stars are first rotated by the (xy and xz) angles required to bring the target to the center of camera view. Optional camera rotation (yz), pan (xy) and tilt (xz) can then be applied in that order. All rotations are combined during initialization into a single rotation quaternion which is used to rotate each selected star in a single rotation operation during processing.

After translation and rotation stars are filtered by field of view and mapped to an image composition buffer pixel by the selected raster projection. A star's linear intensity (adjusted for distance) is multiplied by the r,g,b lookup table for the star's apparent color temperature. If Airy disks are enabled then the pre-computed Airy disk map is used to generate additional pixels up to 'Airy\_disk\_max\_radius' around the central pixel and the star's intensity\*(r,g,b) values are multiplied by the Airy map factor for each Airy disk pixel. Output pixels can optionally be anti-aliased to simulate common consumer/DSLR sensors. Pixels are stored in a double-precision floating-point image composition buffer where they are added to any pixels from previous stars at the same location.

After the image composition buffer is complete post-processing involves several steps all performed in double precision floating point format:
    
  - Normalizing pixel intensity where 'camera\_pixel\_limit\_mag' = 1.0
  - Optionally applying 'camera\_gamma'
  - Limiting maximum values for any r,g,b channel to 1.0 by either saturating to white or preserving color (see 'camera\_pixel\_limit\_mode')
  - Optionally applying Gaussian blur (see 'Gaussian\_blur\_radius');
  - Optionally applying Lanczos2 image resizing (see 'output\_scaling\_factor')
  - Limiting maximum values to 1.0 again
  - Applying appropriate encoding gamma if a color profile is selected

Finally, the resulting image is converted to the selected image format and output to file or stream, with optional color profile information. Note: other than optionally applying encoding gamma, RGB values are not changed (rebalanced) based on a selected color profile. 

## Data Credit

This work has made use of data from the European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.

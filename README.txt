
 Billion Star 3D Rendering Engine
 Kevin M. Loch

  bsrender is a 3D rendering engine for large star databases such as the ESA's Gaia EDR3 data set with over a billion stars.  It generates PNG images from any position inside or outside the galaxy and can be run from the command line or cgi mode.  Direct rendering is used for every star without the use of pre-rendered frames or low-res previews.

  A live demonstration: https://bsrender.io/demo/
  Sample renderings and documentation: https://bsrender.io/
  Sample binary galaxy data files: https://bsrender.io/sample_data/

Key features:

  - 3D translations and rotations of camera position and aiming using ICRS spherical or Euclidian coordinates
  - Customizable camera resolution, field of view, sensitivity, white balance, color saturation, and gamma
  - Camera color is accurately modeled with Planck spectrum and customizable bandpass filters for each color channel
  - Several raster projection modes are supported: lat/lon (equirectangular), Spherical (forward hemisphere centered), Spherical (front/rear hemispheres), Hammer, Mollewide
  - Optional Airy disks dramatically improve the appearance of fields with individual stars and clusters
  - Stars can be selected by parallax quality (parallax over error), distance from camera or target, and effective color temperature
  - Optional Gaussian blur and/or Lanczos2 output scaling.  This allows very high resolution renderings to be smoothed and downsampled on a server before downloading
  - Good performance when the data files can be cached in ram.  On an AWS c6gd.12xlarge instance (48vcpu/96Gram) the full pq000 dataset (1.4B stars) can be rendered in 14 seconds at 2k resolution.  The default pq010 data set (98M stars) renders in just 2 seconds
  - Support for user-supplied stars.  This can be used to add stars that are too bright or dim to have their parallax measured by the Gaia satellite.  A sample external.csv is provided with the Sun and all stars brighter than magnitude 3 that are not included in the Gaia dataset or are not able to be imported into bsrender because they lack parallax and/or G-band flux
  - Effective star temperature (color) is derived from Gaia bp/G and/or rp/G flux ratios
  - 8 or 16 bits per color PNG output, with or without sRGB encoding gamma for either bit depth
  - Sample web interface includes presets for a few camera targets and several common Hubble bandpass filter settings (along with typical LRGB)
  
Installation:

  - Requires gcc and GNU make
  - Requires libpng (static library)
  - Compile with make (make on Linux and Mac w/Xcode, gmake on FreeBSD).  There is no 'configure' script yet
  - To extract and process your own data files, make sure the ./gaia_source contains the Gaia EDR3 compresed csv files
    - run 'gaia-edr3-extract.sh'
    - run 'mkgalaxy' in the same directory as gaia-edr3-extracted.csv to generate the data files for Gaia stars
    - edit 'external/external.csv' for user-supplied stars
    - run 'mkexternal' in the same directory as external.csv to generate the data files for user-supplied stars
  - To download sample data files (43GB total size) run the script 'getgalaxydata.sh'.  It will download the data files to the current directory
  - Move galaxy-pq*.dat files to ./galaxydata or path specified in config file or -d comand line option
  - Edit options in bsrender.cfg for desired rendering settings.  Make sure bsrender.cfg is in current directory or use -c command line option to specify name and location
  - run 'bsrender', it will print progress messages and output an image to 'galaxy.png'
  - To use in cgi mode set 'cgi_mode' in bsrender.cfg.  All status messages are suppressed, and an html header followed by the png image is printed to stdout
  - A sample html page to call bsrender in cgi mode is provided

Operation:

  By default it will render a 360 degree lat/lon (equirectangular) projected panorama of the entire sky from the sun with the following camera settings:
    Resolution: 2000x1000
    Parallax quality: 10 [98M stars]
    White balance: 4300K
    Color saturation: 4.0
    Pixel saturation magnitude: 8.0
    Raster projection: lat/lon (equirectangular)
    And other default settings identical to the sample bsrender.cfg

  Descriptions of configuration options and their function are provided in the sample configuration file.  Options are set in the following sequence:
    - Compiled-in defaults
    - Command line options
    - Configuration file options
    - CGI options
  Some options are privileged and cannot be set by CGI users.

Notes:

  - Due to uncertainty in the parallax data of approximately 20 microarcseconds, things start to look weird as the camera is positioned a short distance away from the sun (more than 10 parsecs).  This is a limitation of the source data and not any bug or problem with the rendering engine.  If override parallax is enabled in mkgalaxy (by setting -p > 0), there will be a spherical shell of residual stars at 1000 / minimum_parallax parsecs from the Sun.  This is of course artificial but is better than having some stars (like LMC and SMC) much farther away from the galaxy than they really are.  The sample data files were generated with a 20 microarcsecond minimum parallax enforced and a 50 kpc artifical shell of distance-limited stars.
  - Rendering with Airy disks enabled can be extremely slow.   Extra processing increases with the square of 'Airy_disk_max_extent' but is also increased by small values of 'Airy_disk_first_null' less than 1.0.  Additional cpus/rendering threads may help somewhat but the extra Airy disk pixels can still saturate memory bandwidth and the main thread.
  - Otherwise, performance strongly depends on whether the binary data files can be completely cached in memory by the operating system. Multiple rendering threads are supported and performance typically limited by memory bandwidth on systems with more than 48 cores that have enough ram to cache all of the data files.

CGI mode:

  when cgi_mode=yes is set in the config file html headers and png data will be output to stdout, with all other output suppressed (unless run with -h).
  Cgi requests should be made with http GET requests using the same key/value pairs as in the config file. Some options (data_file_directory, num_threads, per_thread_buffer, cgi_) cannot be overridden via CGI and some are limited by the cgi_ options in config file.

Methodology:

  Gaia source data is pre-processed with 'gaia-edr3-extract.sh' and then 'mkgalaxy' to tranform the relevant source data feilds into the most efficient form for direct renderng. Spherical ICRS coordinates ('ra', 'dec', and r derived from 'parallax') are transformed into double precision Euclidian x,y,z coordinates. The star's linear intensity (relative to Vega at 1pc) is derived from 'phot_G_mean_flux'.  The effective star color temperature is derived by finding the best match (to the closest integer Kelvin) for bp/G and/or rp/G flux ratios to a Planck spectrum integrated within the Gaia rp, bp, and G passbands.  If reliable bp and rp flux are not available the color wavenumber ('nu_eff_used_in_astrometry' or 'pseudocolor') is treated as the peak wavelength of a Planck spectrum and converted into an effective color temperature.  These five derived fields (x, y, z, color_temperature, linear_1pc_intensity) are encoded into binary data files for use by 'bsrender'.  The Gaia source 'parralax_over_error' value is used to split stars into 10 data files by "parallax quality".

  The rendering engine 'bsrender' uses these data files to generate a PNG file (or stream in CGI mode).  Extensive options for configuring rendering are provided through a configuration file 'bsrender.cfg' and most of those can also be set via CGI GET request.  

  Before rendering begins an rgb table is initialized using the 'camera_wb_temp', 'camera_color_saturation', and camera color channel passband options.  A Planck spectrum is simulated for the white balance temperature to generate white balance factors for each color channel.  Then r,g,b values (normalized to the integrated wide band flux) are calculated for each temperature between 0-32767K.   Later during rendering the r,g,b values for a star's effective temperature is multiplied by the linear star intensity (adjusted for distance) to generate the star's contribution to a pixel or Airy disk map of pixels.  If Airy disks are enabled then a pre-computed map of Airy disk pixel factors is generated. The first null of the Airy disk in the green channel is scaled to match 'Airy_disk_first_null'.  The included Bessel function table supports almost 1000 orders of diffraction and the maximum Airy disk total size is set by 'Airy_disk_max_extent'.

  Stars can be filtered by parallax quality, distance from either camera or target, and effective color temperature.

  3D rotations for the camera target and stars are done with a 'triple-azimuth' method.  Three full 360 degree orthoganal angles (x,y), (x,z), and (y,z) provide redundancy for constant precision and avoiding gimbal lock.  After camera and star translation, the stars are rotated so the camera target is at x=0,z=0,r=r (if pan/tilt options are used then stars are further rotated by these angles). This allows for simple and fast mapping to verious raster projections (especially lat/lon).

  After translation and rotation stars are filtered by field of view and mapped to an image composition buffer pixel by the selected raster projection.  A star's linear intensity (adjusted for distance) is multiplied by the r,g,b lookup table for the star's effective color temperature.   If Airy disks are enabled then the pre-computed Airy disk map is used to generate additional pixels up to 'Airy_disk_max_radius' around the central pixel and the star's intensity*(r,g,b) values are multiplied by the Airy map factor for each Airy disk pixel.  These pixel(s) are added to any existing values for those pixel(s) and stored in the image composition buffer. The image composition buffer uses double precision floating point format so it can accurately handle many stars (or Airy map pixels) at the same location and to provide maximum flexibility for the post-processing steps.

  After the image composition buffer is complete post-processing involves several steps all performed in double precision floating point format:
    - Normalizing pixel intensity where 'camera_pixel_limit_mag' = 1.0
    - Optionally applying 'camera_gamma'
    - Limiting maximum values for any r,g,b channel to 1.0 by either saturating to white or preserving color (see 'camera_pixel_limit_mode')
    - Optionally applying Gaussian blur (see 'Gaussian_blur_radius');
    - Optionally applying Lanczos2 image resizing (see 'output_scaling_factor')
    - Limiting maximum values to 1.0 again
    - Optionally applying sRGB encoding gamma

  Finally, the resulting image is converted to 8 or 16 bits per color and output to a PNG file or stream.

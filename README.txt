
 Billion Star 3D Rendering Engine
 Kevin M. Loch

  bsrender is a 3D rendering engine for the ESA's Gaia EDR3 data set of over a billion stars (with parallax data).  It generates PNG images from any position inside or outside the galaxy.  It can be run on the command line with .png file output or in cgi mode with html headers and png data written to stdout.

  Sample binary galaxy data files are available at https://kevinloch.com/bsrender/sample_data/
  Sample renderings are available at https://kevinloch.com/bsrender/sample_renderings/
  A live demonstration is available at https://kevinloch.com/galaxy/

Key features:

  - 3D translations/rotations of camera position and aiming using ICRS spherical or Euclidian coordinates
  - Customizable camera resolution, field of view, sensitivity, white balance, color saturation, and gamma
  - Color is accurately modeled with Planck spectrum and customizable bandpass filters for each color channel
  - Raster projection modes supported: lat/lon (equirectangular), Spherical (forward hemisphere centered), Spherical (front/rear hemispheres), Hammer, Mollewide
  - Optional accurately modeled and adjustable Airy disks dramatically improve the appearance of fields with individual stars and clusters
  - Stars can be selected by parallax quality (parallax over error), distance from camera or target, and color temperature
  - Optional Gaussian blur and/or Lanczos2 output scaling.  This allows very high resolution renderings to be smoothed and downsampled without using 3rd party image manipulation tools
  - Good performance when the data files can be cached in ram.  On an AWS c6gd.12xlarge instance (48vcpu/96Gram) the full pq0 dataset (1.4B stars) can be rendered in 14 seconds at 2k resolution.  The pq10 set (98M stars) renders in 2 seconds.
  - Support for user-supplied stars.  This can be used to add the Sun and other stars that are too bright or dim to have their parallax measured by the Gaia satellite.
  
Installation:

  - Requires gcc and GNU make
  - Requires libpng (static library)
  - Compile with make (make on Linux and Mac w/Xcode, gmake on FreeBSD).  There is no 'configure' script yet
  - To extract and process your own data files, make sure the ./gaia_source contains the Gaia EDR3 compresed csv files
    - run 'gaia-edr3-extract.sh'
    - run 'mkgalaxy' to generate the data files for Gaia stars
    - edit 'external.csv' for user-supplied stars
    - run 'mkexternal' to generate the data files for user-supplied stars
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

  - Due to uncertainty in the parallax data of approximately 15 microarcseconds, things start to look weird as the camera is positioned a short distance away from the sun (more than 10 parsecs).  This is a limitation of the source data and not any bug or problem with the rendering engine.  If override_parallax_toolow is set to 1 in mkgalaxy.c, there will be a spherical shell of residual stars at 1000 / minimum_parallax parsecs from the Sun.  This is of course artificial but is better than having some stars (like LMC and SMC) much farther away from the galaxy than they really are.  The sample data files were generated with a 15 microarcsecond minimum parallax enforced.
  - Rendering with Airy disks enabled can be extremely slow.   Extra processing increases with the square of 'Airy_disk_max_extent' but is also increased by small values of 'Airy_disk_first_null' less than 1.0.  Additional cpus/rendering threads may help somewhat but the extra Airy disk pixels can still saturate memory bandwidth and the main thread.
  - Rendering with Gaussian blur and/or Lanczos output scaling can be slow for very large image sizes and/or blur radius.  These operations are done by the main thread so they are limited by single cpu performance.
  - Otherwise, performance strongly depends on whether the binary data files can be completely cached in memory by the operating system. Multiple rendering threads are supported and performance typically limited by memory bandwidth on systems with more than 48 cores that have enough ram to cache all of the data files.

CGI mode:

  when cgi_mode=yes is set in the config file html headers and png data will be output to stdout, with all other output suppressed (unless run with -h).
  Cgi requests should be made with http GET requests using the same key/value pairs as in the config file. Some options (data_file_directory, num_threads, per_thread_buffer, cgi_) cannot be overridden via CGI and some are limited by the cgi_ options in config file.

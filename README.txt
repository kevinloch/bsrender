
 Billion Star 3D Rendering Engine
 Kevin M. Loch

  bsrender is a 3d rendering engine for the ESA's Gaia EDR3 data set of over a billion stars (with parallax data).  It generates PNG images from any position inside or outside the galaxy.

Installation:

- Requires libpng (static library)
- Compile with gnu make (make on Linux and Mac w/Xcode, gmake on FreeBSD).
- To extract and process your own data files, make sure the ./gaia_source contains the Gaia EDR3 compresed csv files
  - run 'gaia-edr3-extract.sh'
  - run 'mkgalaxy'
- Move galaxy-pq*.dat files to ./galaxydata or path specified in config file or -d comand line option
- Edit options in bsrender.cfg for desired rendering settings.  Make sure bsrender.cfg is in current directory or use -c command line option to specify name and location
- run './bsrender', it will print a few progress messages and output an image to 'galaxy.png'
- To use in cgi mode, set 'cgi_mode' in bsrender.cfg.  Status messages are suppressed, an html header followed by the png image is printed to stdout

By default it will render a 360 degree lat/lon (equirectangular) projected panorama of the entire sky from the sun with the following camera settings:
 Resolution: 1920x1080
 Parallax quality: 10
 White balance: 4200K
 Color saturation: 4.0
 Pixel saturation magnitude: 7.0
 Camera gamma: 1.0
 sRGB encode gamma: yes
 Raster projection: lat/lon (equirectangular)

Due to uncertainty in the parallax data of approximately 15 microarcseconds, things start to look weird as the camera is positioned a short distance away from the sun (more than 10 parsecs).  This is a limitation of the source data and not any bug or problem with the rendering engine.  If override_parallax_toolow is set to 1 in mkgalaxy.c, there will be a spherical shell of residual stars at 1000 / minimum_parallax parsecs from the Sun.  This is of course artificial but is better than having some stars (like LMC and SMC) much farther away from the galaxy than they really are.

Sample binary data files are available at https://kevinloch.com/bsrender/sample_data/
Sample renderings are available at https://kevinloch.com/bsrender/sample_renderings/

Performance strongly depends on whether the binary data files can be completely cached in memory by the operating system.
Multiple rendering threads are supported and performance typically limited by memory bandwidth on systems with more than 48 cores that have enough ram to cache all of the data files.

The following raster projection modes are supported:

lat/lon (equirectangular)
Spherical - forward hemisphere centered
Spherical - forward hemisphere on left, rear hemisphere on right
Hammer
Mollewide

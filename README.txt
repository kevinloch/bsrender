
 Billion Star 3D Rendering Engine
 Kevin M. Loch

  bsrender is a 3d rendering engine for the ESA's Gaia EDR3 data set of over a billion stars (with parallax data).  It generates PNG images from any position inside or outside the galaxy.

Installation:

- Requires libpng (static)
- Compile with gnu make.
- Edit options in bsrender.cfg for desired camera settings.  Make sure bsrender.cfg is in current directory or use -c command line option to specify name and location
- To extract and process your own data files, make sure the ./gaia_source contains the Gaia EDR3 compresed csv files
  - run 'gaia-edr3-extract.sh'
  - run 'mkgalaxy'
- Move galaxy-pq*.dat files to ./galaxydata or path specified in config file / -d comand line option
- run './bsrender', it will print a few progress messages and output an image to 'galaxy.png'

By default it will render a 360 degree lat/lon (equirectangular) projected panorama of the entire sky from the sun with the following camera settings:
 Resolution: 1920x1080
 White balance: 4200K
 Color saturation: 4.0
 Pixel saturation magnitude: 7.0

Due to uncertainty in the parallax data of approximately 15 microarcseconds, things start to look weird as the camera is positioned a short distance away from the sun (more than 10 parsecs).  This is a limitation of the source data and not any bug or problem with the rendering engine.

Sample binary data files are available at https://kevinloch.com/bsrender/sample_data/
Sample renderings are available at https://kevinloch.com/bsrender/sample_renderings/

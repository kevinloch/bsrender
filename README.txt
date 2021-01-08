
 Billion Star 3D Rendering Engine Proof of Concept
 Kevin M. Loch

 bsrender is a 3d rendering engine for the ESA's Gaia EDR3 data set of over a billion stars (with parallax data).  It generates PPM images (which can be
converted to PNG or other formats with the netpbm collection of utilities) from any position inside or outside the galaxy.   This is an early proof of conecpt
release which is missing many usability features but is able to generate beautiful and scientifically useful images very quickly, especially if the input data
file is cached in memory by the operating system.

Installation:

- Compile with gnu make.
- make sure the gaia_source containing the Gaia EDR3 compresed csv files is in (or symlinked to) the current directory
- run 'gaia-edr3-extract.sh'
- Create a 'galaxy.dat' file by concatenating one or more of the galazy-pqxxx.dat files.  For example to generate a galaxy.dat file with a minimum parallax quality
  of 10: 'cat galaxy-pq100.dat galaxy-pq050.dat galaxy-pq030.dat galaxy-pq020.dat galaxy-pq010.dat > galaxy.dat'
- Edit options in bsrender.c for desired camera settings and recompaile with gnu make
- run 'bsrender > output.ppm', then convert ppm file to png with 'pnmtopng output.ppm > output.png'

By default it will render a 360 degree cylindrically projected panorama of the entire sky from the sun with the following camera settings:
 Resolution: 4096x2048
 White balance: 4100K
 Color saturation: 4.0
 Pixel saturation magnitude: 8.2

Due to uncertainty in the parallax data of approximately 15 microarcseconds, things start to look weird as the camera is positioned a short distance away from the sun (more than 10 parsecs).  This is a limitation of the source data and not any bug or problem with the rendering engine.

Sample binary data files are available from https://kevinloch.com/bsrender/sample_data/poc/

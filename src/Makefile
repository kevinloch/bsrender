CC = gcc
CFLAGS = -I. -I/usr/local/include -Wall -Ofast
# if compile failes with -Ofast, use -O3 instead
#CFLAGS = -I. -I/usr/local/include -Wall -O3

# to compile without support for specific output formats, comment out BSR_USE_<format> in bsrender.h and remove -l<library> from BSR_LIBS below:
# PNG: -lpng
# EXR: -lz
# JPEG: -ljpeg
# AVIF: -lavif
# HEIF: -lheif
BSR_LIBS = -L/usr/local/lib -L/usr/lib -L/usr/lib64 -L/usr/local/lib64 -pthread -lm -lpng -lz -ljpeg -lavif -lheif

LIBS = -L/usr/local/lib -lm
BSR_OBJ = sequence-pixels.o file.o memory.o image-composition.o Gaia-passbands.o Lanczos.o post-process.o Gaussian-blur.o rgb.o diffraction.o cgi.o init-state.o process-stars.o overlay.o icc-profiles.o bsr-png.o bsr-exr.o bsr-jpeg.o bsr-avif.o bsr-heif.o usage.o util.o bsr-config.o bsrender.o
BSR_DEPS = sequence-pixels.h file.h memory.h image-composition.h Gaia-passbands.h Lanczos.h post-process.h Gaussian-blur.h rgb.h diffraction.h cgi.h init-state.h process-stars.h overlay.h icc-profiles.h bsr-png.h bsr-exr.h bsr-jpeg.h bsr-avif.h bsr-heif.h usage.h util.h bsr-config.h bsrender.h Bessel.h Gaia-DR3-transmissivity.h
MKGALAXY_OBJ = util.o Gaia-passbands.o bandpass-ratio.o mkgalaxy.o
MKGALAXY_DEPS = util.h Gaia-passbands.h bandpass-ratio.h Gaia-DR3-transmissivity.h
MKEXTERNAL_OBJ = util.o mkexternal.o
MKEXTERNAL_DEPS = util.h
MKBESSEL_OBJ = mkBessel.o
MKBESSEL_DEPS = Bessel.h

.PHONY: all clean

all: mkBessel mkgalaxy mkexternal bsrender

clean:
	rm -f mkBessel mkgalaxy mkexternal bsrender *.o

$(BSR_OBJ): %.o : %.c $(BSR_DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(MKGALAXY_OBJ): %.o : %.c $(MKGALAXY_DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(MKEXTERNAL_OBJ): %.o : %.c $(MKEXTERNAL_DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(MKBESSEL_OBJ): %.o : %.c $(MKBESSEL_DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

mkBessel: $(MKBESSEL_OBJ)
	$(CC) $(CFLAGS) -o mkBessel $^ $(LIBS)

mkgalaxy: $(MKGALAXY_OBJ)
	$(CC) $(CFLAGS) -o mkgalaxy $^ $(LIBS)

mkexternal: $(MKEXTERNAL_OBJ)
	$(CC) $(CFLAGS) -o mkexternal $^ $(LIBS)

bsrender: $(BSR_OBJ)
	$(CC) $(CFLAGS) $(BSR_LIBS) -o bsrender $^ $(BSR_LIBS)

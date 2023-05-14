CC = gcc
CFLAGS = -I. -I/usr/local/include -Wall -Ofast
# if compile failes with -Ofast, use -O3 instead
#CFLAGS = -I. -I/usr/local/include -Wall -O3
BSR_LIBS = -lm -lz -lpng -L/usr/local/lib -L/usr/lib
LIBS = -lm -L/usr/local/lib
BSR_OBJ = sequence-pixels.o file.o memory.o image-composition.o Gaia-passbands.o Lanczos.o post-process.o Gaussian-blur.o rgb.o diffraction.o cgi.o init-state.o process-stars.o overlay.o icc-profiles.o bsr-png.o bsr-exr.o usage.o util.o bsr-config.o bsrender.o
BSR_DEPS = sequence-pixels.h file.h memory.h image-composition.h Gaia-passbands.h Lanczos.h post-process.h Gaussian-blur.h rgb.h diffraction.h cgi.h init-state.h process-stars.h overlay.h icc-profiles.h bsr-png.h bsr-exr.h usage.h util.h bsr-config.h bsrender.h Bessel.h Gaia-DR3-transmissivity.h
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
	$(CC) $(CFLAGS) $(LIBS) -o mkBessel $^

mkgalaxy: $(MKGALAXY_OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o mkgalaxy $^

mkexternal: $(MKEXTERNAL_OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o mkexternal $^

bsrender: $(BSR_OBJ)
	$(CC) $(CFLAGS) $(BSR_LIBS) -o bsrender $^

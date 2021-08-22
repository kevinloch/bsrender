CC = gcc
CFLAGS = -I. -I/usr/local/include -Wall -Ofast
BSR_LIBS = -lm -lpng -L/usr/local/lib
LIBS = -lm -L/usr/local/lib
BSR_OBJ = quantize.o memory.o image-composition.o Gaia-passbands.o Lanczos.o post-process.o Gaussian-blur.o rgb.o diffraction.o cgi.o init-state.o process-stars.o overlay.o bsr-png.o usage.o util.o bsr-config.o bsrender.o
BSR_DEPS = quantize.h memory.h image-composition.h Gaia-passbands.h Lanczos.h post-process.h Gaussian-blur.h rgb.h diffraction.h cgi.h init-state.h process-stars.h overlay.h bsr-png.h usage.h util.h bsr-config.h bsrender.h Bessel.h Gaia-EDR3-transmissivity.h
MKGALAXY_OBJ = Gaia-passbands.o bandpass-ratio.o mkgalaxy.o
MKGALAXY_DEPS = Gaia-passbands.h bandpass-ratio.h Gaia-EDR3-transmissivity.h
MKEXTERNAL_OBJ = mkexternal.o
MKEXTERNAL_DEPS =
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

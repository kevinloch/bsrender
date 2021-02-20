CC = gcc
CFLAGS = -I. -I/usr/local/include -Wall -Ofast
LIBS = -lm -lpng -L/usr/local/lib
OBJ = rgb.o diffraction.o cgi.o init-state.o process-stars.o overlay.o bsr-png.o usage.o util.o bsr-config.o bsrender.o
DEPS = rgb.h diffraction.h cgi.h init-state.h process-stars.h overlay.h bsr-png.h usage.h util.h bsr-config.h bsrender.h Bessel.h

.PHONY: all clean

all: mkgalaxy bsrender

clean:
	rm -f mkgalaxy bsrender *.o

$(OBJ): %.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

mkgalaxy: mkgalaxy.o
	$(CC) $(CFLAGS) $(LIBS) -o mkgalaxy $^

bsrender: $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o bsrender $^

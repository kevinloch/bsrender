CC = gcc
CFLAGS = -I. -I/usr/local/include -Wall -Ofast
LIBS = -lm -lpng -L/usr/local/lib
OBJ = cgi.o init-state.o process-stars.o overlay.o bsr-png.o rgb.o usage.o util.o bsr-config.o bsrender.o
DEPS = cgi.h init-state.h process-stars.h overlay.h bsr-png.h rgb.h usage.h util.h bsr-config.h bsrender.h

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

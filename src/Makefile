CC = gcc
CFLAGS = -I. -I/usr/local/include -Wall -Ofast
LIBS = -lm -lpng -L/usr/local/lib
OBJ = mkgalaxy.o bsrender.o
DEPS = bsrender.h

.PHONY: all clean

all: mkgalaxy bsrender

clean:
	rm -f mkgalaxy bsrender $(OBJ)

$(OBJ): %.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

mkgalaxy: mkgalaxy.o
	$(CC) $(CFLAGS) $(LIBS) -o mkgalaxy $^

bsrender: bsrender.o
	$(CC) $(CFLAGS) $(LIBS) -o bsrender $^

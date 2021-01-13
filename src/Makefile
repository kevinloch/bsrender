CC = gcc
CFLAGS = -I. -I/usr/local/include -Wall -Ofast
LIBS = -lm -lpng -L/usr/local/lib
OBJ = usage.o bsr-config.o bsrender.o
DEPS = usage.h bsr-config.h bsrender.h

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

#Makefile
TARGET = driver
CC = gcc
CFLAGS = -std=c99 -g -Wall -DDEBUG
INC_DIR = -I../tiff-4.0.2/libtiff
INC_DIR += -I../
INC_DIR += -I../perf
LIB_DIR = ../tiff-4.0.2/libtiff/.libs
LIBS_STATIC = -ltiff
LIBS_DYNAMIC = -lm
LIBS = -Wl,-Bstatic $(LIBS_STATIC) -Wl,-Bdynamic $(LIBS_DYNAMIC)

TIFFLIB = ../tiff-4.0.2/libtiff/.libs/libtiff.a

.PHONY: distclean clean all default

default: $(TARGET)
all: default

OBJ = $(patsubst %.c, %.o, $(wildcard *.c))
DEP = $(wildcard *.h) $(TIFFLIB)

%.o: %.c $(DEP)
	$(CC) -c $(CFLAGS) $(INC_DIR) $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TIFFLIB): ../tiff-4.0.2.tar.gz
	cd ..
	rm -rf tiff-4.0.2
	tar -xzf tiff-4.0.2.tar.gz
	cd tiff-4.0.2 && ./configure --enable-static --disable-shared --disable-jpeg --disable-lzma --disable-zlib --disable-jbig
	cd tiff-4.0.2 && make -j16

$(TARGET): $(OBJ) $(TIFFLIB)
	$(CC) $(OBJ) -L$(LIB_DIR) $(LIBS) -o $@

#clean default
clean:
	rm -f *.o
	rm -f $(TARGET)
	rm -f debug.tif

distclean: clean
	rm -rf ../tiff-4.0.2

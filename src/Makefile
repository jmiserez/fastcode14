#Makefile
TARGET = self

# macros/hardcoded values
N = 1
NB = 1
NBDEF = -DNB=$(NB)
NDEF = -DN=$(N)
CPUMHZDEF = -DCPUMHZ=2192.758
OPCOUNT0 = -DOPCOUNT=$$(( (2*($(NB)*$(NB)*$(NB)))+($(NB)*$(NB)) ))

#compilation options
SHELL = /bin/bash
CC = gcc
AFLAGS = -O3 -std=c99 -m64 -march=native -mno-abm -fno-tree-vectorize
DFLAGS = -g -Wall
PFLAGS = -D_GNU_SOURCE -pthread
CFLAGS = $(AFLAGS) $(PFLAGS) $(DFLAGS) $(NBDEF) $(CPUMHZDEF)
OBJDUMP = objdump -d
INC_DIR = -Itiff-4.0.2/libtiff -Ilibpfm-4.5.0/include -Ilibpfm-4.5.0/perf_examples
LIB_DIR = -Ltiff-4.0.2/libtiff/.libs -Llibpfm-4.5.0/lib
LIBS_STATIC = -ltiff -lpfm
LIBS_DYNAMIC = -lm
LIBS = -Wl,-Bstatic $(LIBS_STATIC) -Wl,-Bdynamic $(LIBS_DYNAMIC)

TIFFLIB = tiff-4.0.2/libtiff/.libs/libtiff.a
PFMLIB = libpfm-4.5.0/lib/libpfm.a

#libpfm Makefile options
PERF_EVENT_HDR=libpfm-4.5.0/include/perfmon/pfmlib_perf_event.h
LPC_UTILS=libpfm-4.5.0/perf_examples/perf_util.o

.PHONY: distclean clean all default

default: $(TARGET)
all: default

#OBJ = $(patsubst %.c, %.o, $(wildcard *.c))
# cannot use automatic resolution, as not all .c files should be used
DEP = $(wildcard *.h) $(TIFFLIB) $(LPC_UTILS) $(PFMLIB) $(PERF_EVENT_HDR)

%.o: %.c $(DEP)
	$(CC) -c $(CFLAGS) $(INC_DIR) $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TIFFLIB): tiff-4.0.2.tar.gz
	rm -rf tiff-4.0.2
	tar -xzf tiff-4.0.2.tar.gz
	cd tiff-4.0.2 && ./configure --enable-static --disable-shared --disable-jpeg --disable-lzma --disable-zlib --disable-jbig
	cd tiff-4.0.2 && make -j16
	
$(PFMLIB): libpfm-4.5.0.tar.gz
	rm -rf libpfm-4.5.0
	tar -xzf libpfm-4.5.0.tar.gz
	cd libpfm-4.5.0 && make -j16
	
$(LPC_UTILS): $(PFMLIB) 

#$(TARGET): $(OBJ) $(TIFFLIB) $(PFMLIB)
#	$(CC) $(OBJ) $(INC_DIR) $(LIB_DIR) $(LIBS) -o $@

main0: main.c code0.o $(DEP)
	$(CC) main.c code0.o $(CFLAGS) $(NBDEF) $(OPCOUNT0) $(INC_DIR) $(LIB_DIR) $(LIBS) -o $@
	
self: self.c $(DEP)
	$(CC) $(CFLAGS) $(INC_DIR) $(LIB_DIR) $(LIBS) -o $@ mak


#clean default
clean:
	rm -f *.o
	rm -f $(TARGET)
	rm -f result.tif
	rm -f debug.tif
	rm -f $(wildcard _*.tif)

distclean: clean
	rm -rf tiff-4.0.2
	rm -rf libpfm-4.5.0
	

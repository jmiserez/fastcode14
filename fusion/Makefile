#Makefile
TARGET = main
LIBS = -ltiff
LIBTIFFSRCDIR = ./libs/tiff-4.0.3
LIBTIFFBUILDDIR = ./libs/tiff-4.0.3-build
LDIRS = -L$(LIBTIFFBUILDDIR)/libtiff/.libs
IDIRS = -I$(LIBTIFFSRCDIR)/libtiff -I$(LIBTIFFBUILDDIR)/libtiff
CC = gcc
CFLAGS = -std=c99 -g -Wall $(IDIRS)

.PHONY: clean cleanlibs cleanall all libs default

#all c/h files in directory, depend on $(LIBTIFFBUILDDIR)
default: libs $(TARGET)

#build libtiff from scratch in separate directory only if it does not yet exist
libs: | $(LIBTIFFBUILDDIR)
$(LIBTIFFBUILDDIR):
	cd ./libs && tar xf tiff-4.0.3.tar.gz
	mkdir ./libs/tiff-4.0.3-build
	cd ./libs/tiff-4.0.3-build && ../tiff-4.0.3/configure
	make -C ./libs/tiff-4.0.3-build

#build everything
all: default
	
OBJ = $(patsubst %.c, %.o, $(wildcard *.c))
DEP = $(wildcard *.h)

%.o: %.c $(DEP)
	$(CC) -c $(CFLAGS) $< -o $@

#prevent removal of generated files
.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LDIRS) $(LIBS) -o $@

#clean default
clean:
	rm -f *.o
	rm -f $(TARGET)
	
cleanlibs:
	rm -rf $(LIBTIFFBUILDDIR)
	rm -rf $(LIBTIFFSRCDIR)

cleanall: clean	cleanlibs
	

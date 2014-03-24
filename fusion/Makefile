#Makefile
TARGET = main
LIBS = -lm -ltiff
CC = gcc
CFLAGS = -std=c99 -g -Wall

.PHONY: default all clean

default: $(TARGET)
all: default

OBJ = $(patsubst %.c, %.o, $(wildcard *.c))
DEP = $(wildcard *.h)
#OBJ = main.o
#DEP = main.h

%.o: %.c $(DEP)
	$(CC) -c $(CFLAGS) $< -o $@

#.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -Wall $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f $(TARGET)

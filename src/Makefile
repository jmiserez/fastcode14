#Makefile
TARGET = fusion
LIBS = -lm -ltiff
CC = gcc
CFLAGS = -std=c99 -g -Wall

.PHONY: clean all default

default: $(TARGET)
all: default
	
OBJ = $(patsubst %.c, %.o, $(wildcard *.c))
DEP = $(wildcard *.h)

%.o: %.c $(DEP)
	$(CC) -c $(CFLAGS) $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LIBS) -o $@

#clean default
clean:
	rm -f *.o
	rm -f $(TARGET)

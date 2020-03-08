CC=gcc
FORTRANC=gfortran

CLASS_DIRECTORY = /home/qvgd89/class_n/class_public
CLASS_INCLUDES = -I$(CLASS_DIRECTORY)/include
CLASS_LIBRARIES = -L$(CLASS_DIRECTORY) -lclass -Wl,-rpath=$(CLASS_DIRECTORY)

INCLUDES=-I/usr/include/hdf5/serial $(CLASS_INCLUDES)
LIBRARIES=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 $(CLASS_LIBRARIES)
INTERNAL=int.o

all:
	$(CC) int.c -c $(INCLUDES)
	$(FORTRANC) interface.f90 -o renderer $(INTERNAL) $(LIBRARIES)

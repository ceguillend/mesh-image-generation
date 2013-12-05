# Makefile. 
#
# Use 'make' for release version.
# Use 'make debug' for debug version.
# Use 'make clean' for clean everything.
#

FLAGS=`pkg-config opencv --cflags`
LIBS=`pkg-config opencv --libs`
CC=g++
	
all: debug 

clean: 
	rm -rf *.o method

debug: FLAGS += -g
release: FLAGS += -s -O3

debug: build
release: build

build: curve.o thinning.o basic_geo.o delaunay.o main.o 
	$(CC) $(FLAGS) main.o curve.o thinning.o basic_geo.o delaunay.o -o method $(LIBS)

main.o: main.cpp
	$(CC) -c $(FLAGS) main.cpp

curve.o: curve.cpp
	$(CC) -c $(FLAGS) curve.cpp

thinning.o: thinning.cpp
	$(CC) -c $(FLAGS) thinning.cpp

basic_geo.o: basic_geo.cpp
	$(CC) -c $(FLAGS) basic_geo.cpp

delaunay.o: delaunay.cpp
	$(CC) -c $(FLAGS) delaunay.cpp

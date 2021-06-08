CXX = g++
CXXFLAGS = -Iinclude
LDFLAGS = -Xcompiler

all: build/lung.o
	$(CXX) $^ $(LDFLAGS) -o lungmodel

debug: build/dbg_lung.o
	$(CXX) $^ $(LDFLAGS) -g -G -o lungmodel

clean:
	rm build/lung.o

build/lung.o: src/lung.cpp build
	$(CXX) -c -o $@ $(CXXFLAGS) src/lung.cpp

build/dbg_lung.o: src/lung.cpp build
	$(CXX) -c -o $@ $(CXXFLAGS) src/lung.cpp

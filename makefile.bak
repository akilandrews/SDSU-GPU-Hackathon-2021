RAND123 = /gpfs/wolf/gen157/proj-shared/simcov/curand/Random123-1.13.2

CXX = g++
NVCC = nvcc
CXXFLAGS = -I$(RAND123)/include -Xcompiler -fopenmp -Iinclude
LDFLAGS = -Xcompiler -fopenmp

all: build/main.o build/system.o
	$(NVCC) $^ $(LDFLAGS) -o simcov_simple


debug: build/dbg_main.o build/dbg_system.o
	$(NVCC) $^ $(LDFLAGS) -g -G -o simcov_simple

clean:
	rm build/main.o build/system.o

build:
	mkdir -p build

build/main.o: src/main.cu build
	$(NVCC) -c -o $@ $(CXXFLAGS) src/main.cu

build/system.o: src/system.cu build
	$(NVCC) -c -o $@ $(CXXFLAGS) src/system.cu

build/dbg_main.o: src/main.cu
	$(NVCC) -c -o $@ $(CXXFLAGS) src/main.cu

build/dbg_system.o: src/system.cu build
	$(NVCC) -c -o $@ $(CXXFLAGS) src/system.cu

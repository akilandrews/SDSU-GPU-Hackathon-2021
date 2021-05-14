#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <omp.h>
#include <fstream>
#include <iostream>

#include "tissue.h"

void dump_state(System& sys, Device dev, std::string output_file, int iteration){    

    const Geom& geo = sys.geo;
    Block* tissue;

    if(dev == Device::CUDA) {
        tissue = (Block *)malloc(geo.num_blocks()*sizeof(Block));
        cudaMemcpy(tissue, sys.tissue, geo.num_blocks()*sizeof(Block), cudaMemcpyDeviceToHost);
    } else if(dev == Device::CPU) {
        tissue = sys.tissue;
    }

    std::ofstream fp;
    fp.open(output_file, std::ios_base::app);
    
    


    idx_t nblk = geo.num_blocks();
    idx_t grids_per_block = geo.wx*geo.wy*geo.wz;

    for(idx_t i = 0; i < nblk; i++) {
        for(idx_t j = 0; j < grids_per_block; j++) {
            int x, y, z;
            geo.to_xyz(i, j, &x, &y, &z);

            if(x < geo.nx) {

                int epitype = tissue[i].epitype[j];
                double virions = (double)tissue[i].virions[j];
            
                std::string output = std::to_string(iteration) + "," +
                                    std::to_string(x) + "," +
                                    std::to_string(y) + "," +
                                    std::to_string(z) + "," +
                                    std::to_string(epitype) + "," +
                                    std::to_string(virions) + "\n";
                fp << output;   
            }        
        }
    }    
    fp.close();
}

// generates a "Tissue object" containing a grid of blocks
void init_tissue_cpu(Block *tissue, const Geom &geo)
{
    idx_t nblk = geo.num_blocks();
    idx_t grids_per_block = geo.wx*geo.wy*geo.wz;
    
    #pragma omp parallel for
    for(idx_t i = 0; i < nblk; i++)
    {
        for(idx_t j = 0; j < grids_per_block; j++)
        {
            geo.to_xyz(i, j, &tissue[i].x[j], &tissue[i].y[j], &tissue[i].z[j]);
            tissue[i].epitype[j] = 0;
            tissue[i].virions[j] = 0.0;
        }
    }
}

Block *allocate_tissue(const idx_t num_blocks, Device d) {
    Block *ret;
    if(d == Device::CPU) {
        ret = (Block *)malloc(num_blocks*sizeof(Block));
    } else {
        CHECK_CUDA( cudaMalloc(&ret, num_blocks*sizeof(Block)) );
    }
    return ret;
}

void free_tissue(Block *tissue, Device d) {
    if(d == Device::CPU) {
        free(tissue);
    } else {
        cudaFree(tissue);
    }
}

/*
    TODO: Break out these parameters so we can parameterize experiments.
    Just keeping them here for testing for now.
*/
int main(int argc, char *argv[]) {
    Device dev = Device::CPU;
    
    std::string output_file;

    if(argc >= 5 && !strcmp(argv[1], "-gpu")) {
        printf("Using CUDA device.\n");
        dev = Device::CUDA;
        for(int i=1; i<argc-1; i++) {
            argv[i] = argv[i+1];
        }
        argc--;
    }
    if(argc < 5) {
        printf("Usage: %s [-gpu] <dim_x> <dim_y> <dim_z> <output_file>\n", argv[0]);
        return -1;
    }
    int dim_x = atoi(argv[1]);
    int dim_y = atoi(argv[2]);
    int dim_z = atoi(argv[3]);
    output_file = argv[4];

    // Create the Grid of GridPoints, currently full of epicells.
    Params param(0.1, 0.1, // diffusion_rate, infectivity
                10, 10, 10, // periods( incb, expr, apop )
                100, // virion_production
                1 // seed
                );
    System sys(dim_x, dim_y, dim_z, dev, param);

    const Geom &geo = sys.geo;

    int initial_infection_x =  0 % geo.nx;
    int initial_infection_y = 49 % geo.ny;
    int initial_infection_z = 49 % geo.nz;
    real initial_virions = 1000.0;
    int maxiters = 100;

    sys.init();
    printf("#Generated grid of dimensions (%d, %d, %d)\n", geo.nx, geo.ny, geo.nz);
 
    // infect initial cells
    printf("Initializing virion at (%d, %d, %d)\n",
                initial_infection_x, initial_infection_y, initial_infection_z);
    sys.infect_point(initial_infection_x, initial_infection_y,
                     initial_infection_z, initial_virions);   

    printf("#Executing Simulation... \n");
    double t0 = omp_get_wtime();

    for(int iter=0; iter < maxiters; iter++){
        sys.step(iter);
    }
    if(dev == Device::CUDA) {
        cudaDeviceSynchronize();
    }
    double t1 = omp_get_wtime();
    printf("Time: %f\n", t1 - t0); 
    dump_state(sys, dev, output_file, maxiters-1);    
    printf("#Freeing memory...\n");

    return 0;
}

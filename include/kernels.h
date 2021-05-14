#pragma once

#include <tissue.h>

GLOBAL_FN inline float crand(unsigned int seed, unsigned int t, uint32_t id) {
    philox4x32_key_t k = {{id, seed}};
    philox4x32_ctr_t c = {{t}};

    return 1.0 - r123::u01<float,uint32_t>( philox4x32(c, k)[0] );
}

GLOBAL_FN inline void spread_virions(const Geom geo, Block* tissue, Block* next,
                                    const int dx[], const int dy[], const int dz[],
                                    const int ND, real diffusion_rate, idx_t i, idx_t j) {
    //diffuse virus
    int ax, ay, az;
    geo.to_xyz(i, j, &ax, &ay, &az);
    //printf("%d,%d,%d,%f\n", iter, i, j, tissue[i].virions[j]);
    double virions = 0.0;
    int nnbr = 0;
    for(int k=0; k<ND; k++) {
        idx_t i2, j2;
        int bx = (ax+dx[k]+geo.nx) % geo.nx;
        int by = (ay+dy[k]+geo.ny) % geo.ny;
        int bz = (az+dz[k]+geo.nz) % geo.nz;
        geo.to_idx(&i2, &j2, bx, by, bz);
        virions += tissue[i2].virions[j2];
        nnbr += 1;
    }
    real self = tissue[i].virions[j];
    virions = self + diffusion_rate*(virions - nnbr*self);
    next[i].virions[j] = virions;
}

GLOBAL_FN inline void check_epicell_state(const Geom geo, Block* tissue, Block* next,
                                          int incb_period, int expr_period, int apop_period,
                                          int virion_production, real infectivity, float r, idx_t i, idx_t j) {

   
    //start with the same epitype and info
    next[i].epitype[j] = tissue[i].epitype[j];
    next[i].was_expressing[j] = tissue[i].was_expressing[j];
    next[i].incb_period[j] = tissue[i].incb_period[j];
    next[i].expr_period[j] = tissue[i].expr_period[j];
    next[i].apop_period[j] = tissue[i].apop_period[j];
    bool produce_virions = false;
    switch(tissue[i].epitype[j]) {
        case 0: //healthy
        if(tissue[i].virions[j] > 0) {
            if(r < infectivity*tissue[i].virions[j]) {
                //become incubating
                next[i].epitype[j] = 1;
                next[i].incb_period[j] = incb_period;
            }
        }
        break;
        case 1: //incubating
            if(tissue[i].incb_period[j] <= 0) {
                //become expressing
                next[i].epitype[j] = 2;
                next[i].was_expressing[j] = true;
                next[i].expr_period[j] = expr_period;
            } else {
                next[i].incb_period[j]--;
            }
        break;
        case 2: //expressing
            if(tissue[i].expr_period[j] <= 0) {
                next[i].epitype[j] = 4; // dies
            } else {
                next[i].expr_period[j]--;
                produce_virions = true;
            }
        break;
        case 3: //apoptotic
            if(tissue[i].apop_period[j] <= 0) {
                next[i].epitype[j] = 4; // dies
            } else {
                next[i].apop_period[j]--;
                if(tissue[i].was_expressing[j]) {
                    produce_virions = true;
                }
            }
        break;
        default: break;
    }
    if(produce_virions) {
        next[i].virions[j] += virion_production;
    }


}


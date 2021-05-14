#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
//#include <curand.h>
//#include <curand_kernel.h>
#include <string>

// for small quantity comparison
//static const float eps = 1e-7;

using idx_t = uint32_t;
using real  = float;
#define WIDTH (64)

#define GLOBAL_FN __host__ __device__
// FIXME
#define CHECK_CUDA(call) { call; }

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if(! (condition)) { \
            fprintf(stderr, "Assertion `%s` failed in %s line %d: %s\n", \
                    #condition, __FILE__, __LINE__, message); \
            exit(1); \
        } \
    } while(false);
#else
#   define ASSERT(condition, message) do {} while(false)
#endif

enum class Device {
  CPU, CUDA
};

struct Geom {
    int nx, ny, nz; // total size
    int wx, wy, wz; // width per Block
    int ux, uy, uz; // blocks in x,y,z directions
    Geom(int nx_, int ny_, int nz_, int wx_, int wy_, int wz_)
            : nx(nx_), ny(ny_), nz(nz_), wx(wx_), wy(wy_), wz(wz_),
              ux( (nx_+wx_-1)/wx_ ),
              uy( (ny_+wy_-1)/wy_ ),
              uz( (nz_+wz_-1)/wz_ ) {
        ASSERT( wx*wy*wz <= WIDTH, "block widths exceed WIDTH" );
    }

    /** Decode Block and thread indices to absolute x,y,z */
    GLOBAL_FN inline void to_xyz(idx_t i, idx_t j, int *ax, int *ay, int *az) const {
        *ax = (i/(uy*uz))*wx + j/(wy*wz);
        *ay = ((i/uz)%uy)*wy + (j/wz)%wy;
        *az = (i%uz)*wz + j%wz;
    }
    /** Encode absolute x,y,z to Block and thread indices */
    GLOBAL_FN inline void to_idx(idx_t *i, idx_t *j, int ax, int ay, int az) const {
        *i = ((ax/wx)*uy + (ay/wy))*uz + (az/wz);
        *j = ((ax%wx)*wy + (ay%wy))*wz + (az%wz);
    }
    GLOBAL_FN inline idx_t num_blocks() const {
        return ux*uy*uz;
    }
};

/* A block contains a list of gridpoint data */
/* This is the array-of-structure-of-arrays pattern. */
struct Block {
    int x[WIDTH];
    int y[WIDTH];
    int z[WIDTH];
    int epitype[WIDTH];
    int incb_period[WIDTH];
    int expr_period[WIDTH];
    int apop_period[WIDTH];
    bool was_expressing[WIDTH];
    real virions[WIDTH];
};

struct Params {
    const real diffusion_rate;
    const real infectivity;
    const int incb_period;
    const int expr_period;
    const int apop_period;
    const int virion_production;
    const long unsigned seed;

    Params(real diffusion_rate_, real infectivity_,
           int incb_period_, int expr_periond_,
           int apop_period_, int virion_production_,
           long unsigned seed_)
        : diffusion_rate(diffusion_rate_),
          infectivity(infectivity_),
          incb_period(incb_period_),
          expr_period(expr_periond_),
          apop_period(apop_period_),
          virion_production(virion_production_),
          seed(seed_) {}
};

/** Note: This class performs actions during
 *        construction / destruction that make
 *        it unsuitable for use directly on the GPU.
 *
 *        Instead, use a System_d class to pass to GPU
 *        functions.  That class is a shallow copy of this one
 *        that performs no memory alloc/free on construction/destruction.
 */
enum GPU_RED{VIRION, EPITYPE};
class System {
  public:
    const Device dev;
    const Geom geo;
    const Params param;
    Block *tissue, *next;

    //curandState* g_state;

    System(int nx, int ny, int nz, Device dev, const Params &param);
    ~System();
    void init();
    void step(unsigned int t);
    void infect_point(int x, int y, int z, real virions);
    template <typename T>
    T gpu_reduce(Block* tissue, GPU_RED red_type, int epi_type);
    template <typename T>    
    T cpu_reduce(Block* tissue, GPU_RED red_type, int epi_type);

};

/** Shallow copy of System, meant to be passed to a device function.
 */
struct System_d {
    const Device dev;
    const Geom geo;
    const Params param;
    Block *tissue, *next;

    System_d(System &s)
        : dev(s.dev), geo(s.geo), tissue(s.tissue),
          next(s.next), param(s.param) {}
};

// try a 9-point stencil
/*
#define ND 9
struct Diffusion {
    int dx[ND];
    int dy[ND];
    int dz[ND];
    real coef[ND];
};

static const real dc = 6.58608256e-04;
static const real db = 2.43461476e-02;
static const real da = 8.99980977e-01;
static const Diffusion diffusion = {
  {0, 0, 0, 0, 0, 0, 0, 0, 0},
  {-1, -1, -1, 0, 0, 0, 1, 1, 1},
  {-1, 0, 1, -1, 0, 1, -1, 0, 1},
  {dc, db, dc, db, da, db, dc, db, dc}
};*/

Block *allocate_tissue(const idx_t num_blocks, Device d);
void free_tissue(Block *tissue, Device d);
void init_tissue_cpu(Block *tissue, const Geom &geo);
void run_step_cpu(System &s);
void update_epicells_cpu(System &s, unsigned int time);
void dump_state(Block *tissue, const Geom &geo, std::string output_file, int iteration);


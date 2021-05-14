#include <Random123/philox.h>
#include <Random123/uniform.hpp>
#include "tissue.h"
#include "kernels.h"
#define FULL_MASK 0xffffffff

void System::infect_point(int x, int y, int z, float virions)
{
    idx_t i, j;
    geo.to_idx(&i, &j, x, y, z);
    if(dev == Device::CPU) {
        tissue[i].virions[j] = virions;
    } else {
        cudaMemcpy(&tissue[i].virions[j], &virions, sizeof(float), cudaMemcpyHostToDevice);
    }
}

/** tissue must reside on device memory.
 */
//from developer.nvidia.com
template <typename T>
__device__ T warpReduceType(T val, GPU_RED red_type, int epi_type) {
	int warpSize = 32;
	for (int offset = warpSize/2; offset > 0; offset /= 2){
    		auto in_val= __shfl_down_sync(FULL_MASK, val, offset);
		if(red_type == GPU_RED::VIRION){
			val += in_val;
		}else{
			val += in_val == epi_type;
		}
		
	}
	return val;
}

template <typename T>
__device__ T warpReduceSum(T val) {
	int warpSize = 32;
	for (int offset = warpSize/2; offset > 0; offset /= 2){
    		val+= __shfl_down_sync(FULL_MASK, val, offset);
		
	}
	return val;
}
//for block (CUDA) wide reduction
template <typename T>
__device__ T blockReduceSum(T val) {

	static __shared__ T shared[32];
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;

	val = warpReduceSum(val);

	if (lane==0) shared[wid]=val;
	__syncthreads();

	val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
	if (wid==0) val = warpReduceSum(val); // reducing partial sums, this will be diff
	return val;
}

//for reducing an array
template <typename T>
__global__ void array_reduce(T *in, T* out, int N) {
	T sum = 0;
	for (idx_t i = blockIdx.x * blockDim.x + threadIdx.x; 
		i < N; 
		i += blockDim.x * gridDim.x) {
		sum += in[i];
	}
	sum = blockReduceSum(sum);
	if (threadIdx.x==0)
		out[blockIdx.x]=sum;
}

//this kernel needs to be launched such that BlockDim = block WIDTH
template <typename T>
__global__ void grid_block_reduce(T* block_sums, Block *tissue, Geom geo, GPU_RED red_type, int epi_type){
	idx_t nblk = geo.num_blocks();
	idx_t grids_per_block = geo.wx*geo.wy*geo.wz;
	idx_t j = threadIdx.x;
	for(idx_t i = blockIdx.x; i < nblk; i +=gridDim.x){
		if(j < grids_per_block){
			//if(red_type == GPU_RED::VIRION)
				auto sum_val = (red_type == GPU_RED::VIRION)?tissue[i].virions[j]:tissue[i].epitype[j];
	//		else
	//			auto sum_val = tissue[i].epitype[j];
			if(red_type == EPITYPE)
				sum_val = sum_val == epi_type;
			sum_val = blockReduceSum(sum_val);
			__syncthreads();
			if(j == 0) block_sums[i] = sum_val;
		}		
	}
}



__global__ void init_tissue_gpu(Block *tissue, Geom geo)
{
    idx_t nblk = geo.num_blocks();
    idx_t grids_per_block = geo.wx*geo.wy*geo.wz;
    
    for(idx_t i = blockIdx.x; i < nblk; i += gridDim.x)
    {
        idx_t j = threadIdx.x;
	printf("%d, %d", &i, &j);
        geo.to_xyz(i, j, &tissue[i].x[j], &tissue[i].y[j], &tissue[i].z[j]);
      	tissue[i].epitype[j] = 0;
        tissue[i].virions[j] = 0.0;
    }
}

/** tissue and next must reside on device memory
 */
__global__ void run_step_gpu(System_d sys) {
    const int dx[] = { 0, 0, 0, 0, 0, 0, 0, 0};
    const int dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const int dz[] = {-1, 0, 1, -1, 1, -1, 0, 1};
    const int ND = sizeof(dx) / sizeof(int);

    const Geom &geo = sys.geo;
    Block *tissue = sys.tissue;
    Block *next = sys.next;
    
    idx_t nblk = geo.num_blocks();
    idx_t grids_per_block = geo.wx*geo.wy*geo.wz;
    
    for(idx_t i = blockIdx.x; i < nblk; i += gridDim.x)
    {
        idx_t j = threadIdx.x;
        spread_virions(geo, tissue, next, dx, dy, dz, ND, sys.diffusion_rate, i, j);
    }
}



__global__ void update_epicells_gpu(System_d sys, unsigned int t) {
    real infectivity = sys.param.infectivity;
    int virion_production = sys.param.virion_production;
    int incb_period = sys.param.incb_period;
    int expr_period = sys.param.expr_period;
    int apop_period = sys.param.apop_period;

    const Geom &geo = sys.geo;
    Block *tissue = sys.tissue;
    Block *next = sys.next;
    
    idx_t nblk = geo.num_blocks();
    idx_t grids_per_block = geo.wx*geo.wy*geo.wz;
    
    idx_t t_id  = threadIdx.x;

    for(idx_t i = blockIdx.x; i < nblk; i += gridDim.x)
    {
	    idx_t j = threadIdx.x;
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
                    float r = 0.0;
                    r = curand_uniform((&g_state)[t_id]);
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
}

System::System(int nx, int ny, int nz, Device dev_, const Params &param_)
        : dev(dev_), geo(nx, ny, nz, 4, 4, 4), param(param_)
{
    idx_t num_blocks = geo.num_blocks();

    tissue = allocate_tissue(num_blocks, dev);
    next   = allocate_tissue(num_blocks, dev);

    //set up curand
    /*if(dev == Device::CUDA) {
        cudaMalloc(&g_state, WIDTH*sizeof(curandState));
        setup_curand<<< 256,WIDTH >>>(g_state, param.seed);
    }*/
}

System::~System() {
    free_tissue(tissue, dev);
    free_tissue(next, dev);

    /*if(dev == Device::CUDA) {
        cudaFree(&g_state);
    }*/
}

void System::init() {
    if(dev == Device::CPU) {
        init_tissue_cpu(tissue, geo);
    } else {
        dim3 threadsPerBlock(4, 4, 4);
	dim3 numBlocks(geo.ux / threadsPerBlock.x, geo.uy / threadsPerBlock.y, geo.uz / threadsPerBlock.z);
        init_tissue_gpu<<< numBlocks , threadsPerBlock >>>(tissue, geo);
    }
}
void System::step(unsigned int t) {
    if(dev == Device::CPU) {
        run_step_cpu(*this);
        update_epicells_cpu(*this, t);
	//auto sum_virions = cpu_reduce<real>(this->tissue, GPU_RED::VIRION, 0);
	//auto sum_epi_one = cpu_reduce<int>(this->tissue, GPU_RED::EPITYPE, 1);
	//auto sum_epi_two = cpu_reduce<int>(this->tissue, GPU_RED::EPITYPE, 2);
	//printf("\nt:%d\tVirions:%f\tEpitype_1:%d\tEpitype_2:%d\n",t, sum_virions,sum_epi_one, sum_epi_two);
    } else {
        System_d sys(*this);
	//auto sum_virions = gpu_reduce<real>(sys.tissue, GPU_RED::VIRION, 0);
	//auto sum_epi_one = gpu_reduce<int>(sys.tissue, GPU_RED::EPITYPE, 1);
	//auto sum_epi_two = gpu_reduce<int>(sys.tissue, GPU_RED::EPITYPE, 2);
	 //printf("\nt:%d\tVirions:%f\tEpitype_1:%d\tEpitype_2:%d\n",t, sum_virions,sum_epi_one, sum_epi_two);
	dim3 threadsPerBlock(4, 4, 4);
	dim3 numBlocks(geo.ux / threadsPerBlock.x, geo.uy / threadsPerBlock.y, geo.uz / threadsPerBlock.z);
       
        run_step_gpu<<< numBlocks, threadsPerBlock >>>(sys);
        update_epicells_gpu<<< numBlocks, threadsPerBlock >>>(sys,
            g_state,
            infectivity, virion_production, incb_period, expr_period, apop_period);
    }
    // swap current and next tissue block
    Block *tmp = next;
    next = tissue;
    tissue = tmp;
}

void run_step_cpu(System &s)
{
    const int dx[] = { 0, 0, 0, 0, 0, 0, 0, 0};
    const int dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const int dz[] = {-1, 0, 1, -1, 1, -1, 0, 1};
    const int ND = sizeof(dx) / sizeof(int);

    const Geom &geo = s.geo;
    Block *tissue = s.tissue;
    Block *next = s.next;

    idx_t nblk = geo.num_blocks();
    idx_t grids_per_block = geo.wx*geo.wy*geo.wz;
    
    #pragma omp parallel for
    for(idx_t i = 0; i < nblk; i++)
    {
        for(idx_t j = 0; j < grids_per_block; j++)
        {
            spread_virions(geo, tissue, next,  dx, dy, dz, ND, s.param.diffusion_rate, i, j);
        }
    }
}

void update_epicells_cpu(System& sys, unsigned int t) {
    real infectivity = sys.param.infectivity;
    int virion_production = sys.param.virion_production;
    int incb_period = sys.param.incb_period;
    int expr_period = sys.param.expr_period;
    int apop_period = sys.param.apop_period;

    const Geom &geo = sys.geo;
    Block *tissue = sys.tissue;
    Block *next = sys.next;
    
    idx_t nblk = geo.num_blocks();
    idx_t grids_per_block = geo.wx*geo.wy*geo.wz;

    #pragma omp parallel for
    for(idx_t i = 0; i < nblk; i++)
    {
        for(idx_t j = 0; j < grids_per_block; j++)
        {
            real r = crand(sys.param.seed, t, i*WIDTH+j);
            check_epicell_state(geo, tissue, next,
                                incb_period, expr_period, apop_period,
                                virion_production, infectivity, r, i, j);
        }
    }
}

template <typename T>
T System::cpu_reduce(Block* tissue, GPU_RED red_type, int epi_type){
	T sum = 0;
	idx_t grids_per_block = geo.wx*geo.wy*geo.wz;

	for(int i = 0; i < geo.num_blocks(); i++){
		for(int j = 0; j < grids_per_block; j++){
			if(red_type == GPU_RED::VIRION)
				sum += tissue[i].virions[j];	
			else
				sum += tissue[i].epitype[j] == epi_type;
		}
	}
	return sum;
}

template <typename T>
T System::gpu_reduce(Block* tissue, GPU_RED red_type, int epi_type){
	System_d sys(*this);

	T *block_sums;
	T *block_sums_partial;
	cudaMalloc(&block_sums, sizeof(T)*geo.num_blocks());
	grid_block_reduce<<<256, WIDTH>>>(block_sums, tissue, geo, red_type, epi_type);

	int threads = 512;
	int N = this->geo.num_blocks();
	int blocks = min((N + threads - 1) / threads, 1024);
	cudaMalloc(&block_sums_partial, sizeof(T)*blocks);

	array_reduce<<<blocks, threads>>>(block_sums, block_sums_partial, N);
	array_reduce<<<1, 1024>>>(block_sums_partial, block_sums_partial, blocks);
	
	T *blocks_partial_cpu = new T[blocks];
	cudaMemcpy(blocks_partial_cpu, block_sums_partial, sizeof(T)*blocks, cudaMemcpyDeviceToHost);
	return blocks_partial_cpu[0];
}

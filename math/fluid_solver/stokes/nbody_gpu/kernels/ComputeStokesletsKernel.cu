
#include <cuda.h>
#include <cuda_runtime.h>

// Macros to simplify shared memory addressing
// the first half of the shared memory allocation is pos
#define sources(i) shared_data[i]
// and the second half is strength
#define forces(i) shared_data[i+blockDim.x]


/**
 * @brief ...
 *
 * @param velocity ...
 * @param target ...
 * @param source ...
 * @param force ...
 * @param delta ...
 **/
__device__ 
float3 ComputeStokeslets(float3 target, float3 velocity, float3 source, float3 force, float delta)
{
    float3 dx;
    dx.x = target.x - source.x;
    dx.y = target.y - source.y;
    dx.z = target.z - source.z;

    float r2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;
    float d2 = delta * delta;
    float H1 = 0.5f / sqrtf(r2 + d2) + d2 / 2.0f / powf( sqrtf(r2 + d2), 3);
    float H2 = 0.5f / powf( sqrtf(r2 + d2), 3);

    float A11 = H1 + dx.x * dx.x * H2;
    float A12 = dx.x * dx.y * H2;
    float A13 = dx.x * dx.z * H2;

    float A21 = dx.y * dx.x * H2;
    float A22 = H1 + dx.y * dx.y * H2;
    float A23 = dx.y * dx.z * H2;

    float A31 = dx.z * dx.x * H2;
    float A32 = dx.z * dx.y * H2;
    float A33 = H1 + dx.z * dx.z * H2;

    velocity.x += (A11 * force.x + A12 * force.y + A13 * force.z) / 12.566370614359172f;
    velocity.y += (A21 * force.x + A22 * force.y + A23 * force.z) / 12.566370614359172f;
    velocity.z += (A31 * force.x + A32 * force.y + A33 * force.z) / 12.566370614359172f;
    return velocity;
}

/**
 * @brief This is the "tile_calculation" function from the GPUG3 article.
 *
 * @param target ...
 * @param velocity ...
 * @param numBody ...
 * @param delta ...
 **/
__device__ 
float3 ComputeStokesletsTile(float3 target, float3 velocity, int numBody, float delta)
{
    extern __shared__ float3 shared_data[];
    // The CUDA 1.1 compiler cannot determine that i is not going to  overflow in the loop below.  
	// Therefore if int is used on 64-bit linux or windows (or long instead of long long on win64), 
	// the compiler generates suboptimal code.  
	// Therefore, we use long long on win64 and  long on everything else. (Workaround for Bug ID 347697)
	#ifdef _Win64
		unsigned long long i = 0;
	#else
		unsigned long i = 0;
	#endif
    // Note that having an unsigned int loop counter and an unsigned long index helps the compiler generate efficient code on 64-bit  OSes.  
	// The compiler can't assume the 64-bit index won't overflow     // so it incurs extra integer operations.  This is a standard issue     
	// in porting 32-bit code to 64-bit OSes.
    for (unsigned int counter = 0; counter < numBody; ++i, ++counter )
        velocity = ComputeStokeslets(target, velocity, sources(i), forces(i), delta);

    return velocity;
}

/**
 * @brief ...
 *
 * @param target ...
 * @param source_array ...
 * @param force_array ...
 * @param num_sources ...
 * @param num_targets ...
 * @param delta ...
 **/
__device__ 
float3 ComputeStokesletsTiles(float3 target, float3* source_array, float3* force_array, int num_sources, int num_targets, float delta)
{
    extern __shared__ float3 shared_data[];

    float3 velocity = {0.0f, 0.0f, 0.0f};

    int p = blockDim.x;
    int q = blockDim.y;

    int n = num_sources;
    int numTiles = n / (p * q);
    int index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    // all blocks must start from index 0.
	if (index >= num_targets)
		index = 0;

    for (int tile = 0; tile < numTiles; ++tile)
	{
        int posIdx = tile * blockDim.x + threadIdx.x;
        sources(threadIdx.x) = source_array[posIdx];
        forces(threadIdx.x) = force_array[posIdx];

        __syncthreads();

        if (index < num_targets)
            velocity = ComputeStokesletsTile(target, velocity, blockDim.x, delta);
        
        __syncthreads();
    }

    // handle if num_targets != 32 * n  (the last tile)
    int bodyLeft = n - numTiles * (p * q);
    if (bodyLeft > 0){
        // only 'bodyLeft' threads have to write the entries in share memory, (32-bodyLeft) threads are idle.
        if (threadIdx.x < bodyLeft){
            int posIdx = numTiles * blockDim.x + threadIdx.x;            
            sources(threadIdx.x) = source_array[posIdx];
            forces(threadIdx.x) = force_array[posIdx];
        }
        __syncthreads();

        if (index < num_targets)
            velocity = ComputeStokesletsTile(target, velocity, bodyLeft, delta);
        
        __syncthreads();
    }
    return velocity;
}



/**
 * @brief ...
 *
 * @param target_array ...
 * @param source_array ...
 * @param force_array ...
 * @param velocity_array ...
 * @param deltas ...
 * @param num_sources ...
 * @param num_targets ...
 * @return void
 **/
__global__ 
void ComputeStokesletsKernel(float3* target_array, float3* velocity_array, float3* source_array, float3* force_array, float deltas, int num_sources, int  num_targets)
{
	int index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
	if (index >= num_targets)
    	index = 0;  

	float3 target = target_array[index];
	float3 velocity = ComputeStokesletsTiles(target, source_array, force_array, num_sources, num_targets, deltas);
	
	if (index < num_targets)
		velocity_array[index] = velocity;
}


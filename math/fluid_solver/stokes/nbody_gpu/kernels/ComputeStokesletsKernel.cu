#ifndef COMPUTE_STOKESLETS_KERNEL_CU
#define COMPUTE_STOKESLETS_KERNEL_CU

#include "ComputeTiles.cuh"


/**
 * @brief Main stokelets computation kernel.
 *
 * @param target_array array containing all targets
 * @param source_array array containing all sources
 * @param force_array array containing all force vectors
 * @param velocity_array array containing all velocity vectors
 * @param deltas regularization parameter
 * @param num_sources total number of sources 
 * @param num_targets total number of targets
 **/
template<typename vector3_type, typename int_type, typename real_type>
__global__ 
void ComputeStokesletsKernel(const vector3_type* target_array, vector3_type* velocity_array, const vector3_type* source_array, const vector3_type* force_array, real_type deltas, int_type num_sources, int_type num_targets, bool with_images = false)
{
	int_type index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
	if (index >= num_targets)
    	index = 0;  

	vector3_type target = target_array[index];
	vector3_type velocity = ComputeTiles(target, source_array, force_array, num_sources, num_targets, deltas, with_images);
	
	if (index < num_targets)
            velocity_array[index] = velocity;
}

#endif

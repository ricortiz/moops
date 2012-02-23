#ifndef COMPUTE_STOKESLETS_NEAR_FIELD_KERNEL_CU
#define COMPUTE_STOKESLETS_NEAR_FIELD_KERNEL_CU

#include "ComputeStokesletsKernel.cu"

template<typename vector3_type, typename int_type, typename real_type>
__global__
void ComputeStokesletsNearFieldKernel(vector3_type* target_array, vector3_type* velocity_array, vector3_type* source_array, vector3_type* force_array, real_type delta, int_type  num_targets, bool with_images)
{
	int_type index_start = __mul24(blockIdx.x,blockDim.x);
    int_type index = index_start + threadIdx.x;
    
    int_type num_sources = blockDim.x;
	vector3_type target = target_array[index];
    vector3_type velocity = ComputeTiles(target, &source_array[index_start], &force_array[index_start], num_sources, num_targets, delta, with_images);

    if (index < num_targets)
    {
        velocity_array[index].x += velocity.x;      //returning addition of far & far_image
        velocity_array[index].y += velocity.y;
        velocity_array[index].z += velocity.z;
    }
}

#endif

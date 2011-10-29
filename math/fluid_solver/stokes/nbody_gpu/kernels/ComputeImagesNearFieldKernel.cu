
#include "ComputeImagesKernel.cu"

/**
 * @brief ...
 *
 * @param target_array ...
 * @param velocity_array ...
 * @param source_array ...
 * @param force_array ...
 * @param deltas ...
 * @param num_stokeslets ...
 * @param num_towers ...
 * @param num_sources ...
 * @param num_targets ...
 * @return void
 **/
__global__ 
void ComputeImagesNearFieldKernel(float3* target_array, float3* velocity_array, float3* source_array, float3* force_array, float delta,  int num_stokeslets, int  num_targets)
{
	int index_start = __mul24(blockIdx.x,blockDim.x);
    int index = index_start + threadIdx.x;
	if (index >= num_targets)
    	index = 0;  

	float3 target = target_array[index];
	float3 velocity = ComputeImagesTiles(target, &source_array[index_start], &force_array[index_start], blockDim.x, num_targets, delta);

	if (index < num_targets)
	{
		velocity_array[index].x += velocity.x;		//returning addition of far & far_image
		velocity_array[index].y += velocity.y;
		velocity_array[index].z += velocity.z;
	}

}


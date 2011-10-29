
#include "ComputeStokesletsKernel.cu"

__global__
void ComputeStokesletsNearFieldKernel(float3* target_array, float3* velocity_array, float3* source_array, float3* force_array, float delta, int  num_targets)
{

	int index_start = __mul24(blockIdx.x,blockDim.x);
    int index = index_start + threadIdx.x;

	float3 target = target_array[index];
    velocity_array[index] = ComputeStokesletsTiles(target, &source_array[index_start], &force_array[index_start], blockDim.x, num_targets, delta);
}


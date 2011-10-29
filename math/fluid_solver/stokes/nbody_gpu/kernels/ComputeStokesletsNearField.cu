#include "cuda.h"
#include "ComputeStokesletsNearFieldKernel.cu"


	/**
	 * @brief ...
	 *
	 * @param h_targets ...
	 * @param h_velocities ...
	 * @param h_sources ...
	 * @param h_forces ...
	 * @param deltas ...
	 * @param num_sources ...
	 * @param num_targets ...
	 * @param num_clusters ...
	 * 
	 **/
	void ComputeStokesletsNearField( const float *h_targets, float *h_velocities, const float * h_sources, const float * h_forces, float delta, int num_sources, int num_targets, int num_clusters ) 
	{
        int size_targets = 3 * num_targets;
        int size_sources = 3 * num_sources;
        float *d_velocities, *d_targets, *d_sources, *d_strengths;
        cudaMalloc((void**)&d_velocities, sizeof( float) * size_targets);
        cudaMalloc((void**)&d_targets,    sizeof( float) * size_targets);
        cudaMalloc((void**)&d_sources,    sizeof( float) * size_sources);
        cudaMalloc((void**)&d_strengths,  sizeof( float) * size_sources);
        cudaMemcpy(d_velocities, h_velocities, size_targets*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_targets,    h_targets,    size_targets*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sources,    h_sources,    size_sources*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_strengths,  h_forces,  size_sources*sizeof(float), cudaMemcpyHostToDevice);

        int block_threads = num_sources/num_clusters;
        int num_blocks = num_sources / block_threads + (num_sources % block_threads == 0 ? 0 : 1);
        int sharedMemSize = 2 * block_threads * sizeof(float3);
        dim3 threads(block_threads, 1, 1);
        dim3 grid(num_blocks, 1, 1);
        ComputeStokesletsNearFieldKernel<<<grid,threads,sharedMemSize>>>(
			(float3 *)d_targets, 
			(float3 *)d_velocities, 
			(float3 *)d_sources, 
			(float3 *)d_strengths, 
			delta, 
			num_targets
		);
		//Copying data from host to device
        cudaMemcpy(h_velocities,  d_velocities,  size_targets*sizeof(float), cudaMemcpyDeviceToHost);
        cudaFree((void**)d_velocities);
        cudaFree((void**)d_targets);
        cudaFree((void**)d_sources);
        cudaFree((void**)d_strengths);
	}


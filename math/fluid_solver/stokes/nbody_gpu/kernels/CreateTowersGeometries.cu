
#include <iostream>
#include <iterator>

#include <vector>
#include "cuda.h"
#include "math_constants.h"
#include "CreateTowersGeometriesKernel.cu"
 
extern "C" {

template<int num_particles>
void CreateTowersGeometries(float *h_sources, int num_particles, int num_blocks_x, int num_blocks_y, int block_threads_x, int block_threads_y)
{

	// First, allocate space on device
	std::vector<int> h_faces(4*block_threads_x*(block_threads_y-1),0);
    int size_sources = num_particles*3;
    float *d_sources;
	int *d_faces;

    dim3 blocks(block_threads_x, block_threads_y, 1);
    dim3 grid(num_blocks_x,num_blocks_y, 1);
	int2 tower_dimensions = {block_threads_x,block_threads_y};
	float2 cylynder_dimensions = {20.0f,.25f};
	int2 grid_dimensions = {num_blocks_x,num_blocks_y};
	TowersGeometries<float,int,num_particles> geometries(tower_dimensions,cylynder_dimensions,grid_dimensions);
	CreateTowersGeometriesKernel<<<grid,blocks>>>(&geometries);

	blocks.y--;
	SetMeshKernel<<<1,blocks>>>((int4*)&h_faces[0]);

	cudaMemcpy(h_sources,  d_sources,  size_sources*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&h_faces[0],  d_faces,  h_faces.size()*sizeof(int), cudaMemcpyDeviceToHost);

// 	std::copy(h_sources,h_sources+size_sources,std::ostream_iterator<float>(std::cout ," "));
// 	std::cout << std::endl;
	std::copy(h_faces.begin(),h_faces.end(),std::ostream_iterator<int>(std::cout ," "));
	std::cout << std::endl;
	cudaFree((void**)d_sources);
	cudaFree((void**)d_faces);
}


}
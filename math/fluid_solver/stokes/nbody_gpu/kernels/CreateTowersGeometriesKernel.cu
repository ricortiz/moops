#include "math_constants.h"

__device__
float distance(const float3 &a, const float3 &b)
{
	float dx = b.x - a.x;
	float dy = b.y - a.y;
	float dz = b.z - a.z;

	return sqrt(dx*dx+dy*dy+dz*dz);
}


__device__
void connection_count(const int2 *springs_idx, int index)
{
	extern __shared__ int spring_counter[];
	int i = 0;
	while (springs_idx[i].x != index) 
		++i;

	int counter = 0;
	while(springs_idx[i++].x == index)
		++counter;

	spring_counter[index] = counter;

	__syncthreads();
}

__device__
int4 set_face(int index)
{
    int4 face;
    face.x = index;
    if (threadIdx.x == blockDim.x - 1)
    {
        face.y = face.x + 1 - blockDim.x;
        face.z = index + 1;
    }
    else
    {
        face.y = face.x + 1;
        face.z = index + blockDim.x + 1;
    }
    face.w = index + blockDim.x;
    return face;
}



__device__
float3 point_on_surface(float L, float r)
{
    float3 position;
    float dx = 2.0f;
    float dy = 2.0f;
    float dz = L / blockDim.y;
    float dt = 2.0f * CUDART_PI_F / blockDim.x;

    position.x = r * cosf(dt * threadIdx.x) + blockIdx.x * dx + r;
    position.y = r * sinf(dt * threadIdx.x) + blockIdx.y * dy + r;
    position.z = dz * threadIdx.y;

    return position;
}

__global__
void compute_resting_lengths(float3 *sources, int *idx_ptr, int *col_idx, float *resting_lengths)
{
	int index = threadIdx.y*blockDim.x + threadIdx.x;
	for(int i = idx_ptr[index]; i < idx_ptr[index+1]; ++i)
		resting_lengths[i] = distance(sources[index],sources[col_idx[i]]);
}


__global__
void get_adjacency_matrix(const int2 *springs_idx, int *col_ptr, int *col_idx)
{
	extern __shared__ int spring_counter[];
	int index = threadIdx.y*blockDim.x + threadIdx.x;
	
	connection_count(springs_idx,index);

	if(index == 0)
		col_ptr[index] = 0;

	int counter = 0;
	for(int i = 0; i <= index; ++i)
		counter += spring_counter[i];
	col_ptr[index+1] = counter;
	
	for(int i = col_ptr[index]; i < col_ptr[index+1]; ++i)
		col_idx[i] = springs_idx[i].y;
}

__global__
void create_connections(int2 *springs_idx, int size)
{
    int total_pairs = size/4;
    int index = __mul24(threadIdx.y, blockDim.x) + threadIdx.x;

    int2 springs[4];

    springs[0].x = index;
    if (threadIdx.x == blockDim.x - 1)
        springs[0].y = index + 1 - blockDim.x;
    else
        springs[0].y = index + 1;

    springs[1].x = index;
    springs[1].y = index + blockDim.x;

    springs[2].x = index;
    if (threadIdx.x == blockDim.x - 1)
        springs[2].y = index + blockDim.x + 1 - blockDim.x;
    else
        springs[2].y = index + blockDim.x + 1;

    if (threadIdx.x == blockDim.x - 1)
        springs[3].x = index + 1 - blockDim.x;
    else
        springs[3].x = index + 1;
    springs[3].y = index + blockDim.x;

    int idx = 4 * index;
    for (int i = 0; i < 4; ++i)
    {
        springs_idx[idx+i]   = springs[i];
        springs_idx[idx+total_pairs+i].x   = springs[i].y;
        springs_idx[idx+total_pairs+i].y   = springs[i].x;
    }

    __syncthreads();

    int  offset = __mul24(__mul24(4, blockDim.x), blockDim.y);
    int2 spring;
    if (threadIdx.y == blockDim.y - 1)
    {
        spring.x = index + blockDim.x;
        if (threadIdx.x == blockDim.x - 1)
            spring.y = spring.x + 1 - blockDim.x;
        else
            spring.y = spring.x + 1;
        idx = __mul24(threadIdx.y, blockDim.x) + threadIdx.x - (blockDim.y - 1) * blockDim.x;
        offset = __mul24(__mul24(4, blockDim.x), blockDim.y);
        springs_idx[idx+offset] = spring;
		springs_idx[idx+total_pairs+offset].x = spring.y;
		springs_idx[idx+total_pairs+offset].y = spring.x;
    }
    if (threadIdx.x < blockDim.x / 2)
    {
        spring.x = index;
        spring.y = index + 3;
        idx = __mul24(threadIdx.y, blockDim.x / 2) + threadIdx.x + blockDim.x;
        springs_idx[idx+offset] = spring;
		springs_idx[idx+total_pairs+offset].x = spring.y;
		springs_idx[idx+total_pairs+offset].y = spring.x;
    }
    if (threadIdx.y == blockDim.y - 1 && threadIdx.x < blockDim.x / 2)
    {
        spring.x = index + blockDim.x;
        spring.y = spring.x + 3;
        idx = __mul24(threadIdx.y, blockDim.x / 2) + threadIdx.x + blockDim.x + blockDim.x / 2;
        springs_idx[idx+offset] = spring;
		springs_idx[idx+total_pairs+offset].x = spring.y;
		springs_idx[idx+total_pairs+offset].y = spring.x;
    }

}

__global__
void create_faces(int4* faces)
{
    int index = __mul24(threadIdx.y, blockDim.x) + threadIdx.x;
    faces[index] = set_face(index);
}

__global__
void create_geometries(float3 *positions, float L, float r)
{
    int num_threads = __mul24(blockDim.x, blockDim.y);
    int block_index = __mul24(blockIdx.y, gridDim.x) + blockIdx.x;

    int tx = threadIdx.x + __mul24(block_index, num_threads);
    int ty = threadIdx.y;

    int thread_index = __mul24(ty, blockDim.x) + tx;
    positions[thread_index] = point_on_surface(L, r);
}


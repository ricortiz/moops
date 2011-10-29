
// Macros to simplify shared memory addressing
// the first half of the shared memory allocation is pos
#define sources(i) shared_data[i]
// and the second half is strength
#define forces(i) shared_data[i+blockDim.x*blockDim.y]

__device__
void compute_force(float3 &force, const float3 &sourceA, const float3 &sourceB, const float &resting_length, float spring_constant)
{
  float3 dx = {sourceB.x - sourceA.x, 
			   sourceB.y - sourceA.y, 
			   sourceB.z - sourceA.z};
  float magnitude = sqrt(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);
  
  float L = spring_constant * (magnitude/resting_length - 1.0f) / magnitude;
  force.x += L * dx.x;
  force.y += L * dx.y;
  force.z += L * dx.z;
}

__global__
void kernel_compute_force(float3 *force_array, const float3 *source_array, const int *col_ptr, const int *col_idx, const float *resting_lengths, float k)
{
  extern __shared__ float3 shared_data[];
  int grid_idx = blockIdx.y*gridDim.x + blockIdx.x;
  int local_idx = threadIdx.y*blockDim.x + threadIdx.x;
  
  int block_dim = blockDim.x*blockDim.y;
  int offset = grid_idx*block_dim;
  int global_idx = local_idx + offset;

  sources(local_idx)  = source_array[global_idx];
  forces(local_idx).x = 0.0f;
  forces(local_idx).y = 0.0f;
  forces(local_idx).z = 0.0f;

  __syncthreads();

  int start = col_ptr[local_idx];
  int end 	= col_ptr[local_idx+1];
  for(int i = start; i < end; ++i)
      compute_force(forces(local_idx),sources(local_idx),sources(col_idx[i]),resting_lengths[i], k);

  force_array[global_idx] = forces(local_idx);
}



#ifndef SYSTEM_CUH
#define SYSTEM_CUH

#define position_array(i) shared_data[i]
#define velocity_array(i) shared_data[i+3*num_particles_block]
#define forces_array(i) shared_data[i+6*num_particles_block]

template<class T>
struct SharedMemory
{
    __device__ 
	inline operator T*()
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }

    __device__ 
	inline operator const T*() const
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }

};

template<typename Types>
struct Particle 
{

  SharedMemory<typename Types::typle_type> shared_memory;

  __host__ __device__
  Particle(unsigned int i, typename Types::typle_type *position, typename Types::value_typle *velocities, typename Types::value_typle *forces) 
  {    
	
    position_array(i) = position[i];
    velocity_array(i) = velocity[i];
    forces_array(i)   = force[i];

    __syncthreads();

  }

  __host__ __device__
  typename Types::tuple_type position(Types::int_type i) {return position_array(i); }
  __host__ __device__
  typename Types::tuple_type velocity(Types::int_type i) {return velocity_array(i); }
  __host__ __device__
  typename Types::tuple_type force(Types::int_type i) 	  {return forces_array(i); }
};

template<typename Types, unsigned int num_particles, int dim = 3, int data_size = dim*num_particles>
class ParticleSystem
{
  __host__ __device__ 
  typename Types::value_type *m_positions;
  __host__ __device__ 
  typename Types::value_type *m_forces;
  __host__ __device__ 
  typename Types::value_type *m_velocities;
  __host__ __device__
  typename Types::particle_type m_particles[num_particles];

  public:
  
    __host__ __device__
    typename Types::particle_type *particle_begin() { return &m_particles[0];} 
    __host__ __device__
    typename Types::particle_type *particle_end() { return &m_particles[0]+m_num_particles; } 

    ParticleSystem(typename Types::value_type *positions,
		   typename Types::value_type *velocities,
		   typename Types::value_type *forces, ) : m_positions(positions),
							   m_velocities(velocities),
							   m_forces(forces)
    {
      

    }

    ParticleSystem() : m_position(0), m_velocities(0), m_forces(0)
    {
      cudaMalloc((void**)&m_velocities,sizeof(typename Types::value_type)*data_size);
      cudaMalloc((void**)&m_forcess,sizeof(typename Types::value_type)*data_size);
    }
   
};





#endif
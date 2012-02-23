#ifndef GPU_COMPUTE_VELOCITY_HPP
#define GPU_COMPUTE_VELOCITY_HPP

template<typename value_type>
void ComputeStokeslets(const value_type *, value_type *, const value_type *, const value_type *, value_type, size_t, size_t, bool with_images = false);

#endif

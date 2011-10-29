#ifndef GPU_COMPUTE_VELOCITY_HPP
#define GPU_COMPUTE_VELOCITY_HPP

// template<typename float>
void ComputeStokeslets(const float *, float *, const float *, const float *, float, size_t, size_t);

// template<typename float>
void ComputeImages(const float *, float *, const float *, const float *, float, size_t, size_t);

#endif

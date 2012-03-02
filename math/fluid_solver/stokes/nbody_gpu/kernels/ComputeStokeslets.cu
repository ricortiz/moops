/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
#include <cuda_runtime.h>
#include "../gpu_compute_velocity.hpp"
#include "ComputeStokesletsKernel.cu"

/**
 * @brief ...
 *
 * @param h_targets ...
 * @param h_sources ...
 * @param h_strengths ...
 * @param h_velocities ...
 * @param deltas ...
 **/
template<>
void ComputeStokeslets<float>(const float *h_targets, float *h_velocities, const float * h_sources, const float * h_strengths, float delta, size_t num_sources, size_t num_targets, bool with_images)
{
    cudaSetDeviceFlags(cudaDeviceMapHost);
    // First, allocate space on device
    size_t size_targets = 3 * num_targets;
    size_t size_sources = 3 * num_sources;
    float *d_velocities, *d_targets, *d_sources, *d_strengths;
    cudaHostRegister(const_cast<float*>(h_targets), sizeof(float) * size_targets, cudaHostRegisterMapped);
    cudaHostRegister(h_velocities, sizeof(float) * size_targets, cudaHostRegisterMapped);
    cudaHostRegister(const_cast<float*>(h_sources), sizeof(float) * size_sources, cudaHostRegisterMapped);
    cudaHostRegister(const_cast<float*>(h_strengths), sizeof(float) * size_sources, cudaHostRegisterMapped);

    cudaHostGetDevicePointer((void **) & d_targets, (void *)h_targets, 0);
    cudaHostGetDevicePointer((void **) & d_velocities, (void *)h_velocities, 0);
    cudaHostGetDevicePointer((void **) & d_sources, (void *)h_sources, 0);
    cudaHostGetDevicePointer((void **) & d_strengths, (void *)h_strengths, 0);

    size_t block_threads = 32;
    size_t num_blocks = num_targets / block_threads + (num_targets % block_threads == 0 ? 0 : 1);
    size_t sharedMemSize = 2 * block_threads * sizeof(float3);
    dim3 threads(block_threads, 1, 1);
    dim3 grid(num_blocks, 1, 1);
    ComputeStokesletsKernel<<<grid,threads,sharedMemSize>>>(
        (float3 *)d_targets,
        (float3 *)d_velocities,
        (float3 *)d_sources,
        (float3 *)d_strengths,
        delta,
        num_sources,
        num_targets, with_images
    );

    //Copying data from host to device
    cudaHostUnregister(const_cast<float*>(h_targets));
    cudaHostUnregister(h_velocities);
    cudaHostUnregister(const_cast<float*>(h_sources));
    cudaHostUnregister(const_cast<float*>(h_strengths));

}







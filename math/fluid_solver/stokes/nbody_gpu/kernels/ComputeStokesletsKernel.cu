#ifndef COMPUTE_STOKESLETS_KERNEL_CU
#define COMPUTE_STOKESLETS_KERNEL_CU
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
#include "ComputeTiles.cuh"


/**
 * @brief Main stokelets computation kernel.
 *
 * @param target_array array containing all targets
 * @param source_array array containing all sources
 * @param force_array array containing all force vectors
 * @param velocity_array array containing all velocity vectors
 * @param deltas regularization parameter
 * @param num_sources total number of sources 
 * @param num_targets total number of targets
 **/
template<typename vector3_type, typename int_type, typename real_type>
__global__ 
void ComputeStokesletsKernel(const vector3_type* target_array, vector3_type* velocity_array, const vector3_type* source_array, const vector3_type* force_array, real_type deltas, int_type num_sources, int_type num_targets, bool with_images = false)
{
	int_type index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
	if (index >= num_targets)
    	index = 0;  

	vector3_type target = target_array[index];
	vector3_type velocity = ComputeTiles(target, source_array, force_array, deltas, num_sources, num_targets, with_images);
	
	if (index < num_targets)
            velocity_array[index] = velocity;
}

#endif

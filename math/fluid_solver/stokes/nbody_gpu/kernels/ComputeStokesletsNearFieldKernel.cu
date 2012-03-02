#ifndef COMPUTE_STOKESLETS_NEAR_FIELD_KERNEL_CU
#define COMPUTE_STOKESLETS_NEAR_FIELD_KERNEL_CU
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
#include "ComputeStokesletsKernel.cu"

template<typename vector3_type, typename int_type, typename real_type>
__global__
void ComputeStokesletsNearFieldKernel(vector3_type* target_array, vector3_type* velocity_array, vector3_type* source_array, vector3_type* force_array, real_type delta, int_type  num_targets, bool with_images)
{
	int_type index_start = __mul24(blockIdx.x,blockDim.x);
    int_type index = index_start + threadIdx.x;
    
    int_type num_sources = blockDim.x;
	vector3_type target = target_array[index];
    vector3_type velocity = ComputeTiles(target, &source_array[index_start], &force_array[index_start], num_sources, num_targets, delta, with_images);

    if (index < num_targets)
    {
        velocity_array[index].x += velocity.x;      //returning addition of far & far_image
        velocity_array[index].y += velocity.y;
        velocity_array[index].z += velocity.z;
    }
}

#endif

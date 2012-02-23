#ifndef COMPUTE_TILES_CUH
#define COMPUTE_TILES_CUH

/**
 * @brief This function is the fine grained stokeslet computation.
 *      It computes the velocity at the target point due to a delta-regularized force at source point
 *
 * @param velocity velocity vector
 * @param target target point
 * @param source source point
 * @param force force vector
 * @param delta regularization parameter
 * @return updated velocity vector
 **/
template<typename vector3_type, typename real_type> __device__
void ComputeStokeslets(const vector3_type &target, vector3_type &velocity, const vector3_type &source, const vector3_type &force, real_type delta)
{
    vector3_type dx = {target.x - source.x,target.y - source.y,target.z - source.z};

    real_type r2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;
    real_type d2 = delta * delta;
    real_type R1 = r2 + d2;
    real_type R2 = R1 + d2;
    real_type invR = 1.0f/R1;
    real_type H = 0.5f * sqrtf(invR)*invR * 7.957747154594767e-02f;

    real_type fdx = (force.x*dx.x+force.y*dx.y+force.z*dx.z);

    velocity.x += H*(force.x*R2+fdx*dx.x);
    velocity.y += H*(force.y*R2+fdx*dx.y);
    velocity.z += H*(force.z*R2+fdx*dx.z);
}

/**
 * @brief
 *
 * @param velocity velocity vector
 * @param target target point
 * @param source source point
 * @param force force vector
 * @param delta regularization parameter
 * @return updated velocity vector
 **/
template<typename vector3_type, typename real_type> __device__
void ComputeImages(const vector3_type &target, vector3_type &velocity, const vector3_type &source, const vector3_type &force, real_type delta)
{
    real_type h0 = source.z;

    vector3_type dx;
    dx.x = target.x - source.x;
    dx.y = target.y - source.y;
    dx.z = target.z + h0; // image part

    real_type r2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;

    real_type d2 = delta * delta;

    real_type H1 = 0.5f / sqrtf(r2 + d2) + d2 / 2 / powf( sqrtf(r2 + d2), 3);
    real_type H2 = 0.5f / powf( sqrtf(r2 + d2), 3);

    real_type Hdp1 = -h0*h0*( 1.0f/powf( sqrtf(r2 + d2), 3) - 3*d2/powf( sqrtf(r2 + d2), 5) );
    real_type Hdp2 = -h0*h0*( -3.0f/powf( sqrtf(r2 + d2), 5) );

    real_type Hdb1 = 2*h0*( -3.0f / 2 * d2 / powf( sqrtf(r2 + d2), 5) - 1.0f / 2 / powf( sqrtf(r2 + d2), 3) );
    real_type Hdb2 = 2*h0*( 1.0f / 2 / powf( sqrtf(r2 + d2), 3) );
    real_type Hdb3 = 2*h0*( -3.0f / 2 / powf( sqrtf(r2 + d2), 5) );

    real_type Hr1 = h0 * 3 * d2 / powf( sqrtf(r2 + d2), 5);

    // 21 flops
    real_type A11 = -(H1 + dx.x * dx.x * H2);
    real_type A12 = -(dx.x * dx.y * H2);
    real_type A13 = -(dx.x * dx.z * H2);

    real_type A21 = -(dx.y * dx.x * H2);
    real_type A22 = -(H1 + dx.y * dx.y * H2);
    real_type A23 = -(dx.y * dx.z * H2);

    real_type A31 = -(dx.z * dx.x * H2);
    real_type A32 = -(dx.z * dx.y * H2);
    real_type A33 = -(H1 + dx.z * dx.z * H2);

    //% A pos Dipole
    A11 = A11-( Hdp1 + dx.x*dx.x*Hdp2 ); // %i=1 j=1
    A12 = A12-( dx.x*dx.y*Hdp2 );          //  %i=1 j=2
    A13 = A13+( dx.x*dx.z*Hdp2 );          //  %i=1 j=3

    A21 = A21-( dx.y*dx.x*Hdp2 );          // %i=2 j=1
    A22 = A22-( Hdp1 + dx.y*dx.y*Hdp2 ); //%i=2 j=2
    A23 = A23+ ( dx.y*dx.z*Hdp2 );        //  %i=2 j=3

    A31 = A31-( dx.z*dx.x*Hdp2 );         //  %i=3 j=1
    A32 = A32-( dx.z*dx.y*Hdp2 );         //  %i=3 j=2
    A33 = A33+( Hdp1 + dx.z*dx.z*Hdp2 ); //%i=3 j=3

    //% A pos Doublet
    A11 = A11-( dx.z*(Hdb2 + dx.x*dx.x*Hdb3) ) ;
    A12 = A12-( dx.x*dx.y*dx.z*Hdb3 );     //% i=1 j=2
    A13 = A13+dx.x*(Hdb2 + dx.z*dx.z*Hdb3 );  // % i=1 j=3

    A21 = A21-( dx.y*dx.x*dx.z*Hdb3 ) ;          // % i=2 j=1
    A22 = A22-( dx.z*(Hdb2 + dx.y*dx.y*Hdb3) ) ;
    A23 = A23+dx.y*(Hdb2+ dx.z*dx.z*Hdb3 ) ; //% i=2 j=3

    A31 = A31-( dx.x*(Hdb1 + dx.z*dx.z*Hdb3 ) );
    A32 = A32-( dx.y*(Hdb1 + dx.z*dx.z*Hdb3 ) );
    A33 = A33+ dx.z*(Hdb1 + 2*Hdb2 + dx.z*dx.z*Hdb3);

    //% A pos Rotlet
    A11 = A11 + dx.z*Hr1;
    A22 = A22 + dx.z*Hr1;
    A31 = A31-dx.x*Hr1;
    A32 = A32-dx.y*Hr1;
// TODO: check the constant term correctness!
    velocity.x += (A11 * force.x + A12 * force.y + A13 * force.z) * 7.957747154594767e-02f;
    velocity.y += (A21 * force.x + A22 * force.y + A23 * force.z) * 7.957747154594767e-02f;
    velocity.z += (A31 * force.x + A32 * force.y + A33 * force.z) * 7.957747154594767e-02f;
}

/**
 * @brief This is the "tile_calculation" function from the GPUG3 article.
 *
 * @param target target point
 * @param velocity velocity vector to update
 * @param num_sources total number of sources in the tile
 * @param delta regularization parameter
 * @return updated velocity vector
 **/
template<typename vector3_type, typename int_type, typename real_type> __device__
void ComputeTile(const vector3_type &target, vector3_type &velocity, int_type num_sources, real_type delta, bool with_image)
{
    extern __shared__ vector3_type shared_data[];
    vector3_type *sources = &shared_data[0];
    vector3_type *forces  = &shared_data[blockDim.x];
    // The CUDA 1.1 compiler cannot determine that i is not going to  overflow in the loop below.
        // Therefore if int is used on 64-bit linux or windows (or long instead of long long on win64),
        // the compiler generates suboptimal code.
        // Therefore, we use long long on win64 and  long on everything else. (Workaround for Bug ID 347697)
        #ifdef _Win64
                unsigned long long i = 0;
        #else
                unsigned long i = 0;
        #endif
    // Note that having an unsigned int loop counter and an unsigned long index helps the compiler generate efficient code on 64-bit  OSes.
        // The compiler can't assume the 64-bit index won't overflow     // so it incurs extra integer operations.  This is a standard issue
        // in porting 32-bit code to 64-bit OSes.
    for (unsigned int counter = 0; counter < num_sources; ++i, ++counter )
        ComputeStokeslets(target, velocity, sources[i], forces[i], delta);
    if(with_image)
        for (unsigned int counter = 0; counter < num_sources; ++i, ++counter )
            ComputeImages(target, velocity, sources[i], forces[i], delta);


}

/**
 * @brief This function sets up computation tiles to update velocities.
 *
 * @param target target where verlocity is going to be updated
 * @param source_array array of point where forces are exceted
 * @param force_array array of force vectors
 * @param num_sources total number of sources
 * @param num_targets total number of targets
 * @param delta regularization parameter
 * @return updated velocity vector
 **/
template<typename vector3_type, typename int_type, typename real_type> __device__
vector3_type ComputeTiles(const vector3_type &target, const vector3_type* source_array, const vector3_type* force_array, int_type num_sources, int_type num_targets, real_type delta, bool with_image)
{
    extern __shared__ vector3_type shared_data[];
    vector3_type *sources = &shared_data[0];
    vector3_type *forces  = &shared_data[blockDim.x];

    vector3_type velocity = {0.0f, 0.0f, 0.0f};

    int_type p = blockDim.x;
    int_type q = blockDim.y;

    int_type n = num_targets;
    int_type numTiles = n / (p * q);
    int_type index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    // all blocks must start from index 0.
        if (index >= num_targets)
                index = 0;

    for (int_type tile = 0; tile < numTiles; ++tile)
        {
        int_type posIdx = tile * blockDim.x + threadIdx.x;
        sources[threadIdx.x] = source_array[posIdx];
        forces[threadIdx.x] = force_array[posIdx];

        __syncthreads();

        if (index < num_targets)
            ComputeTile(target, velocity, blockDim.x, delta, with_image);

        __syncthreads();
    }

    // handle if num_targets != 32 * n  (the last tile)
    int_type bodyLeft = n - numTiles * (p * q);
    if (bodyLeft > 0){
        // only 'bodyLeft' threads have to write the entries in share memory, (32-bodyLeft) threads are idle.
        if (threadIdx.x < bodyLeft)
        {
            int_type posIdx = numTiles * blockDim.x + threadIdx.x;
            sources[threadIdx.x] = source_array[posIdx];
            forces[threadIdx.x] = force_array[posIdx];
        }
        __syncthreads();

        if (index < num_targets)
            ComputeTile(target, velocity, bodyLeft, delta, with_image);

        __syncthreads();
    }
    return velocity;
}


#endif

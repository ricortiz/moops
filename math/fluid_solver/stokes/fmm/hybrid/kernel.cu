#include <stdio.h>
#include <stdlib.h>
#include<cuda_runtime.h>
#include <cuda.h>
#include <sys/time.h>

//*****************************USE THESE PARAMETERS TO CHANGE THE BEHAVIOR********************************

// #define IGNORE_FIRST_DEVICE             //ignore first device, (on systems where first device might be hooked to display) or might not be a faster GPU
// #define LOAD_BALANCING                  //to enable or disable load balancing, if disabled each GPU executes same amount of CUDA blocks, which might change actual load

// #define CUDA_LOG                        //enable logs (minimal)
// #define GPU_INFO                        //print which GPUs are being used

// #define ENABLE_BLOCKING                 //CPU thread that invokes cuda gets blocked
// #define CPU_EXECUTE                     //execute on CPU as well, just for validating GPU results with CPU
// #define CUDA_DEBUG                      //detailed GPU logging

//globals
FILE *fp;
float *GL_positions, *GL_gpuVelocities, *cpu_velocities;
unsigned int GL_NUM_BODIES, *source_list_temp, *source_start_indices_for_blocks_temp, *dnum_interaction_pairs;
unsigned int *source_list, *source_start_indices_for_blocks, *GL_num_interaction_pairs, *GL_target_list;
double startTime;


__device__ float3 evaluate(float3 a, float3 b, float3 force, float3 velocity, float delta)
{
    float3 dx = {b.x - a.x, b.y - a.y, b.z - a.z};

    float r2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;
    float d2 = delta * delta;
    float R1 = r2 + d2;
    float R2 = R1 + d2;
    float invR = 1.0f / R1;
    float H = sqrtf(invR) * invR * 0.039788735772974;

    float fdx = (force.x * dx.x + force.y * dx.y + force.z * dx.z);

    velocity.x += H * (force.x * R2 + fdx * dx.x);
    velocity.y += H * (force.y * R2 + fdx * dx.y);
    velocity.z += H * (force.z * R2 + fdx * dx.z);

    return velocity;
}

__device__ float3 TileCalculate(float3 pos, float3 velocity, int numBody, int NUM_THREADS, float delta)
{
    extern __shared__ float3 sharedPos[];
#ifdef _Win64
    unsigned long long i = 0;
#else
    unsigned long i = 0;
#endif
    for (unsigned int counter = 0; counter < numBody; i++, counter++)
    {
        velocity = evaluate(pos, sharedPos[i], sharedPos[i+NUM_THREADS], velocity, delta);
    }
    return velocity;
}


__device__ float3 compute_velocity(float3 current_target_pos,
                                   float3 * positions,
                                   unsigned int * source_list,
                                   int num_source_pairs_to_interact_with,
                                   int source_start,
                                   int total_num_sources_for_this_loop,
                                   int num_targets_for_this_loop,
                                   int NUM_THREADS, float delta)
{
    int num_threads_in_this_block = blockDim.x * blockDim.y;
    extern __shared__ float3 sharedPos[];

    float3 vel = {0.0f, 0.0f, 0.0f};
    int num_horiz_tiles = total_num_sources_for_this_loop / num_threads_in_this_block;
    int num_sources_completed_for_this_loop = 0;
    int index_to_load_source_from = 0;
    int subtract_offset = 0;

    for (int current_tile = 0; current_tile < num_horiz_tiles; current_tile++)
    {

        int current_source = current_tile * num_threads_in_this_block + threadIdx.x;

        //calculate the index from where sources will be loaded
        num_sources_completed_for_this_loop = 0;
        subtract_offset = 0;
        for (int i = 0; i < num_source_pairs_to_interact_with; i++)
        {
            num_sources_completed_for_this_loop += source_list[source_start + (2*i+1)];
            if (current_source < num_sources_completed_for_this_loop)
            {
                index_to_load_source_from = source_list[source_start + (2*i)] + current_source - subtract_offset;
                break;
            }
            subtract_offset += source_list[source_start + (2*i+1)];
        }

        sharedPos[threadIdx.x] = positions[2*index_to_load_source_from];
        sharedPos[threadIdx.x+NUM_THREADS] = positions[2*index_to_load_source_from+1];

        __syncthreads();

        if (threadIdx.x < num_targets_for_this_loop)
            vel = TileCalculate(current_target_pos, vel, num_threads_in_this_block, NUM_THREADS, delta);

        __syncthreads();
    }

    int num_sources_left = total_num_sources_for_this_loop - num_horiz_tiles * num_threads_in_this_block;
    if (num_sources_left > 0)
    {
        if (threadIdx.x < num_sources_left)
        {
            int current_source = num_horiz_tiles * num_threads_in_this_block + threadIdx.x;

            //calculate the index from where sources will be loaded
            num_sources_completed_for_this_loop = 0;
            subtract_offset = 0;
            for (int i = 0; i < num_source_pairs_to_interact_with; i++)
            {
                num_sources_completed_for_this_loop += source_list[source_start + (2*i+1)];
                if (current_source < num_sources_completed_for_this_loop)
                {
                    index_to_load_source_from = source_list[source_start + (2*i)] + current_source - subtract_offset;
                    break;
                }
                subtract_offset += source_list[source_start + (2*i+1)];
            }

            sharedPos[threadIdx.x] = positions[2*index_to_load_source_from];
            sharedPos[threadIdx.x+NUM_THREADS] = positions[2*index_to_load_source_from+1];

        }
        __syncthreads();

        if (threadIdx.x < num_targets_for_this_loop)
            vel = TileCalculate(current_target_pos, vel, num_sources_left, NUM_THREADS, delta);
        __syncthreads();
    }
    return vel;
}

__global__ void kernel(int block_offset,
                       float3 * positions,
                       float3 * velocities,
                       unsigned int * target_list,
                       unsigned int * source_list,
                       unsigned int * num_interaction_pairs,
                       unsigned int * source_start_array,
                       int NUM_THREADS, float delta)
{
    int blockId = block_offset + blockIdx.x;
    int target_start = target_list[blockId * 2];
    int num_targets = target_list[blockId * 2 + 1];
    int source_start_for_this_block = source_start_array[blockId] * 2;
    int num_source_pairs_to_interact_with = num_interaction_pairs[blockId];
    if (!num_source_pairs_to_interact_with)
        return;
    int total_num_sources = 0;


    for (int i = 0; i < num_source_pairs_to_interact_with; i++)
        total_num_sources += source_list[source_start_for_this_block + (2*i+1)];

    int num_vertical_tiles = num_targets / (blockDim.x  * blockDim.y);
    float3 current_target;
    float3 vel;

    for (int i = 0; i < num_vertical_tiles; i++)
    {
        current_target = positions[ 2*(target_start + i*blockDim.x + threadIdx.x)];

        velocities[ target_start + i*blockDim.x + threadIdx.x ] = compute_velocity(current_target,
                positions,
                source_list,
                num_source_pairs_to_interact_with,
                source_start_for_this_block,
                total_num_sources,
                blockDim.x * blockDim.y,
                NUM_THREADS, delta);
    }

    int targetsLeft = num_targets - num_vertical_tiles * (blockDim.x * blockDim.y);
    if (targetsLeft)
    {
        if (threadIdx.x < targetsLeft)
            current_target = positions[ 2*(target_start + num_vertical_tiles*blockDim.x + threadIdx.x)];

        //call for all threads as they have to load data in shared memory
        vel = compute_velocity(current_target,
                               positions,
                               source_list,
                               num_source_pairs_to_interact_with,
                               source_start_for_this_block,
                               total_num_sources,
                               targetsLeft,
                               NUM_THREADS, delta);

        if (threadIdx.x < targetsLeft)
            velocities[ target_start + num_vertical_tiles*blockDim.x + threadIdx.x ] = vel;
    }
}

// inline unsigned int getAlignedSize(unsigned int size)
// {
//     return (size + MEMORY_ALIGNMENT - (size % MEMORY_ALIGNMENT));
// }

inline double wcTime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

inline void checkError(cudaError_t err)
{
#ifdef CUDA_DEBUG
    if (err != cudaSuccess)
        printf("CUDA_DEBUG:: cudaError = %s\n", cudaGetErrorString(err));
#endif
}



extern "C"
{

    void gpuInit()
    {
        cudaSetDevice(0);
        cudaSetDeviceFlags(cudaDeviceMapHost);
    }

    double gpuGetVelocities()
    {

        //just for making it async
        cudaMemcpy(GL_num_interaction_pairs, (unsigned int *)dnum_interaction_pairs, 4*sizeof(unsigned int), cudaMemcpyDeviceToHost);

        //unregister
        cudaError_t err;
        err = cudaHostUnregister(GL_positions);
        checkError(err);
        err = cudaHostUnregister(GL_gpuVelocities);
        checkError(err);
        err = cudaHostUnregister(GL_target_list);
        checkError(err);
        err = cudaHostUnregister(GL_num_interaction_pairs);
        checkError(err);
        err = cudaHostUnregister(source_list);
        checkError(err);
        err = cudaHostUnregister(source_start_indices_for_blocks);
        checkError(err);

        //free CPU memory
        free(source_list_temp);
        free(source_start_indices_for_blocks_temp);

#ifdef CUDA_DEBUG
        fprintf(fp, "\n\nGPU Velocities:\n");
        for (int i = 0; i < (3*GL_NUM_BODIES);i += 3)
        {
            fprintf(fp, "Velocity %d: %.10f %.10f %.10f\n", i / 3, GL_gpuVelocities[i], GL_gpuVelocities[i+1], GL_gpuVelocities[i+2]);
        }
#ifdef CPU_EXECUTE
        float total_deviation = 0.0f , max_deviation = 0.0f, deviation = 0.0f, cpuvelocity = 0.0f, gpuvelocity = 0.0f;
        int body = 0;
        for (int i = 0; i < (3*GL_NUM_BODIES);i += 3)
        {
            fprintf(fp, "Diff Velocity %d: %.10f %.10f %.10f\n", i / 3, GL_gpuVelocities[i] - cpu_velocities[i], GL_gpuVelocities[i+1] - cpu_velocities[i+1], GL_gpuVelocities[i+2] - cpu_velocities[i+2]);
            deviation = fabs(cpu_velocities[i] - GL_gpuVelocities[i]);
            total_deviation += deviation;
            if (deviation > max_deviation)
            {
                max_deviation = deviation;
                cpuvelocity = cpu_velocities[i];
                gpuvelocity = GL_gpuVelocities[i];
                body = i / 3;
            }
        }
        free(cpu_velocities);
        fprintf(stdout, "CUDA_DEBUG::Total Deviation: %f\n", total_deviation);
        fprintf(stdout, "CUDA_DEBUG::Max Deviation in body %d, (Deviation) = %f, cpu velocity = %f, gpu velocity = %f \n", body, max_deviation, cpuvelocity, gpuvelocity);
#endif
        fclose(fp);
#endif
        double t = wcTime() - startTime;
#ifdef CUDA_LOG
        fprintf(stdout, "CUDA_LOG::Total GPU time: %f\n", t);
        fprintf(stdout, "--------------------------------------------------\n");
#endif
        return t;
    }


    void gpuVelocitiesEval(const unsigned int NUM_BODIES,
                           const unsigned int NUM_LEAF_NODES,
                           const unsigned int TOTAL_NUM_SOURCES,
                           float * positions,
                           float * gpuVelocities,
                           unsigned int * target_list,
                           unsigned int * num_interaction_pairs,
                           unsigned int * interaction_pairs, float delta)
    {

        startTime = wcTime();
#ifdef CUDA_LOG
        fprintf(stdout, "--------------------------------------------------\n");
        fprintf(stdout, "CUDA_LOG::GPU Start time: %f\n", wcTime() - startTime);
#endif

        cudaSetDevice(0);
        cudaSetDeviceFlags(cudaDeviceMapHost);

#ifdef CUDA_LOG
        fprintf(stdout, "CUDA_LOG::Actual GPU start time: %f\n", wcTime() - startTime);
        startTime = wcTime();
#endif

#ifdef CUDA_DEBUG
        fp = fopen("LogGPUhost.txt", "w");

        //quit if there is an overlap in target_list
        bool overlap_check[NUM_BODIES];
        for (int i = 0; i < NUM_BODIES; i++)
            overlap_check[i] = false;
        for (int i = 0; i < NUM_LEAF_NODES; i++)
        {
            for (int j = target_list[2*i]; j < (target_list[2*i] + target_list[2*i+1]) ; j++)
            {
                if (overlap_check[j])
                {
                    printf("CUDA_DEBUG::Overlap found.. Terminating\n");
                    return;
                }
                overlap_check[j] = true;
            }
        }
        printf("CUDA_DEBUG::No Overlap found.\n");
#endif


#ifdef CPU_EXECUTE
        FILE *cfp;

        cpu_velocities  = (float *) malloc(sizeof(float) * 3 * NUM_BODIES);
        memset(cpu_velocities, 0, sizeof(float)*3*NUM_BODIES);

        int counter = 0 ;
        for (int i = 0; i < NUM_LEAF_NODES; i++)
        {

            fprintf(cfp, "\n\nLeaf node: %d", i);
            int target_index_start = target_list[2*i];
            int target_size  = target_list[2*i + 1];

            int j;
            for (j = 0; j < num_interaction_pairs[i]; j++)
            {

                counter = 0;
                for (int c = 0; c < i; c++)
                    counter += num_interaction_pairs[c];

                int source_body = interaction_pairs[counter + j];
                int source_index_start = target_list[2*source_body];
                int source_size = target_list[2*source_body + 1];

                fprintf(cfp, "\nInteracting target (%d): %d %d  with source (%d) %d %d", i, target_index_start, target_size, source_body, source_index_start, source_size);

                for (int k = target_index_start; k < (target_index_start + target_size) ; k++)
                {

                    for (int l = source_index_start; l < (source_index_start + source_size) ; l++)
                    {

                        float rx = positions[6*l]   - positions[6*k];
                        float ry = positions[6*l+1] - positions[6*k+1];
                        float rz = positions[6*l+2] - positions[6*k+2];

                        float r2 = rx * rx + ry * ry + rz * rz ;
                        float d2 = delta * delta;
                        float R1 = r2 + d2;
                        float R2 = R1 + d2;
                        float invR = 1.0 / R1;
                        float H = sqrt(invR) * invR * 0.039788735772974;

                        float fdx =  positions[6*l+3] * rx + positions[6*l+4] * ry + positions[6*l+5] * rz;

                        cpu_velocities[3*k]   += H * (positions[6*l+3] * R2 + fdx * rx);
                        cpu_velocities[3*k+1] += H * (positions[6*l+4] * R2 + fdx * ry);
                        cpu_velocities[3*k+2] += H * (positions[6*l+5] * R2 + fdx * rz);

                    }
                }
            }
        }

        printf("cpu_velocities = [");
        for (int i = 0; i < (3*NUM_BODIES);++i)
        {
            printf("%f ", cpu_velocities[i]);
        }
        printf("]\n");

#endif

        //create source lists
        unsigned long src_list_count = 0;      //get size for source_list array
        unsigned long total_work = 0;
        unsigned num_sources_per_target = 0;

#ifdef CUDA_DEBUG
        printf("CUDA_DEBUG::Malloc on CPU for source_list : %f MB and source_start_indices_for_blocks = %f MB\n",
               (float)sizeof(unsigned int)*2.0f* TOTAL_NUM_SOURCES / (1024.0f*1024.0f),
               (float)sizeof(unsigned int)* NUM_LEAF_NODES / (1024.0f*1024.0f));
#endif

        source_start_indices_for_blocks_temp = (unsigned int *) malloc(sizeof(unsigned int) * (NUM_LEAF_NODES + 1) /*+ MEMORY_ALIGNMENT*/);
        source_list_temp = (unsigned int *) malloc(sizeof(unsigned int) * 2 * TOTAL_NUM_SOURCES/* + MEMORY_ALIGNMENT*/);
        if (source_list_temp == NULL || source_start_indices_for_blocks_temp == NULL)
        {
            printf("CUDA ::Not enough memory on CPU.. Exiting\n");
            return;
        }

        //align
        source_list = (unsigned int *) /*ALIGN_UP(*/ source_list_temp/*, MEMORY_ALIGNMENT )*/;
        source_start_indices_for_blocks = (unsigned int *) /*ALIGN_UP( */source_start_indices_for_blocks_temp/*, MEMORY_ALIGNMENT )*/;

        //copy source_list and source_indices array
        unsigned int j = 0, k = 0, num_sources_to_interact_with = 0;
        unsigned int nodeID_to_interact_with = 0;
        source_start_indices_for_blocks[0] = 0;
        for (int i = 0; i < NUM_LEAF_NODES; ++i)
        {
            num_sources_to_interact_with = num_interaction_pairs[i];
            source_start_indices_for_blocks[i+1] = source_start_indices_for_blocks[i] + num_interaction_pairs[i];
            num_sources_per_target = 0;
            for (j = 0; j < num_sources_to_interact_with ; ++j)
            {
                nodeID_to_interact_with = interaction_pairs[j+k];
                source_list[src_list_count] = target_list[2*nodeID_to_interact_with];
                source_list[src_list_count+1] = target_list[2*nodeID_to_interact_with+1];
                num_sources_per_target += source_list[src_list_count+1];
                src_list_count += 2;
            }
            total_work += num_sources_per_target * target_list[2*i+1];
            k += j;
        }

        //print source lists created
#ifdef CUDA_DEBUG
        int max_count = 0;
        printf("Total Work = %lu\n", total_work);
        for (int i = 0; i < NUM_LEAF_NODES; i++)
        {
            printf(/*fp,*/ "\n\nLeaf node: %d", i);
            int source_pair_offset = 2 * source_start_indices_for_blocks[i];
            for (int j = 0; j < num_interaction_pairs[i] ; j++)
            {
                printf(/*fp,*/ "\nInteracting target (%d): %d %d  with source (%d) %d %d",
                               i, target_list[2*i], target_list[2*i+1], interaction_pairs[max_count] , source_list[source_pair_offset], source_list[source_pair_offset+1]);
                source_pair_offset += 2;
                max_count++;
            }
        }
#endif
        if (src_list_count != 2*TOTAL_NUM_SOURCES)
        {
            printf("CUDA ::Incorrect copy in Sources...Exiting\n");
            return;
        }
        else
        {
#ifdef CUDA_DEBUG
            printf("CUDA_DEBUG::Copied Sources Correctly...\n");
#endif
        }


        //GPU Execution begin - create GPU pointers
        float *dpositions, *dvelocities;
        unsigned int * dtarget_list , *dsource_list, *dsource_start_indices_for_blocks;

        //check all memory is aligned
#ifdef CUDA_DEBUG
//         if ((ulong) positions % MEMORY_ALIGNMENT != 0)
//             printf("CUDA_DEBUG::Positions not aligned correctly \n");
//         if ((ulong)gpuVelocities % MEMORY_ALIGNMENT != 0)
//             printf("CUDA_DEBUG::gpuVelocities not aligned correctly \n");
//         if ((ulong)target_list % MEMORY_ALIGNMENT != 0)
//             printf("CUDA_DEBUG::target_list not aligned correctly \n");
//         if ((ulong)num_interaction_pairs % MEMORY_ALIGNMENT != 0)
//             printf("CUDA_DEBUG::num_interaction_pairs not aligned correctly \n");
//         if ((ulong)source_list % MEMORY_ALIGNMENT != 0)
//             printf("CUDA_DEBUG::source_list not aligned correctly \n");
//         if ((ulong)source_start_indices_for_blocks % MEMORY_ALIGNMENT != 0)
//             printf("CUDA_DEBUG::source_start_indices_for_blocks not aligned correctly \n");
#endif

        //register memory
        cudaError_t err;
        err =   cudaHostRegister(positions, sizeof(float) * 6 * NUM_BODIES, cudaHostRegisterMapped);
        checkError(err);
        err =   cudaHostRegister(gpuVelocities, sizeof(float) * 3 * NUM_BODIES, cudaHostRegisterMapped);
        checkError(err);
        err =   cudaHostRegister(target_list, sizeof(unsigned int) * 2 * NUM_LEAF_NODES, cudaHostRegisterMapped);
        checkError(err);
        err =   cudaHostRegister(num_interaction_pairs, sizeof(unsigned int) * NUM_LEAF_NODES, cudaHostRegisterMapped);
        checkError(err);
        err =   cudaHostRegister(source_list, sizeof(unsigned int) * 2 * TOTAL_NUM_SOURCES, cudaHostRegisterMapped);
        checkError(err);
        err =   cudaHostRegister(source_start_indices_for_blocks, sizeof(unsigned int) * (NUM_LEAF_NODES + 1), cudaHostRegisterMapped);
        checkError(err);

        err =   cudaHostGetDevicePointer((void **) & dpositions, (void *)positions, 0);
        checkError(err);
        err =   cudaHostGetDevicePointer((void **) & dvelocities, (void *)gpuVelocities, 0);
        checkError(err);
        err =   cudaHostGetDevicePointer((void **) & dtarget_list, (void *)target_list, 0);
        checkError(err);
        err =   cudaHostGetDevicePointer((void **) & dnum_interaction_pairs, (void *)num_interaction_pairs, 0);
        checkError(err);
        err =   cudaHostGetDevicePointer((void **) & dsource_list, (void *)source_list, 0);
        checkError(err);
        err =   cudaHostGetDevicePointer((void **) & dsource_start_indices_for_blocks, (void *)source_start_indices_for_blocks, 0);
        checkError(err);

        //CUDA BLOCK params
        int NUM_THREADS = 64;
        int MAX_BLOCK_SIZE = 65535;
        int SHARED_MEM_SIZE = 2 * NUM_THREADS * sizeof(float3);
        if (SHARED_MEM_SIZE > 16384)
        {
            printf("CUDA ::Exceeded shared memory limit: Shared Mem Size = %d KB. Exiting..\n", SHARED_MEM_SIZE / 1024);
            return;
        }
        dim3 threads(NUM_THREADS, 1, 1);

#ifdef CUDA_DEBUG
        printf("CUDA_DEBUG::NUM_LEAF_NODES %d \n", NUM_LEAF_NODES);
        printf("CUDA_DEBUG::Num threads = %d, Max Block size = %d, Shared Mem Size = %d KB\n", NUM_THREADS, MAX_BLOCK_SIZE, SHARED_MEM_SIZE / 1024);
#endif

        int num_gpus = 1;
        cudaGetDeviceCount(&num_gpus);

#ifdef IGNORE_FIRST_DEVICE
        if (num_gpus > 1) // WARNING: rortiz - Make sure num_gpus > 1 before doing this
            num_gpus--;         //ignore the default (0th) device that is hooked to display
#endif

        total_work /= num_gpus;
#ifdef CUDA_LOG
        printf("CUDA_LOG::NUM GPUs : %d \n", num_gpus);
#endif

        //load balancing of work on GPUs
        // WARNING: rortiz - Initialize to zero
        int *load_balance_offset = (int *)malloc(sizeof(int) * num_gpus);
        int *load_balance_length = (int *)malloc(sizeof(int) * num_gpus);
        for (int i = 0; i < num_gpus; ++i)
        {
            load_balance_offset[i] = 0;
            load_balance_length[i] = 0;
        }
#ifdef LOAD_BALANCING
        j = 0; k = 0;
        unsigned long work = 0;
        int load_counter = 0;
        for (int i = 0; i < NUM_LEAF_NODES; ++i)
        {
            num_sources_per_target = 0;
            for (j = 0; j < num_interaction_pairs[i] ; ++j)
                num_sources_per_target += target_list[2*interaction_pairs[j+k] + 1];
            work = work + target_list[2*i+1] * num_sources_per_target;
            if (work > total_work)
            {
                load_balance_offset[load_counter++] = i + 1;
                work = 0;
                if ((load_counter + 1) == num_gpus)
                {
                    load_balance_offset[load_counter] = NUM_LEAF_NODES;
                    break;
                }
            }
            k += j;
        }
#else
        int constant_offset = NUM_LEAF_NODES / num_gpus + 1;
        for (int i = 0; i < num_gpus - 1; i++)
            load_balance_offset[i] = (i + 1) * constant_offset;
        load_balance_offset[num_gpus-1] = NUM_LEAF_NODES;
#endif

        load_balance_length[0] = load_balance_offset[0];
        for (int i = 1; i < num_gpus; i++)
            load_balance_length[i] = load_balance_offset[i] - load_balance_offset[i-1];

        for (int i = num_gpus - 1; i > 0; i--)
            load_balance_offset[i] = load_balance_offset[i-1];
        load_balance_offset[0] = 0;


        //invoke kernels
        for (int current_gpu = 0; current_gpu < num_gpus; current_gpu++)
        {

#ifdef IGNORE_FIRST_DEVICE
            cudaSetDevice(current_gpu + 1);
#else
            cudaSetDevice(current_gpu);
#endif

#ifdef GPU_INFO
            const int num_gpu = num_gpus;
            cudaDeviceProp p[num_gpu];
#ifdef IGNORE_FIRST_DEVICE
            err = cudaGetDeviceProperties(&p[current_gpu], current_gpu + 1);
#else
            err = cudaGetDeviceProperties(&p[current_gpu], current_gpu);
#endif

            checkError(err);
            /*
            #ifdef IGNORE_FIRST_DEVICE
                printf("\nGPU_INFO:: Device : %d \n", current_gpu+1);
            #else
                printf("\nGPU_INFO:: Device : %d \n", current_gpu);
            #endif
            printf("GPU_INFO:: Name : %s \n", p[current_gpu].name);
            printf("GPU_INFO:: totalGlobalMem : %u \n", p[current_gpu].totalGlobalMem);
            printf("GPU_INFO:: sharedMemPerBlock : %u \n", p[current_gpu].sharedMemPerBlock);
            printf("GPU_INFO:: regsPerBlock : %d \n", p[current_gpu].regsPerBlock);
            printf("GPU_INFO:: warpSize : %d \n", p[current_gpu].warpSize);
            printf("GPU_INFO:: maxThreadsPerBlock : %d \n", p[current_gpu].maxThreadsPerBlock);
            printf("GPU_INFO:: maxThreadsDim : %d  %d  %d \n", p[current_gpu].maxThreadsDim[0],p[current_gpu].maxThreadsDim[1],p[current_gpu].maxThreadsDim[2]);
            printf("GPU_INFO:: maxGridSize : %d  %d  %d\n", p[current_gpu].maxGridSize[0],p[current_gpu].maxGridSize[1],p[current_gpu].maxGridSize[2]);
            printf("GPU_INFO:: totalConstMem : %u \n", p[current_gpu].totalConstMem);
            printf("GPU_INFO:: deviceOverlap : %d \n", p[current_gpu].deviceOverlap);
            printf("GPU_INFO:: multiProcessorCount : %d \n", p[current_gpu].multiProcessorCount);
            printf("GPU_INFO:: kernelExecTimeoutEnabled : %d \n", p[current_gpu].kernelExecTimeoutEnabled);
            printf("GPU_INFO:: canMapHostMemory : %d \n", p[current_gpu].canMapHostMemory);
            printf("GPU_INFO:: computeMode : %d \n", p[current_gpu].computeMode);*/
#endif
            int j = 0;
            while (load_balance_length[current_gpu] > MAX_BLOCK_SIZE)
            {
                dim3 grid(MAX_BLOCK_SIZE, 1, 1);
#ifdef CUDA_LOG

#ifdef IGNORE_FIRST_DEVICE
                printf("CUDA_LOG::Calling %d CUDA blocks starting from offset %d, on GPU %d ", MAX_BLOCK_SIZE, j*MAX_BLOCK_SIZE + load_balance_offset[current_gpu], current_gpu + 1);
#else
                printf("CUDA_LOG::Calling %d CUDA blocks starting from offset %d, on GPU %d ", MAX_BLOCK_SIZE, j*MAX_BLOCK_SIZE + load_balance_offset[current_gpu], current_gpu);
#endif

#ifdef GPU_INFO
                printf("(%s)", p[current_gpu].name);
#endif
                printf("\n");
#endif
                kernel <<< grid, threads, SHARED_MEM_SIZE>>>(j*MAX_BLOCK_SIZE + load_balance_offset[current_gpu],
                        (float3 *)dpositions,
                        (float3 *)dvelocities,
                        (unsigned int *)dtarget_list,
                        (unsigned int *)dsource_list,
                        (unsigned int *)dnum_interaction_pairs,
                        (unsigned int *)dsource_start_indices_for_blocks,
                        NUM_THREADS, delta);
                j++;
                load_balance_length[current_gpu] -= MAX_BLOCK_SIZE;
            }
            dim3 grid(load_balance_length[current_gpu], 1, 1);
#ifdef CUDA_LOG

#ifdef IGNORE_FIRST_DEVICE
            printf("CUDA_LOG::Calling %d CUDA blocks starting from offset %d, on GPU %d ", load_balance_length[current_gpu], j*MAX_BLOCK_SIZE + load_balance_offset[current_gpu], current_gpu + 1);
#else
            printf("CUDA_LOG::Calling %d CUDA blocks starting from offset %d, on GPU %d ", load_balance_length[current_gpu], j*MAX_BLOCK_SIZE + load_balance_offset[current_gpu], current_gpu);
#endif

#ifdef GPU_INFO
            printf("(%s)", p[current_gpu].name);
#endif

            printf("\n");
#endif

            kernel <<< grid, threads, SHARED_MEM_SIZE>>>(j*MAX_BLOCK_SIZE + load_balance_offset[current_gpu],
                    (float3 *)dpositions,
                    (float3 *)dvelocities,
                    (unsigned int *)dtarget_list,
                    (unsigned int *)dsource_list,
                    (unsigned int *)dnum_interaction_pairs,
                    (unsigned int *)dsource_start_indices_for_blocks,
                    NUM_THREADS, delta);
        }

#ifdef ENABLE_BLOCKING
        cudaDeviceSynchronize();
        fprintf(stdout, "CUDA_LOG::Total GPU Time(Blocking call): %f\n,", wcTime() - startTime);
#endif
//         err = cudaDeviceSynchronize();
//         checkError(err);
        //assign globals
        GL_positions              = positions;
        GL_gpuVelocities          = gpuVelocities;
        GL_target_list            = target_list;
        GL_num_interaction_pairs  = num_interaction_pairs;
        GL_NUM_BODIES             = NUM_BODIES;

    }
}

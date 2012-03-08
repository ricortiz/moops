/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#define KERNEL
#include "kernel.hpp"
#undef KERNEL
#include <cutil.h>

__device__ __constant__ gpureal constDevc[514];                 // Constants on device

namespace                                                       // Prevent overlap of definitions among equations
{
__device__ void cart2sph(gpureal& r, gpureal& theta, gpureal& phi,// Get r,theta,phi from x,y,z on GPU
                         gpureal dx, gpureal dy, gpureal dz)
{
    r = sqrtf(dx * dx + dy * dy + dz * dz) + EPS;                 // r = sqrt(x^2 + y^2 + z^2) + eps
    theta = acosf(dz / r);                                        // theta = acos(z / r)
    if ( fabs(dx) + fabs(dy) < EPS )                              // If |x| < eps & |y| < eps
    {
        phi = 0;                                                    //  phi can be anything so we set it to 0
    }
    else if ( fabs(dx) < EPS )                                  // If |x| < eps
    {
        phi = dy / fabs(dy) * M_PI * 0.5;                           //  phi = sign(y) * pi / 2
    }
    else if ( dx > 0 )                                          // If x > 0
    {
        phi = atanf(dy / dx);                                       //  phi = atan(y / x)
    }
    else                                                        // If x < 0
    {
        phi = atanf(dy / dx) + M_PI;                                //  phi = atan(y / x) + pi
    }                                                             // End if for x,y cases
}

__device__ void sph2cart(gpureal r, gpureal theta, gpureal phi, // Spherical to cartesian coordinates on GPU
                         gpureal *spherical, gpureal *cartesian)
{
    cartesian[0] = sinf(theta) * cosf(phi) * spherical[0]         // x component (not x itself)
                   + cosf(theta) * cosf(phi) / r * spherical[1]
                   - sinf(phi) / r / sinf(theta) * spherical[2];
    cartesian[1] = sinf(theta) * sinf(phi) * spherical[0]         // y component (not y itself)
                   + cosf(theta) * sinf(phi) / r * spherical[1]
                   + cosf(phi) / r / sinf(theta) * spherical[2];
    cartesian[2] = cosf(theta) * spherical[0]                     // z component (not z itself)
                   - sinf(theta) / r * spherical[1];
}

__device__ void evalMultipole(gpureal *YnmShrd, gpureal rho,    // Evaluate solid harmonics r^n * Ynm on GPU
                              gpureal alpha, gpureal *factShrd)
{
    gpureal x = cosf(alpha);                                      // x = cos(alpha)
    gpureal y = sinf(alpha);                                      // y = sin(alpha)
    gpureal fact = 1;                                             // Initialize 2 * m + 1
    gpureal pn = 1;                                               // Initialize Legendre polynomial Pn
    gpureal rhom = 1;                                             // Initialize rho^m
    for ( int m = 0; m < P; ++m )                                 // Loop over m in Ynm
    {
        gpureal p = pn;                                             //  Associate Legendre polynomial Pnm
        int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
        int nmn = m * m;                                            //  Index of Ynm for m < 0
        YnmShrd[npn] = rhom * p / factShrd[2*m];                    //  rho^m * Ynm for m > 0
        YnmShrd[nmn] = YnmShrd[npn];                                //  Use conjugate relation for m < 0
        gpureal p1 = p;                                             //  Pnm-1
        p = x * (2 * m + 1) * p;                                    //  Pnm using recurrence relation
        rhom *= -rho;                                               //  rho^m
        gpureal rhon = rhom;                                        //  rho^n
        for ( int n = m + 1; n < P; ++n )                           //  Loop over n in Yn
        {
            int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
            int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
            YnmShrd[npm] = rhon * p / factShrd[n+m];                  //   rho^n * Ynm
            YnmShrd[nmm] = YnmShrd[npm];                              //   Use conjugate relation for m < 0
            gpureal p2 = p1;                                          //   Pnm-2
            p1 = p;                                                   //   Pnm-1
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
            rhon *= -rho;                                             //   Update rho^n
        }                                                           //  End loop over n in Ynm
        pn = -pn * fact * y;                                        //  Pn
        fact += 2;                                                  //  2 * m + 1
    }                                                             // End loop over m in Ynm
}

__device__ void evalLocal(gpureal *YnmShrd, gpureal rho,        // Evaluate singular harmonics r^(-n-1) * Ynm
                          gpureal alpha, gpureal *factShrd)
{
    gpureal x = cosf(alpha);                                      // x = cos(alpha)
    gpureal y = sinf(alpha);                                      // y = sin(alpha)
    gpureal rho_1 = 1 / rho;                                      // 1 / rho
    for ( int l = threadIdx.x; l < (2*P + 1)*P; l += THREADS )    // Loop over coefficients in Ynm
    {
        gpureal fact = 1;                                           //  Initialize 2 * m + 1
        gpureal pn = 1;                                             //  Initialize Legendre polynomial Pn
        gpureal rhom = rho_1;                                       //  Initialize rho^(-m-1)
        int nn = floor(sqrtf(2 * l + 0.25) - 0.5);                  //  Calculate index n of Ynm
        int mm = 0;                                                 //  Initialize index m of Ynm
        gpureal Ynm;                                                //  Define temporary Ynm
        for ( int i = 0; i <= nn; ++i ) mm += i;                    //  Offset of m
        mm = l - mm;                                                //  Calculate index m of Ynm
        int n;                                                      //  Define temporary n
        for ( int m = 0; m < mm; ++m )                              //  Loop up to m
        {
            rhom *= rho_1;                                            //   rho^(-m-1)
            pn = -pn * fact * y;                                      //   Pn
            fact += 2;                                                //   2 * m + 1
        }                                                           //  End loop up to m
        int m = mm;                                                 //  Define temporary m
        gpureal p = pn;                                             //  Associated Legendre polynomial Pnm
        if ( mm == nn ) Ynm = rhom * p * EPS;                       //  Ynm for n == m
        gpureal p1 = p;                                             //  Pnm-1
        p = x * (2 * m + 1) * p;                                    //  Pnm
        rhom *= rho_1;                                              //  rho^(-m-1)
        gpureal rhon = rhom;                                        //  rho^(-n-1)
        for ( n = m + 1; n < nn; ++n )                              //  Loop up to n
        {
            gpureal p2 = p1;                                          //   Pnm-2
            p1 = p;                                                   //   Pnm-1
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm
            rhon *= rho_1;                                            //   rho^(-n-1)
        }                                                           //  End loop up to n
        if ( n <= nn ) Ynm = rhon * p * factShrd[n-m];              //  rho^(-n-1) * Ynm
        YnmShrd[l] = Ynm;                                           //  Put Ynm in shared memory
    }                                                             // End loop over coefficients in Ynm
    __syncthreads();                                              // Syncronize threads
}
}                                                               // End anonymous namespace

void Kernel<Stokes>::initialize()
{
    startTimer("Init GPU     ");                                  // Start timer
    cudaThreadExit();                                             // Exit GPU thread
    cudaSetDevice(DEVICE);                                        // Set GPU device
    cudaThreadSynchronize();                                      // Sync GPU threads
#ifdef CUDA_4_1
    cudaSetDeviceFlags(cudaDeviceMapHost);
#endif
    stopTimer("Init GPU     ", MPIRANK == 0);                     // Stop timer & print
    eraseTimer("Init GPU     ");                                  // Erase timer
}

void Kernel<Stokes>::finalize()
{
}

void Kernel<Stokes>::allocate()
{
    cudaThreadSynchronize();
    startTimer("cudaMalloc   ");
#ifdef CUDA_4_1
    if ( keysHost.size() > keysDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&keysHost[0], sizeof(int) * keysHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &keysDevc, (void *)&keysHost[0], 0));
        keysDevcSize = keysHost.size();
    }
    if ( rangeHost.size() > rangeDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&rangeHost[0], sizeof(int) * rangeHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &rangeDevc, (void *)&rangeHost[0], 0));
        rangeDevcSize = rangeHost.size();
    }
    if ( sourceHost.size() > sourceDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&sourceDevc[0], sizeof(gpureal) * sourceHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &sourceDevc, (void *)&sourceDevc[0], 0));
        sourceDevcSize = sourceHost.size();
    }
    if ( targetHost.size() > targetDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&targetHost[0], sizeof(gpureal) * targetHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &targetDevc, (void *)&targetHost[0], 0));
        targetDevcSize = targetHost.size();
    }
#else
    if ( keysHost.size() > keysDevcSize )
    {
        if ( keysDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(keysDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &keysDevc, keysHost.size()*sizeof(int) ));
        keysDevcSize = keysHost.size();
    }
    if ( rangeHost.size() > rangeDevcSize )
    {
        if ( rangeDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(rangeDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &rangeDevc, rangeHost.size()*sizeof(int) ));
        rangeDevcSize = rangeHost.size();
    }
    if ( sourceHost.size() > sourceDevcSize )
    {
        if ( sourceDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(sourceDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(gpureal) ));
        sourceDevcSize = sourceHost.size();
    }
    if ( targetHost.size() > targetDevcSize )
    {
        if ( targetDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(targetDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(gpureal) ));
        targetDevcSize = targetHost.size();
    }
#endif
    cudaThreadSynchronize();
    stopTimer("cudaMalloc   ");
}


void Kernel<Stokes>::hostToDevice()
{
    cudaThreadSynchronize();
    startTimer("cudaMemcpy   ");
#ifndef CUDA_4_1
    CUDA_SAFE_CALL(cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sourceDevc, &sourceHost[0], sourceHost.size()*sizeof(gpureal), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(targetDevc, &targetHost[0], targetHost.size()*sizeof(gpureal), cudaMemcpyHostToDevice));
#endif
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(constDevc, &constHost[0], constHost.size()*sizeof(gpureal)));
    cudaThreadSynchronize();
    stopTimer("cudaMemcpy   ");
}

void Kernel<Stokes>::deviceToHost()
{
    cudaThreadSynchronize();
    startTimer("cudaMemcpy   ");
#ifdef CUDA_4_1
    CUDA_SAFE_CALL(cudaHostUnregister(&keysHost[0]));
    CUDA_SAFE_CALL(cudaHostUnregister(&rangeHost[0]));
    CUDA_SAFE_CALL(cudaHostUnregister(&sourceHost[0]));
    CUDA_SAFE_CALL(cudaHostUnregister(&targetHost[0]));
#else
    CUDA_SAFE_CALL(cudaMemcpy(&targetHost[0], targetDevc, targetHost.size()*sizeof(gpureal), cudaMemcpyDeviceToHost));
#endif
    cudaThreadSynchronize();
    stopTimer("cudaMemcpy   ");
}

__device__ void StokesP2M_core(gpureal *target, gpureal rho, gpureal alpha, gpureal beta, gpureal source)
{
    __shared__ gpureal factShrd[2*P];
    gpureal Ynm;
    gpureal fact = 1;
    for ( int i = 0; i < 2*P; ++i )
    {
        factShrd[i] = fact;
        fact *= i + 1;
    }
    __syncthreads();
    int nn = floorf(sqrtf(2 * threadIdx.x + 0.25) - 0.5);
    int mm = 0;
    for ( int i = 0; i <= nn; ++i ) mm += i;
    mm = threadIdx.x - mm;
    if ( threadIdx.x >= NTERM ) nn = mm = 0;
    gpureal x = cosf(alpha);
    gpureal s = sqrtf(1 - x * x);
    fact = 1;
    gpureal pn = 1;
    gpureal rhom = 1;
    for ( int m = 0; m < mm; ++m )
    {
        rhom *= rho;
        pn = -pn * fact * s;
        fact += 2;
    }
    int m = mm;
    gpureal p = pn;
    if (mm == nn) Ynm = rhom * p * rsqrtf(factShrd[2*m]);
    gpureal p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= rho;
    gpureal rhon = rhom;
    for ( int n = m + 1; n <= nn; ++n )
    {
        if (n == nn)
        {
            Ynm = rhon * p * rsqrtf(factShrd[n+m] / factShrd[n-m]);
        }
        gpureal p2 = p1;
        p1 = p;
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
        rhon *= rho;
    }
    gpureal ere = cosf(-mm * beta);
    gpureal eim = sinf(-mm * beta);
    target[0] += source * Ynm * ere;
    target[1] += source * Ynm * eim;
}

__global__ void StokesP2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal target[2] = {0, 0};
    __shared__ gpureal targetShrd[3];
    __shared__ gpureal sourceShrd[4*THREADS];
    int itarget = blockIdx.x * THREADS;
    targetShrd[0] = targetGlob[2*itarget+0];
    targetShrd[1] = targetGlob[2*itarget+1];
    targetShrd[2] = targetGlob[2*itarget+2];
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin = rangeGlob[keys+3*ilist+1];
        int size  = rangeGlob[keys+3*ilist+2];
        for ( int iblok = 0; iblok < (size - 1) / THREADS; ++iblok )
        {
            int isource = begin + iblok * THREADS + threadIdx.x;
            __syncthreads();
            sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
            sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
            sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
            sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
            __syncthreads();
            for ( int i = 0; i < THREADS; ++i )
            {
                float3 d;
                d.x = sourceShrd[4*i+0] - targetShrd[0];
                d.y = sourceShrd[4*i+1] - targetShrd[1];
                d.z = sourceShrd[4*i+2] - targetShrd[2];
                gpureal rho, alpha, beta;
                cart2sph(rho, alpha, beta, d.x, d.y, d.z);
                StokesP2M_core(target, rho, alpha, beta, sourceShrd[4*i+3]);
            }
        }
        int iblok = (size - 1) / THREADS;
        int isource = begin + iblok * THREADS + threadIdx.x;
        __syncthreads();
        if ( threadIdx.x < size - iblok * THREADS )
        {
            sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
            sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
            sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
            sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
        }
        __syncthreads();
        for ( int i = 0; i < size - iblok*THREADS; ++i )
        {
            float3 d;
            d.x = sourceShrd[4*i+0] - targetShrd[0];
            d.y = sourceShrd[4*i+1] - targetShrd[1];
            d.z = sourceShrd[4*i+2] - targetShrd[2];
            gpureal rho, alpha, beta;
            cart2sph(rho, alpha, beta, d.x, d.y, d.z);
            StokesP2M_core(target, rho, alpha, beta, sourceShrd[4*i+3]);
        }
    }
    itarget = blockIdx.x * THREADS + threadIdx.x;
    targetGlob[2*itarget+0] = target[0];
    targetGlob[2*itarget+1] = target[1];
}

void Kernel<Stokes>::P2M()
{
    cudaThreadSynchronize();
    startTimer("P2M GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
        StokesP2M_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc);

    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("P2M GPUkernel");
}

__device__ void StokesM2M_core(gpureal *target, gpureal beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd)
{
    int j = floorf(sqrtf(2 * threadIdx.x + 0.25) - 0.5);
    int k = 0;
    for ( int i = 0; i <= j; ++i ) k += i;
    k = threadIdx.x - k;
    if ( threadIdx.x >= NTERM ) j = k = 0;
    gpureal ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
    for ( int n = 0; n <= j; ++n )
    {
        for ( int m = -n; m <= min(k - 1, n); ++m )
        {
            if ( j - n >= k - m )
            {
                int nm = n * n + n + m;
                int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                gpureal ere = cosf(-m * beta);
                gpureal eim = sinf(-m * beta);
                gpureal ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
                gpureal cnm = ODDEVEN((m - abs(m)) / 2 + j);
                cnm *= ajnkm / ajk * YnmShrd[nm];
                gpureal CnmReal = cnm * ere;
                gpureal CnmImag = cnm * eim;
                target[0] += sourceShrd[2*jnkms+0] * CnmReal;
                target[0] -= sourceShrd[2*jnkms+1] * CnmImag;
                target[1] += sourceShrd[2*jnkms+0] * CnmImag;
                target[1] += sourceShrd[2*jnkms+1] * CnmReal;
            }
        }
        for ( int m = k; m <= n; ++m )
        {
            if ( j - n >= m - k )
            {
                int nm = n * n + n + m;
                int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                gpureal ere = cosf(-m * beta);
                gpureal eim = sinf(-m * beta);
                gpureal ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
                gpureal cnm = ODDEVEN(k + j + m);
                cnm *= ajnkm / ajk * YnmShrd[nm];
                gpureal CnmReal = cnm * ere;
                gpureal CnmImag = cnm * eim;
                target[0] += sourceShrd[2*jnkms+0] * CnmReal;
                target[0] += sourceShrd[2*jnkms+1] * CnmImag;
                target[1] += sourceShrd[2*jnkms+0] * CnmImag;
                target[1] -= sourceShrd[2*jnkms+1] * CnmReal;
            }
        }
    }
}

__global__ void StokesM2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal target[2] = {0, 0};
    __shared__ gpureal sourceShrd[2*THREADS];
    __shared__ gpureal factShrd[2*P];
    __shared__ gpureal YnmShrd[P*P];
    gpureal fact = 1;
    for ( int i = 0; i < 2*P; ++i )
    {
        factShrd[i] = fact;
        fact *= i + 1;
    }
    __syncthreads();
    int itarget = blockIdx.x * THREADS;
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin = rangeGlob[keys+3*ilist+1];
        float3 d;
        d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
        d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
        d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
        __syncthreads();
        if ( threadIdx.x < NTERM )
        {
            sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
            sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
        }
        __syncthreads();
        gpureal rho, alpha, beta;
        cart2sph(rho, alpha, beta, d.x, d.y, d.z);
        evalMultipole(YnmShrd, rho, alpha, factShrd);
        StokesM2M_core(target, beta, factShrd, YnmShrd, sourceShrd);
    }
    itarget = blockIdx.x * THREADS + threadIdx.x;
    targetGlob[2*itarget+0] = target[0];
    targetGlob[2*itarget+1] = target[1];
}

void Kernel<Stokes>::M2M()
{
    cudaThreadSynchronize();
    startTimer("M2M GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
        StokesM2M_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc);

    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("M2M GPUkernel");
}

__device__ void StokesM2L_core(gpureal *target, gpureal  beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd)
{
    int j = floorf(sqrtf(2 * threadIdx.x + 0.25) - 0.5);
    int k = 0;
    for ( int i = 0; i <= j; ++i ) k += i;
    k = threadIdx.x - k;
    if ( threadIdx.x >= NTERM ) j = k = 0;
    gpureal ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
    for ( int n = 0; n < P; ++n )
    {
        for ( int m = -n; m < 0; ++m )
        {
            int jnkm = (j + n) * (j + n + 1) / 2 - m + k;
            gpureal ere = cosf((m - k) * beta);
            gpureal eim = sinf((m - k) * beta);
            gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
            gpureal cnm = anm * ajk * YnmShrd[jnkm];
            gpureal CnmReal = cnm * ere;
            gpureal CnmImag = cnm * eim;
            int i = n * (n + 1) / 2 - m;
            target[0] += sourceShrd[2*i+0] * CnmReal;
            target[0] += sourceShrd[2*i+1] * CnmImag;
            target[1] += sourceShrd[2*i+0] * CnmImag;
            target[1] -= sourceShrd[2*i+1] * CnmReal;
        }
        for ( int m = 0; m <= n; ++m )
        {
            int jnkm = (j + n) * (j + n + 1) / 2 + abs(m - k);
            gpureal ere = cosf((m - k) * beta);
            gpureal eim = sinf((m - k) * beta);
            gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
            gpureal cnm = ODDEVEN((abs(k - m) - k - m) / 2);
            cnm *= anm * ajk * YnmShrd[jnkm];
            gpureal CnmReal = cnm * ere;
            gpureal CnmImag = cnm * eim;
            int i = n * (n + 1) / 2 + m;
            target[0] += sourceShrd[2*i+0] * CnmReal;
            target[0] -= sourceShrd[2*i+1] * CnmImag;
            target[1] += sourceShrd[2*i+0] * CnmImag;
            target[1] += sourceShrd[2*i+1] * CnmReal;
        }
    }
}

__global__ void StokesM2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal D0 = -constDevc[0];
    gpureal target[2] = {0, 0};
    __shared__ gpureal sourceShrd[2*THREADS];
    __shared__ gpureal factShrd[2*P];
    __shared__ gpureal YnmShrd[4*NTERM];
    gpureal fact = EPS;
    for ( int i = 0; i < 2*P; ++i )
    {
        factShrd[i] = fact;
        fact *= i + 1;
    }
    __syncthreads();
    int itarget = blockIdx.x * THREADS;
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin     = rangeGlob[keys+3*ilist+1];
        int Iperiodic = rangeGlob[keys+3*ilist+3];
        __syncthreads();
        if ( threadIdx.x < NTERM )
        {
            sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
            sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
        }
        __syncthreads();
        int I = 0;
        for ( int ix = -1; ix <= 1; ++ix )
        {
            for ( int iy = -1; iy <= 1; ++iy )
            {
                for ( int iz = -1; iz <= 1; ++iz, ++I )
                {
                    if ( Iperiodic & (1 << I) )
                    {
                        float3 d;
                        d.x = ix * D0;
                        d.y = iy * D0;
                        d.z = iz * D0;
                        d.x += targetGlob[2*itarget+0] - sourceGlob[begin+0];
                        d.y += targetGlob[2*itarget+1] - sourceGlob[begin+1];
                        d.z += targetGlob[2*itarget+2] - sourceGlob[begin+2];
                        gpureal rho, alpha, beta;
                        cart2sph(rho, alpha, beta, d.x, d.y, d.z);
                        evalLocal(YnmShrd, rho, alpha, factShrd);
                        StokesM2L_core(target, beta, factShrd, YnmShrd, sourceShrd);
                    }
                }
            }
        }
    }
    itarget = blockIdx.x * THREADS + threadIdx.x;
    targetGlob[2*itarget+0] = target[0] * EPS;
    targetGlob[2*itarget+1] = target[1] * EPS;
}

void Kernel<Stokes>::M2L()
{
    cudaThreadSynchronize();
    startTimer("M2L GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
        StokesM2L_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc);

    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("M2L GPUkernel");
}

__device__ void StokesM2P_core(gpureal *target, gpureal r, gpureal theta, gpureal phi, gpureal *factShrd, gpureal *sourceShrd)
{
    gpureal x = cosf(theta);
    gpureal y = sinf(theta);
    if ( fabsf(y) < EPS ) y = 1 / EPS;
    gpureal s = sqrtf(1 - x * x);
    gpureal spherical[3] = {0, 0, 0};
    gpureal cartesian[3] = {0, 0, 0};
    gpureal fact = 1;
    gpureal pn = 1;
    gpureal rhom = 1.0 / r;
    for ( int m = 0; m < P; ++m )
    {
        gpureal p = pn;
        int i = m * (m + 1) / 2 + m;
        gpureal ere = cosf(m * phi);
        if ( m == 0 ) ere = 0.5;
        gpureal eim = sinf(m * phi);
        gpureal anm = rhom * rsqrtf(factShrd[2*m]);
        gpureal Ynm = anm * p;
        gpureal p1 = p;
        p = x * (2 * m + 1) * p;
        gpureal YnmTheta = anm * (p - (m + 1) * x * p1) / y;
        gpureal realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
        gpureal imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
        target[0] += 2 * Ynm * realj;
        spherical[0] -= 2 * (m + 1) / r * Ynm * realj;
        spherical[1] += 2 * YnmTheta * realj;
        spherical[2] -= 2 * m * Ynm * imagj;
        rhom /= r;
        gpureal rhon = rhom;
        for ( int n = m + 1; n < P; ++n )
        {
            i = n * (n + 1) / 2 + m;
            anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
            Ynm = anm * p;
            gpureal p2 = p1;
            p1 = p;
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
            YnmTheta = anm * ((n - m + 1) * p - (n + 1) * x * p1) / y;
            realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
            imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
            target[0] += 2 * Ynm * realj;
            spherical[0] -= 2 * (n + 1) / r * Ynm * realj;
            spherical[1] += 2 * YnmTheta * realj;
            spherical[2] -= 2 * m * Ynm * imagj;
            rhon /= r;
        }
        pn = -pn * fact * s;
        fact += 2;
    }
    sph2cart(r, theta, phi, spherical, cartesian);
    target[1] += cartesian[0];
    target[2] += cartesian[1];
    target[3] += cartesian[2];
}

__global__ void StokesM2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal D0 = -constDevc[0];
    gpureal targetX[3];
    gpureal target[4] = {0, 0, 0, 0};
    __shared__ gpureal sourceShrd[2*THREADS];
    __shared__ gpureal factShrd[2*P];
    gpureal fact = 1;
    for ( int i = 0; i < 2*P; ++i )
    {
        factShrd[i] = fact;
        fact *= i + 1;
    }
    __syncthreads();
    int itarget = blockIdx.x * THREADS + threadIdx.x;
    targetX[0] = targetGlob[4*itarget+0];
    targetX[1] = targetGlob[4*itarget+1];
    targetX[2] = targetGlob[4*itarget+2];
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin     = rangeGlob[keys+3*ilist+1];
        int Iperiodic = rangeGlob[keys+3*ilist+3];
        __syncthreads();
        if ( threadIdx.x < NTERM )
        {
            sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
            sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
        }
        __syncthreads();
        int I = 0;
        for ( int ix = -1; ix <= 1; ++ix )
        {
            for ( int iy = -1; iy <= 1; ++iy )
            {
                for ( int iz = -1; iz <= 1; ++iz, ++I )
                {
                    if ( Iperiodic & (1 << I) )
                    {
                        float3 d;
                        d.x = ix * D0;
                        d.y = iy * D0;
                        d.z = iz * D0;
                        d.x += targetX[0] - sourceGlob[begin+0];
                        d.y += targetX[1] - sourceGlob[begin+1];
                        d.z += targetX[2] - sourceGlob[begin+2];
                        gpureal r, theta, phi;
                        cart2sph(r, theta, phi, d.x, d.y, d.z);
                        StokesM2P_core(target, r, theta, phi, factShrd, sourceShrd);
                    }
                }
            }
        }
    }
    targetGlob[4*itarget+0] = target[0];
    targetGlob[4*itarget+1] = target[1];
    targetGlob[4*itarget+2] = target[2];
    targetGlob[4*itarget+3] = target[3];
}

void Kernel<Stokes>::M2P()
{
    cudaThreadSynchronize();
    startTimer("M2P GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
    {
        StokesM2P_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc);
    }
    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("M2P GPUkernel");
}

__device__ void StokesP2P_core(gpureal *target, gpureal *targetX, gpureal *sourceShrd, float3 d, int i, float delta)
{
    d.x += targetX[0];
    d.x -= sourceShrd[6*i+0];
    d.y += targetX[1];
    d.y -= sourceShrd[6*i+1];
    d.z += targetX[2];
    d.z -= sourceShrd[6*i+2];

    float3 force = {sourceShrd[6*i+3], sourceShrd[6*i+4], sourceShrd[6*i+5]};

    float r2 = d.x * d.x + d.y * d.y + d.z * d.z;
    float d2 = delta * delta;
    float R1 = r2 + d2;
    float R2 = R1 + d2;
    float invR = 1.0 / R1;
    float H = sqrt(invR) * invR;

    float fdx =  force.x * d.x + force.y * d.y + force.z * d.z;

    target[0] += H * (force.x * R2 + fdx * d.x);
    target[1] += H * (force.y * R2 + fdx * d.y);
    target[2] += H * (force.z * R2 + fdx * d.z);

}

__global__ void StokesP2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob, float delta)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal D0 = -constDevc[0];
    gpureal targetX[3];
    gpureal target[4] = {0, 0, 0, 0};
    __shared__ gpureal sourceShrd[6*THREADS];
    int itarget = blockIdx.x * THREADS + threadIdx.x;
    targetX[0] = targetGlob[4*itarget+0];
    targetX[1] = targetGlob[4*itarget+1];
    targetX[2] = targetGlob[4*itarget+2];
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin     = rangeGlob[keys+3*ilist+1];
        int size      = rangeGlob[keys+3*ilist+2];
        int Iperiodic = rangeGlob[keys+3*ilist+3];
        for ( int iblok = 0; iblok < (size - 1) / THREADS; ++iblok )
        {
            int isource = begin + iblok * THREADS + threadIdx.x;
            __syncthreads();
            sourceShrd[6*threadIdx.x+0] = sourceGlob[6*isource+0];
            sourceShrd[6*threadIdx.x+1] = sourceGlob[6*isource+1];
            sourceShrd[6*threadIdx.x+2] = sourceGlob[6*isource+2];
            sourceShrd[6*threadIdx.x+3] = sourceGlob[6*isource+3];
            sourceShrd[6*threadIdx.x+4] = sourceGlob[6*isource+4];
            sourceShrd[6*threadIdx.x+5] = sourceGlob[6*isource+5];
            __syncthreads();
            int I = 0;
            for ( int ix = -1; ix <= 1; ++ix )
            {
                for ( int iy = -1; iy <= 1; ++iy )
                {
                    for ( int iz = -1; iz <= 1; ++iz, ++I )
                    {
                        if ( Iperiodic & (1 << I) )
                        {
                            float3 d;
                            d.x = ix * D0;
                            d.y = iy * D0;
                            d.z = iz * D0;
#pragma unroll 64
                            for ( int i = 0; i < THREADS; ++i )
                            {
                                StokesP2P_core(target, targetX, sourceShrd, d, i, delta);
                            }
                        }
                    }
                }
            }
        }
        int iblok = (size - 1) / THREADS;
        int isource = begin + iblok * THREADS + threadIdx.x;
        __syncthreads();
        if ( threadIdx.x < size - iblok * THREADS )
        {
            sourceShrd[6*threadIdx.x+0] = sourceGlob[6*isource+0];
            sourceShrd[6*threadIdx.x+1] = sourceGlob[6*isource+1];
            sourceShrd[6*threadIdx.x+2] = sourceGlob[6*isource+2];
            sourceShrd[6*threadIdx.x+3] = sourceGlob[6*isource+3];
            sourceShrd[6*threadIdx.x+4] = sourceGlob[6*isource+4];
            sourceShrd[6*threadIdx.x+5] = sourceGlob[6*isource+5];
        }
        __syncthreads();
        int I = 0;
        int icounter = 0;
        for ( int ix = -1; ix <= 1; ++ix )
        {
            for ( int iy = -1; iy <= 1; ++iy )
            {
                for ( int iz = -1; iz <= 1; ++iz, ++I )
                {
                    if ( Iperiodic & (1 << I) )
                    {
                        icounter++;
                        float3 d;
                        d.x = ix * D0;
                        d.y = iy * D0;
                        d.z = iz * D0;
                        for ( int i = 0; i < size - iblok*THREADS; ++i )
                        {
                            StokesP2P_core(target, targetX, sourceShrd, d, i, delta);
                        }
                    }
                }
            }
        }
    }
    targetGlob[4*itarget+0] = target[0];
    targetGlob[4*itarget+1] = target[1];
    targetGlob[4*itarget+2] = target[2];
    targetGlob[4*itarget+3] = target[3];
}

void Kernel<Stokes>::P2P()
{
    cudaThreadSynchronize();
    startTimer("P2P GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
    {
        StokesP2P_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc, delta);
    }
    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("P2P GPUkernel");
}

__device__ void StokesL2L_core(gpureal *target, gpureal beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd)
{
    int j = floorf(sqrtf(2 * threadIdx.x + 0.25) - 0.5);
    int k = 0;
    for ( int i = 0; i <= j; ++i ) k += i;
    k = threadIdx.x - k;
    if ( threadIdx.x >= NTERM ) j = k = 0;
    gpureal ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
    for ( int n = 0; n < P; ++n )
    {
        for ( int m = j + k - n; m < 0; ++m )
        {
            int nms = n * (n + 1) / 2 - m;
            int jnkm = (n - j) * (n - j) + n - j + m - k;
            gpureal ere = cosf((m - k) * beta);
            gpureal eim = sinf((m - k) * beta);
            gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
            gpureal cnm = ODDEVEN(k - n) * ajk / anm * YnmShrd[jnkm];
            gpureal CnmReal = cnm * ere;
            gpureal CnmImag = cnm * eim;
            target[0] += sourceShrd[2*nms+0] * CnmReal;
            target[0] += sourceShrd[2*nms+1] * CnmImag;
            target[1] += sourceShrd[2*nms+0] * CnmImag;
            target[1] -= sourceShrd[2*nms+1] * CnmReal;
        }
        for ( int m = 0; m <= n; ++m )
        {
            if ( n - j >= abs(m - k) )
            {
                int nms = n * (n + 1) / 2 + m;
                int jnkm = (n - j) * (n - j) + n - j + m - k;
                gpureal ere = cosf((m - k) * beta);
                gpureal eim = sinf((m - k) * beta);
                gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
                gpureal cnm = ODDEVEN((m - k - abs(m - k)) / 2 - n);
                cnm *= ajk / anm * YnmShrd[jnkm];
                gpureal CnmReal = cnm * ere;
                gpureal CnmImag = cnm * eim;
                target[0] += sourceShrd[2*nms+0] * CnmReal;
                target[0] -= sourceShrd[2*nms+1] * CnmImag;
                target[1] += sourceShrd[2*nms+0] * CnmImag;
                target[1] += sourceShrd[2*nms+1] * CnmReal;
            }
        }
    }
}

__global__ void StokesL2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal target[2] = {0, 0};
    __shared__ gpureal sourceShrd[2*THREADS];
    __shared__ gpureal factShrd[2*P];
    __shared__ gpureal YnmShrd[P*P];
    gpureal fact = 1;
    for ( int i = 0; i < 2*P; ++i )
    {
        factShrd[i] = fact;
        fact *= i + 1;
    }
    __syncthreads();
    int itarget = blockIdx.x * THREADS;
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin = rangeGlob[keys+3*ilist+1];
        float3 d;
        d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
        d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
        d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
        __syncthreads();
        if ( threadIdx.x < NTERM )
        {
            sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
            sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
        }
        __syncthreads();
        gpureal rho, alpha, beta;
        cart2sph(rho, alpha, beta, d.x, d.y, d.z);
        evalMultipole(YnmShrd, rho, alpha, factShrd);
        StokesL2L_core(target, beta, factShrd, YnmShrd, sourceShrd);
    }
    itarget = blockIdx.x * THREADS + threadIdx.x;
    targetGlob[2*itarget+0] = target[0];
    targetGlob[2*itarget+1] = target[1];
}

void Kernel<Stokes>::L2L()
{
    cudaThreadSynchronize();
    startTimer("L2L GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
        StokesL2L_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc);
    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("L2L GPUkernel");
}

__device__ void StokesL2P_core(gpureal *target, gpureal r, gpureal theta, gpureal phi, gpureal *factShrd, gpureal *sourceShrd)
{
    gpureal x = cosf(theta);
    gpureal y = sinf(theta);
    if ( fabsf(y) < EPS ) y = 1 / EPS;
    gpureal s = sqrtf(1 - x * x);
    gpureal spherical[3] = {0, 0, 0};
    gpureal cartesian[3] = {0, 0, 0};
    gpureal fact = 1;
    gpureal pn = 1;
    gpureal rhom = 1;
    for ( int m = 0; m < P; ++m )
    {
        gpureal p = pn;
        int i = m * (m + 1) / 2 + m;
        gpureal ere = cosf(m * phi);
        if ( m == 0 ) ere = 0.5;
        gpureal eim = sinf(m * phi);
        gpureal anm = rhom * rsqrtf(factShrd[2*m]);
        gpureal Ynm = anm * p;
        gpureal p1 = p;
        p = x * (2 * m + 1) * p;
        gpureal YnmTheta = anm * (p - (m + 1) * x * p1) / y;
        gpureal realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
        gpureal imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
        target[0] += 2 * Ynm * realj;
        spherical[0] += 2 * m / r * Ynm * realj;
        spherical[1] += 2 * YnmTheta * realj;
        spherical[2] -= 2 * m * Ynm * imagj;
        rhom *= r;
        gpureal rhon = rhom;
        for ( int n = m + 1; n < P; ++n )
        {
            i = n * (n + 1) / 2 + m;
            anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
            Ynm = anm * p;
            gpureal p2 = p1;
            p1 = p;
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
            YnmTheta = anm * ((n - m + 1) * p - (n + 1) * x * p1) / y;
            realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
            imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
            target[0] += 2 * Ynm * realj;
            spherical[0] += 2 * n / r * Ynm * realj;
            spherical[1] += 2 * YnmTheta * realj;
            spherical[2] -= 2 * m * Ynm * imagj;
            rhon *= r;
        }
        pn = -pn * fact * s;
        fact += 2;
    }
    sph2cart(r, theta, phi, spherical, cartesian);
    target[1] += cartesian[0];
    target[2] += cartesian[1];
    target[3] += cartesian[2];
}

__global__ void StokesL2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal targetX[3];
    gpureal target[4] = {0, 0, 0, 0};
    __shared__ gpureal sourceShrd[2*THREADS];
    __shared__ gpureal factShrd[2*P];
    gpureal fact = 1;
    for ( int i = 0; i < 2*P; ++i )
    {
        factShrd[i] = fact;
        fact *= i + 1;
    }
    __syncthreads();
    int itarget = blockIdx.x * THREADS + threadIdx.x;
    targetX[0] = targetGlob[4*itarget+0];
    targetX[1] = targetGlob[4*itarget+1];
    targetX[2] = targetGlob[4*itarget+2];
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin = rangeGlob[keys+3*ilist+1];
        float3 d;
        d.x = targetX[0] - sourceGlob[begin+0];
        d.y = targetX[1] - sourceGlob[begin+1];
        d.z = targetX[2] - sourceGlob[begin+2];
        __syncthreads();
        if ( threadIdx.x < NTERM )
        {
            sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
            sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
        }
        __syncthreads();
        gpureal r, theta, phi;
        cart2sph(r, theta, phi, d.x, d.y, d.z);
        StokesL2P_core(target, r, theta, phi, factShrd, sourceShrd);
    }
    targetGlob[4*itarget+0] = target[0];
    targetGlob[4*itarget+1] = target[1];
    targetGlob[4*itarget+2] = target[2];
    targetGlob[4*itarget+3] = target[3];
}

void Kernel<Stokes>::L2P()
{
    cudaThreadSynchronize();
    startTimer("L2P GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
        StokesL2P_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc);
    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("L2P GPUkernel");
}


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

void Kernel<Stokes>::P2P(C_iter Ci, C_iter Cj) const           // Stokes P2P kernel on CPU
{
    for ( B_iter Bi = Ci->LEAF; Bi != Ci->LEAF + Ci->NDLEAF; ++Bi )      // Loop over target bodies
    {
        for ( B_iter Bj = Cj->LEAF; Bj != Cj->LEAF + Cj->NDLEAF; ++Bj )    //  Loop over source bodies
        {
            vect dx = Bi->X - Bj->X;
            real r2 = norm(dx);
            real d2 = delta * delta;
            real R1 = r2 + d2;
            real R2 = R1 + d2;
            real invR = 1.0 / R1;
            real H = std::sqrt(invR) * invR;
            
            real fdx = (Bj->FORCE[0] * dx[0] + Bj->FORCE[1] * dx[1] + Bj->FORCE[2] * dx[2]);
            
            Bi->TRG[0] += H * (Bj->FORCE[0] * R2 + fdx * dx[0]);
            Bi->TRG[1] += H * (Bj->FORCE[1] * R2 + fdx * dx[1]);
            Bi->TRG[2] += H * (Bj->FORCE[2] * R2 + fdx * dx[2]);
        }                                                           //  End loop over source bodies
    }                                                             // End loop over target bodies
}

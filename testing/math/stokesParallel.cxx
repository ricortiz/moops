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
#include "math/fluid_solver/stokes/fmm/hybrid/include/parallelfmm.hpp"

int stokesParallel(int, char**)
{
    int numBodies = 10000;                                        // Number of bodies
    int numTarget = 100;                                          // Number of target points to be used for error eval.
    IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic FMM)
    THETA = 1 / sqrtf(4);                                         // Multipole acceptace criteria
    Bodies bodies, jbodies;                                       // Define vector of bodies
    Cells cells, jcells;                                          // Define vector of cells
    ParallelFMM<Stokes> FMM;                                     // Instantiate ParallelFMM class
    FMM.initialize();                                             // Initialize FMM
    FMM.setDelta(.01);
    bool printNow = MPIRANK == 0;                                 // Print only if MPIRANK == 0

    for ( int it=0; it!=25; ++it ) {                              // Loop from N = 10^4 to 10^7
        numBodies = int(pow(10,(it+32)/8.0));                       //  Exponentially increase N
        if (printNow) std::cout << "N             : " << numBodies << std::endl;// Print N
        bodies.resize(numBodies);                                   //  Resize bodies vector
        FMM.random(bodies,MPIRANK+1);                               //  Initialize bodies with random coordinates
        FMM.startTimer("FMM          ");                            //  Start timer
        FMM.setGlobDomain(bodies);                                  //  Set global domain size of FMM
        FMM.octsection(bodies);                                     //  Partition domain and redistribute bodies
        cells.clear();                                              //  Make sure cells vector is empty
#ifdef TOPDOWN
        FMM.topdown(bodies,cells);                                  //  Tree construction (top down) & upward sweep
#else
        FMM.bottomup(bodies,cells);                                 //  Tree construction (bottom up) & upward sweep
#endif
        FMM.commBodies(cells);                                      //  Send bodies (not receiving yet)
        jbodies = bodies;                                           //  Vector of source bodies
        jcells = cells;                                             //  Vector of source cells
        FMM.commCells(jbodies,jcells);                              //  Communicate cells (receive bodies here)

        FMM.downward(cells,jcells);                                 //  Downward sweep
        FMM.stopTimer("FMM          ",printNow);                    //  Stop timer
        FMM.eraseTimer("FMM          ");                            //  Erase entry from timer to avoid timer overlap

        FMM.startTimer("Direct sum   ");                            //  Start timer
        Bodies bodies2 = bodies;                                    //  Define new bodies vector for direct sum
#if 1
        FMM.initTarget(bodies2);                                    //  Reset target values to 0
        if ( IMAGES != 0 ) {                                        //  For periodic boundary condition
            jbodies = FMM.periodicBodies(bodies2);                    //   Copy source bodies for all periodic images
        } else {                                                    //  For free field boundary condition
            jbodies = bodies2;                                        //   Source bodies = target bodies
        }                                                           //  End if for periodic boundary condition
        bodies2.resize(numTarget);                                  //  Shrink target bodies vector to save time
        for ( int i=0; i!=MPISIZE; ++i ) {                          //  Loop over all MPI processes
            FMM.shiftBodies(jbodies);                                 //   Communicate bodies round-robin
            FMM.evalP2P(bodies2,jbodies,1);                             //   Direct summation between bodies2 and jbodies
            if (FMM.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;// Print loop counter
        }                                                           //  End loop over all MPI processes
        FMM.writeTarget(bodies2);                                   //  Write direct summation results to file
#else
        FMM.readTarget(bodies2);                                    //  Read direct summation results from file
#endif
        FMM.stopTimer("Direct sum   ",printNow);                    //  Stop timer
        FMM.eraseTimer("Direct sum   ");                            //  Erase entry from timer to avoid timer overlap
        if (printNow) FMM.writeTime();                              //  Write timings of all events to file
        if (printNow) FMM.resetTimer();                             //  Erase all events in timer

        real diff = 0, norm = 0, diff2 = 0, norm2 = 0;
        bodies.resize(numTarget);                                   //  Shrink bodies to match bodies2
        FMM.evalError(bodies,bodies2,diff,norm);      //  Evaluate error on the reduced set of bodies
        MPI_Datatype MPI_TYPE = FMM.getType(diff);                 //  Get MPI datatype
        MPI_Reduce(&diff,&diff2,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);// Reduce difference in potential
        MPI_Reduce(&norm,&norm2,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);// Reduce norm of potential
        if (printNow) FMM.printError(diff2,norm2);      //  Print the L2 norm error of potential & force
    }                                                             // End loop over N
    FMM.finalize();                                               //
    std::cout << "v1 = [";
    for (B_iter B = bodies.begin(), end = bodies.end(); B != end; ++B)
        std::cout << B->TRG[0] << "," << B->TRG[1] << "," << B->TRG[2] << ";";
    std::cout << "];" << std::endl;
    std::cout << "v2 = [";
    for (B_iter B = FMM.buffer.begin(), end = FMM.buffer.end(); B != end; ++B)
        std::cout << B->TRG[0] << "," << B->TRG[1] << "," << B->TRG[2] << ";";
    std::cout << "];" << std::endl;
    return 0;
}

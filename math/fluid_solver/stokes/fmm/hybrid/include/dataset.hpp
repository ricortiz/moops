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
#ifndef STOKES_DATASET_HPP
#define STOKES_DATASET_HPP

#include "external/exafmm/include/dataset.h"

template<>
class Dataset<Stokes> : public Kernel<Stokes>
{
    private:
        long filePosition;                                            //!< Position of file stream

    public:
        //! Constructor
        Dataset() : filePosition(0) {}
        //! Destructor
        ~Dataset() {}

        //! Initialize source values
        void initSource(Bodies &bodies)
        {
            for ( B_iter B = bodies.begin(); B != bodies.end(); ++B )   // Loop over bodies
            {
                B->IBODY = B - bodies.begin();                            //  Tag body with initial index
                B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
                B->FORCE[0] = drand48();                    //  Initialize force
                B->FORCE[1] = drand48();                    //  Initialize force
                B->FORCE[2] = drand48();                    //  Initialize force
            }                                                           // End loop over bodies
        }

        //! Initialize target values
        void initTarget(Bodies &bodies)
        {
            srand48(0);                                                 // Set seed for random number generator
            for ( B_iter B = bodies.begin(); B != bodies.end(); ++B )   // Loop over bodies
            {
                B->IBODY = B - bodies.begin();                            //  Tag body with initial index
                B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
                B->TRG = 0;                                               //  Clear previous target values 
            }                                                           // End loop over bodies
        }

        //! Read target values from file
        void readTarget(Bodies &bodies)
        {
            char fname[256];                                            // File name for saving direct calculation values
            sprintf(fname, "direct%4.4d", MPIRANK);                     // Set file name
            std::ifstream file(fname, std::ios::in | std::ios::binary); // Open file
            file.seekg(filePosition);                                   // Set position in file
            for ( B_iter B = bodies.begin(); B != bodies.end(); ++B )   // Loop over bodies
            {
                file >> B->TRG[0];                                        //  Read data for x velocity
                file >> B->TRG[1];                                        //  Read data for y velocity
                file >> B->TRG[2];                                        //  Read data for z velocity
            }                                                           // End loop over bodies
            filePosition = file.tellg();                                // Get position in file
            file.close();                                               // Close file
        }

        //! Write target values to file
        void writeTarget(Bodies &bodies)
        {
            char fname[256];                                            // File name for saving direct calculation values
            sprintf(fname, "direct%4.4d", MPIRANK);                     // Set file name
            std::ofstream file(fname, std::ios::out | std::ios::app | std::ios::binary);// Open file
            for ( B_iter B = bodies.begin(); B != bodies.end(); ++B )   // Loop over bodies
            {
                file << B->TRG[0];                                        //  Write data for x velocity
                file << B->TRG[1];                                        //  Write data for y velocity
                file << B->TRG[2];                                        //  Write data for z velocity
            }                                                           // End loop over bodies
            file.close();                                               // Close file
        }

        //! Evaluate relative L2 norm error
        //! Evaluate relative L2 norm error
        void evalError(Bodies &bodies, Bodies &bodies2, real &diff, real &norm)
        {
            B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
            for ( B_iter B = bodies.begin(); B != bodies.end(); ++B, ++B2 )
            {                                                           // Loop over bodies & bodies2
#ifdef DEBUG
                std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
                diff += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of x velocity
                diff += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of y velocity
                diff += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of z velocity
                norm += B2->TRG[0] * B2->TRG[0];                         //  Value of x velocity
                norm += B2->TRG[1] * B2->TRG[1];                         //  Value of y velocity
                norm += B2->TRG[2] * B2->TRG[2];                         //  Value of z velocity
            }                                                           //  End loop over bodies & bodies2
        }

        //! Print relative L2 norm error
        void printError(real diff, real norm)
        {
            std::cout << "Error (vel)   : " << std::sqrt(diff / norm) << std::endl;
        }
};

#endif

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
#ifndef STOKES_KERNEL_HPP
#define STOKES_KERNEL_HPP

#include "external/exafmm/include/types.h"

template<Equation e>
class Kernel;

#include "external/exafmm/include/kernel.h"

template<>
class Kernel<Stokes> : public KernelBase {
public:
  void initialize();                                            //!< Initialize kernels
  void P2M(C_iter Ci) const;                                    //!< Evaluate P2M kernel on CPU
  void M2M(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2M kernel on CPU
  void M2L(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2L kernel on CPU
  void M2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2P kernel on CPU
  void P2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate P2P kernel on CPU
  void L2L(C_iter Ci, C_iter Cj) const;                         //!< Evaluate L2L kernel on CPU
  void L2P(C_iter Ci) const;                                    //!< Evaluate L2P kernel on CPU
  void D2M(C_iter Ci) const;                                    //!< Evaluate D2P kernel on CPU
  void P2M();                                                   //!< Evaluate P2M kernel on GPU
  void M2M();                                                   //!< Evaluate M2M kernel on GPU
  void M2L();                                                   //!< Evaluate M2L kernel on GPU
  void M2P();                                                   //!< Evaluate M2P kernel on GPU
  void P2P();                                                   //!< Evalaute P2P kernel on GPU
  void L2L();                                                   //!< Evaluate L2L kernel on GPU
  void L2P();                                                   //!< Evaluate L2P kernel on GPU
  void finalize();                                              //!< Finalize kernels

  void allocate();                                              //!< Allocate GPU variables
  void hostToDevice();                                          //!< Copy from host to device
  void deviceToHost();
  void setDelta(real _delta) { delta = _delta; }
  real delta;
};

#endif

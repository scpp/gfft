/***************************************************************************
 *   Copyright (C) 2006-2015 by Vladimir Mirnyy                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 ***************************************************************************/

#ifndef __gfftstdalg_h
#define __gfftstdalg_h

/** \file
    \brief Recursive FFT algorithms on "complex" data types like std::complex
    
*/

#include "gfftstdspec.h"
#include "gfftfactor.h"
#include "gfftswap.h"

#include "metacomplex.h"
#include "metaroot.h"

namespace GFFT {

using namespace MF;


/// In-place scaled FFT algorithm
/**
\tparam K first factor
\tparam M second factor (N=K*M) 
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam W compile-time root of unity
\tparam doStaticLoop rely on template instantiation loop IterateInTime (only for short loops)

The notation for this template class follows SPIRAL. 
The class performs DFT(k) with the Kronecker product by the mxm identity matrix Im
and twiddle factors (T).
\sa InTime
*/

template<long_t K, long_t LastK, long_t M, long_t Step, typename VType, int S, class W1, long_t SimpleSpec>
class DFTk_x_Im_T<K,LastK,M,Step,VType,S,W1,SimpleSpec,false>
{
   typedef typename VType::ValueType CT;
   static const long_t N = K*M;
   DFTk_inp<K,M,VType,S> spec_inp;
public:
   void apply(CT* data) 
   {
      spec_inp.apply(data);

      ComputeRootsStd<K,VType,W1> roots;
      
      spec_inp.apply(data+Step, roots.get());
      for (long_t i=Step+Step; i<M; i+=Step) {
        roots.step();
        spec_inp.apply(data+i, roots.get());
      }
   }
  
};


// Specialization for radix 2
template<long_t M, long_t LastK, long_t Step, typename VType, int S, class W1, long_t SimpleSpec>
class DFTk_x_Im_T<2,LastK,M,Step,VType,S,W1,SimpleSpec,false>
{
   typedef typename VType::ValueType CT;
//   typedef typename VType::ValueType T;
   //static const long_t N = 2*M;
   //static const long_t S2 = 2*Step;
   //static const long_t NR = (4*PrecomputeRoots > LastK*M) ? LastK*M/4 : PrecomputeRoots;
//   static const long_t NR = (PrecomputeRoots>=N) ? 2 : PrecomputeRoots/2;
   //static const long_t K = N/PrecomputeRoots;
//   static const long_t K2 = 2*K;
   //typedef typename GetFirstRoot<N,S,VType::Accuracy>::Result W1;
   typedef Compute<typename W1::Re,VType::Accuracy> WR;
   typedef Compute<typename W1::Im,VType::Accuracy> WI;
   DFTk_inp<2,M,VType,S> spec_inp;
public:
   void apply(CT* data)
   {
      spec_inp.apply(data);
      if (M%2 == 0)
        spec_inp.apply_1(data + M/2);

      const CT wp(WR::value(), WI::value());
      CT t, w(wp);

      spec_inp.apply(data+Step, &w);
      t = CT(-w.real(), w.imag());
      spec_inp.apply(data+M-Step, &t);
      for (long_t i=Step+Step; i<M/2; i+=Step) {
          w *= wp;
          spec_inp.apply(data+i, &w);
          t = CT(-w.real(), w.imag());
          spec_inp.apply(data+M-i, &t);
      }
   }
};

// Specialization for radix 2
template<long_t LastK, long_t M, long_t Step, typename VType, int S, class W1>
class DFTk_x_Im_T<2,LastK,M,Step,VType,S,W1,2,false>
{
   typedef typename VType::ValueType CT;
   DFTk_inp<2,M,VType,S> spec_inp;
public:
   void apply(CT* data)
   {
      spec_inp.apply(data);
      spec_inp.apply_1(data+Step);
   }
};

}  //namespace DFT

#endif /*__gfftstdalg_h*/

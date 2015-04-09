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

template<int_t K, int_t LastK, int_t M, int_t Step, typename VType, int S, class W1, int_t SimpleSpec>
class DFTk_x_Im_T<K,LastK,M,Step,VType,S,W1,SimpleSpec,false>
{
   typedef typename VType::ValueType CT;
   static const int_t N = K*M;
   DFTk_inp<K,M,VType,S> spec_inp;
public:
   void apply(CT* data) 
   {
      spec_inp.apply(data);

      ComputeRootsStd<K,VType,W1> roots;
      
      spec_inp.apply(data+Step, roots.get());
      for (int_t i=Step+Step; i<M; i+=Step) {
	roots.step();
	spec_inp.apply(data+i, roots.get());
      }
   }
  
};

  
}  //namespace DFT

#endif /*__gfftstdalg_h*/

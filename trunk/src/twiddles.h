/***************************************************************************
 *   Copyright (C) 2006-2014 by Vladimir Mirnyy                            *
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

#ifndef __twiddles_h
#define __twiddles_h

/** \file
    \brief Short-radix in-place FFT specifications 
*/

#include "metafunc.h"

namespace GFFT {

using namespace MF;

/// Twiddle factors computation
/*!
\tparam T value type
\tparam N length of the data
\tparam S sign of the transform (-1 for inverse)

Computes twiddle factors, 
e.g. cos(2*pi/N),...,cos(2*K*pi/N) and sin(2*pi/N),...,sin(2*K*pi/N).
It is used for prime factors greater than 3.
*/
template<typename T, int_t N, int S, int_t K>
struct ComputeTwiddles
{
  static const int_t K2 = K*2;
  
  static void apply(T* c, T* s)
  {
    ComputeTwiddles<T,N,S,K-1>::apply(c, s);
    c[K-1] = Cos<N,K2,T>::value();
    s[K-1] = S*Sin<N,K2,T>::value();
  }
};

template<typename T, int_t N, int S>
struct ComputeTwiddles<T,N,S,1>
{
  static void apply(T* c, T* s)
  {
    c[0] = Cos<N,2,T>::value();
    s[0] = S*Sin<N,2,T>::value();
  }
};

}  //namespace DFT

#endif /*__twiddles_h*/

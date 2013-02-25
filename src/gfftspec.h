/***************************************************************************
 *   Copyright (C) 2006-2013 by Volodymyr Myrnyy                           *
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

#ifndef __gfftspec_h
#define __gfftspec_h

/** \file
    \brief Recursive algorithms and short-radix FFT specifications
*/

#include "metafunc.h"


namespace GFFT {

using namespace MF;

/// Specialization for complex-valued radix 2 FFT in-place
/// \tparam T is value type
/// \param data is the array of length 4, containing two complex numbers (real,imag,real,imag).
template<typename T>
inline void _spec2(T* data) 
{
      const T tr = data[2];
      const T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
}

template<typename T>
inline void _spec2(const T* src, T* dst) 
{ 
    const T *a(src + 2), *a2(src + 1), *a3(src + 3);
    *(dst) = (*(src) + *(a));
    *((dst + 1)) = (*(a2) + *(a3));
    *((dst + 2)) = (*(src) - *(a));
    *((dst + 3)) = (*(a2) - *(a3));
// test against above version (from spiral)
//     const T v1(src[1]), v2(src[2]), v3(src[3]);
//     dst[0] = (*src + v2);
//     dst[1] = (v1 + v3);
//     dst[2] = (*src - v2);
//     dst[3] = (v1 - v3);
}

template<typename T>
inline void _spec3_fwd(T* data) 
{ // 5 mult, 12 add
  const T c = Sqrt<3, T>::value() * 0.5;
  const T sr = data[2] + data[4];
  const T dr = c*(data[2] - data[4]);
  const T si = data[3] + data[5];
  const T di = c*(data[3] - data[5]);
  const T tr = data[0] - 0.5*sr;
  const T ti = data[1] - 0.5*si;
  data[0] += sr;
  data[1] += si;
  data[4] = tr + di;
  data[5] = ti - dr;
  data[2] = tr - di;
  data[3] = ti + dr;
}

template<typename T>
inline void _spec3_fwd(const T* src, T* dst) 
{ // 5 mult, 12 add
  const T c = Sqrt<3, T>::value() * 0.5;
  const T sr = src[2] + src[4];
  const T dr = c*(src[2] - src[4]);
  const T si = src[3] + src[5];
  const T di = c*(src[3] - src[5]);
  const T tr = src[0] - 0.5*sr;
  const T ti = src[1] - 0.5*si;
  dst[0] = src[0] + sr;
  dst[1] = src[1] + si;
  dst[2] = tr + di;
  dst[3] = ti - dr;
  dst[4] = tr - di;
  dst[5] = ti + dr;
}

template<typename T>
inline void _spec4_fwd(const T* src, T* dst) 
{
  const double  *a, *a2, *a3, *a4, *a5, *a6, *a7;
    double t, t2, t3, t4, t5, t6, t7, t8;
    a = (src + 4);
    t = (*(src) + *(a));
    a2 = (src + 1);
    a3 = (src + 5);
    t2 = (*(a2) + *(a3));
    t3 = (*(src) - *(a));
    t4 = (*(a2) - *(a3));
    a4 = (src + 2);
    a5 = (src + 6);
    t5 = (*(a4) + *(a5));
    a6 = (src + 3);
    a7 = (src + 7);
    t6 = (*(a6) + *(a7));
    t7 = (*(a4) - *(a5));
    t8 = (*(a6) - *(a7));
    *(dst) = (t + t5);
    *((dst + 1)) = (t2 + t6);
    *((dst + 4)) = (t - t5);
    *((dst + 5)) = (t2 - t6);
    *((dst + 2)) = (t3 + t8);
    *((dst + 3)) = (t4 - t7);
    *((dst + 6)) = (t3 - t8);
    *((dst + 7)) = (t4 + t7);
}



}  //namespace DFT

#endif /*__gfftspec_h*/

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

typedef unsigned long int_t;

namespace GFFT {

using namespace MF;

template<int_t N, int_t M, typename T, int S>
class DFTk_inplace;

template<int_t M, typename T, int S>
class DFTk_inplace<3,M,T,S> 
{
  T m_coef;
  
public:
  DFTk_inplace() : m_coef(S * Sqrt<3, T>::value() * 0.5) { }
  
  void apply(T* data) 
  { 
      const int_t i10 = M;
      const int_t i11 = i10+1;
      const int_t i20 = i10+M;
      const int_t i21 = i20+1;

      const T sr = data[i10] + data[i20];
      const T dr = m_coef * (data[i10] - data[i20]);
      const T si = data[i11] + data[i21];
      const T di = m_coef * (data[i11] - data[i21]);
      const T tr = data[0] - 0.5*sr;
      const T ti = data[1] - 0.5*si;
      data[0] += sr;
      data[1] += si;
      data[i10] = tr + di;
      data[i11] = ti - dr;
      data[i20] = tr - di;
      data[i21] = ti + dr;
  }
  template<class LT>
  void apply(T* data, const LT* wr, const LT* wi) 
  { 
	const int_t i10 = M;
	const int_t i11 = i10+1;
	const int_t i20 = i10+M;
	const int_t i21 = i20+1;
        const T tr1 = data[i10]*wr[0] - data[i11]*wi[0];
        const T ti1 = data[i10]*wi[0] + data[i11]*wr[0];
        const T tr2 = data[i20]*wr[1] - data[i21]*wi[1];
        const T ti2 = data[i20]*wi[1] + data[i21]*wr[1];

	const T sr = tr1 + tr2;
	const T dr = m_coef * (tr1 - tr2);
	const T si = ti1 + ti2;
	const T di = m_coef * (ti1 - ti2);
	const T tr = data[0] - 0.5*sr;
	const T ti = data[1] - 0.5*si;
	data[0] += sr;
	data[1] += si;
	data[i10] = tr + di;
	data[i11] = ti - dr;
	data[i20] = tr - di;
	data[i21] = ti + dr;
  }
};

template<int_t M, typename T, int S>
class DFTk_inplace<2,M,T,S> 
{
public:
  void apply(T* data) 
  { 
      const T tr = data[M];
      const T ti = data[M+1];
      data[M] = data[0]-tr;
      data[M+1] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
// To test
//       const T tr = data[0] - data[M];
//       const T ti = data[1] - data[M+1];
//       data[0] += data[M];
//       data[1] += data[M+1];
//       data[M] = tr;
//       data[M+1] = ti;
  }
  template<class LT>
  void apply(T* data, const LT* wr, const LT* wi) 
  { 
        const T tr = data[M] * (*wr) - data[M+1] * (*wi);
        const T ti = data[M] * (*wi) + data[M+1] * (*wr);
        data[M] = data[0]-tr;
        data[M+1] = data[1]-ti;
        data[0] += tr;
        data[1] += ti;
  }
  template<class LT>
  void apply(const LT* wr, const LT* wi, T* data) 
  { 
        const T tr = data[0] - data[M];
        const T ti = data[1] - data[M+1];
        data[0] += data[M];
        data[1] += data[M+1];
        data[M]   = tr * (*wr) - ti * (*wi);
        data[M+1] = ti * (*wr) + tr * (*wi);
  }  
};



template<int_t N, int_t M, typename T, int S>
class DFTk;

template<int_t M, typename T, int S>
class DFTk<3,M,T,S> 
{
  static const int_t M2 = M+M;
  T m_coef;
  
public:
  DFTk() : m_coef(S * Sqrt<3, T>::value() * 0.5) { } // sqrt(3)/2 = sin(2pi/3)
  
  void apply(const T* src, T* dst) 
  { 
    // 4 mult, 12 add
      const T sr = src[M] + src[M2];
      const T dr = m_coef * (src[M] - src[M2]);
      const T si = src[M+1] + src[M2+1];
      const T di = m_coef * (src[M+1] - src[M2+1]);
      const T tr = src[0] - 0.5*sr;
      const T ti = src[1] - 0.5*si;
      dst[0] = src[0] + sr;
      dst[1] = src[1] + si;
      dst[2] = tr + di;
      dst[3] = ti - dr;
      dst[4] = tr - di;
      dst[5] = ti + dr;     
  }
};

template<int_t M, typename T, int S>
class DFTk<2,M,T,S> 
{
public:
  void apply(const T* src, T* dst) 
  { 
    const T v1(src[1]), v2(src[M]), v3(src[M+1]);
    dst[0] = (*src + v2);
    dst[1] = (v1 + v3);
    dst[2] = (*src - v2);
    dst[3] = (v1 - v3);
  }
  template<class LT>
  void apply(const LT* wr, const LT* wi, const T* src, T* dst) 
  { 
        const T tr = src[0] - src[M];
        const T ti = src[1] - src[M+1];
        dst[0] = src[0] + src[M];
        dst[1] = src[1] + src[M+1];
        dst[M]   = tr * (*wr) - ti * (*wi);
        dst[M+1] = ti * (*wr) + tr * (*wi);
  }
};

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

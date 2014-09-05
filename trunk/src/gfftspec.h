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

#ifndef __gfftspec_h
#define __gfftspec_h

/** \file
    \brief Short-radix FFT specifications
*/

#include "twiddles.h"

namespace GFFT {

/// Out-of-place DFT
/*!
\tparam N length of the data
\tparam SI step in the source data
\tparam DI step in the result data
\tparam T value type
\tparam S sign of the transform (-1 for inverse)

Non-recursive out-of-place DFT for a general (odd) length with 
short-radix specializations for N=2,3
*/
template<int_t N, int_t SI, int_t DI, typename VType, int S,
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
class DFTk;

template<int_t N, int_t SI, int_t DI, typename VType, int S>
class DFTk<N,SI,DI,VType,S,true>
{
  // N is assumed odd, otherwise compiler would not come here
  
  typedef typename VType::ValueType T;
  typedef typename VType::TempType LocalVType;
  static const int_t K = (N-1)/2; 
  static const int_t NSI = N*SI; 
  static const int_t NDI = N*DI; 
   
  LocalVType m_c[K], m_s[K];
  
public:
  DFTk() 
  { 
    ComputeTwiddles<LocalVType, N, S, K>::apply(m_c, m_s);
  }
  
  void apply(const T* src, T* dst) 
  { 
    T sr[K], si[K], dr[K], di[K];
    for (int_t i=0; i<K; ++i) {
      const int_t k = (i+1)*SI;
      sr[i] = src[k]   + src[NSI-k];
      si[i] = src[k+1] + src[NSI-k+1];
      dr[i] = src[k]   - src[NSI-k];
      di[i] = src[k+1] - src[NSI-k+1];
    }
    
    for (int_t i=1; i<K+1; ++i) {
      T re1(0), re2(0), im1(0), im2(0);
      for (int_t j=0; j<K; ++j) {
	const bool sign_change = (i*(j+1) % N) > K;
	const int_t kk = (i+j*i)%N;
	const int_t k = (kk>K) ? N-kk-1 : kk-1;
	const T s1 = m_s[k]*di[j];
	const T s2 = m_s[k]*dr[j];
	re1 += m_c[k]*sr[j];
	im1 += m_c[k]*si[j];
	re2 += sign_change ? -s1 : s1;
	im2 -= sign_change ? -s2 : s2;
      }
      const int_t k = i*DI;
      dst[k] = src[0] + re1 + re2;
      dst[k+1] = src[1] + im1 + im2;
      dst[NDI-k] = src[0] + re1 - re2;
      dst[NDI-k+1] = src[1] + im1 - im2;
    }

    dst[0] = src[0];
    dst[1] = src[1];
    for (int_t i=0; i<K; ++i) {
      dst[0] += sr[i];
      dst[1] += si[i];
    }
  }
};

template<int_t SI, int_t DI, typename VType, int S>
class DFTk<3,SI,DI,VType,S,true> 
{
  typedef typename VType::ValueType T;
  static const int_t SI2 = SI+SI;
  static const int_t DI2 = DI+DI;
  T m_coef;
  
public:
  DFTk() : m_coef(S * MF::Sqrt<3, T>::value() * 0.5) { } // sqrt(3)/2 = sin(2pi/3)
  
  void apply(const T* src, T* dst) 
  { 
    // 4 mult, 12 add
      const T sr = src[SI] + src[SI2];
      const T dr = m_coef * (src[SI] - src[SI2]);
      const T si = src[SI+1] + src[SI2+1];
      const T di = m_coef * (src[SI+1] - src[SI2+1]);
      const T tr = src[0] - 0.5*sr;
      const T ti = src[1] - 0.5*si;
      dst[0]     = src[0] + sr;
      dst[1]     = src[1] + si;
      dst[DI]    = tr + di;
      dst[DI+1]  = ti - dr;
      dst[DI2]   = tr - di;
      dst[DI2+1] = ti + dr;     
  }
/*
  template<class LT>
  void apply(const LT* wr, const LT* wi, const T* src, T* dst) 
  { 
    // 4 mult, 12 add
      const T sr = src[SI] + src[SI2];
      const T dr = m_coef * (src[SI] - src[SI2]);
      const T si = src[SI+1] + src[SI2+1];
      const T di = m_coef * (src[SI+1] - src[SI2+1]);
      const T tr = src[0] - 0.5*sr;
      const T ti = src[1] - 0.5*si;
      const T trpdi = tr + di;
      const T trmdi = tr - di;
      const T tipdr = ti + dr;
      const T timdr = ti - dr;
      dst[0]     = src[0] + sr;
      dst[1]     = src[1] + si;
      dst[DI]    = trpdi*wr[0] - timdr*wi[0];
      dst[DI+1]  = trpdi*wi[0] + timdr*wr[0];
      dst[DI2]   = trmdi*wr[1] - tipdr*wi[1];
      dst[DI2+1] = trmdi*wi[1] + tipdr*wr[1];     
  }
  */
};

template<int_t SI, int_t DI, typename VType, int S>
class DFTk<2,SI,DI,VType,S,true> 
{
  typedef typename VType::ValueType T;
public:
  void apply(const T* src, T* dst) 
  { 
    // the temporaries tr, ti are necessary, because may happen src == dst
        const T tr = src[0] - src[SI];
        const T ti = src[1] - src[SI+1];
        dst[0] = src[0] + src[SI];
        dst[1] = src[1] + src[SI+1];
        dst[DI]   = tr;
        dst[DI+1] = ti;
  }
  /*
  template<class LT>
  void apply(const LT* wr, const LT* wi, const T* src, T* dst) 
  { 
        const T tr = src[0] - src[SI];
        const T ti = src[1] - src[SI+1];
        dst[0] = src[0] + src[SI];
        dst[1] = src[1] + src[SI+1];
        dst[DI]   = tr * (*wr) - ti * (*wi);
        dst[DI+1] = ti * (*wr) + tr * (*wi);
  }
  */
};

/// Out-of-place specialization for complex-valued radix 2 FFT 
/// \tparam T is value type
/// \param data is the array of length 4, containing two complex numbers (real,imag,real,imag).
template<typename T>
inline void _spec2(const T* src, T* dst) 
{ 
    const T v1(src[1]), v2(src[2]), v3(src[3]);
    dst[0] = (*src + v2);
    dst[1] = (v1 + v3);
    dst[2] = (*src - v2);
    dst[3] = (v1 - v3);
}



/// Out-of-place DCT-2
/*!
\tparam N length of the data
\tparam SI step in the source data
\tparam DI step in the result data
\tparam T value type
\tparam S sign of the transform (-1 for inverse)

Non-recursive out-of-place DFT for a general (odd) length with 
short-radix specializations for N=2,3
*/
template<int_t N, int_t SI, int_t DI, typename VType, int S,
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
class DCT2k;
/*
template<int_t N, int_t SI, int_t DI, typename VType, int S>
class DCT2k<N,SI,DI,VType,S,true>
{
  // N is assumed odd, otherwise compiler would not come here
  
  typedef typename VType::ValueType T;
  typedef typename VType::TempType LocalVType;
  static const int_t K = (N-1)/2; 
  static const int_t NSI = N*SI; 
  static const int_t NDI = N*DI; 
   
  LocalVType m_c[K], m_s[K];
  
public:
  DCT2k() 
  { 
    ComputeTwiddles<LocalVType, N, S, K>::apply(m_c, m_s);
  }

  void apply(const T* src, T* dst) 
  { 
    T sr[K], si[K], dr[K], di[K];
    for (int_t i=0; i<K; ++i) {
      const int_t k = (i+1)*SI;
      sr[i] = src[k]   + src[NSI-k];
      si[i] = src[k+1] + src[NSI-k+1];
      dr[i] = src[k]   - src[NSI-k];
      di[i] = src[k+1] - src[NSI-k+1];
    }
    
    for (int_t i=1; i<K+1; ++i) {
      T re1(0), re2(0), im1(0), im2(0);
      for (int_t j=0; j<K; ++j) {
	const bool sign_change = (i*(j+1) % N) > K;
	const int_t kk = (i+j*i)%N;
	const int_t k = (kk>K) ? N-kk-1 : kk-1;
	const T s1 = m_s[k]*di[j];
	const T s2 = m_s[k]*dr[j];
	re1 += m_c[k]*sr[j];
	im1 += m_c[k]*si[j];
	re2 += sign_change ? -s1 : s1;
	im2 -= sign_change ? -s2 : s2;
      }
      const int_t k = i*DI;
      dst[k] = src[0] + re1 + re2;
      dst[k+1] = src[1] + im1 + im2;
      dst[NDI-k] = src[0] + re1 - re2;
      dst[NDI-k+1] = src[1] + im1 - im2;
    }

    dst[0] = src[0];
    dst[1] = src[1];
    for (int_t i=0; i<K; ++i) {
      dst[0] += sr[i];
      dst[1] += si[i];
    }
  }
};

template<int_t SI, int_t DI, typename VType, int S>
class DCT2k<3,SI,DI,VType,S,true> 
{
  typedef typename VType::ValueType T;
  static const int_t SI2 = SI+SI;
  static const int_t DI2 = DI+DI;
  T m_coef;
  
public:
  DCT2k() : m_coef(S * Sqrt<3, T>::value() * 0.5) { } // sqrt(3)/2 = sin(2pi/3)
  
  void apply(const T* src, T* dst) 
  { 
    // 4 mult, 12 add
      const T sr = src[SI] + src[SI2];
      const T dr = m_coef * (src[SI] - src[SI2]);
      const T si = src[SI+1] + src[SI2+1];
      const T di = m_coef * (src[SI+1] - src[SI2+1]);
      const T tr = src[0] - 0.5*sr;
      const T ti = src[1] - 0.5*si;
      dst[0]     = src[0] + sr;
      dst[1]     = src[1] + si;
      dst[DI]    = tr + di;
      dst[DI+1]  = ti - dr;
      dst[DI2]   = tr - di;
      dst[DI2+1] = ti + dr;     
  }
};
*/
template<int_t SI, int_t DI, typename VType, int S>
class DCT2k<2,SI,DI,VType,S,true> 
{
  typedef typename VType::ValueType T;
  T m_coef;
public:
  DCT2k() : m_coef(S * MF::Sqrt<2, T>::value() * 0.5) {} // sqrt(2)/2 = cos(pi/4)
  
  void apply(const T* src, T* dst) 
  { 
    // the temporary dif is necessary, because may happen src == dst
        const T dif = src[0] - src[SI];
        dst[0] = src[0] + src[SI];
        dst[DI] = m_coef*dif;
  }
};
  
  
}  //namespace DFT

#endif /*__gfftspec_h*/

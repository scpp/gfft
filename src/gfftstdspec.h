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

#ifndef __gfftstdspec_h
#define __gfftstdspec_h

/** \file
    \brief Short-radix FFT specifications for "complex" types like std::complex
*/

#include "metafunc.h"

namespace GFFT {

using namespace MF;


/// In-place DFT for "complex" types like std::complex
/*!
\tparam N length of the data
\tparam M step in the data
\tparam T value type
\tparam S sign of the transform (-1 for inverse)

Non-recursive in-place DFT for a general (odd) length with 
short-radix specializations for N=2,3
*/
template<long_t N, long_t M, typename VType, int S>
class DFTk_inp<N,M,VType,S,false>
{
  // N is assumed odd, otherwise compiler would not come here
  
  typedef typename VType::ValueType CT;
  typedef typename VType::TempType LocalVType;
  typedef typename CT::value_type T;
  static const long_t K = (N-1)/2;
  static const long_t NM = N*M;
   
  T m_c[K], m_s[K];
  
  void _transform(CT* data, const CT* s, const CT* d) 
  {
    for (long_t i=1; i<K+1; ++i) {
      CT t1(0,0), t2(0,0);
      for (long_t j=0; j<K; ++j) {
	const bool sign_change = (i*(j+1) % N) > K;
    const long_t kk = (i+j*i)%N;
    const long_t k = (kk>K) ? N-kk-1 : kk-1;
	const T s1 = m_s[k]*d[j].imag();
	const T s2 = m_s[k]*d[j].real();
	t1 += m_c[k]*s[j];
	CT tt(sign_change ? -s1 : s1, sign_change ? s2 : -s2);
	t2 += tt;
      }
      const long_t k = i*M;
      data[k] = data[0] + t1 + t2;
      data[NM-k] = data[0] + t1 - t2;
    }
    
    for (long_t i=0; i<K; ++i)
      data[0] += s[i];
  }
  
public:
  DFTk_inp() 
  { 
    ComputeTwiddles<T, N, S, K>::apply(m_c, m_s);
  }
  
  void apply(CT* data) 
  { 
    CT s[K], d[K];
    for (long_t i=0; i<K; ++i) {
      const long_t k = (i+1)*M;
      s[i] = data[k] + data[NM-k];
      d[i] = data[k] - data[NM-k];
    }
    _transform(data, s, d);
  }

  void apply(CT* data, const CT* w) 
  { 
    CT s[K], d[K];
    for (long_t i=0; i<K; ++i) {
      const long_t k = (i+1)*M;
      CT t1(data[k]*w[i]);
      CT t2(data[NM-k]*w[N-i-2]);
      s[i] = t1 + t2;
      d[i] = t1 - t2;
    }
    _transform(data, s, d);
  }

  void apply_m(CT* data, const CT* w) 
  { 
    data[0] *= w[0];
    apply(data, w+1);
  }
};

template<long_t M, typename VType, int S>
class DFTk_inp<3,M,VType,S,false> 
{
  typedef typename VType::ValueType CT;
  typedef typename VType::TempType LocalVType;
  typedef typename CT::value_type T;

  static const long_t I10 = M;
  static const long_t I20 = M+M;
  
  T m_coef;
  
public:
  DFTk_inp() : m_coef(S * Sqrt<3, T>::value() * 0.5) { }
  
  void apply(CT* data) 
  { 
      CT sum(data[I10] + data[I20]);
      CT dif(m_coef * (data[I10] - data[I20]));
      CT t(data[0] - static_cast<T>(0.5)*sum);
      data[0] += sum;
      data[I10] = CT(t.real() + dif.imag(), t.imag() - dif.real());
      data[I20] = CT(t.real() - dif.imag(), t.imag() + dif.real());
  }
  void apply(CT* data, const CT* w) 
  { 
      CT t1(data[I10]*w[0]);
      CT t2(data[I20]*w[1]);

      CT sum(t1 + t2);
      CT dif(m_coef * (t1 - t2));
      CT t(data[0] - static_cast<T>(0.5)*sum);
      data[0] += sum;
      data[I10] = CT(t.real() + dif.imag(), t.imag() - dif.real());
      data[I20] = CT(t.real() - dif.imag(), t.imag() + dif.real());
  }
  void apply_m(CT* data, const CT* w) 
  { 
      data[0] *= w[0];
      apply(data, w+1);
  }
};

template<long_t M, typename VType, int S>
class DFTk_inp<2,M,VType,S,false> 
{
  typedef typename VType::ValueType CT;
  typedef typename VType::TempType LocalVType;
public:
  void apply(CT* data) 
  { 
     const CT t(data[M]);
     data[M] = data[0] - t;
     data[0] += t;
  }
  // For decimation-in-time
  void apply(CT* data, const CT* w) 
  { 
     const CT t(data[M] * (*w));
     data[M] = data[0] - t;
     data[0] += t;
  }
  void apply_m(CT* data, const CT* w) 
  { 
     const CT t0(data[0] * w[0]);
     const CT t1(data[M] * w[1]);
     data[M] = t0-t1;
     data[0] = t0+t1;
  }
  // For decimation-in-frequency
//   void apply(const Complex<T>* w, Complex<T>* data) 
//   { 
//         Complex<T> t(data[0] - data[M]);
//         data[0] += data[M];
//         data[M]  = t * (*w);
//   }  
};

/// Out-of-place DFT for "complex" types like std::complex
/*!
\tparam N length of the data
\tparam SI step in the source data
\tparam DI step in the result data
\tparam T value type
\tparam S sign of the transform (-1 for inverse)

Non-recursive out-of-place DFT for a general (odd) length with 
short-radix specializations for N=2,3
*/
template<long_t N, long_t SI, long_t DI, typename VType, int S>
class DFTk<N,SI,DI,VType,S,false>
{
  // N is assumed odd, otherwise compiler would not come here
  
  typedef typename VType::ValueType CT;
  typedef typename VType::TempType LocalVType;
  typedef typename CT::value_type T;

  static const long_t K = (N-1)/2;
  static const long_t NSI = N*SI;
  static const long_t NDI = N*DI;
   
  T m_c[K], m_s[K];
  
public:
  DFTk() 
  { 
    ComputeTwiddles<T, N, S, K>::apply(m_c, m_s);
  }
  
  void apply(const CT* src, CT* dst) 
  { 
    CT s[K], d[K];
    for (long_t i=0; i<K; ++i) {
      const long_t k = (i+1)*SI;
      s[i] = src[k] + src[NSI-k];
      d[i] = src[k] - src[NSI-k];
    }
    
    for (long_t i=1; i<K+1; ++i) {
      CT t1(0,0), t2(0,0);
      for (long_t j=0; j<K; ++j) {
	const bool sign_change = (i*(j+1) % N) > K;
    const long_t kk = (i+j*i)%N;
    const long_t k = (kk>K) ? N-kk-1 : kk-1;
	const T s1 = m_s[k]*d[j].imag();
	const T s2 = m_s[k]*d[j].real();
	t1 += m_c[k]*s[j];
	CT tt(sign_change ? -s1 : s1, sign_change ? s2 : -s2);
	t2 += tt;
      }
      const long_t k = i*DI;
      dst[k] = src[0] + t1 + t2;
      dst[NDI-k] = src[0] + t1 - t2;
    }

    dst[0] = src[0];
    for (long_t i=0; i<K; ++i)
      dst[0] += s[i];
  }
};

template<long_t SI, long_t DI, typename VType, int S>
class DFTk<3,SI,DI,VType,S,false> 
{
  typedef typename VType::ValueType CT;
  typedef typename VType::TempType LocalVType;
  typedef typename CT::value_type T;

  static const long_t SI2 = SI+SI;
  static const long_t DI2 = DI+DI;
  T m_coef;
  
public:
  DFTk() : m_coef(S * Sqrt<3, T>::value() * 0.5) { } // sqrt(3)/2 = sin(2pi/3)
  
  void apply(const CT* src, CT* dst) 
  { 
      CT s(src[SI] + src[SI2]);
      CT d(m_coef * (src[SI] - src[SI2]));
      CT t(src[0] - static_cast<T>(0.5)*s);
      dst[0]   = src[0] + s;
      dst[DI]  = CT(t.real() + d.imag(), t.imag() - d.real());
      dst[DI2] = CT(t.real() - d.imag(), t.imag() + d.real());     
  }
};

template<long_t SI, long_t DI, typename VType, int S>
class DFTk<2,SI,DI,VType,S,false> 
{
  typedef typename VType::ValueType CT;
  typedef typename VType::TempType LocalVType;

public:
  void apply(const CT* src, CT* dst) 
  { 
    // the temporary t is necessary, because may happen src == dst
        CT t(src[0] - src[SI]);
        dst[0]  = src[0] + src[SI];
        dst[DI] = t;
  }
};

/// Specialization for complex-valued radix 2 FFT in-place
/// \tparam T is value type
/// \tparam Complex<T> is a generic type representing complex numbers (like std::complex)
/// \param data is the array containing two complex numbers of type Complex<T>.
template<typename T,
template<typename> class Complex>
inline void _spec2(Complex<T>* data) {
      Complex<T> t(data[1]);
      data[1] = data[0]-t;
      data[0] += t;
}

template<typename T,
template<typename> class Complex>
inline void _spec2(const Complex<T>* src, Complex<T>* dst) 
{ 
    dst[0] = src[0] + src[1];
    dst[1] = src[0] - src[1];
}


}  //namespace DFT

#endif /*__gfftspec_h*/

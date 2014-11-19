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

#ifndef __gfftspec_inp_h
#define __gfftspec_inp_h

/** \file
    \brief Short-radix in-place FFT specifications 
*/

#include "twiddles.h"

namespace GFFT {

/// In-place DFT
/*!
\tparam N length of the data
\tparam M step in the data
\tparam T value type
\tparam S sign of the transform (-1 for inverse)

Non-recursive in-place DFT for a general (odd) length with 
short-radix specializations for N=2,3
*/
template<int_t N, int_t M, typename VType, int S, typename W,
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
class DFTk_inp;

template<int_t N, int_t M, typename VType, int S, typename W>
class DFTk_inp<N,M,VType,S,W,true>
{
  // N is assumed odd, otherwise compiler would not come here
  
  typedef typename VType::ValueType T;
  typedef typename VType::TempType LocalVType;
  static const int_t K = (N-1)/2; 
  static const int_t NM = N*M; 
   
  LocalVType m_c[K], m_s[K];
  
public:
  DFTk_inp() 
  { 
    ComputeTwiddles<LocalVType, N, S, K>::apply(m_c, m_s);
  }
  
  void apply(T* data) 
  { 
    T sr[K], si[K], dr[K], di[K];
    for (int_t i=0; i<K; ++i) {
      const int_t k = (i+1)*M;
      sr[i] = data[k]   + data[NM-k];
      si[i] = data[k+1] + data[NM-k+1];
      dr[i] = data[k]   - data[NM-k];
      di[i] = data[k+1] - data[NM-k+1];
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
      const int_t k = i*M;
      data[k] = data[0] + re1 + re2;
      data[k+1] = data[1] + im1 + im2;
      data[NM-k] = data[0] + re1 - re2;
      data[NM-k+1] = data[1] + im1 - im2;
    }
    
    for (int_t i=0; i<K; ++i) {
      data[0] += sr[i];
      data[1] += si[i];
    }
  }

  template<class LT>
  void apply(T* data, const LT* wr, const LT* wi) 
  { 
    T sr[K], si[K], dr[K], di[K];
    for (int_t i=0; i<K; ++i) {
      const int_t k = (i+1)*M;
      const T tr1 = data[k]*wr[i] - data[k+1]*wi[i];
      const T ti1 = data[k]*wi[i] + data[k+1]*wr[i];
      const T tr2 = data[NM-k]*wr[N-i-2] - data[NM-k+1]*wi[N-i-2];
      const T ti2 = data[NM-k]*wi[N-i-2] + data[NM-k+1]*wr[N-i-2];
      sr[i] = tr1 + tr2;
      si[i] = ti1 + ti2;
      dr[i] = tr1 - tr2;
      di[i] = ti1 - ti2;
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
      const int_t k = i*M;
      data[k] = data[0] + re1 + re2;
      data[k+1] = data[1] + im1 + im2;
      data[NM-k] = data[0] + re1 - re2;
      data[NM-k+1] = data[1] + im1 - im2;
    }
    
    for (int_t i=0; i<K; ++i) {
      data[0] += sr[i];
      data[1] += si[i];
    }
  }
  
  template<class LT>
  void apply_m(T* data, const LT* wr, const LT* wi) 
  { 
    const T t = data[0];
    data[0] = t*wr[0] - data[1]*wi[0];
    data[1] = t*wi[0] + data[1]*wr[0];
    T sr[K], si[K], dr[K], di[K];
    for (int_t i=1; i<K+1; ++i) {
      const int_t k = i*M;
      const T tr1 = data[k]*wr[i] - data[k+1]*wi[i];
      const T ti1 = data[k]*wi[i] + data[k+1]*wr[i];
      const T tr2 = data[NM-k]*wr[N-i] - data[NM-k+1]*wi[N-i];
      const T ti2 = data[NM-k]*wi[N-i] + data[NM-k+1]*wr[N-i];
      sr[i-1] = tr1 + tr2;
      si[i-1] = ti1 + ti2;
      dr[i-1] = tr1 - tr2;
      di[i-1] = ti1 - ti2;
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
      const int_t k = i*M;
      data[k] = data[0] + re1 + re2;
      data[k+1] = data[1] + im1 + im2;
      data[NM-k] = data[0] + re1 - re2;
      data[NM-k+1] = data[1] + im1 - im2;
    }
    
    for (int_t i=0; i<K; ++i) {
      data[0] += sr[i];
      data[1] += si[i];
    }
  }

  template<class LT>
  void apply(const LT* wr, const LT* wi, T* data) 
  { 
    T sr[K], si[K], dr[K], di[K];
    for (int_t i=0; i<K; ++i) {
      const int_t k = (i+1)*M;
      sr[i] = data[k]   + data[NM-k];
      si[i] = data[k+1] + data[NM-k+1];
      dr[i] = data[k]   - data[NM-k];
      di[i] = data[k+1] - data[NM-k+1];
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
      const int_t k = i*M;
      const T tr1 = data[0] + re1 + re2;
      const T ti1 = data[1] + im1 + im2;
      const T tr2 = data[0] + re1 - re2;
      const T ti2 = data[1] + im1 - im2;
      data[k] = tr1*wr[i-1] - ti1*wi[i-1];
      data[k+1] = tr1*wi[i-1] + ti1*wr[i-1];
      data[NM-k] = tr2*wr[K-i+1] - ti2*wi[K-i+1];
      data[NM-k+1] = tr2*wi[K-i+1] + ti2*wr[K-i+1];
    }
    
    for (int_t i=0; i<K; ++i) {
      data[0] += sr[i];
      data[1] += si[i];
    }
  }
};
/*
template<int_t M, typename VType, int S>
class DFTk_inp<4,M,VType,S,true> 
{
  typedef typename VType::ValueType T;

  static const int_t I10 = M;
  static const int_t I11 = M+1;
  static const int_t I20 = M+M;
  static const int_t I21 = I20+1;
  static const int_t I30 = I20+M;
  static const int_t I31 = I30+1;
  
public:
  DFTk_inp() { }
  
//   void apply(T* data) 
//   { 
//       T tr = data[I20];
//       T ti = data[I21];
//       data[I20] = data[0]-tr;
//       data[I21] = data[1]-ti;
//       data[0] += tr;
//       data[1] += ti;
//       tr = data[I30];
//       ti = data[I31];
//       data[I30] = S*(data[I11]-ti);
//       data[I31] = S*(tr-data[I10]);
//       data[I10] += tr;
//       data[I11] += ti;
// 
//       tr = data[I10];
//       ti = data[I11];
//       data[I10] = data[0]-tr;
//       data[I11] = data[1]-ti;
//       data[0] += tr;
//       data[1] += ti;
//       tr = data[I30];
//       ti = data[I31];
//       data[I30] = data[I20]-tr;
//       data[I31] = data[I21]-ti;
//       data[I20] += tr;
//       data[I21] += ti;
//   }

  void apply(T* data) 
  {
      const T sr1 = data[0] + data[I20];
      const T dr1 = data[0] - data[I20];
      const T sr2 = data[I10] + data[I30];
      const T dr2 = data[I10] - data[I30];
      const T si1 = data[1] + data[I21];
      const T di1 = data[1] - data[I21];
      const T si2 = data[I11] + data[I31];
      const T di2 = data[I11] - data[I31];
      data[0]   = sr1 + sr2;
      data[1]   = si1 + si2;
      data[I10] = dr1 + di2;
      data[I11] = di1 - dr2;
      data[I20] = sr1 - sr2;
      data[I21] = si1 - si2;
      data[I30] = dr1 - di2;
      data[I31] = di1 + dr2;
  }
  template<class LT>
  void apply(T* data, const LT* wr, const LT* wi) 
  { 
      const T tr1 = data[I10]*wr[0] - data[I11]*wi[0];
      const T ti1 = data[I10]*wi[0] + data[I11]*wr[0];
      const T tr2 = data[I20]*wr[1] - data[I21]*wi[1];
      const T ti2 = data[I20]*wi[1] + data[I21]*wr[1];
      const T tr3 = data[I30]*wr[2] - data[I31]*wi[2];
      const T ti3 = data[I30]*wi[2] + data[I31]*wr[2];
      
      const T sr1 = data[0] + tr2;
      const T dr1 = data[0] - tr2;
      const T sr2 = tr1 + tr3;
      const T dr2 = tr1 - tr3;
      const T si1 = data[1] + ti2;
      const T di1 = data[1] - ti2;
      const T si2 = ti1 + ti3;
      const T di2 = ti1 - ti3;
      
      data[0]   = sr1 + sr2;
      data[1]   = si1 + si2;
      data[I10] = dr1 + di2;
      data[I11] = di1 - dr2;
      data[I20] = sr1 - sr2;
      data[I21] = si1 - si2;
      data[I30] = dr1 - di2;
      data[I31] = di1 + dr2;
  }
  
  template<class LT>
  void apply(const LT* wr, const LT* wi, T* data) 
  { 
  }
};
*/

// Specialization for N=3
template<int_t M, typename VType, int S, typename W>
class DFTk_inp<3,M,VType,S,W,true> 
{
  static const int_t I10 = M;
  static const int_t I11 = M+1;
  static const int_t I20 = M+M;
  static const int_t I21 = I20+1;
  
  typedef typename VType::ValueType T;
  static const int Acc = VType::Accuracy;
  typedef Compute<typename MF::SqrtDecAcc<SInt<3>,Acc>::Result,Acc,T> CSqrt3;

  T m_coef;
  
public:
  DFTk_inp() : m_coef(S * CSqrt3::value() * 0.5) { }
  
  void apply(T* data) 
  { 
      const T sum_r = data[I10] + data[I20];
      const T dif_r = m_coef * (data[I10] - data[I20]);
      const T sum_i = data[I11] + data[I21];
      const T dif_i = m_coef * (data[I11] - data[I21]);
      const T tr = data[0] - 0.5*sum_r;
      const T ti = data[1] - 0.5*sum_i;
      data[0] += sum_r;
      data[1] += sum_i;
      data[I10] = tr + dif_i;
      data[I11] = ti - dif_r;
      data[I20] = tr - dif_i;
      data[I21] = ti + dif_r;
  }
  template<class LT>
  void apply(T* data, const LT* wr, const LT* wi) 
  { 
        const T tr1 = data[I10]*wr[0] - data[I11]*wi[0];
        const T ti1 = data[I10]*wi[0] + data[I11]*wr[0];
        const T tr2 = data[I20]*wr[1] - data[I21]*wi[1];
        const T ti2 = data[I20]*wi[1] + data[I21]*wr[1];

	const T sum_r = tr1 + tr2;
	const T dif_r = m_coef * (tr1 - tr2);
	const T sum_i = ti1 + ti2;
	const T dif_i = m_coef * (ti1 - ti2);
	const T tr = data[0] - 0.5*sum_r;
	const T ti = data[1] - 0.5*sum_i;
	data[0] += sum_r;
	data[1] += sum_i;
	data[I10] = tr + dif_i;
	data[I11] = ti - dif_r;
	data[I20] = tr - dif_i;
	data[I21] = ti + dif_r;
  }
  template<class LT>
  void apply_m(T* data, const LT* wr, const LT* wi) 
  { 
        const T tr0 = data[0]*wr[0] - data[1]*wi[0];
        const T ti0 = data[0]*wi[0] + data[1]*wr[0];
        const T tr1 = data[I10]*wr[1] - data[I11]*wi[1];
        const T ti1 = data[I10]*wi[1] + data[I11]*wr[1];
        const T tr2 = data[I20]*wr[2] - data[I21]*wi[2];
        const T ti2 = data[I20]*wi[2] + data[I21]*wr[2];

	const T sum_r = tr1 + tr2;
	const T dif_r = m_coef * (tr1 - tr2);
	const T sum_i = ti1 + ti2;
	const T dif_i = m_coef * (ti1 - ti2);
	const T tr = tr0 - 0.5*sum_r;
	const T ti = ti0 - 0.5*sum_i;
	data[0] = tr0 + sum_r;
	data[1] = ti0 + sum_i;
	data[I10] = tr + dif_i;
	data[I11] = ti - dif_r;
	data[I20] = tr - dif_i;
	data[I21] = ti + dif_r;
  }
  template<class LT>
  void apply(const LT* wr, const LT* wi, T* data) 
  { 
      const T sum_r = data[I10] + data[I20];
      const T dif_r = m_coef * (data[I10] - data[I20]);
      const T sum_i = data[I11] + data[I21];
      const T dif_i = m_coef * (data[I11] - data[I21]);
      const T tr = data[0] - 0.5*sum_r;
      const T ti = data[1] - 0.5*sum_i;
      const T trpdi = tr + dif_i;
      const T trmdi = tr - dif_i;
      const T tipdr = ti + dif_r;
      const T timdr = ti - dif_r;
      data[0] += sum_r;
      data[1] += sum_i;
      data[I10] = trpdi*wr[0] - timdr*wi[0];
      data[I11] = trpdi*wi[0] + timdr*wr[0];
      data[I20] = trmdi*wr[1] - tipdr*wi[1];
      data[I21] = trmdi*wi[1] + tipdr*wr[1];
  }
};

// Specialization for N=2
template<int_t M, typename VType, int S, typename W>
class DFTk_inp<2,M,VType,S,W,true> 
{
  typedef typename VType::ValueType T;
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
  // For decimation-in-time
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
  // as one above with wr = 0, wi = -1
  void apply_1(T* data) 
  { 
        const T tr = data[M+1];
        const T ti = -data[M];
        data[M] = data[0]-tr;
        data[M+1] = data[1]-ti;
        data[0] += tr;
        data[1] += ti;

//         data[0] += data[M+1];
//         data[1] -= data[M];
// 	const T tr = ;
//         const T ti = ;
//         data[M] = data[0]-tr;
//         data[M+1] = data[1]-ti;
  }
  template<class LT>
  void apply_m(T* data, const LT* wr, const LT* wi) 
  { 
        const T tr0 = data[0] * wr[0] - data[1] * wi[0];
        const T ti0 = data[0] * wi[0] + data[1] * wr[0];
        const T tr1 = data[M] * wr[1] - data[M+1] * wi[1];
        const T ti1 = data[M] * wi[1] + data[M+1] * wr[1];
        data[M] = tr0-tr1;
        data[M+1] = ti0-ti1;
        data[0] = tr0+tr1;
        data[1] = ti0+ti1;
  }
  // For decimation-in-frequency
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

template<int_t M, typename VType, int S, typename W>
class DFTk_inp<1,M,VType,S,W,true> 
{
};

/// In-place specialization for complex-valued radix 2 FFT 
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


}  //namespace DFT

#endif /*__gfftspec_inp_h*/

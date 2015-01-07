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

#include "metaroot.h"

namespace GFFT {

using namespace MF;

static const int uninitialized_flag = 1000; 

template<int_t N, typename NList, typename VType, typename W1, int S, int_t LastK=1, bool C = (N>=4)>
class _RootsCompute;

template<int_t N, typename Head, typename Tail, typename VType, typename W1, int S, int_t LastK>
class _RootsCompute<N,Loki::Typelist<Head,Tail>,VType,W1,S,LastK, true>
{
   typedef typename VType::ValueType T;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   static const int_t NN = N*LastK-2;
   static const int_t Step = 2*LastK;
   
   typedef Compute<typename W1::Re,VType::Accuracy> WR;
   typedef Compute<typename W1::Im,VType::Accuracy> WI;

   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   
   _RootsCompute<M,NFactNext,VType,WK,S,K*LastK> next;
public:
  void init(std::vector<T>& w)
  {
    next.init(w);
    
    const T wpr = WR::value();
    const T wpi = WI::value();
    w[Step-2] = wpr;
    w[Step-1] = wpi;

    for (int_t i=Step+Step-2; i<NN; i+=Step) {
      if (w[i] == uninitialized_flag) {
	w[i]   = w[i-Step]*wpr - w[i-Step+1]*wpi;
	w[i+1] = w[i-Step+1]*wpr + w[i-Step]*wpi;
      }
    }
  }
};

template<int_t N, int_t K, typename Tail, typename VType, class W, int S, int_t LastK>
class _RootsCompute<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, VType, W, S, LastK, true>
: public _RootsCompute<N, Tail, VType, W, S, LastK, true> {};

template<int_t N, int_t K, typename Tail, typename VType, class W, int S, int_t LastK>
class _RootsCompute<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, VType, W, S, LastK, false>
: public _RootsCompute<N, Tail, VType, W, S, LastK, false> {};

// Specialization for a prime N
template<int_t N, typename Head, typename VType, class W, int S, int_t LastK>
class _RootsCompute<N,Loki::Typelist<Head, Loki::NullType>, VType, W, S, LastK, false> 
{
   typedef typename VType::ValueType T;
   typedef Compute<typename W::Re,VType::Accuracy> WR;
   typedef Compute<typename W::Im,VType::Accuracy> WI;
   static const int_t NN = N*LastK-2;
   static const int_t Step = 2*LastK;
public:
  void init(std::vector<T>& w)
  {
    const T wpr = WR::value();
    const T wpi = WI::value();
    w[Step-2] = wpr;
    w[Step-1] = wpi;
    for (int_t i=Step+Step-2; i<NN; i+=Step) {
      w[i]   = w[i-Step]*wpr - w[i-Step+1]*wpi;
      w[i+1] = w[i-Step+1]*wpr + w[i-Step]*wpi;
    }
  }  
};
  
// Specialization for N=4
template<typename Head, typename VType, class W, int S, int_t LastK>
class _RootsCompute<4,Loki::Typelist<Head, Loki::NullType>, VType, W, S, LastK, false> 
{
   typedef typename VType::ValueType T;
   static const int_t Step = 2*LastK;
public:
  void init(std::vector<T>& w)
  {
    w[Step-2] = 0;
    w[Step-1] = -1;
  }  
};

template<int_t N, typename NList, typename VType, int Sign, bool C = (N>=4)>
class RootsHolder;

template<int_t N, typename NList, typename VType, int Sign>
class RootsHolder<N,NList,VType,Sign,true>
{
  typedef typename VType::ValueType T;
  //typedef typename VType::TempType LT;

  typedef typename GetFirstRoot<N,Sign,VType::Accuracy>::Result W1;
  typedef Compute<typename W1::Re,VType::Accuracy> WR;
  typedef Compute<typename W1::Im,VType::Accuracy> WI;
  
  std::vector<T> m_data;
  
//  _RootsCompute<N,NList,VType,W1,Sign> comp;
public:
  static const int_t Length = (N-1)/2;
  
  RootsHolder() { init(); }
  
  T* getData() { return &(m_data[0]); }
  
  void init()
  {
      m_data.resize(2*Length);
      std::fill(m_data.begin(), m_data.end(), uninitialized_flag);
//      comp.init(m_data);

      const T wpr = WR::value();
      const T wpi = WI::value();
      m_data[0] = wpr;
      m_data[1] = wpi;

      for (int_t i=2; i<2*Length; i+=2) {
        m_data[i]   = m_data[i-2]*wpr - m_data[i-1]*wpi;
        m_data[i+1] = m_data[i-1]*wpr + m_data[i-2]*wpi;
      }
  }
};
  
template<int_t N, typename NList, typename VType, int Sign>
class RootsHolder<N,NList,VType,Sign,false>
{
  typedef typename VType::ValueType T;
public:
  RootsHolder() {}
  T* getData() { return NULL; }
  void init() {}
};


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

///////////////////////////////////////////////////////////////////

template<int_t K, typename VType>
class RootsContainer
{
protected:
  typedef typename VType::TempType T;
  T wr[K-1], wi[K-1], wpr[K-1], wpi[K-1];
public:
  const T* get_real() const { return wr; }
  const T* get_imag() const { return wi; }
  
  void step()
  {
    T t;
    for (int_t i=0; i<K-1; ++i) {
      t = wr[i];
      wr[i] = t*wpr[i] - wi[i]*wpi[i];
      wi[i] = wi[i]*wpr[i] + t*wpi[i];
    }
  }
};


template<int_t K, typename VType, typename W1, typename Permut = Loki::NullType>
class ComputeRoots : public RootsContainer<K,VType> 
{
  typedef typename VType::TempType T;
  typedef RootsContainer<K,VType> Base;
  
  typedef Compute<typename W1::Re,VType::Accuracy> WR;
  typedef Compute<typename W1::Im,VType::Accuracy> WI;

  using Base::wpr;
  using Base::wpi;
  using Base::wr;
  using Base::wi;
  
public:
  ComputeRoots() { init(); }
  ComputeRoots(const int nthreads, const int n) 
  { 
    init(); 
    T tr, ti, t;
    for (int_t i=0; i<K-1; ++i) {
      tr = wr[i];
      ti = wi[i];
      for (int_t j=1; j<(n==0 ? nthreads : n); ++j) {
	t = wr[i];
	wr[i] = t*tr - wi[i]*ti;
	wi[i] = t*ti + tr*wi[i];
      }
    }
    for (int_t i=0; i<K-1; ++i) {
      tr = wpr[i];
      ti = wpi[i];
      for (int_t j=1; j<nthreads; ++j) {
	t = wpr[i];
	wpr[i] = t*tr - wpi[i]*ti;
	wpi[i] = t*ti + tr*wpi[i];
      }
    }
  }
  
  void init() 
  {
	// W = (wpr[0], wpi[0])
	wpr[0] = WR::value();
	wpi[0] = WI::value();
	
	// W^i = (wpr[i], wpi[i])
	for (int_t i=0; i<K-2; ++i) {
	  wpr[i+1] = wpr[i]*wpr[0] - wpi[i]*wpi[0];
	  wpi[i+1] = wpr[i]*wpi[0] + wpr[0]*wpi[i];
	}
	
	for (int_t i=0; i<K-1; ++i) {
	  int_t ii = Permut::value(i+1) - 1;
	  wr[ii] = wpr[i];
	  wi[ii] = wpi[i];
	}
	
	for (int_t i=0; i<K-1; ++i) {
	  wpr[i] = wr[i];
	  wpi[i] = wi[i];
	}
  }  
};

// Specialization without permutation
template<int_t K, typename VType, typename W1>
class ComputeRoots<K,VType,W1,Loki::NullType> : public RootsContainer<K,VType> 
{
  typedef typename VType::TempType T;
  typedef RootsContainer<K,VType> Base;
  
  typedef Compute<typename W1::Re,VType::Accuracy> WR;
  typedef Compute<typename W1::Im,VType::Accuracy> WI;

  using Base::wpr;
  using Base::wpi;
  using Base::wr;
  using Base::wi;
  
public:
  ComputeRoots() { init(); }
  
  const T* get_real() const { return wr; }
  const T* get_imag() const { return wi; }
  
  void init() 
  {
      wpr[0] = WR::value();
      wpi[0] = WI::value();
      //t = Sin<N,1,LocalVType>::value();
//       wpr[0] = 1 - 2.0*t*t;
//       wpi[0] = -S*Sin<N,2,LocalVType>::value();
      
      // W^i = (wpr[i], wpi[i])
      for (int_t i=0; i<K-2; ++i) {
	wpr[i+1] = wpr[i]*wpr[0] - wpi[i]*wpi[0];
	wpi[i+1] = wpr[i]*wpi[0] + wpr[0]*wpi[i];
      }
      
      for (int_t i=0; i<K-1; ++i) {
	wr[i] = wpr[i];
	wi[i] = wpi[i];
      }
  }
};



template<int_t K, typename VType, typename W1, typename Permut = Loki::NullType>
class ComputeRootsStd
{
  
};

// Specialization without permutation
template<int_t K, typename VType, typename W1>
class ComputeRootsStd<K,VType,W1,Loki::NullType> 
{
  typedef typename VType::ValueType CT;
  
  typedef Compute<typename W1::Re,VType::Accuracy> WR;
  typedef Compute<typename W1::Im,VType::Accuracy> WI;

  CT w[K-1], wp[K-1];
  
public:
  ComputeRootsStd() { init(); }
  
  const CT* get() const { return w; }
  
  void init() 
  {
    wp[0] = CT(WR::value(), WI::value());
    //LocalVType t = Sin<N,1,LocalVType>::value();
//       wp[0] = Complex<LocalVType>(1 - 2.0*t*t, -S*Sin<N,2,LocalVType>::value());
      
    // W^i = (wpr2, wpi2)
    for (int_t i=0; i<K-2; ++i) 
      wp[i+1] = wp[i]*wp[0];
      
    for (int_t i=0; i<K-1; ++i) 
      w[i] = wp[i];
  }
  
  void step()
  {
    for (int_t i=0; i<K-1; ++i) 
      w[i] = w[i]*wp[i];
  }
};

}  //namespace DFT

#endif /*__twiddles_h*/

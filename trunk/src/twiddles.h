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

#include "metasincos.h"
#include "metasqrt.h"

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
  ComputeRoots(const int n) 
  { 
    init(); 
    for (int_t j=1; j<n; ++j) 
      Base::step();
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

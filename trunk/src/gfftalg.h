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

#ifndef __gfftalg_h
#define __gfftalg_h

/** \file
    \brief Recursive FFT algorithms 
*/

#include "gfftspec.h"
#include "gfftspec_inp.h"
#include "gfftfactor.h"
#include "gfftswap.h"

#include "metacomplex.h"
#include "metaroot.h"

namespace GFFT {

using namespace MF;


static const int_t StaticLoopLimit = 8;


template<int_t K, int_t M, typename VType, int S, class W1, int NIter = 1, class W = W1>
class IterateInTime
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   typedef Compute<typename W::Re,VType::Accuracy> WR;
   typedef Compute<typename W::Im,VType::Accuracy> WI;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N = K*M;

   typedef typename GetNextRoot<NIter+1,N,W1,W,VType::Accuracy>::Result Wnext;
   IterateInTime<K,M,VType,S,W1,NIter+1,Wnext> next;
   DFTk_inp<K,M2,VType,S> spec_inp;
   
public:
   void apply(T* data) 
   {
      const LocalVType wr = WR::value();
      const LocalVType wi = WI::value();

      spec_inp.apply(data + (NIter-1)*C, &wr, &wi);

      next.apply(data);
   }
};

// Last step of the loop
template<int_t K, int_t M, typename VType, int S, class W1, class W>
class IterateInTime<K,M,VType,S,W1,M,W> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   typedef Compute<typename W::Re,VType::Accuracy> WR;
   typedef Compute<typename W::Im,VType::Accuracy> WI;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N = K*M;
   DFTk_inp<K,M2,VType,S> spec_inp;
public:
   void apply(T* data) 
   {
      const LocalVType wr = WR::value();
      const LocalVType wi = WI::value();

      spec_inp.apply(data + (M-1)*C, &wr, &wi);
   }
};

// First step in the loop
template<int_t K, int_t M, typename VType, int S, class W1, class W>
class IterateInTime<K,M,VType,S,W1,1,W> 
{
   typedef typename VType::ValueType T;
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   DFTk_inp<K,M2,VType,S> spec_inp;
   IterateInTime<K,M,VType,S,W1,2,W> next;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);
      next.apply(data);
   }
};


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
\sa InTime, IterateInTime
*/
template<int_t K, typename KFact, int_t M, int_t Step, typename VType, int S, class W, bool doStaticLoop,
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
class DFTk_x_Im_T;

// Rely on the static template loop
template<int_t K, int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<K,Loki::Typelist<Pair<SInt<K>,SInt<1> >,Loki::NullType>,M,Step,VType,S,W,true,true> 
: public IterateInTime<K,M,VType,S,W> {};

// General implementation
template<int_t K, int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<K,Loki::Typelist<Pair<SInt<K>,SInt<1> >,Loki::NullType>,M,Step,VType,S,W,false,true>
{
   typedef typename VType::ValueType T;
   static const int_t N = K*M;
   static const int_t M2 = M*2;
   static const int_t S2 = 2*Step;
   DFTk_inp<K,M2,VType,S> spec_inp;
   
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      ComputeRoots<K,VType,W> roots;

      spec_inp.apply(data+S2, roots.get_real(), roots.get_imag());
      for (int_t j=S2+S2; j<M2; j+=S2) {
	roots.step();
	spec_inp.apply(data+j, roots.get_real(), roots.get_imag());
      }
   }
  
};

template<int_t K, int_t KK, int_t M, int_t Step, typename Tail, typename VType, int S, class W>
class DFTk_x_Im_T<K,Loki::Typelist<Pair<SInt<KK>,SInt<0> >,Tail>,M,Step,VType,S,W,false,true>
: public DFTk_x_Im_T<K,Tail,M,Step,VType,S,W,false,true> {};


/// In-place decimation-in-time FFT version
/**
\tparam N current transform length
\tparam NFact factorization list
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam W1 compile-time root of unity

This is the core of decimation in-time FFT algorithm:
Strided DFT runs K times recursively, where the next 
factor K is taken from the compile-time list.
The scaled DFT is performed afterwards.
\sa InFreq, DFTk_x_Im_T
*/
template<int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1>
class InTime;

template<int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTime<N, Loki::Typelist<Head,Tail>, VType, S, W1, LastK>
{
   typedef typename VType::ValueType T;
//   // Not implemented, because not allowed
   void apply(T* data) 
   {
//#error Transforms in-place are allowed for powers of primes only!!!
   }
};

template<int_t N, typename Head, typename VType, int S, class W1, int_t LastK>
class InTime<N, Loki::Typelist<Head,Loki::NullType>, VType, S, W1, LastK>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Loki::NullType> NFactNext;
   typedef Loki::Typelist<Pair<SInt<K>,SInt<1> >, Loki::NullType> KFact;
   InTime<M,NFactNext,VType,S,WK,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T<K,KFact,M,1,VType,S,W1,false> dft_scaled;
public:
   void apply(T* data) 
   {
     // run strided DFT recursively K times
      for (int_t m=0; m < N2; m+=M2) 
	dft_str.apply(data + m);

      dft_scaled.apply(data);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTime<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, VType, S, W1, LastK>
: public InTime<N, Tail, VType, S, W1, LastK> {};


// Specialization for a prime N
template<int_t N, typename VType, int S, class W1, int_t LastK>
class InTime<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,VType,S,W1,LastK> 
{
  typedef typename VType::ValueType T;
  static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
  DFTk_inp<N, C, VType, S> spec_inp;
public:
  void apply(T* data) 
  { 
    spec_inp.apply(data);
  }
};


/// Out-of-place decimation-in-time FFT version
/**
\tparam N current transform length
\tparam NFact factorization list
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam W1 compile-time root of unity

This is the core of decimation in-time FFT algorithm:
Strided DFT runs K times recursively, where the next 
factor K is taken from the compile-time list.
The scaled DFT is performed afterwards.
\sa DFTk_x_Im_T
*/
template<int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1>
class InTimeOOP;

template<int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTimeOOP<N, Loki::Typelist<Head,Tail>, VType, S, W1, LastK>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;
   static const int_t LastK2 = LastK*C;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   typedef Loki::Typelist<Pair<SInt<K>,SInt<1> >, Loki::NullType> KFact;
   InTimeOOP<M,NFactNext,VType,S,WK,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T<K,KFact,M,1,VType,S,W1,false> dft_scaled;
public:

   void apply(const T* src, T* dst) 
   {
     // run strided DFT recursively K times
      int_t lk = 0;
      for (int_t m = 0; m < N2; m+=M2, lk+=LastK2)
        dft_str.apply(src + lk, dst + m);

      dft_scaled.apply(dst);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTimeOOP<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, VType, S, W1, LastK>
: public InTimeOOP<N, Tail, VType, S, W1, LastK> {};


// Specialization for prime N
template<int_t N, typename VType, int S, class W1, int_t LastK>
class InTimeOOP<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,VType,S,W1,LastK> 
{
   typedef typename VType::ValueType T;
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   DFTk<N, LastK*C, C, VType, S> spec;
public:
   void apply(const T* src, T* dst) { spec.apply(src, dst); }
};

  

/// Out-of-place DCT-2
/**
\tparam N current transform length
\tparam NFact factorization list
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam W1 compile-time root of unity

This is the core of decimation in-time FFT algorithm:
Strided DFT runs K times recursively, where the next 
factor K is taken from the compile-time list.
The scaled DFT is performed afterwards.
\sa DFTk_x_Im_T
*/
template<int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1>
class DCT2_impl;

template<int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class DCT2_impl<N, Loki::Typelist<Head,Tail>, VType, S, W1, LastK>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;
   static const int_t LastK2 = LastK*C;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   DCT2_impl<M,NFactNext,VType,S,WK,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,1,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
//   DCTk_x_Im_T<K,M,1,VType,S,W1,false> dft_scaled;
public:

   void apply(const T* src, T* dst) 
   {
     // run strided DCT recursively K times
      int_t lk = 0;
      for (int_t m = 0; m < N2; m+=M2, lk+=LastK2)
        dft_str.apply(src + lk, dst + m);

      //dft_scaled.apply(dst);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename VType, int S, class W1, int_t LastK>
class DCT2_impl<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, VType, S, W1, LastK>
: public DCT2_impl<N, Tail, VType, S, W1, LastK> {};


// Specialization for a prime N
template<int_t N, typename VType, int S, class W1, int_t LastK>
class DCT2_impl<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,VType,S,W1,LastK> 
{
   typedef typename VType::ValueType T;
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   DCT2k<N, LastK*C, C, VType, S> spec;
public:
   void apply(const T* src, T* dst) { spec.apply(src, dst); }
};


}  //namespace DFT

#endif /*__gfftalg_h*/

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

template<int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1>
class InTime;

template<int_t K, int_t M, int_t Step, typename VType, int S, class W, bool doStaticLoop,
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
class DFTk_x_Im_T;

template<int_t K, int_t M, typename VType, int S>
struct DFTk_inp_adapter;

template<int_t M, typename VType, int S>
struct DFTk_inp_adapter<2,M,VType,S> : public DFTk_inp<2,M,VType,S> { };

template<int_t M, typename VType, int S>
struct DFTk_inp_adapter<3,M,VType,S> : public DFTk_inp<3,M,VType,S> { };

template<int_t M, typename VType, int S>
struct DFTk_inp_adapter<4,M,VType,S> 
{ 
   typedef typename VType::ValueType T;

   typedef Pair<SInt<2>,SInt<2> > T1;
   typedef TYPELIST_1(T1) NFact;

   typedef typename GetFirstRoot<4,S,VType::Accuracy>::Result W1;

   //InTime<4,NFact,VType,S,W1,2> dft_str;
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*2;
   DFTk_inp<2, M2, VType, S> spec_inp;
   DFTk_x_Im_T<2,M2,M,VType,S,W1,false> dft_scaled;
public:

   void apply(T* data) 
   {
      std::swap(data[M2],data[M2+M2]);
      std::swap(data[M2+1],data[M2+M2+1]);
      std::swap(data[M2+2],data[M2+M2+2]);
      std::swap(data[M2+3],data[M2+M2+3]);

      spec_inp.apply(data);
      spec_inp.apply(data + M2+M2);

      dft_scaled.apply(data);
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
// template<int_t K, int_t M, typename VType, int S, class W, bool doStaticLoop,
// bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
// class DFTk_x_Im_T;

// Rely on the static template loop
template<int_t K, int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<K,M,Step,VType,S,W,true,true> : public IterateInTime<K,M,VType,S,W> {};

// General implementation
template<int_t K, int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<K,M,Step,VType,S,W,false,true>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   typedef Compute<typename W::Re,VType::Accuracy> WR;
   typedef Compute<typename W::Im,VType::Accuracy> WI;
   static const int_t N = K*M;
   static const int_t M2 = M*2;
   static const int_t S2 = 2*Step;
   DFTk_inp<K,M2,VType,S> spec_inp;
   //DFTk_inp_adapter<K,M,VType,S> spec_inp_a;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[K-1], wi[K-1], wpr[K-1], wpi[K-1], t;

      // W = (wpr[0], wpi[0])
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
      
      for (int_t i=S2; i<M2; i+=S2) {
	spec_inp.apply(data+i, wr, wi);

	for (int_t i=0; i<K-1; ++i) {
	  t = wr[i];
	  wr[i] = t*wpr[i] - wi[i]*wpi[i];
	  wi[i] = wi[i]*wpr[i] + t*wpi[i];
	}
      }
   }
  
};

// Specialization for radix 3
template<int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<3,M,Step,VType,S,W,false,true> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   typedef Compute<typename W::Re,VType::Accuracy> WR;
   typedef Compute<typename W::Im,VType::Accuracy> WI;
   static const int_t N = 3*M;
   static const int_t M2 = M*2;
   DFTk_inp<3,M2,VType,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[2],wi[2],t;

      // W = (wpr1, wpi1)
//       t = Sin<N,1,LocalVType>::value();
//       const LocalVType wpr1 = 1 - 2.0*t*t;
//       const LocalVType wpi1 = -S*Sin<N,2,LocalVType>::value();
      const LocalVType wpr1 = WR::value();
      const LocalVType wpi1 = WI::value();
      
      // W^2 = (wpr2, wpi2)
      const LocalVType wpr2 = wpr1*wpr1 - wpi1*wpi1;
      const LocalVType wpi2 = 2*wpr1*wpi1;
      
      wr[0] = wpr1;
      wi[0] = wpi1;
      wr[1] = wpr2;
      wi[1] = wpi2;
      for (int_t i=2; i<M2; i+=2) {
	spec_inp.apply(data+i, wr, wi);

        t = wr[0];
        wr[0] = t*wpr1 - wi[0]*wpi1;
        wi[0] = wi[0]*wpr1 + t*wpi1;
        t = wr[1];
        wr[1] = t*wpr2 - wi[1]*wpi2;
        wi[1] = wi[1]*wpr2 + t*wpi2;
      }
   }
};

// Specialization for radix 2
template<int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<2,M,Step,VType,S,W,false,true> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   typedef Compute<typename W::Re,VType::Accuracy> WR;
   typedef Compute<typename W::Im,VType::Accuracy> WI;
   static const int_t N = 2*M;
   static const int_t S2 = 2*Step;
   DFTk_inp<2,N,VType,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr,wi,t;
//       t = Sin<N,1,LocalVType>::value();
//       const LocalVType wpr = 1-2.0*t*t;
//       const LocalVType wpi = -S*Sin<N,2,LocalVType>::value();
      const LocalVType wpr = WR::value();
      const LocalVType wpi = WI::value();
      wr = wpr;
      wi = wpi;
      for (int_t i=S2; i<N; i+=S2) {
	spec_inp.apply(data+i, &wr, &wi);

        t = wr;
        wr = wr*wpr - wi*wpi;
        wi = wi*wpr + t*wpi;
      }
   }
};

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
// template<int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1>
// class InTime;

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
   InTime<M,NFactNext,VType,S,WK,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T<K,M,1,VType,S,W1,false> dft_scaled;
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
   InTimeOOP<M,NFactNext,VType,S,WK,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T<K,M,1,VType,S,W1,false> dft_scaled;
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

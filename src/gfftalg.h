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

#ifndef __gfftalg_h
#define __gfftalg_h

/** \file
    \brief Recursive algorithms and short-radix FFT specifications
*/

#include "gfftspec.h"
#include "gfftfactor.h"
#include "gfftswap.h"

#include "metacomplex.h"
#include "metaroot.h"

namespace GFFT {

using namespace MF;


template<bool C>
struct StaticAssert;

template<>
struct StaticAssert<true> { };

#define GFFT_STATIC_ASSERT(c) StaticAssert<c> static__assert; 


static const int_t StaticLoopLimit = 8;


template<int_t K, int_t M, typename T, int S, class W1, int NIter = 1, class W = W1>
class IterateInTime
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t M2 = M*2;
   static const int_t N = K*M;

   typedef typename GetNextRoot<NIter+1,N,W1,W,2>::Result Wnext;
   IterateInTime<K,M,T,S,W1,NIter+1,Wnext> next;
   DFTk_inp<K,M2,T,S> spec_inp;
   
public:
   void apply(T* data) 
   {
      const LocalVType wr = WR::value();
      const LocalVType wi = WI::value();
//std::cout << NIter-1 << "/" << N << ": (" << wr << ", " << wi << ")" << std::endl;

      spec_inp.apply(data + (NIter-1)*2, &wr, &wi);

      next.apply(data);
   }
};

// Last step of the loop
template<int_t K, int_t M, typename T, int S, class W1, class W>
class IterateInTime<K,M,T,S,W1,M,W> 
{
//    typedef typename RList::Head H;
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t M2 = M*2;
   static const int_t N = K*M;
   DFTk_inp<K,M2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      const LocalVType wr = WR::value();
      const LocalVType wi = WI::value();

      spec_inp.apply(data + (M-1)*2, &wr, &wi);
   }
};

// First step in the loop
template<int_t K, int_t M, typename T, int S, class W1, class W>
class IterateInTime<K,M,T,S,W1,1,W> {
   static const int_t M2 = M*2;
   DFTk_inp<K,M2,T,S> spec_inp;
   IterateInTime<K,M,T,S,W1,2,W> next;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);
      next.apply(data);
   }
};


/////////////////////////////////////////////////////////

template<int_t K, int_t M, typename T, int S, class W, bool doStaticLoop>
class DFTk_x_Im_T;

template<int_t K, int_t M, typename T, int S, class W>
class DFTk_x_Im_T<K,M,T,S,W,true>
{
   IterateInTime<K,M,T,S,W> iterate;
public:
   void apply(T* data) 
   {
      iterate.apply(data);
   }
};

template<int_t K, int_t M, typename T, int S, class W>
class DFTk_x_Im_T<K,M,T,S,W,false>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t N = K*M;
   static const int_t M2 = M*2;
   DFTk_inp<K,M2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[K-1], wi[K-1], wpr[K-1], wpi[K-1], t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr[0], wpi[0])
      wpr[0] = WR::value();
      wpi[0] = WI::value();
//       wpr[0] = 1 - 2.0*t*t;
//       wpi[0] = -S*Sin<N,2,LocalVType>::value();
      
      // W^i = (wpr2, wpi2)
      for (int_t i=0; i<K-2; ++i) {
	wpr[i+1] = wpr[i]*wpr[0] - wpi[i]*wpi[0];
	wpi[i+1] = wpr[i]*wpi[0] + wpr[0]*wpi[i];
      }
      
      for (int_t i=0; i<K-1; ++i) {
	wr[i] = wpr[i];
	wi[i] = wpi[i];
      }
      
      for (int_t i=2; i<M2; i+=2) {
	spec_inp.apply(data+i, wr, wi);

	for (int_t i=0; i<K-1; ++i) {
	  t = wr[i];
	  wr[i] = t*wpr[i] - wi[i]*wpi[i];
	  wi[i] = wi[i]*wpr[i] + t*wpi[i];
	}
      }
   }
  
};

template<int_t M, typename T, int S, class W>
class DFTk_x_Im_T<3,M,T,S,W,false> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t N = 3*M;
   static const int_t M2 = M*2;
   DFTk_inp<3,M2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[2],wi[2],t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr1, wpi1)
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

template<int_t M, typename T, int S, class W>
class DFTk_x_Im_T<2,M,T,S,W,false> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t N = 2*M;
   DFTk_inp<2,N,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr,wi,t;
      t = Sin<N,1,LocalVType>::value();
//       const LocalVType wpr = 1-2.0*t*t;
//       const LocalVType wpi = -S*Sin<N,2,LocalVType>::value();
// std::cout << wpr << "  " << wpi << std::endl;
      const LocalVType wpr = WR::value();
      const LocalVType wpi = WI::value();
// std::cout << wpr1 << "  " << wpi1 << std::endl;
      wr = wpr;
      wi = wpi;
      for (int_t i=2; i<N; i+=2) {
	spec_inp.apply(data+i, &wr, &wi);

        t = wr;
        wr = wr*wpr - wi*wpi;
        wi = wi*wpr + t*wpi;
      }
   }
};

/// Danielson-Lanczos section of the decimation-in-time FFT version
/**
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward

This template class implements resursive metaprogram, which
runs funciton apply() twice recursively at the beginning of the function apply()
with the half of the transform length N
until the simplest case N=2 has been reached. Then function \a _spec2 is called.
Therefore, it has two specializations for N=2 and N=1 (the trivial and empty case).
\sa InFreq
*/
template<int_t N, typename NFact, typename T, int S, class W1, int_t LastK = 1>
class InTime;

// template<int_t N, typename Head, typename Tail, typename T, int S, class W1, int_t LastK>
// class InTime<N, Loki::Typelist<Head,Tail>, T, S, W1, LastK>
// {
//   // Not implemented, because not allowed
//   // Transforms in-place are allowed for powers of primes only!!!
// };

template<int_t N, typename Head, typename T, int S, class W1, int_t LastK>
class InTime<N, Loki::Typelist<Head,Loki::NullType>, T, S, W1, LastK>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   static const int_t M2 = M*2;
   static const int_t N2 = N*2;
   static const int_t LastK2 = LastK*2;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Loki::NullType> NFactNext;
   InTime<M,NFactNext,T,S,WK,K*LastK> dft_str;
   DFTk_x_Im_T<K,M,T,S,W1,(N<=StaticLoopLimit)> dft_scaled;
//   DFTk_x_Im_T<K,M,T,S,W1,false> dft_scaled;
public:
   void apply(T* data) 
   {
      for (int_t m=0; m < N2; m+=M2) 
	dft_str.apply(data + m);

      dft_scaled.apply(data);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename T, int S, class W1, int_t LastK>
class InTime<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, T, S, W1, LastK>
: public InTime<N, Tail, T, S, W1, LastK> {};


// Specialization for prime N
template<int_t N, typename T, int S, class W1, int_t LastK>
class InTime<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,T,S,W1,LastK> {
  DFTk_inp<N, 2, T, S> spec_inp;
public:
  void apply(T* data) 
  { 
    spec_inp.apply(data);
  }
};


template<int_t N, typename NFact, typename T, int S, class W1, int_t LastK = 1>
class InTimeOOP;

template<int_t N, typename Head, typename Tail, typename T, int S, class W1, int_t LastK>
class InTimeOOP<N, Loki::Typelist<Head,Tail>, T, S, W1, LastK>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   static const int_t M2 = M*2;
   static const int_t N2 = N*2;
   static const int_t LastK2 = LastK*2;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   InTimeOOP<M,NFactNext,T,S,WK,K*LastK> dft_str;
   DFTk_x_Im_T<K,M,T,S,W1,(N<=StaticLoopLimit)> dft_scaled;
//   DFTk_x_Im_T<K,M,T,S,W1,false> dft_scaled;
public:

   void apply(const T* src, T* dst) 
   {
      int_t lk = 0;
      for (int_t m = 0; m < N2; m+=M2, lk+=LastK2)
        dft_str.apply(src + lk, dst + m);

      dft_scaled.apply(dst);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename T, int S, class W1, int_t LastK>
class InTimeOOP<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, T, S, W1, LastK>
: public InTimeOOP<N, Tail, T, S, W1, LastK> {};


// Specialization for prime N
template<int_t N, typename T, int S, class W1, int_t LastK>
class InTimeOOP<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,T,S,W1,LastK> {
  DFTk<N, LastK*2, 2, T, S> spec;
public:
  void apply(const T* src, T* dst) 
  { 
    spec.apply(src, dst);
  }
};

  
}  //namespace DFT

#endif /*__gfftalg_h*/

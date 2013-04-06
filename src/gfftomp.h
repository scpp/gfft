/***************************************************************************
 *   Copyright (C) 2009-2013 by Volodymyr Myrnyy                           *
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

#ifndef __gfftomp_h
#define __gfftomp_h

/** \file
    \brief Parallelization template classes using %OpenMP standard
*/

#include "gfftalg.h"

#include <omp.h>

namespace GFFT {

/** \var static const uint SwitchToOMP
This static constant defines FFT length for that 
OpenMP parallelization is switched on. The overhead is
too large for transforms with smaller sizes.
*/
static const unsigned int SwitchToOMP = (1<<12);

/** \class {GFFT::InTimeOMP}
\brief %OpenMP parallelized Danielson-Lanczos section of the decimation-in-time FFT version.
\tparam NThreads is number of threads
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam C condition to ensure that (N>NThreads) and (N>=SwitchToOMP), otherwise 
        parallelization is meaningless and sequential implementation InTime
        is inherited.

Comparing to sequential implementation in template class InTime, this class
runs apply() function of both instances of the half length (N/2) in the separated
threads and so on until NThreads has become equal 1. Then the sequential version
in template class InTime is inherited.
\sa InFreqOMP, InTime, InFreq
*/
template<short_t NThreads, int_t N, typename NFact, typename T, int S, int_t LastK = 1, 
bool C = ((N>NThreads) && (N>=SwitchToOMP))>
class InTimeOMP;

template<short_t NThreads, int_t N, typename Head, typename Tail, typename T, int S, int_t LastK>
class InTimeOMP<NThreads,N,Loki::Typelist<Head,Tail>,T,S,LastK,true> 
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   static const int_t M2 = M*2;
   static const int_t N2 = N*2;
   static const int_t LastK2 = LastK*2;
   static const short_t NThreadsNext = (NThreads > K) ? NThreads/K : 1;
   
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   //IterateInTime<M,T,LastK,K-1> iter;
   InTimeOMP<NThreadsNext,M,NFactNext,T,S,K*LastK> dft_str;
   DFTk_x_Im_T<K,M,T,S> dft_scaled;
public:
   void apply(T* data) 
   {
      #pragma omp parallel shared(data) 
      {
	#pragma omp for schedule(static)
	for (int_t m = 0; m < N2; m+=M2)
	  dft_str.apply(data + m);
      }
      
      dft_scaled.apply(data);
   }

   void apply(const T* src, T* dst, T* buf) 
   {
      #pragma omp parallel shared(src,dst,buf) 
      {
	int_t lk = 0;
	#pragma omp for schedule(static) private(lk)
	for (int_t m = 0; m < N2; m+=M2) {
	  dft_str.apply(src + lk, dst + m, buf);
//        iter.apply(step, src, dst);
	  lk+=LastK2;
	}
      }
      
      dft_scaled.apply(dst);
   }
};

template<int_t N, typename NFact, typename T, int S, int_t LastK>
class InTimeOMP<1,N,NFact,T,S,LastK,true> : public InTime<N,NFact,T,S,LastK> { };

template<short_t NThreads, int_t N, typename NFact, typename T, int S, int_t LastK>
class InTimeOMP<NThreads,N,NFact,T,S,LastK,false> : public InTime<N,NFact,T,S,LastK> { };


/** \class {GFFT::InFreqOMP}
\brief %OpenMP parallelized Danielson-Lanczos section of the decimation-in-time FFT version.
\tparam NThreads is number of threads
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam C condition to ensure that (N>NThreads) and (N>=SwitchToOMP), otherwise 
        parallelization is meaningless and sequential implementation InFreq
        is inherited.

Comparing to sequential implementation in template class InFreq, this class
runs apply() function of both instances of the half length (N/2) in the separated
threads and so on until NThreads has become equal 1. Then the sequential version
in template class InTime is inherited.
\sa InFreqOMP, InTime, InFreq
*/
template<short_t NThreads, int_t N, typename NFact, typename T, int S, int_t LastK, 
bool C=((N>NThreads) && (N>=SwitchToOMP))>
class InFreqOMP;

template<unsigned int NThreads, int_t N, typename Head, typename Tail, typename T, int S, int_t LastK>
class InFreqOMP<NThreads,N,Loki::Typelist<Head,Tail>,T,S,LastK,true> 
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   static const int_t M2 = M*2;
   static const int_t N2 = N*2;
   static const int_t LastK2 = LastK*2;
   static const short_t NThreadsNext = (NThreads > K) ? NThreads/K : 1;
   
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   //IterateInFreq<M,T,LastK,K-1> iter;
   InFreqOMP<NThreadsNext,M,NFactNext,T,S,K*LastK> dft_str;
   T_DFTk_x_Im<K,M,T,S> dft_scaled;

public:
   void apply(T* data) 
   {
      dft_scaled.apply(data);

      #pragma omp parallel shared(data) 
      {
      // K times call to dft_str.apply()
	#pragma omp for schedule(static)
	for (int_t m = 0; m < N2; m+=M2)
	  dft_str.apply(data + m);
      }
   }

   void apply(const T* src, T* dst, T* buf) 
   { 
      dft_scaled.apply(src, buf);

      #pragma omp parallel shared(src,dst,buf) 
      {
	int_t lk = 0;
	#pragma omp for schedule(static) private(lk)
	for (int_t m = 0; m < N2; m+=M2) {
	  dft_str.apply(buf + m, dst + lk, buf + m);
	  lk+=LastK2;
	}
      }
   }
};

template<int_t N, typename NFact, typename T, int S, int_t LastK>
class InFreqOMP<1,N,NFact,T,S,LastK,true> : public InFreq<N,NFact,T,S,LastK> { };

template<short_t NThreads, int_t N, typename NFact, typename T, int S, int_t LastK>
class InFreqOMP<NThreads,N,NFact,T,S,LastK,false> : public InFreq<N,NFact,T,S,LastK> { };




// doesn't work
template<short_t NThreads, int_t N, typename T>
class GFFTswapOMP {
public:
   void apply(T* data) {
     unsigned long i,m,j=0;
     unsigned long n = N;
     #pragma omp parallel firstprivate(data) private(i,j,m)
     {
       #pragma omp for ordered
         for (i=0; i<2*n-1; i+=2) {
           if (j>i) {
             std::swap(data[j], data[i]);
             std::swap(data[j+1], data[i+1]);
           }
           m = n;
           while (m>=2 && j>=m) {
             j -= m;
             m >>= 1;
           }
           j += m;
         }
     }
   }
};

/** \class {GFFT::GFFTswap2OMP}
\brief Binary reordering parallelized by %OpenMP
\tparam NThreads is number of threads
\tparam P current transform length as power of two
\tparam T value type of the data array
\tparam I is internal counter
\tparam C condition to ensure that (N>NThreads) and (N>=SwitchToOMP), otherwise 
        parallelization is meaningless and the sequential implementation GFFTswap2
        is inherited.
*/
template<short_t NThreads, int_t P, typename T,
unsigned int I=0, bool C=(((1<<P)>NThreads) && ((1<<P)>=SwitchToOMP))>
class GFFTswap2OMP;

template<short_t NThreads, int_t P, typename T, int_t I>
class GFFTswap2OMP<NThreads,P,T,I,true> {
   static const int_t BN = 1<<(I+1);
   static const int_t BR = 1<<(P-I);
   GFFTswap2OMP<NThreads/2,P,T,I+1> next;
public:
   void apply(T* data, const int_t n=0, const int_t r=0) {
     #pragma omp parallel shared(data)
     {
       #pragma omp sections
       {
         #pragma omp section
         next.apply(data, n, r);

         #pragma omp section
         next.apply(data, n|BN, r|BR);
       }
     }
   }
};

template<short_t NThreads, int_t P, typename T>
class GFFTswap2OMP<NThreads,P,T,P,true> {
public:
   void apply(T* data, const int_t n, const int_t r) {
      if (n>r) {
        swap(data[n],data[r]);
        swap(data[n+1],data[r+1]);
      }
   }
};

template<int_t P, typename T, int_t I>
class GFFTswap2OMP<1,P,T,I,true> : public GFFTswap2<P,T,I> { };

template<int_t P, typename T>
class GFFTswap2OMP<1,P,T,P,true> : public GFFTswap2<P,T,P> { };

template<short_t NThreads, int_t P, typename T, int_t I>
class GFFTswap2OMP<NThreads,P,T,I,false> : public GFFTswap2<P,T,I> { };



} //namespace

#endif

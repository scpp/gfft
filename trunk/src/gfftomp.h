/***************************************************************************
 *   Copyright (C) 2009-2014 by Vladimir Mirnyy                            *
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
//#include "gfftstdalg.h"
#include "gfftalgfreq.h"
#include "gfftswap.h"

#include <omp.h>

namespace GFFT {

/** \var static const uint SwitchToOMP
This static constant defines FFT length for that 
OpenMP parallelization is switched on. The overhead is
too large for transforms with smaller sizes.
*/
//static const int_t SwitchToOMP = (1<<6);
static const int_t SwitchToOMP = (1<<1);


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
template<int_t NThreads, int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1, 
bool C = ((N>NThreads) && (N>=SwitchToOMP))>
class InTimeOMP;

template<int_t NThreads, int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTimeOMP<NThreads,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK,true> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;

   static const int_t NThreadsCreate = (NThreads > K) ? K : NThreads;
   static const int_t NThreadsNext = (NThreads != NThreadsCreate) ? NThreads-NThreadsCreate : 1;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   InTimeOMP<NThreadsNext,M,NFactNext,VType,S,WK,K*LastK> dft_str;
   DFTk_x_Im_T<K,M,1,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
public:
   void apply(T* data) 
   {
      #pragma omp parallel for shared(data) schedule(static) num_threads(NThreadsCreate)
      for (int_t m = 0; m < N2; m+=M2)
	dft_str.apply(data + m);
      
      dft_scaled.apply(data);
   }
};

template<int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTimeOMP<1,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK,true> 
: public InTime<N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> { };

template<int_t NThreads, int_t N, typename NFact, typename VType, int S, class W1, int_t LastK>
class InTimeOMP<NThreads,N,NFact,VType,S,W1,LastK,false> : public InTime<N,NFact,VType,S,W1,LastK> { };

///////////////////////

// Assume: K >= NThreadsCreate
template<int_t N2, int_t M2, int_t LastK2, int_t NThreads, int K, 
int I = 0, bool C = (K>NThreads)>
struct ParallLoopOOP;

template<int_t N2, int_t M2, int_t LastK2, int_t NThreads, int K, int I>
struct ParallLoopOOP<N2,M2,LastK2,NThreads,K,I,true>
{
  static const int_t NThreadsAlreadyCreated = I*NThreads;
  
  ParallLoopOOP<N2,M2,LastK2,NThreads,K-NThreads,I+1> NextStep;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, const T* src, T* dst)
  {
      #pragma omp parallel num_threads(NThreads)
      {
	int tid = omp_get_thread_num() + NThreadsAlreadyCreated;
	dft_str.apply(src + tid*LastK2, dst + tid*M2);
	//#pragma omp barrier
      }
      
      NextStep.apply(dft_str, src, dst);
  }
};

template<int_t N2, int_t M2, int_t LastK2, int_t NThreads, int K, int I>
struct ParallLoopOOP<N2,M2,LastK2,NThreads,K,I,false>
{
  static const int_t NThreadsAlreadyCreated = I*NThreads;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, const T* src, T* dst)
  {
      #pragma omp parallel num_threads(K)
      {
	int tid = omp_get_thread_num() + NThreadsAlreadyCreated;
	dft_str.apply(src + tid*LastK2, dst + tid*M2);
	#pragma omp barrier
      }
  }
};
/*
template<int_t N2, int_t LastK2>
struct ParallSpecOOP<N2,2,LastK2,2>
{
  static const int_t M2 = N2/2;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, const T* src, T* dst)
  {
      // Two times call to dft_str.apply()
      #pragma omp parallel sections
      {
	#pragma omp section
	dft_str.apply(src, dst);

	#pragma omp section
	dft_str.apply(src + LastK2, dst + M2);
      }
  }
};

template<int_t N2, int_t LastK2>
struct ParallSpecOOP<N2,4,LastK2,4>
{
  static const int_t M2 = N2/4;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, const T* src, T* dst)
  {
      // Two times call to dft_str.apply()
      #pragma omp parallel sections
      {
	#pragma omp section
	dft_str.apply(src, dst);

	#pragma omp section
	dft_str.apply(src + LastK2, dst + M2);

	#pragma omp section
	dft_str.apply(src + 2*LastK2, dst + 2*M2);

	#pragma omp section
	dft_str.apply(src + 3*LastK2, dst + 3*M2);
      }
  }
};
*/

template<int_t NThreads, int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1, 
bool C = ((N>NThreads) && (N>=SwitchToOMP))>
class InTimeOOP_omp;

template<int_t NThreads, int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTimeOOP_omp<NThreads,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK,true> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;
   static const int_t LastK2 = LastK*C;
   static const int_t NThreadsCreate = (NThreads > K) ? K : NThreads;
   static const int_t NThreadsNext = (NThreads != NThreadsCreate) ? NThreads-NThreadsCreate : 1;
   
//   typedef typename Factorization<SIntID<K>, SInt>::Result KFact;

   typedef typename IPowBig<W1,K>::Result WK;
//   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
//    InTimeOOP_omp<NThreadsNext,M,NFactNext,VType,S,WK,K*LastK> dft_str;
   InTimeOOP<M,Tail,VType,S,WK,K*LastK> dft_str;
//    DFTk_x_Im_T<K,M,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T<K,M,1,VType,S,W1,false> dft_scaled;

   ParallLoopOOP<N2,M2,LastK2,NThreadsCreate,K> parall;
public:

   void apply(const T* src, T* dst) 
   {
      parall.apply(dft_str, src, dst);
      // K times call to dft_str.apply()
      //int_t m, lk;
      //#pragma omp parallel for shared(src,dst) private(m,lk) schedule(static) num_threads(NThreadsCreate)
//       for (m = lk = 0; m < N2; m+=M2) {
// 	dft_str.apply(src + lk, dst + m);
// 	lk += LastK2;
//       }
      
      dft_scaled.apply(dst);
   }
};

template<int_t N, typename Head, typename Tail, typename T, int S, class W1, int_t LastK>
class InTimeOOP_omp<1,N,Loki::Typelist<Head,Tail>,T,S,W1,LastK,true> 
: public InTimeOOP<N,Loki::Typelist<Head,Tail>,T,S,W1,LastK> { };

template<int_t NThreads, int_t N, typename NFact, typename T, int S, class W1, int_t LastK>
class InTimeOOP_omp<NThreads,N,NFact,T,S,W1,LastK,false> : public InTimeOOP<N,NFact,T,S,W1,LastK> { };


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
template<short_t NThreads, int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1, 
bool C=((N>NThreads) && (N>=SwitchToOMP))>
class InFreqOMP;

template<unsigned int NThreads, int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InFreqOMP<NThreads,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK,true> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;

   static const short_t NThreadsCreate = (NThreads > K) ? K : NThreads;
   static const short_t NThreadsNext = (NThreads != NThreadsCreate) ? NThreads-NThreadsCreate : 1;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   InFreqOMP<NThreadsNext,M,NFactNext,VType,S,WK,K*LastK> dft_str;
   T_DFTk_x_Im<K,M,VType,S,W1,true> dft_scaled;

public:
   void apply(T* data) 
   {
      dft_scaled.apply(data);

      // K times call to dft_str.apply()
      #pragma omp parallel for shared(data) schedule(static) num_threads(NThreadsCreate)
      for (int_t m = 0; m < N2; m+=M2)
	dft_str.apply(data + m);
   }
};

template<int_t N, typename Head, typename Tail, typename T, int S, class W1, int_t LastK>
class InFreqOMP<1,N,Loki::Typelist<Head,Tail>,T,S,W1,LastK,true> 
: public InFreq<N,Loki::Typelist<Head,Tail>,T,S,W1,LastK> { };

template<short_t NThreads, int_t N, typename NFact, typename T, int S, class W1, int_t LastK>
class InFreqOMP<NThreads,N,NFact,T,S,W1,LastK,false> : public InFreq<N,NFact,T,S,W1,LastK> { };



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
/*
template<short_t NThreads, uint_t M, uint_t P, typename T,
unsigned int I=0, bool C=(((1<<P)>NThreads) && ((1<<P)>=SwitchToOMP))>
class GFFTswap2OMP;

template<short_t NThreads, uint_t P, typename T, int_t I>
class GFFTswap2OMP<NThreads,2,P,T,I,true> {
   static const int_t BN = 1<<(I+1);
   static const int_t BR = 1<<(P-I);
   GFFTswap2OMP<NThreads/2,2,P,T,I+1> next;
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

template<short_t NThreads, uint_t P, typename T>
class GFFTswap2OMP<NThreads,2,P,T,P,true> {
public:
   void apply(T* data, const int_t n, const int_t r) {
      if (n>r) {
        swap(data[n],data[r]);
        swap(data[n+1],data[r+1]);
      }
   }
};

template<int_t P, typename T, int_t I>
class GFFTswap2OMP<1,2,P,T,I,true> : public GFFTswap2<2,P,T,I> { };

template<int_t P, typename T>
class GFFTswap2OMP<1,2,P,T,P,true> : public GFFTswap2<2,P,T,P> { };

template<short_t NThreads, int_t P, typename T, int_t I>
class GFFTswap2OMP<NThreads,2,P,T,I,false> : public GFFTswap2<2,P,T,I> { };


template<unsigned int NThreads, uint_t M, uint_t P, typename T,
template<typename> class Complex, unsigned int I>
class GFFTswap2OMP<NThreads,M,P,Complex<T>,I,true> {
   static const int_t BN = 1<<I;
   static const int_t BR = 1<<(P-I-1);
   GFFTswap2OMP<NThreads/2,M,P,Complex<T>,I+1> next;
public:
   void apply(Complex<T>* data, const int_t n=0, const int_t r=0) {
     #pragma omp parallel shared(data)
     {
       #pragma omp sections
       {
         #pragma omp section
         next.apply(data,n,r);

         #pragma omp section
         next.apply(data,n|BN,r|BR);
       }
     }
   }
};

template<unsigned int NThreads, uint_t M, uint_t P, typename T,
template<typename> class Complex>
class GFFTswap2OMP<NThreads,M,P,Complex<T>,true> {
public:
   void apply(Complex<T>* data, const int_t n, const int_t r) {
      if (n>r)
        swap(data[n],data[r]);
   }
};

template<uint_t M, uint_t P, typename T, unsigned int I,
template<typename> class Complex>
class GFFTswap2OMP<1,M,P,Complex<T>,I,true> : public GFFTswap2<M,P,Complex<T>,I> { };

template<uint_t M, uint_t P, typename T,
template<typename> class Complex>
class GFFTswap2OMP<1,M,P,Complex<T>,P,true> : public GFFTswap2<M,P,Complex<T>,P> { };

template<unsigned int NThreads, uint_t M, uint_t P, typename T, unsigned int I,
template<typename> class Complex>
class GFFTswap2OMP<NThreads,M,P,Complex<T>,I,false> : public GFFTswap2<M,P,Complex<T>,I> { };
*/

} //namespace

#endif

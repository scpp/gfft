/***************************************************************************
 *   Copyright (C) 2009-2015 by Vladimir Mirnyy                            *
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
#include "gfftstdalg.h"
#include "gfftalgfreq.h"
#include "gfftswap.h"

#include <omp.h>

namespace GFFT {

/** \var static const uint SwitchToOMP
This static constant defines FFT length for that 
OpenMP parallelization is switched on. The overhead is
too large for transforms with smaller sizes.
*/
//static const long_t SwitchToOMP = (1<<6);



// Assume: K >= NThreads
template<long_t M2, long_t NThreads, int K, int I = 0, bool C = (K>NThreads)>
struct ParallLoop;

template<long_t M2, long_t NThreads, int K, int I>
struct ParallLoop<M2,NThreads,K,I,true>
{
  static const long_t NThreadsAlreadyCreated = I*NThreads;
  
  ParallLoop<M2,NThreads,K-NThreads,I+1> NextStep;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, T* data)
  {
      #pragma omp parallel num_threads(NThreads)
      {
	int tid = omp_get_thread_num() + NThreadsAlreadyCreated;
	dft_str.apply(data + tid*M2);
	//#pragma omp barrier
      }
      
      NextStep.apply(dft_str, data);
  }
};

template<long_t M2, long_t NThreads, int K, int I>
struct ParallLoop<M2,NThreads,K,I,false>
{
  static const long_t NThreadsAlreadyCreated = I*NThreads;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, T* data)
  {
      #pragma omp parallel num_threads(K)
      {
	int tid = omp_get_thread_num() + NThreadsAlreadyCreated;
	dft_str.apply(data + tid*M2);
	#pragma omp barrier
      }
  }
};


template<long_t NThreads, long_t K, typename KFact, long_t M, long_t Step, typename VType, int S, class W1,
long_t SimpleSpec = (M / Step),
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
class DFTk_x_Im_T_omp;


/** \class {GFFT::InTime_omp}
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
template<long_t NThreads, long_t N, typename NFact, typename VType, int S, class W1, long_t LastK = 1>
class InTime_omp;

template<long_t NThreads, long_t N, typename Head, typename Tail, typename VType, int S, class W1, long_t LastK>
class InTime_omp<NThreads,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const long_t K = Head::first::value;
   static const long_t M = N/K;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const long_t M2 = M*C;
   static const long_t N2 = N*C;
   static const long_t NThreadsCreate = (NThreads > K) ? K : NThreads;
   
   typedef typename Factorize<ulong_<K> >::Result KFact;
   
   typedef typename IPowBig<W1,K>::Result WK;

   InTime<M,Tail,VType,S,WK,K*LastK> dft_str;
//    DFTk_x_Im_T<K,KFact,M,1,VType,S,W1> dft_scaled;
   DFTk_x_Im_T_omp<NThreads,K,KFact,M,1,VType,S,W1> dft_scaled;

   ParallLoop<M2,NThreadsCreate,K> parall;
public:
   void apply(T* data) 
   {
      parall.apply(dft_str, data);
//      #pragma omp parallel for shared(data) schedule(static) num_threads(NThreadsCreate)
//       for (long_t m = 0; m < N2; m+=M2)
// 	dft_str.apply(data + m);
      
      dft_scaled.apply(data);
   }
};

template<long_t N, typename Head, typename Tail, typename VType, int S, class W1, long_t LastK>
class InTime_omp<1,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
: public InTime<N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> {};

///////////////////////

// Assume: K >= NThreadsCreate
template<typename Perm, long_t N2, long_t M2, long_t LastK2, long_t NThreads, int K,
int I = 0, bool C = (K>NThreads)>
struct ParallLoopOOP;

template<typename Perm, long_t N2, long_t M2, long_t LastK2, long_t NThreads, int K, int I>
struct ParallLoopOOP<Perm,N2,M2,LastK2,NThreads,K,I,true>
{
  static const long_t NThreadsAlreadyCreated = I*NThreads;
  
  ParallLoopOOP<Perm,N2,M2,LastK2,NThreads,K-NThreads,I+1> NextStep;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, const T* src, T* dst)
  {
      #pragma omp parallel num_threads(NThreads)
      {
	int tid = omp_get_thread_num() + NThreadsAlreadyCreated;
	dft_str.apply(src + Perm::value(tid)*LastK2, dst + tid*M2);
	//#pragma omp barrier
      }
      
      NextStep.apply(dft_str, src, dst);
  }
};

template<typename Perm, long_t N2, long_t M2, long_t LastK2, long_t NThreads, int K, int I>
struct ParallLoopOOP<Perm,N2,M2,LastK2,NThreads,K,I,false>
{
  static const long_t NThreadsAlreadyCreated = I*NThreads;
  
  template<class DftStr, class T>
  void apply(DftStr& dft_str, const T* src, T* dst)
  {
      #pragma omp parallel num_threads(K)
      {
	int tid = omp_get_thread_num() + NThreadsAlreadyCreated;
	dft_str.apply(src + Perm::value(tid)*LastK2, dst + tid*M2);
	#pragma omp barrier
      }
  }
};


template<long_t K, typename KFact>
struct Permutation;
  
template<long_t K, long_t N, long_t P, typename Tail>
struct Permutation<K, Loki::Typelist<pair_<long_<N>,long_<P> >,Tail> >
{
  static const long_t M = K/N;   // K = M*N
  typedef Permutation<K/N, Loki::Typelist<pair_<long_<N>,long_<P-1> >,Tail> > Next;
  static long_t value(const long_t ii)
  {
    assert(ii >= 0 && ii < K);
    return (ii%N)*M + Next::value(ii/N);
  }
};

template<long_t K, long_t N, typename Tail>
struct Permutation<K, Loki::Typelist<pair_<long_<N>,long_<0> >,Tail> >
: public Permutation<K, Tail> {};

// template<>
// struct Permutation<4, Loki::Typelist<Pair<long_<2>,long_<2> >,Loki::NullType> >
// {
//   static long_t value(const long_t ii) { return ii; }
// };

template<long_t K>
struct Permutation<K, Loki::NullType>
{
  static long_t value(const long_t ii) { return ii; }
};


template<long_t K, typename KFact, long_t M, typename VType, int S, typename W1,
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
struct DFTk_inp_adapter;

template<long_t K, typename Head, typename Tail, long_t M, typename VType, int S, typename W1>
struct DFTk_inp_adapter<K, Loki::Typelist<Head, Tail>, M, VType, S, W1, true>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const long_t KF = Head::first::value;
   static const long_t KNext = K/KF;
   
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const long_t M2 = M*C;
   static const long_t MKF = M2*KF;
   //static const long_t N2 = K*M2;
   
   typedef typename IPowBig<W1,KF>::Result WK;
   typedef Loki::Typelist<pair_<typename Head::first, long_<Head::second::value-1> >, Tail> KFactNext;

//    DFTk_inp<KF,M2,VType,S> dft_str;
//    DFTk_x_Im_T_omp<1,KNext,KFactNext,KF*M,M,VType,S,W1> dft_scaled;

   DFTk_inp_adapter<KNext,KFactNext,M,VType,S,WK> dft_str;
   DFTk_x_Im_T<KF,KF*KNext,KNext*M,M,VType,S,W1> dft_scaled;
public:

   void apply(T* data) 
   {
     // run strided DFT recursively KF times
//       for (long_t i=0; i < KNext; ++i)
// 	dft_str.apply(data + i*MKF);
      for (long_t i=0; i < KF; ++i)
	dft_str.apply(data + i*M2*KNext);

      dft_scaled.apply(data);
   }

   template<class LT>
   void apply(T* data, LT* wr, LT* wi) 
   {
      dft_str.apply(data, wr, wi);

//       for (long_t i=1; i < KNext; ++i)
// 	dft_str.apply_m(data + i*MKF, wr+i*KF-1, wi+i*KF-1);
      for (long_t i=1; i < KF; ++i)
	dft_str.apply_m(data + i*M2*KNext, wr+i*KNext-1, wi+i*KNext-1);

      dft_scaled.apply(data);
   }
   template<class LT>
   void apply_m(T* data, const LT* wr, const LT* wi) 
   {
//       for (long_t i=0; i < KNext; ++i)
// 	dft_str.apply_m(data + i*MKF, wr+i*KF, wi+i*KF);
      for (long_t i=0; i < KF; ++i)
	dft_str.apply_m(data + i*M2*KNext, wr+i*KNext, wi+i*KNext);

      dft_scaled.apply(data);
   }
};

template<long_t K, long_t KF, typename Tail, long_t M, typename VType, int S, typename W1>
struct DFTk_inp_adapter<K, Loki::Typelist<pair_<long_<KF>, long_<0> >, Tail>, M, VType, S, W1, true>
: public DFTk_inp_adapter<K, Tail, M, VType, S, W1> { };

// Specialization for prime K
template<long_t K, long_t M, typename VType, int S, class W1>
struct DFTk_inp_adapter<K,Loki::Typelist<pair_<long_<K>, long_<1> >, Loki::NullType>,M,VType,S,W1,true>
: public DFTk_inp<K, M*2, VType, S> { };

// Specialization for K=4
// template<long_t M, typename VType, int S, class W1>
// struct DFTk_inp_adapter<4,Loki::Typelist<Pair<long_<2>, long_<2> >, Loki::NullType>,M,VType,S,W1> 
// : public DFTk_inp<4, M*(Loki::TypeTraits<typename VType::ValueType>::isStdFundamental ? 2 : 1), VType, S> { };

///////////////////////////////////////////////////

template<long_t K, typename Head, typename Tail, long_t M, typename VType, int S, typename W1>
struct DFTk_inp_adapter<K, Loki::Typelist<Head, Tail>, M, VType, S, W1, false>
{
   typedef typename VType::ValueType CT;
   static const long_t KF = Head::first::value;
   static const long_t KNext = K/KF;
   
   static const int C = Loki::TypeTraits<CT>::isStdFundamental ? 2 : 1;
   static const long_t M2 = M*C;
   static const long_t MKF = M2*KF;
   
   typedef typename IPowBig<W1,KF>::Result WK;
   typedef Loki::Typelist<pair_<typename Head::first, long_<Head::second::value-1> >, Tail> KFactNext;

//    DFTk_inp<KF,M2,VType,S> dft_str;
//    DFTk_x_Im_T_omp<1,KNext,KFactNext,KF*M,M,VType,S,W1> dft_scaled;

   DFTk_inp_adapter<KNext,KFactNext,M,VType,S,WK> dft_str;
   DFTk_x_Im_T<KF,KF*KNext,KNext*M,M,VType,S,W1> dft_scaled;
public:

   void apply(CT* data) 
   {
     // run strided DFT recursively KF times
      for (long_t i=0; i < KF; ++i)
	dft_str.apply(data + i*M2*KNext);

      dft_scaled.apply(data);
   }

   void apply(CT* data, const CT* w) 
   {
      dft_str.apply(data, w);

      for (long_t i=1; i < KF; ++i)
	dft_str.apply_m(data + i*M2*KNext, w+i*KNext-1);

      dft_scaled.apply(data);
   }

   void apply_m(CT* data, const CT* w) 
   {
      for (long_t i=0; i < KF; ++i)
	dft_str.apply_m(data + i*M2*KNext, w+i*KNext);

      dft_scaled.apply(data);
   }
};

template<long_t K, long_t KF, typename Tail, long_t M, typename VType, int S, typename W1>
struct DFTk_inp_adapter<K, Loki::Typelist<pair_<long_<KF>, long_<0> >, Tail>, M, VType, S, W1, false>
: public DFTk_inp_adapter<K, Tail, M, VType, S, W1> { };

// Specialization for prime K
template<long_t K, long_t M, typename VType, int S, class W1>
struct DFTk_inp_adapter<K,Loki::Typelist<pair_<long_<K>, long_<1> >, Loki::NullType>,M,VType,S,W1,false>
: public DFTk_inp<K, M, VType, S> { };

///////////////////////////////////////////////////


// General implementation
template<long_t NThreads, long_t K, typename KFact, long_t M, long_t Step, typename VType, int S, class W1, long_t SimpleSpec>
class DFTk_x_Im_T_omp<NThreads,K,KFact,M,Step,VType,S,W1,SimpleSpec,true>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const long_t N = K*M;
   static const long_t M2 = M*2;
   static const long_t S2 = 2*Step;
   
   typedef typename GetFirstRoot<K,S,VType::Accuracy>::Result W;
   DFTk_inp_adapter<K,KFact,M,VType,S,W> spec_inp_a;

//   typedef Permutation<K,typename Loki::TL::Reverse<KFact>::Result> Perm;
   typedef Permutation<K,KFact> Perm;

public:
  /*  // This has low performance benefit and doesn't work for powers other than 2
   void apply(T* data) 
   {
      #pragma omp parallel num_threads(K) shared(data)
      {
	int tid = omp_get_thread_num();// + NThreadsAlreadyCreated;
	if (tid == 0) {
	  spec_inp_a.apply(data);
	  ComputeRoots<K,VType,W1,Perm> roots(NThreads,tid);
	  spec_inp_a.apply(data+S2*NThreads, roots.get_real(), roots.get_imag());
      for (long_t j=2*NThreads*S2; j<M2; j+=S2*NThreads) {
	    roots.step();
	    spec_inp_a.apply(data+j, roots.get_real(), roots.get_imag());
	  }
	}
	else {
	  ComputeRoots<K,VType,W1,Perm> roots(NThreads,tid);
	  spec_inp_a.apply(data+S2*tid, roots.get_real(), roots.get_imag());
      for (long_t j=(NThreads+tid)*S2; j<M2; j+=S2*NThreads) {
	    roots.step();
	    spec_inp_a.apply(data+j, roots.get_real(), roots.get_imag());
	  }
	}
	#pragma omp barrier
      }
   }
   */
  
   // Sequential version
   void apply(T* data) 
   {
     // M times call to spec_inp_a.apply()
      spec_inp_a.apply(data);
      ComputeRoots<K,VType,W1,Perm> roots;

      spec_inp_a.apply(data+S2, roots.get_real(), roots.get_imag());
      for (long_t j=S2+S2; j<M2; j+=S2) {
	roots.step();
	spec_inp_a.apply(data+j, roots.get_real(), roots.get_imag());
      }
   }
    
};

/*
template<long_t NThreads, typename KFact, long_t M, long_t Step, typename VType, int S, class W>
class DFTk_x_Im_T_omp<NThreads,3,KFact,M,Step,VType,S,W,false,true>
: public DFTk_x_Im_T<3,KFact,M,Step,VType,S,W,false,true> {};

template<long_t NThreads, typename KFact, long_t M, long_t Step, typename VType, int S, class W>
class DFTk_x_Im_T_omp<NThreads,2,KFact,M,Step,VType,S,W,false,true>
: public DFTk_x_Im_T<2,KFact,M,Step,VType,S,W,false,true> {};
*/

template<long_t NThreads, long_t K, typename KFact, long_t M, long_t Step, typename VType, int S, class W1, long_t SimpleSpec>
class DFTk_x_Im_T_omp<NThreads,K,KFact,M,Step,VType,S,W1,SimpleSpec,false>
{
   typedef typename VType::ValueType CT;
   static const long_t N = K*M;
   
   typedef typename GetFirstRoot<K,S,VType::Accuracy>::Result W;
   DFTk_inp_adapter<K,KFact,M,VType,S,W> spec_inp_a;

   typedef Permutation<K,KFact> Perm;

public:
  
   /*  // This has low performance benefit and doesn't work for powers other than 2
  void apply(CT* data) 
   {
      #pragma omp parallel num_threads(K) shared(data)
      {
	int tid = omp_get_thread_num();
	if (tid == 0) {
	  spec_inp_a.apply(data);
	  ComputeRootsStd<K,VType,W1,Perm> roots(NThreads,tid);
	  spec_inp_a.apply(data+Step*NThreads, roots.get());
      for (long_t j=2*NThreads*Step; j<M; j+=Step*NThreads) {
	    roots.step();
	    spec_inp_a.apply(data+j, roots.get());
	  }
	}
	else {
	  ComputeRootsStd<K,VType,W1,Perm> roots(NThreads,tid);
	  spec_inp_a.apply(data+Step*tid, roots.get());
      for (long_t j=(NThreads+tid)*Step; j<M; j+=Step*NThreads) {
	    roots.step();
	    spec_inp_a.apply(data+j, roots.get());
	  }
	}
	#pragma omp barrier
      }
   }
   */
  
   // Sequential version
   void apply(CT* data) 
   {
     // M times call to spec_inp_a.apply()
      spec_inp_a.apply(data);
      ComputeRootsStd<K,VType,W1,Perm> roots;

      spec_inp_a.apply(data+Step, roots.get());
      for (long_t j=Step+Step; j<M; j+=Step) {
	roots.step();
	spec_inp_a.apply(data+j, roots.get());
      }
   }
    
};

///////////////////////////////////////////////////////////

template<long_t NThreads, long_t N, typename NFact, typename VType, int S, class W1, long_t LastK = 1>
class InTimeOOP_omp;

template<long_t NThreads, long_t N, typename Head, typename Tail, typename VType, int S, class W1, long_t LastK>
class InTimeOOP_omp<NThreads,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const long_t K = Head::first::value;
   static const long_t M = N/K;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const long_t M2 = M*C;
   static const long_t N2 = N*C;
   static const long_t LastK2 = LastK*C;
   static const long_t NThreadsCreate = (NThreads > K) ? K : NThreads;
   
   typedef typename Factorize<ulong_<K> >::Result KFact;
   //typedef Permutation<K,KFact> Perm;
   typedef Permutation<K,typename Loki::TL::Reverse<KFact>::Result> Perm;

   typedef typename IPowBig<W1,K>::Result WK;
   InTimeOOP<M,Tail,VType,S,WK,K*LastK> dft_str;
//    DFTk_x_Im_T<K,M,VType,S,W1> dft_scaled;
   DFTk_x_Im_T_omp<NThreadsCreate,K,KFact,M,1,VType,S,W1> dft_scaled;

   ParallLoopOOP<Perm,N2,M2,LastK2,NThreadsCreate,K> parall;
public:

   void apply(const T* src, T* dst) 
   {
      parall.apply(dft_str, src, dst);

//       #pragma omp parallel for shared(src,dst) schedule(static) num_threads(NThreadsCreate)
//       for (long_t i = 0; i < K; ++i) {
// 	long_t ii = Perm::value(i);
// 	dft_str.apply(src + ii*LastK2, dst + i*M2);
//       }
      
      dft_scaled.apply(dst);
   }
};

template<long_t N, typename Head, typename Tail, typename VType, int S, class W1, long_t LastK>
class InTimeOOP_omp<1,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
: public InTimeOOP<N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> {};


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
template<short_t NThreads, ulong_t M, ulong_t P, typename T,
unsigned int I=0, bool C=(((1<<P)>NThreads) && ((1<<P)>=SwitchToOMP))>
class GFFTswap2OMP;

template<short_t NThreads, ulong_t P, typename T, long_t I>
class GFFTswap2OMP<NThreads,2,P,T,I,true> {
   static const long_t BN = 1<<(I+1);
   static const long_t BR = 1<<(P-I);
   GFFTswap2OMP<NThreads/2,2,P,T,I+1> next;
public:
   void apply(T* data, const long_t n=0, const long_t r=0) {
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

template<short_t NThreads, ulong_t P, typename T>
class GFFTswap2OMP<NThreads,2,P,T,P,true> {
public:
   void apply(T* data, const long_t n, const long_t r) {
      if (n>r) {
        swap(data[n],data[r]);
        swap(data[n+1],data[r+1]);
      }
   }
};

template<long_t P, typename T, long_t I>
class GFFTswap2OMP<1,2,P,T,I,true> : public GFFTswap2<2,P,T,I> { };

template<long_t P, typename T>
class GFFTswap2OMP<1,2,P,T,P,true> : public GFFTswap2<2,P,T,P> { };

template<short_t NThreads, long_t P, typename T, long_t I>
class GFFTswap2OMP<NThreads,2,P,T,I,false> : public GFFTswap2<2,P,T,I> { };


template<unsigned int NThreads, ulong_t M, ulong_t P, typename T,
template<typename> class Complex, unsigned int I>
class GFFTswap2OMP<NThreads,M,P,Complex<T>,I,true> {
   static const long_t BN = 1<<I;
   static const long_t BR = 1<<(P-I-1);
   GFFTswap2OMP<NThreads/2,M,P,Complex<T>,I+1> next;
public:
   void apply(Complex<T>* data, const long_t n=0, const long_t r=0) {
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

template<unsigned int NThreads, ulong_t M, ulong_t P, typename T,
template<typename> class Complex>
class GFFTswap2OMP<NThreads,M,P,Complex<T>,true> {
public:
   void apply(Complex<T>* data, const long_t n, const long_t r) {
      if (n>r)
        swap(data[n],data[r]);
   }
};

template<ulong_t M, ulong_t P, typename T, unsigned int I,
template<typename> class Complex>
class GFFTswap2OMP<1,M,P,Complex<T>,I,true> : public GFFTswap2<M,P,Complex<T>,I> { };

template<ulong_t M, ulong_t P, typename T,
template<typename> class Complex>
class GFFTswap2OMP<1,M,P,Complex<T>,P,true> : public GFFTswap2<M,P,Complex<T>,P> { };

template<unsigned int NThreads, ulong_t M, ulong_t P, typename T, unsigned int I,
template<typename> class Complex>
class GFFTswap2OMP<NThreads,M,P,Complex<T>,I,false> : public GFFTswap2<M,P,Complex<T>,I> { };
*/

} //namespace

#endif

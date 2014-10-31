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
//static const int_t SwitchToOMP = (1<<6);



// Assume: K >= NThreads
template<int_t M2, int_t NThreads, int K, int I = 0, bool C = (K>NThreads)>
struct ParallLoop;

template<int_t M2, int_t NThreads, int K, int I>
struct ParallLoop<M2,NThreads,K,I,true>
{
  static const int_t NThreadsAlreadyCreated = I*NThreads;
  
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

template<int_t M2, int_t NThreads, int K, int I>
struct ParallLoop<M2,NThreads,K,I,false>
{
  static const int_t NThreadsAlreadyCreated = I*NThreads;
  
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


template<int_t NThreads, int_t K, typename KFact, int_t M, int_t Step, typename VType, int S, class W, bool doStaticLoop,
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
template<int_t NThreads, int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1>
class InTime_omp;

template<int_t NThreads, int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTime_omp<NThreads,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;

   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;
   static const int_t NThreadsCreate = (NThreads > K) ? K : NThreads;
   
   typedef typename Factorization<SIntID<K>, SInt>::Result KFact;
   
   typedef typename IPowBig<W1,K>::Result WK;

   InTime<M,Tail,VType,S,WK,K*LastK> dft_str;
//    DFTk_x_Im_T<K,KFact,M,1,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T_omp<NThreads,K,KFact,M,1,VType,S,W1,false> dft_scaled;

   ParallLoop<M2,NThreadsCreate,K> parall;
public:
   void apply(T* data) 
   {
      parall.apply(dft_str, data);
//      #pragma omp parallel for shared(data) schedule(static) num_threads(NThreadsCreate)
//       for (int_t m = 0; m < N2; m+=M2) 
// 	dft_str.apply(data + m);
      
      dft_scaled.apply(data);
   }
};

template<int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTime_omp<1,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
: public InTime<N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> {};

///////////////////////

// Assume: K >= NThreadsCreate
template<typename Perm, int_t N2, int_t M2, int_t LastK2, int_t NThreads, int K, 
int I = 0, bool C = (K>NThreads)>
struct ParallLoopOOP;

template<typename Perm, int_t N2, int_t M2, int_t LastK2, int_t NThreads, int K, int I>
struct ParallLoopOOP<Perm,N2,M2,LastK2,NThreads,K,I,true>
{
  static const int_t NThreadsAlreadyCreated = I*NThreads;
  
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

template<typename Perm, int_t N2, int_t M2, int_t LastK2, int_t NThreads, int K, int I>
struct ParallLoopOOP<Perm,N2,M2,LastK2,NThreads,K,I,false>
{
  static const int_t NThreadsAlreadyCreated = I*NThreads;
  
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


template<int_t K, typename KFact>
struct Permutation;
  
template<int_t K, int_t N, int_t P, typename Tail>
struct Permutation<K, Loki::Typelist<Pair<SInt<N>,SInt<P> >,Tail> >
{
  static const int_t M = K/N;   // K = M*N
  typedef Permutation<K/N, Loki::Typelist<Pair<SInt<N>,SInt<P-1> >,Tail> > Next;
  static int_t value(const int_t ii)
  {
    assert(ii >= 0 && ii < K);
    return (ii%N)*M + Next::value(ii/N);
  }
};

template<int_t K, int_t N, typename Tail>
struct Permutation<K, Loki::Typelist<Pair<SInt<N>,SInt<0> >,Tail> >
: public Permutation<K, Tail> {};

// template<>
// struct Permutation<4, Loki::Typelist<Pair<SInt<2>,SInt<2> >,Loki::NullType> >
// {
//   static int_t value(const int_t ii) { return ii; }
// };

template<int_t K>
struct Permutation<K, Loki::NullType>
{
  static int_t value(const int_t ii) { return ii; }
};


template<int_t K, typename KFact, int_t M, typename VType, int S, typename W>
struct DFTk_inp_adapter;

template<int_t K, typename Head, typename Tail, int_t M, typename VType, int S, typename W1>
struct DFTk_inp_adapter<K, Loki::Typelist<Head, Tail>, M, VType, S, W1>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t KF = Head::first::value;
   static const int_t KNext = K/KF;
   
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t MKF = M2*KF;
   //static const int_t N2 = K*M2;
   
   typedef typename IPowBig<W1,KF>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> KFactNext;

//    DFTk_inp<KF,M2,VType,S> dft_str;
//    DFTk_x_Im_T_omp<1,KNext,KFactNext,KF*M,M,VType,S,W1,false> dft_scaled;

   typedef Loki::Typelist<Pair<SInt<KF>,SInt<1> >, Loki::NullType> KF_fact;
   DFTk_inp_adapter<KNext,KFactNext,M,VType,S,WK> dft_str;
   DFTk_x_Im_T<KF,KF_fact,KNext*M,M,VType,S,W1,false> dft_scaled;
public:

   void apply(T* data) 
   {
     // run strided DFT recursively KF times
//       for (int_t i=0; i < KNext; ++i) 
// 	dft_str.apply(data + i*MKF);
      for (int_t i=0; i < KF; ++i) 
	dft_str.apply(data + i*M2*KNext);

      dft_scaled.apply(data);
   }

   template<class LT>
   void apply(T* data, LT* wr, LT* wi) 
   {
      dft_str.apply(data, wr, wi);

//       for (int_t i=1; i < KNext; ++i) 
// 	dft_str.apply_m(data + i*MKF, wr+i*KF-1, wi+i*KF-1);
      for (int_t i=1; i < KF; ++i) 
	dft_str.apply_m(data + i*M2*KNext, wr+i*KNext-1, wi+i*KNext-1);

      dft_scaled.apply(data);
   }
   template<class LT>
   void apply_m(T* data, const LT* wr, const LT* wi) 
   {
//       for (int_t i=0; i < KNext; ++i) 
// 	dft_str.apply_m(data + i*MKF, wr+i*KF, wi+i*KF);
      for (int_t i=0; i < KF; ++i) 
	dft_str.apply_m(data + i*M2*KNext, wr+i*KNext, wi+i*KNext);

      dft_scaled.apply(data);
   }
};

template<int_t K, int_t KF, typename Tail, int_t M, typename VType, int S, typename W1>
struct DFTk_inp_adapter<K, Loki::Typelist<Pair<SInt<KF>, SInt<0> >, Tail>, M, VType, S, W1>
: public DFTk_inp_adapter<K, Tail, M, VType, S, W1> { };

// Specialization for prime K
template<int_t K, int_t M, typename VType, int S, class W1>
struct DFTk_inp_adapter<K,Loki::Typelist<Pair<SInt<K>, SInt<1> >, Loki::NullType>,M,VType,S,W1> 
: public DFTk_inp<K, M*(Loki::TypeTraits<typename VType::ValueType>::isStdFundamental ? 2 : 1), VType, S> { };

// Specialization for K=4
// template<int_t M, typename VType, int S, class W1>
// struct DFTk_inp_adapter<4,Loki::Typelist<Pair<SInt<2>, SInt<2> >, Loki::NullType>,M,VType,S,W1> 
// : public DFTk_inp<4, M*(Loki::TypeTraits<typename VType::ValueType>::isStdFundamental ? 2 : 1), VType, S> { };


// General implementation
template<int_t NThreads, int_t K, typename KFact, int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T_omp<NThreads,K,KFact,M,Step,VType,S,W,false,true>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t N = K*M;
   static const int_t M2 = M*2;
   static const int_t S2 = 2*Step;
   
   typedef typename GetFirstRoot<K,S,VType::Accuracy>::Result W1;
   DFTk_inp_adapter<K,KFact,M,VType,S,W1> spec_inp_a;

//   typedef Permutation<K,typename Loki::TL::Reverse<KFact>::Result> Perm;
   typedef Permutation<K,KFact> Perm;

public:
  /*
   void apply(T* data) 
   {
      #pragma omp parallel num_threads(K) shared(data)
      {
	int tid = omp_get_thread_num();// + NThreadsAlreadyCreated;
	if (tid == 0) {
	  spec_inp_a.apply(data);
	  ComputeRoots<K,VType,W,Perm> roots(NThreads,tid);
	  spec_inp_a.apply(data+S2*NThreads, roots.get_real(), roots.get_imag());
	  for (int_t j=2*NThreads*S2; j<M2; j+=S2*NThreads) {
	    roots.step();
	    spec_inp_a.apply(data+j, roots.get_real(), roots.get_imag());
	  }
	}
	else {
	  ComputeRoots<K,VType,W,Perm> roots(NThreads,tid);
	  spec_inp_a.apply(data+S2*tid, roots.get_real(), roots.get_imag());
	  for (int_t j=(NThreads+tid)*S2; j<M2; j+=S2*NThreads) {
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
      ComputeRoots<K,VType,W,Perm> roots;

      spec_inp_a.apply(data+S2, roots.get_real(), roots.get_imag());
      for (int_t j=S2+S2; j<M2; j+=S2) {
	roots.step();
	spec_inp_a.apply(data+j, roots.get_real(), roots.get_imag());
      }
   }
    
};

/*
template<int_t NThreads, typename KFact, int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T_omp<NThreads,3,KFact,M,Step,VType,S,W,false,true>
: public DFTk_x_Im_T<3,KFact,M,Step,VType,S,W,false,true> {};

template<int_t NThreads, typename KFact, int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T_omp<NThreads,2,KFact,M,Step,VType,S,W,false,true>
: public DFTk_x_Im_T<2,KFact,M,Step,VType,S,W,false,true> {};
*/

template<int_t NThreads, int_t N, typename NFact, typename VType, int S, class W1, int_t LastK = 1>
class InTimeOOP_omp;

template<int_t NThreads, int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTimeOOP_omp<NThreads,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
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
   
   typedef typename Factorization<SIntID<K>, SInt>::Result KFact;
   //typedef Permutation<K,KFact> Perm;
   typedef Permutation<K,typename Loki::TL::Reverse<KFact>::Result> Perm;

   typedef typename IPowBig<W1,K>::Result WK;
   InTimeOOP<M,Tail,VType,S,WK,K*LastK> dft_str;
//    DFTk_x_Im_T<K,M,VType,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T_omp<NThreadsCreate,K,KFact,M,1,VType,S,W1,false> dft_scaled;

   ParallLoopOOP<Perm,N2,M2,LastK2,NThreadsCreate,K> parall;
public:

   void apply(const T* src, T* dst) 
   {
      parall.apply(dft_str, src, dst);

//       #pragma omp parallel for shared(src,dst) schedule(static) num_threads(NThreadsCreate)
//       for (int_t i = 0; i < K; ++i) {
// 	int_t ii = Perm::value(i);
// 	dft_str.apply(src + ii*LastK2, dst + i*M2);
//       }
      
      dft_scaled.apply(dst);
   }
};

template<int_t N, typename Head, typename Tail, typename VType, int S, class W1, int_t LastK>
class InTimeOOP_omp<1,N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> 
: public InTimeOOP<N,Loki::Typelist<Head,Tail>,VType,S,W1,LastK> { };


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
/*
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

*/

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

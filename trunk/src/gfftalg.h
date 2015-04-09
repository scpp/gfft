/***************************************************************************
 *   Copyright (C) 2006-2015 by Vladimir Mirnyy                            *
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

static const int_t StaticLoopLimit = (1<<10);
static const int_t PrecomputeRoots = StaticLoopLimit;


///////////////////////////////////////////////////
// This class works, but takes too much compile-time
template<class DFTk, int_t M, int_t N, typename VType, int S, class W1, class WK = W1, 
int NIter = M-2, int C = N/M>
struct SmartIterate
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t M2 = 2*M;
   static const int_t Step = PrecomputeRoots/M2;
   static const int_t I = Step*(M-NIter) - 2;
   
   //static const int_t N = M2;
//    typedef typename Simplify<SRational<SInt<I>,SInt<N> > >::Result R;
//    static const int_t I1 = R::Numer::value;
//    static const int_t N1 = R::Denom::value;

//    typedef typename CosPiDecimal<I,N,VType::Accuracy>::Result Re;
//    typedef typename SinPiDecimal<I,N,VType::Accuracy>::Result Sin2;
//    typedef typename Loki::Select<(S<0),Sin2,
//           typename Negate<Sin2>::Result>::Result Im;

//   typedef typename GetFirstRoot<M2,S,VType::Accuracy>::Result W;
//    typedef Compute<Re,VType::Accuracy> WR;
//    typedef Compute<Im,VType::Accuracy> WI;
   typedef Compute<typename WK::Re,VType::Accuracy> WR;
   typedef Compute<typename WK::Im,VType::Accuracy> WI;

   typedef typename DFTk::RootsHolder SW;
   
   //typedef typename IPowBig<W1,2>::Result WK;
   typedef typename Mult<W1,WK>::Result WKnext;
   SmartIterate<DFTk,M,N,VType,S,W1,WKnext,NIter-2> next;
   DFTk spec_inp;

   void apply(T* data) 
   {
     T* roots = SW::Instance().getData(); 
     if (roots[I] == uninitialized_flag) {
       roots[I]   = WR::value();
       roots[I+1] = WI::value();
     }
//      LocalVType wr = WR::value();
//      LocalVType wi = WI::value();
//      spec_inp.apply(data+(M-NIter),&wr,&wi);
     spec_inp.apply(data+(M-NIter),roots+I,roots+I+1);

     T wr = -roots[I];
//      wr = -wr;
//      spec_inp.apply(data+M2-(M-NIter),&wr,&wi);
     spec_inp.apply(data+M2-(M-NIter),&wr,roots+I+1);
     
     next.apply(data);
   }
};
/*
template<class DFTk, int_t M, int_t N, typename VType, int S, class W1, class WK, int NIter>
struct SmartIterate<DFTk, M, N, VType, S, W1, WK, NIter, 2>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t M2 = 2*M;
   //static const int_t Step = N/M2/2;
   //static const int_t I = Step*(M-NIter) - 2;
   
   typedef Compute<typename WK::Re,VType::Accuracy> WR;
   typedef Compute<typename WK::Im,VType::Accuracy> WI;
    
   typedef typename DFTk::RootsHolder SW;
   
   //typedef typename IPowBig<W1,2>::Result WK;
   typedef typename Mult<W1,WK>::Result WKnext;
   SmartIterate<DFTk,M,2*M,VType,S,W1,WKnext,NIter-2> next;
   DFTk spec_inp;

   void apply(T* data) 
   {
     T wr = WR::value();
     T wi = WI::value();
     spec_inp.apply(data+(M-NIter),&wr,&wi);

     wr = -wr;
     spec_inp.apply(data+M2-(M-NIter),&wr,&wi);
     
     next.apply(data);
   }
};
*/
template<class DFTk, int_t M, int_t N, typename VType, int S, class W1, class WK, int C>
struct SmartIterate<DFTk,M,N,VType,S,W1,WK,0,C>
: public SmartIterate<DFTk,M,N,VType,S,W1,WK,0,2> {};

template<class DFTk, int_t M, int_t N, typename VType, int S, class W1, class WK>
struct SmartIterate<DFTk,M,N,VType,S,W1,WK,0,2>
{
   typedef typename VType::ValueType T;
   DFTk spec_inp;

   void apply(T* data) 
   {
     spec_inp.apply(data);
     spec_inp.apply_1(data+M);
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
template<int_t K, int_t LastK, int_t M, int_t Step, typename VType, int S, class W1, 
int_t SimpleSpec = (M / Step),
bool isStd = Loki::TypeTraits<typename VType::ValueType>::isStdFundamental>
class DFTk_x_Im_T;

// Rely on the static template loop
// template<int_t K, int_t LastK, int_t M, int_t Step, typename VType, int S, class W1>
// class DFTk_x_Im_T<K,LastK,M,Step,VType,S,W1,true,true> 
// : public DFTk_x_Im_T<K,LastK,M,Step,VType,S,W1,false,true> {};
//: public SmartIterate<DFTk_inp<K,M*2,VType,S>,M,LastK*M,VType,S,W1> {};

// General implementation
template<int_t K, int_t LastK, int_t M, int_t Step, typename VType, int S, class W1, int_t SimpleSpec>
class DFTk_x_Im_T<K,LastK,M,Step,VType,S,W1,SimpleSpec,true>
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

      ComputeRoots<K,VType,W1> roots;

      spec_inp.apply(data+S2, roots.get_real(), roots.get_imag());
      for (int_t j=S2+S2; j<M2; j+=S2) {
	roots.step();
	spec_inp.apply(data+j, roots.get_real(), roots.get_imag());
      }
   }
  
};
/*
// Specialization for radix 2
template<int_t M, int_t LastK, int_t Step, typename VType, int S, class W1, int_t SimpleSpec>
class DFTk_x_Im_T<2,LastK,M,Step,VType,S,W1,SimpleSpec,true> 
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   typedef Compute<typename W1::Re,VType::Accuracy> WR;
   typedef Compute<typename W1::Im,VType::Accuracy> WI;
   static const int_t N = 2*M;
   static const int_t S2 = 2*Step;
   DFTk_inp<2,N,VType,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr,wi,t;
      const LocalVType wpr = WR::value();
      const LocalVType wpi = WI::value();
      wr = wpr;
      wi = wpi;

      spec_inp.apply(data+S2, &wr, &wi);
      for (int_t i=S2+S2; i<N; i+=S2) {
        t = wr;
        wr = wr*wpr - wi*wpi;
        wi = wi*wpr + t*wpi;
	spec_inp.apply(data+i, &wr, &wi);
      }
   }
};
*/

// Specialization for radix 2
template<int_t M, int_t LastK, int_t Step, typename VType, int S, class W1, int_t SimpleSpec>
class DFTk_x_Im_T<2,LastK,M,Step,VType,S,W1,SimpleSpec,true> 
{
   typedef typename VType::ValueType T;
   static const int_t N = 2*M;
   static const int_t S2 = 2*Step;
   //static const int_t NR = (4*PrecomputeRoots > LastK*M) ? LastK*M/4 : PrecomputeRoots;
   static const int_t NR = (PrecomputeRoots>=N) ? 2 : PrecomputeRoots/2;
   static const int_t K = N/PrecomputeRoots;
   static const int_t K2 = 2*K;
   //typedef typename GetFirstRoot<N,S,VType::Accuracy>::Result W1;
   typedef Compute<typename W1::Re,VType::Accuracy> WR;
   typedef Compute<typename W1::Im,VType::Accuracy> WI;
   DFTk_inp<2,N,VType,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);
      if (M%2 == 0) 
	spec_inp.apply_1(data+M);
      
      T wr,wi,t;
      const T wpr = WR::value();
      const T wpi = WI::value();
      wr = wpr;
      wi = wpi;

      spec_inp.apply(data+S2, &wr, &wi);
      t = -wr;
      spec_inp.apply(data+N-S2, &t, &wi);
      for (int_t i=S2+S2; i<M; i+=S2) {
	  t = wr;
	  wr = wr*wpr - wi*wpi;
	  wi = wi*wpr + t*wpi;
	  spec_inp.apply(data+i, &wr, &wi);
	  t = -wr;
	  spec_inp.apply(data+N-i, &t, &wi);
      }
   }  
};

/*
template<int_t M, int_t LastK, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<2,LastK,M,Step,VType,S,W,false,true> 
{
   typedef typename VType::ValueType T;
//    typedef typename VType::TempType LocalVType;
//    typedef typename GetFirstRoot<2*M,S,VType::Accuracy>::Result W1;
//    typedef Compute<typename W1::Re,VType::Accuracy> WR;
//    typedef Compute<typename W1::Im,VType::Accuracy> WI;
   static const int_t N = 2*M;
   static const int_t S2 = 2*Step;
   static const int_t NR = (2*PrecomputeRoots > LastK*M) ? LastK*M/2 : PrecomputeRoots;
   static const int_t Inc = 2*NR/N;
   DFTk_inp<2,N,VType,S,W> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

//       T wr,wi,t;
//       const T wpr = WR::value();
//       const T wpi = WI::value();
//       wr = wpr;
//       wi = wpi;
// 
//       spec_inp.apply(data+S2, &wr, &wi);
//       for (int_t i=S2+S2; i<N; i+=S2) {
//         t = wr;
//         wr = wr*wpr - wi*wpi;
//         wi = wi*wpr + t*wpi;
// 	spec_inp.apply(data+i, &wr, &wi);
//       }
      
      const T* roots = W::Instance().getData(); 
      int_t lk = Inc-S2;
      for (int_t i=S2; i<N; i+=S2) {
	spec_inp.apply(data+i, roots+lk, roots+lk+1);
	lk += Inc;
      }

//       const T* roots = Loki::SingletonHolder<RootsHolder<N,Loki::NullType,VType,S> >::Instance().getData(); 
//       for (int_t i=S2; i<N; i+=S2) {
// 	spec_inp.apply(data+i, roots+i-S2, roots+i-S2+1);
//       }
   }
};
*/
/*
template<int_t M, int_t Step, typename VType, int S, class W>
class DFTk_x_Im_T<2,2,M,Step,VType,S,W,false,true> 
{
   typedef typename VType::ValueType T;
//    typedef typename VType::TempType LocalVType;
   static const int_t N = 2*M;
   static const int_t S2 = 2*Step;
   typedef typename GetFirstRoot<N,S,VType::Accuracy>::Result W1;
   typedef Compute<typename W1::Re,VType::Accuracy> WR;
   typedef Compute<typename W1::Im,VType::Accuracy> WI;
   DFTk_inp<2,N,VType,S,W> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      T wr,wi,t;
      const T wpr = WR::value();
      const T wpi = WI::value();
      wr = wpr;
      wi = wpi;

      spec_inp.apply(data+S2, &wr, &wi);
      
      const T* roots = W::Instance().getData(); 
      int_t lk = 0;
      for (int_t i=S2+S2; i<N; i+=S2) {
	spec_inp.apply(data+i, roots+lk, roots+lk+1);
	i += S2;
        wr = roots[lk]*wpr - roots[lk+1]*wpi;
        wi = roots[lk+1]*wpr + roots[lk]*wpi;
	spec_inp.apply(data+i, &wr, &wi);
	lk += 2;
      }
   }
};
*/
// Specialization for radix 2
template<int_t LastK, int_t M, int_t Step, typename VType, int S, class W1>
class DFTk_x_Im_T<2,LastK,M,Step,VType,S,W1,2,true> 
{
   typedef typename VType::ValueType T;
   static const int_t S2 = 2*Step;
   DFTk_inp<2,2*M,VType,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);
      spec_inp.apply_1(data+S2);
   }
};

/// In-place decimation-in-time FFT version
/**
\tparam N current transform length
\tparam NFact factorization list
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam W singleton for roots of unity

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
   InTime<M,NFactNext,VType,S,WK,K*LastK> dft_str;
   DFTk_x_Im_T<K,K*LastK,M,1,VType,S,W1> dft_scaled;
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
\tparam W singleton for roots of unity

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
   DFTk_x_Im_T<K,K*LastK,M,1,VType,S,W1> dft_scaled;
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
template<int_t N, typename NFact, typename VType, int S, int_t LastK = 1>
class DCT2_impl;

template<int_t N, typename Head, typename Tail, typename VType, int S, int_t LastK>
class DCT2_impl<N, Loki::Typelist<Head,Tail>, VType, S, LastK>
{
   typedef typename VType::ValueType T;
   typedef typename VType::TempType LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   static const int_t M2 = M*C;
   static const int_t N2 = N*C;
   static const int_t LastK2 = LastK*C;
   
   //typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   DCT2_impl<M,NFactNext,VType,S,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,1,VType,S,W1> dft_scaled;
//   DCTk_x_Im_T<K,M,1,VType,S,W1> dft_scaled;
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
template<int_t N, int_t K, typename Tail, typename VType, int S, int_t LastK>
class DCT2_impl<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, VType, S, LastK>
: public DCT2_impl<N, Tail, VType, S, LastK> {};


// Specialization for a prime N
template<int_t N, typename VType, int S, int_t LastK>
class DCT2_impl<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,VType,S,LastK> 
{
   typedef typename VType::ValueType T;
   static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;
   DCT2k<N, LastK*C, C, VType, S> spec;
public:
   void apply(const T* src, T* dst) { spec.apply(src, dst); }
};


}  //namespace DFT

#endif /*__gfftalg_h*/

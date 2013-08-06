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

#include "metacomplex.h"

namespace GFFT {

using namespace MF;


template<bool C>
struct StaticAssert;

template<>
struct StaticAssert<true> { };

#define GFFT_STATIC_ASSERT(c) StaticAssert<c> static__assert; 




// TODO: compare this with the simple loop
template<int_t M, typename T, int LastK, int NIter>
class IterateInTime {
   static const unsigned int M2 = M*2;
   IterateInTime<M,T,LastK,NIter-1> next;
public:
   template<class InTimeType>
   void apply(InTimeType& step, T* data) 
   {
      next.apply(step, data);
      step.apply(data + NIter*M2);
   }
   template<class InTimeType>
   void apply(InTimeType& step, const T* src, T* dst) 
   {
      next.apply(step, src, dst);
      step.apply(src + NIter*2*LastK, dst + NIter*M2);
   }
};

template<int_t M, typename T, int LastK>
class IterateInTime<M,T,LastK,0> {
public:
   template<class InTimeType>
   void apply(InTimeType& step, T* data) 
   {
      step.apply(data);
   }
   template<class InTimeType>
   void apply(InTimeType& step, const T* src, T* dst) 
   {
      step.apply(src, dst);
   }
};


template<int_t K, int_t M, typename T, int S>
class DFTk_x_Im_T 
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = K*M;
   static const int_t M2 = M*2;
   DFTk_inplace<K,M2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[K-1], wi[K-1], wpr[K-1], wpi[K-1], t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr[0], wpi[0])
      wpr[0] = 1 - 2.0*t*t;
      wpi[0] = -S*Sin<N,2,LocalVType>::value();
      
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

template<int_t M, typename T, int S>
class DFTk_x_Im_T<3,M,T,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = 3*M;
   static const int_t M2 = M*2;
   DFTk_inplace<3,M2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[2],wi[2],t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr1, wpi1)
      const LocalVType wpr1 = 1 - 2.0*t*t;
      const LocalVType wpi1 = -S*Sin<N,2,LocalVType>::value();
      
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

template<int_t M, typename T, int S>
class DFTk_x_Im_T<2,M,T,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = 2*M;
   DFTk_inplace<2,N,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr,wi,t;
      t = Sin<N,1,LocalVType>::value();
      const LocalVType wpr = -2.0*t*t;
      const LocalVType wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1+wpr;
      wi = wpi;
      for (int_t i=2; i<N; i+=2) {
	spec_inp.apply(data+i, &wr, &wi);

        t = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + t*wpi;
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
template<int_t N, typename NFact, typename T, int S, int_t LastK = 1>
class InTime;

template<int_t N, typename Head, typename Tail, typename T, int S, int_t LastK>
class InTime<N, Loki::Typelist<Head,Tail>, T, S, LastK>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   static const int_t M2 = M*2;
   static const int_t N2 = N*2;
   static const int_t LastK2 = LastK*2;
   
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   //IterateInTime<M,T,LastK,K-1> iter;
   InTime<M,NFactNext,T,S,K*LastK> dft_str;
   DFTk_x_Im_T<K,M,T,S> dft_scaled;
public:
   void apply(T* data) 
   {
      for (int_t m = 0; m < N2; m+=M2)
	dft_str.apply(data + m);

      dft_scaled.apply(data);
   }

   void apply(const T* src, T* dst, T* buf) 
   {
      int_t lk = 0;
      for (int_t m = 0; m < N2; m+=M2, lk+=LastK2)
        dft_str.apply(src + lk, dst + m, buf);
//      iter.apply(step, src, dst);

      dft_scaled.apply(dst);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename T, int S, int_t LastK>
class InTime<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, T, S, LastK>
: public InTime<N, Tail, T, S, LastK> {};


// Specialization for prime N
template<int_t N, typename T, int S, int_t LastK>
class InTime<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,T,S,LastK> {
  DFTk<N, LastK*2, 2, T, S> spec;
  DFTk_inplace<N, 2, T, S> spec_inp;
public:
  void apply(T* data) 
  { 
    spec_inp.apply(data);
  }
  void apply(const T* src, T* dst, T*) 
  { 
    spec.apply(src, dst);
  }
};

/*
// Specialization for N=1, decimation-in-time
template<typename T, int S, int_t LastK>
class InTime<1,T,S,LastK> {
public:
   void apply(T* data) { }
   void apply(const T* src, T* dst, T*) { *dst = *src; }
};
*/

/////////////////////////////////////////////////////////

// template<class W, class W1, int_t N, int_t I, int S>
// struct GetNextRoot;
// 
// template<class H, class T, class H1, class T1, int_t N, int_t I, int S>
// struct GetNextRoot<Loki::Typelist<H,T>,Loki::Typelist<H1,T1>,N,I,S>
// {
//    typedef typename Mult<H,H1>::Result W;
//    typedef typename GetNextRoot<T,T1,N,I,S>::Result Next;
//    typedef Loki::Typelist<W,Next> Result;
// };
// 
// template<int_t N, int_t I, int S>
// struct GetNextRoot<Loki::NullType,Loki::NullType,N,I,S>
// {
//    typedef Loki::NullType Result;
// };

template<int_t I, int_t M, class W1, class W, int Accuracy>
struct GetNextRoot {
  typedef typename Mult<W1,W>::Result Result;
};

template<class W1, class W, int_t I, int Accuracy>
struct GetNextRoot<I,2,W1,W,Accuracy> {
  typedef typename EX::CosPiFrac<I,2,Accuracy>::Result Re;
  typedef typename EX::SinPiFrac<I,2,Accuracy>::Result Im;
  typedef typename EX::FractionToDecimal<Re,Accuracy>::Result ReDec;
  typedef typename EX::FractionToDecimal<Im,Accuracy>::Result ImDec;
  typedef MComplex<ReDec,ImDec> Result;
};

template<class W1, class W, int_t I, int Accuracy>
struct GetNextRoot<I,1,W1,W,Accuracy> {
  typedef typename EX::CosPiFrac<I,1,Accuracy>::Result Re;
  typedef typename EX::SinPiFrac<I,1,Accuracy>::Result Im;
  typedef typename EX::FractionToDecimal<Re,Accuracy>::Result ReDec;
  typedef typename EX::FractionToDecimal<Im,Accuracy>::Result ImDec;
  typedef MComplex<ReDec,ImDec> Result;
};


// TODO: compare this with the simple loop
template<int_t K, int_t M, typename T, int S, class W1, int NIter = 1, class W = W1>
class IterateInFreq
{
// template<int_t K, int_t M, typename T, int S, class H, class Tail, int NIter = 1>
// class IterateInFreq<K,M,T,S,Loki::Typelist<H,Tail>,NIter> {
   //typedef typename RList::Head H;
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef EX::Compute<typename W::Re,2> WR;
   typedef EX::Compute<typename W::Im,2> WI;
   static const int_t M2 = M*2;
   static const int_t N = K*M;

   typedef typename GetNextRoot<NIter+1,N,W1,W,2>::Result Wnext;
   IterateInFreq<K,M,T,S,W1,NIter+1,Wnext> next;
   DFTk_inplace<K,M2,T,S> spec_inp;
   
public:
   void apply(T* data) 
   {
//      LocalVType wr[K-1], wi[K-1];
//       const LocalVType t = Sin<N,NIter-1,LocalVType>::value();
//       const LocalVType wr = 1 - 2.0*t*t;
//       const LocalVType wi = -S*Sin<N,2*(NIter-1),LocalVType>::value();
      const LocalVType wr = WR::value();
      const LocalVType wi = WI::value();
//       const LocalVType wr = H::first::value();
//       const LocalVType wi = H::second::value();
std::cout << NIter-1 << "/" << N << ": (" << wr << ", " << wi << ")" << std::endl;

      spec_inp.apply(&wr, &wi, data + (NIter-1)*2);

      next.apply(data);
   }
};

template<int_t K, int_t M, typename T, int S, class W1, class W>
class IterateInFreq<K,M,T,S,W1,M,W> 
{
//    typedef typename RList::Head H;
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef EX::Compute<typename W::Re,2> WR;
   typedef EX::Compute<typename W::Im,2> WI;
   static const int_t N = K*M;
   DFTk_inplace<K,M*2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
//       const LocalVType t = Sin<N,M-1,LocalVType>::value();
//       const LocalVType wr = 1 - 2.0*t*t;
//       const LocalVType wi = -S*Sin<N,2*(M-1),LocalVType>::value();
      const LocalVType wr = WR::value();
      const LocalVType wi = WI::value();
//       const LocalVType wr = H::first::value();
//       const LocalVType wi = H::second::value();
std::cout << M-1 << "/" << N << ": (" << wr << ", " << wi << ")" << std::endl;

      spec_inp.apply(&wr, &wi, data + (M-1)*2);
   }
};

template<int_t K, int_t M, typename T, int S, class W1, class W>
class IterateInFreq<K,M,T,S,W1,1,W> {
   DFTk_inplace<K,M*2,T,S> spec_inp;
   IterateInFreq<K,M,T,S,W1,2,W> next;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);
      next.apply(data);
   }
};

/////////////////////////////////////////////////////////

template<int_t K, int_t M, typename T, int S, class W, bool doStaticLoop>
class T_DFTk_x_Im;

template<int_t K, int_t M, typename T, int S, class W>
class T_DFTk_x_Im<K,M,T,S,W,true>
{
   IterateInFreq<K,M,T,S,W> iterate;
public:
   void apply(T* data) 
   {
      iterate.apply(data);
   }
   void apply(const T* src, T* dst) 
   {
   }
};

template<int_t K, int_t M, typename T, int S, class W>
class T_DFTk_x_Im<K,M,T,S,W,false>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = K*M;
   static const int_t M2 = M*2;
   DFTk<K,M2,M2,T,S> spec;
   DFTk_inplace<K,M2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[K-1], wi[K-1], wpr[K-1], wpi[K-1], t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr[0], wpi[0])
      wpr[0] = 1 - 2.0*t*t;
      wpi[0] = -S*Sin<N,2,LocalVType>::value();
      
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
	spec_inp.apply(wr, wi, data+i);

	for (int_t i=0; i<K-1; ++i) {
	  t = wr[i];
	  wr[i] = t*wpr[i] - wi[i]*wpi[i];
	  wi[i] = wi[i]*wpr[i] + t*wpi[i];
	}
      }
   }
   void apply(const T* src, T* dst) 
   {
      spec.apply(src, dst);

      LocalVType wr[K-1], wi[K-1], wpr[K-1], wpi[K-1], t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr[0], wpi[0])
      wpr[0] = 1 - 2.0*t*t;
      wpi[0] = -S*Sin<N,2,LocalVType>::value();
      
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
	spec.apply(wr, wi, src+i, dst+i);

	for (int_t i=0; i<K-1; ++i) {
	  t = wr[i];
	  wr[i] = t*wpr[i] - wi[i]*wpi[i];
	  wi[i] = wi[i]*wpr[i] + t*wpi[i];
	}
      }
   }
};

template<int_t M, typename T, int S, class W>
class T_DFTk_x_Im<3,M,T,S,W,false> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = 3*M;
   static const int_t M2 = M*2;
   DFTk<3,M2,M2,T,S> spec;
   DFTk_inplace<3,M2,T,S> spec_inp;
public:
   void apply(T* data) 
   {
      spec_inp.apply(data);

      LocalVType wr[2],wi[2],t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr1, wpi1)
      const LocalVType wpr1 = 1 - 2.0*t*t;
      const LocalVType wpi1 = -S*Sin<N,2,LocalVType>::value();
      
      // W^2 = (wpr2, wpi2)
      const LocalVType wpr2 = wpr1*wpr1 - wpi1*wpi1;
      const LocalVType wpi2 = 2*wpr1*wpi1;
      
      wr[0] = wpr1;
      wi[0] = wpi1;
      wr[1] = wpr2;
      wi[1] = wpi2;
      for (int_t i=2; i<M2; i+=2) {
	spec_inp.apply(wr, wi, data+i);

        t = wr[0];
        wr[0] = t*wpr1 - wi[0]*wpi1;
        wi[0] = wi[0]*wpr1 + t*wpi1;
        t = wr[1];
        wr[1] = t*wpr2 - wi[1]*wpi2;
        wi[1] = wi[1]*wpr2 + t*wpi2;
      }
   }
   void apply(const T* src, T* dst) 
   {
      spec.apply(src, dst);

      LocalVType wr[2],wi[2],t;
      t = Sin<N,1,LocalVType>::value();

      // W = (wpr1, wpi1)
      const LocalVType wpr1 = 1 - 2.0*t*t;
      const LocalVType wpi1 = -S*Sin<N,2,LocalVType>::value();
      
      // W^2 = (wpr2, wpi2)
      const LocalVType wpr2 = wpr1*wpr1 - wpi1*wpi1;
      const LocalVType wpi2 = 2*wpr1*wpi1;
      
      wr[0] = wpr1;
      wi[0] = wpi1;
      wr[1] = wpr2;
      wi[1] = wpi2;
      for (int_t i=2; i<M2; i+=2) {
	spec.apply(wr, wi, src+i, dst+i);

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
class T_DFTk_x_Im<2,M,T,S,W,false>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = 2*M;
   DFTk<2,N,N,T,S> spec;
   DFTk_inplace<2,N,T,S> spec_inp;

   IterateInFreq<2,M,T,S,W> iterate;
public:
   void apply(T* data) 
   {
      iterate.apply(data);
//       spec_inp.apply(data);
// 
//       LocalVType t,wr,wi;
//       t = Sin<N,1,LocalVType>::value();
//       const LocalVType wpr = -2.0*t*t;
//       const LocalVType wpi = -S*Sin<N,2,LocalVType>::value();
//       wr = 1+wpr;
//       wi = wpi;
//       for (int_t i=2; i<N; i+=2) {
// std::cout << i << "/" << N << ": " << t << " --- (" << wr << ", " << wi << ")" << std::endl;
// 	spec_inp.apply(&wr, &wi, data+i);
// 
//         t = wr;
//         wr += t*wpr - wi*wpi;
//         wi += wi*wpr + t*wpi;
//       }
   }
   
   void apply(const T* src, T* dst) 
   { 
      spec.apply(src, dst);
     
      LocalVType t,wr,wi;
      t = Sin<N,1,LocalVType>::value();
      const LocalVType wpr = -2.0*t*t;
      const LocalVType wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1+wpr;
      wi = wpi;
      for (int_t i=2; i<N; i+=2) {
	spec.apply(&wr, &wi, src+i, dst+i);

        t = wr;
        wr += t*wpr - wi*wpi;
        wi += wi*wpr + t*wpi;
      }
   }
};

/// Danielson-Lanczos section of the decimation-in-frequency FFT version
/**
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward

This template class implements resursive metaprogram, which
runs funciton apply() twice recursively at the end of the function apply()
with the half of the transform length N
until the simplest case N=2 has been reached. Then function \a _spec2 is called.
Therefore, it has two specializations for N=2 and N=1 (the trivial and empty case).
\sa InTime
*/
template<int_t N, typename NFact, typename T, int S, class W1, int_t LastK = 1>
class InFreq;

template<int_t N, typename Head, typename Tail, typename T, int S, class W1, int_t LastK>
class InFreq<N, Loki::Typelist<Head,Tail>, T, S, W1, LastK>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   static const int_t M2 = M*2;
   static const int_t N2 = N*2;
   static const int_t LastK2 = LastK*2;
   
//    typedef typename Loki::TL::Next<RList,K-1>::Result RListK;
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   InFreq<M,NFactNext,T,S,WK,K*LastK> dft_str;
   T_DFTk_x_Im<K,M,T,S,W1,true> dft_scaled;

public:
   void apply(T* data) 
   {
      dft_scaled.apply(data);

      // K times call to dft_str.apply()
      for (int_t m = 0; m < N2; m+=M2)
        dft_str.apply(data + m);
   }

   void apply(const T* src, T* dst, T* buf) 
   { 
      dft_scaled.apply(src, buf);

      int_t lk = 0;
      for (int_t m = 0; m < N2; m+=M2, lk+=LastK2)
	dft_str.apply(buf + m, dst + lk, buf + m);
   }
};


// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename T, int S, class W1, int_t LastK>
class InFreq<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, T, S, W1, LastK>
: public InFreq<N, Tail, T, S, W1, LastK> {};


// Specialization for prime N
template<int_t N, typename T, int S, class W1, int_t LastK>
class InFreq<N, Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,T,S,W1,LastK> {
  DFTk<N, 2, LastK*2, T, S> spec;
  DFTk_inplace<N, 2, T, S> spec_inp;
public:
  void apply(T* data) 
  { 
    spec_inp.apply(data);
  }
  void apply(const T* src, T* dst, T*) 
  { 
    spec.apply(src, dst);
  }
};


/// Binary reordering of array elements
/*!
\tparam N length of the data
\tparam T value type

This is C-like implementation. It has been written very 
similar to the one presented in the book 
"Numerical recipes in C++".
\sa GFFTswap2
*/
template<int_t N, typename T>
class GFFTswap {
public:
   void apply(T* data) {
     int_t m, j = 0;
     for (int_t i=0; i<2*N-1; i+=2) {
        if (j>i) {
            std::swap(data[j], data[i]);
            std::swap(data[j+1], data[i+1]);
        }
        m = N;
        while (m>=2 && j>=m) {
            j -= m;
            m >>= 1;
        }
        j += m;
     }
   }

   void apply(const T*, T*) { }
   void apply(const T*, T* dst, T*) { }
};


/// Binary reordering of array elements
/*!
\tparam N length of the data
\tparam T value type

This is second version of binary reordering.
It is based on template class recursion
similar to InTime and InFreq template classes,
where member function apply() is called twice recursively
building the parameters n and r, which are at last the
indexes of exchanged data values.

This version is slightly slower than GFFTswap, but 
allows parallelization of this algorithm, which is
implemented in template class GFFTswap2OMP.
\sa GFFTswap, GFFTswap2OMP
*/
template<unsigned int P, typename T,
unsigned int I=0>
class GFFTswap2 {
   static const int_t BN = 1<<(I+1);
   static const int_t BR = 1<<(P-I);
   GFFTswap2<P,T,I+1> next;
public:
   void apply(T* data, int_t n=0, int_t r=0) {
      next.apply(data,n,r);
      next.apply(data,n|BN,r|BR);
   }
};

template<unsigned int P, typename T>
class GFFTswap2<P,T,P> {
public:
   void apply(T* data, int_t n=0, int_t r=0) {
      if (n>r) {
        swap(data[n],data[r]);
        swap(data[n+1],data[r+1]);
      }
   }
};


/// Reordering of data for real-valued transforms
/*!
\tparam N length of the data
\tparam T value type
\tparam S sign of the transform: 1 - forward, -1 - backward
*/
template<int_t N, typename T, int S>
class Separate {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int M = (S==1) ? 2 : 1;
public:
   void apply(T* data) {
      int_t i,i1,i2,i3,i4;
      LocalVType wtemp,wr,wi,wpr,wpi;
      LocalVType h1r,h1i,h2r,h2i,h3r,h3i;
      wtemp = Sin<2*N,1,LocalVType>::value();
      wpr = -2.*wtemp*wtemp;
      wpi = -S*Sin<N,1,LocalVType>::value();
      wr = 1.+wpr;
      wi = wpi;
      for (i=1; i<N/2; ++i) {
        i1 = i+i;
        i2 = i1+1;
        i3 = 2*N-i1;
        i4 = i3+1;
        h1r = 0.5*(data[i1]+data[i3]);
        h1i = 0.5*(data[i2]-data[i4]);
        h2r = S*0.5*(data[i2]+data[i4]);
        h2i =-S*0.5*(data[i1]-data[i3]);
        h3r = wr*h2r - wi*h2i;
        h3i = wr*h2i + wi*h2r;
        data[i1] = h1r + h3r;
        data[i2] = h1i + h3i;
        data[i3] = h1r - h3r;
        data[i4] =-h1i + h3i;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }
      h1r = data[0];
      data[0] = M*0.5*(h1r + data[1]);
      data[1] = M*0.5*(h1r - data[1]);

      if (N>1) data[N+1] = -data[N+1];
   }
   
   void apply(const T*, T*) { }
   void apply(const T*, T*, T*) { }
};

// Policy for a definition of forward FFT
template<int_t N, typename T>
struct Forward {
   enum { Sign = 1 };
   void apply(T*) { }
   void apply(const T*, T*) { }
   void apply(const T*, T*, T*) { }
};

template<int_t N, typename T,
template<typename> class Complex>
struct Forward<N,Complex<T> > {
   enum { Sign = 1 };
   void apply(Complex<T>*) { }
   void apply(const Complex<T>*, Complex<T>*) { }
   void apply(const Complex<T>*, Complex<T>*, Complex<T>*) { }
};

// Policy for a definition of backward FFT
template<int_t N, typename T>
struct Backward {
   enum { Sign = -1 };
   void apply(T* data) {
      for (T* i=data; i<data+2*N; ++i) *i/=N;
   }
   void apply(const T*, T* dst) { apply(dst); }
   void apply(const T*, T* dst, T*) { apply(dst); }
};

template<int_t N, typename T,
template<typename> class Complex>
struct Backward<N,Complex<T> > {
   enum { Sign = -1 };
   void apply(Complex<T>* data) {
      for (int_t i=0; i<N; ++i) {
        data[i]/=N;
      }
   }
   void apply(const Complex<T>*, Complex<T>* dst) { apply(dst); }
   void apply(const Complex<T>*, Complex<T>* dst, Complex<T>*) { apply(dst); }
};


}  //namespace DFT

#endif /*__gfftalg_h*/

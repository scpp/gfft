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


namespace GFFT {

using namespace MF;

typedef unsigned long int_t;

template<typename T>
struct TempTypeTrait;

template<>
struct TempTypeTrait<float> {
   typedef double Result;
};

template<>
struct TempTypeTrait<double> {
   typedef long double Result;
};

template<typename T,
template<typename> class Complex>
struct TempTypeTrait<Complex<T> > {
   typedef typename TempTypeTrait<T>::Result Result;
};

template<typename T, typename A,
template<typename,typename> class Complex>
struct TempTypeTrait<Complex<T,A> > {
   typedef typename TempTypeTrait<T>::Result Result;
};

// template<typename T, typename A,
// template<typename,typename> class Complex>
// struct TempTypeTrait<Complex<T,A> > {
//    typedef T Result;
// };

// Look here for a small prime factor using 6k+1, 6k+5 algorithm
// until some relative small limit (e.g. 100)
// then rely on some prime factor algorithm like Rader (only for primes), Winograd or Bluestein (for any n)
// then come back to factoring
template<int_t N>
struct GetNextFactor {
  // dummy trial: assume multiple of 2 and 3 only
  static const int value = (N%2 == 0) ? 2 : 3; 
};


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
class DFTk_x_Im_T;

template<int_t M, typename T, int S>
class DFTk_x_Im_T<3,M,T,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = 3*M;
   static const int_t M2 = M*2;
  
public:
   void apply(T* data) 
   {
      LocalVType tr1,ti1,tr2,ti2,wr1,wi1,wr2,wi2,t;

      const T c = S * Sqrt<3, T>::value() * 0.5;

      const int_t i10 = M2;
      const int_t i11 = i10+1;
      const int_t i20 = i10+M2;
      const int_t i21 = i20+1;

      const T sr = data[i10] + data[i20];
      const T dr = c*(data[i10] - data[i20]);
      const T si = data[i11] + data[i21];
      const T di = c*(data[i11] - data[i21]);
      const T tr = data[0] - 0.5*sr;
      const T ti = data[1] - 0.5*si;
      data[0] += sr;
      data[1] += si;
      data[i10] = tr + di;
      data[i11] = ti - dr;
      data[i20] = tr - di;
      data[i21] = ti + dr;

      t = Sin<N,1,LocalVType>::value();

      // W = (wpr1, wpi1)
      const LocalVType wpr1 = 1 - 2.0*t*t;
      const LocalVType wpi1 = -S*Sin<N,2,LocalVType>::value();
      
      // W^2 = (wpr2, wpi2)
      const LocalVType wpr2 = wpr1*wpr1 - wpi1*wpi1;
      const LocalVType wpi2 = 2*wpr1*wpi1;
      
      wr1 = wpr1;
      wi1 = wpi1;
      wr2 = wpr2;
      wi2 = wpi2;
      for (int_t i=2; i<M2; i+=2) {
	const int_t i10 = i+M2;
	const int_t i11 = i10+1;
	const int_t i20 = i10+M2;
	const int_t i21 = i20+1;
        tr1 = data[i10]*wr1 - data[i11]*wi1;
        ti1 = data[i10]*wi1 + data[i11]*wr1;
        tr2 = data[i20]*wr2 - data[i21]*wi2;
        ti2 = data[i20]*wi2 + data[i21]*wr2;

	const T sr = tr1 + tr2;
	const T dr = c*(tr1 - tr2);
	const T si = ti1 + ti2;
	const T di = c*(ti1 - ti2);
	const T tr = data[i] - 0.5*sr;
	const T ti = data[i+1] - 0.5*si;
	data[i] += sr;
	data[i+1] += si;
	data[i10] = tr + di;
	data[i11] = ti - dr;
	data[i20] = tr - di;
	data[i21] = ti + dr;

        t = wr1;
        wr1 = t*wpr1 - wi1*wpi1;
        wi1 = wi1*wpr1 + t*wpi1;
        t = wr2;
        wr2 = t*wpr2 - wi2*wpi2;
        wi2 = wi2*wpr2 + t*wpi2;
      }
   }
};

template<int_t M, typename T, int S>
class DFTk_x_Im_T<2,M,T,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t N = 2*M;
   static const int_t M2 = M*2;
  
public:
   void apply(T* data) 
   {
      LocalVType tr,ti,wr,wi,t;

      const int_t i10 = M2;
      const int_t i11 = i10+1;
      tr = data[i10];
      ti = data[i11];
      data[i10] = data[0]-tr;
      data[i11] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;

      t = Sin<N,1,LocalVType>::value();
      const LocalVType wpr = -2.0*t*t;
      const LocalVType wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1+wpr;
      wi = wpi;
      for (int_t i=2; i<M2; i+=2) {
	const int_t i10 = i+M2;
	const int_t i11 = i10+1;
        tr = data[i10]*wr - data[i11]*wi;
        ti = data[i10]*wi + data[i11]*wr;
        data[i10] = data[i]-tr;
        data[i11] = data[i+1]-ti;
        data[i] += tr;
        data[i+1] += ti;

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
template<int_t N, typename T, int S, int LastK = 1>
class InTime {
   static const unsigned int K = GetNextFactor<N>::value;
   static const unsigned int M = N/K;
   static const unsigned int M2 = M*2;
   //IterateInTime<M,T,LastK,K-1> iter;
   InTime<M,T,S,K*LastK> dft_str;
   DFTk_x_Im_T<K,M,T,S> dft_scaled;
public:
   void apply(T* data) 
   {
      for (unsigned int i = 0; i < K; ++i)
	dft_str.apply(data + i*M2);

      dft_scaled.apply(data);
   }

   void apply(const T* src, T* dst) 
   {
      for (unsigned int i = 0; i < K; ++i)
        dft_str.apply(src + i*2*LastK, dst + i*M2);
//      iter.apply(step, src, dst);

      dft_scaled.apply(dst);
   }
};

/*
/// Specialization for N=4, decimation-in-time
template<typename T, int S>
class InTime<4,T,S> {
public:
   void apply(T* data) {
      T tr = data[2];
      T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = S*(data[5]-ti);
      data[7] = S*(tr-data[4]);
      data[4] += tr;
      data[5] += ti;

      tr = data[4];
      ti = data[5];
      data[4] = data[0]-tr;
      data[5] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = data[2]-tr;
      data[7] = data[3]-ti;
      data[2] += tr;
      data[3] += ti;
   }
};
*/

// Specialization for N=3, decimation-in-time
template<typename T, int S, int LastK>
class InTime<3,T,S,LastK> {
  static const int K2 = LastK*2;
  static const int K4 = 2*K2;
public:
   void apply(T* data) 
   { 
    // _spec3_fwd(data); 
      const T c = S * Sqrt<3, T>::value() * 0.5;
      const T sr = data[2] + data[4];
      const T dr = c*(data[2] - data[4]);
      const T si = data[3] + data[5];
      const T di = c*(data[3] - data[5]);
      const T tr = data[0] - 0.5*sr;
      const T ti = data[1] - 0.5*si;
      data[0] += sr;
      data[1] += si;
      data[2] = tr + di;
      data[3] = ti - dr;
      data[4] = tr - di;
      data[5] = ti + dr;
   }
   void apply(const T* src, T* dst) 
   { 
     //_spec3_fwd(src, dst); 
    // 5 mult, 12 add
      const T c = S * Sqrt<3, T>::value() * 0.5;
      const T sr = src[K2] + src[K4];
      const T dr = c*(src[K2] - src[K4]);
      const T si = src[K2+1] + src[K4+1];
      const T di = c*(src[K2+1] - src[K4+1]);
      const T tr = src[0] - 0.5*sr;
      const T ti = src[1] - 0.5*si;
      dst[0] = src[0] + sr;
      dst[1] = src[1] + si;
      dst[2] = tr + di;
      dst[3] = ti - dr;
      dst[4] = tr - di;
      dst[5] = ti + dr;     
   }
};

// Specialization for N=2, decimation-in-time
template<typename T, int S, int LastK>
class InTime<2,T,S,LastK> {
  static const int K2 = LastK*2;
public:
  void apply(T* data) 
  { 
      const T tr = data[2];
      const T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
  }
  void apply(const T* src, T* dst) 
  { 
    const T v1(src[1]), v2(src[K2]), v3(src[K2+1]);
    dst[0] = (*src + v2);
    dst[1] = (v1 + v3);
    dst[2] = (*src - v2);
    dst[3] = (v1 - v3);
  }
};

// Specialization for N=1, decimation-in-time
template<typename T, int S, int LastK>
class InTime<1,T,S,LastK> {
public:
   void apply(T* data) { }
   void apply(const T* src, T* dst) { }
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
template<int_t N, typename T, int S>
class InFreq {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   InFreq<N/2,T,S> next;
public:
   void apply(T* data) {

      LocalVType wtemp,tempr,tempi,wr,wi,wpr,wpi;
//    Change dynamic calculation to the static one
//      wtemp = sin(M_PI/N);
      wtemp = Sin<N,1,LocalVType>::value();
      wpr = -2.0*wtemp*wtemp;
//      wpi = -sin(2*M_PI/N);
      wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1.0;
      wi = 0.0;
      for (int_t i=0; i<N; i+=2) {
        tempr = data[i] - data[i+N];
        tempi = data[i+1] - data[i+N+1];
        data[i] += data[i+N];
        data[i+1] += data[i+N+1];
        data[i+N] = tempr*wr - tempi*wi;
        data[i+N+1] = tempi*wr + tempr*wi;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }

      next.apply(data);
      next.apply(data+N);
   }

   void apply(const T*, T*) { }
};

// Specialization for N=4, decimation-in-frequency
// template<typename T, int S>
// class InFreq<4,T,S> {
// public:
//    void apply(T* data) {
//       T tr = data[4];
//       T ti = data[5];
//       data[4] = data[0]-tr;
//       data[5] = data[1]-ti;
//       data[0] += tr;
//       data[1] += ti;
//       tr = data[6];
//       ti = data[7];
//       data[6] = S*(data[3]-ti);
//       data[7] = S*(tr-data[2]);
//       data[2] += tr;
//       data[3] += ti;
// 
//       tr = data[2];
//       ti = data[3];
//       data[2] = data[0]-tr;
//       data[3] = data[1]-ti;
//       data[0] += tr;
//       data[1] += ti;
//       tr = data[6];
//       ti = data[7];
//       data[6] = data[4]-tr;
//       data[7] = data[5]-ti;
//       data[4] += tr;
//       data[5] += ti;
//    }
// };

// Specialization for N=2, decimation-in-frequency
template<typename T, int S>
class InFreq<2,T,S> {
public:
   void apply(T* data) { _spec2(data); }
   void apply(const T* src, T* dst) { _spec2(src, dst); }
};

// Specialization for N=1, decimation-in-frequency
template<typename T, int S>
class InFreq<1,T,S> {
public:
   void apply(T* data) { }
   void apply(const T*, T*) { }
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
};

// Policy for a definition of forward FFT
template<int_t N, typename T>
struct Forward {
   enum { Sign = 1 };
   void apply(T*) { }
   void apply(const T*, T*) { }
};

template<int_t N, typename T,
template<typename> class Complex>
struct Forward<N,Complex<T> > {
   enum { Sign = 1 };
   void apply(Complex<T>*) { }
   void apply(const Complex<T>*, Complex<T>*) { }
};

// Policy for a definition of backward FFT
template<int_t N, typename T>
struct Backward {
   enum { Sign = -1 };
   void apply(T* data) {
      for (T* i=data; i<data+2*N; ++i) *i/=N;
   }
   void apply(const T*, T* dst) {
      for (T* i=dst; i<dst+2*N; ++i) *i/=N;
   }
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
   void apply(const Complex<T>*, Complex<T>* dst) {
      for (int_t i=0; i<N; ++i) {
        dst[i]/=N;
      }
   }
};


}  //namespace DFT

#endif /*__gfftalg_h*/

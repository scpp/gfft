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

#ifndef __gfftstdalg_h
#define __gfftstdalg_h

/** \file
    \brief Recursive FFT algorithms
*/

#include "gfftstdspec.h"
#include "gfftfactor.h"
#include "gfftswap.h"

#include "metacomplex.h"
#include "metaroot.h"

namespace GFFT {

using namespace MF;


template<int_t K, int_t M, typename T, int S, class W1, int NIter, class W,
template<typename> class Complex>
class IterateInTime<K,M,Complex<T>,S,W1,NIter,W>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t M2 = M*2;
   static const int_t N = K*M;

   typedef typename GetNextRoot<NIter+1,N,W1,W,2>::Result Wnext;
   IterateInTime<K,M,Complex<T>,S,W1,NIter+1,Wnext> next;
   DFTk_inp<K,M2,Complex<T>,S> spec_inp;
   
public:
   void apply(Complex<T>* data) 
   {
      const Complex<LocalVType> w(WR::value(), WI::value());

      spec_inp.apply(data + (NIter-1)*2, &w);

      next.apply(data);
   }
};

// Last step of the loop
template<int_t K, int_t M, typename T, int S, class W1, class W,
template<typename> class Complex>
class IterateInTime<K,M,Complex<T>,S,W1,M,W> 
{
//    typedef typename RList::Head H;
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t M2 = M*2;
   static const int_t N = K*M;
   DFTk_inp<K,M2,Complex<T>,S> spec_inp;
public:
   void apply(Complex<T>* data) 
   {
      const Complex<LocalVType> w(WR::value(), WI::value());

      spec_inp.apply(data + (M-1)*2, &w);
   }
};


// First step in the loop
template<int_t K, int_t M, typename T, int S, class W1, class W,
template<typename> class Complex>
class IterateInTime<K,M,Complex<T>,S,W1,1,W> {
   static const int_t M2 = M*2;
   DFTk_inp<K,M2,Complex<T>,S> spec_inp;
   IterateInTime<K,M,Complex<T>,S,W1,2,W> next;
public:
   void apply(Complex<T>* data) 
   {
      spec_inp.apply(data);
      next.apply(data);
   }
};


/////////////////////////////////////////////////////////

template<int_t K, int_t M, typename T, int S, class W, bool doStaticLoop>
class DFTk_x_Im_T;

template<int_t K, int_t M, typename T, int S, class W,
template<typename> class Complex>
class DFTk_x_Im_T<K,M,Complex<T>,S,W,true>
{
   IterateInTime<K,M,Complex<T>,S,W> iterate;
public:
   void apply(Complex<T>* data) 
   {
      iterate.apply(data);
   }
};

template<int_t K, int_t M, typename T, int S, class W,
template<typename> class Complex>
class DFTk_x_Im_T<K,M,Complex<T>,S,W,false>
{
   //typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t N = K*M;
   //static const int_t M2 = M*2;
   DFTk_inp<K,M,Complex<T>,S> spec_inp;
public:
   void apply(Complex<T>* data) 
   {
      spec_inp.apply(data);

      Complex<T> w[K-1], wp[K-1];

      // W = (wpr[0], wpi[0])
      wp[0] = Complex<T>(WR::value(), WI::value());
      //LocalVType t = Sin<N,1,LocalVType>::value();
//       wp[0] = Complex<LocalVType>(1 - 2.0*t*t, -S*Sin<N,2,LocalVType>::value());
      
      // W^i = (wpr2, wpi2)
      for (int_t i=0; i<K-2; ++i) 
	wp[i+1] = wp[i]*wp[0];
      
      for (int_t i=0; i<K-1; ++i) 
	w[i] = wp[i];
      
      for (int_t i=1; i<M; i++) {
	spec_inp.apply(data+i, w);

	for (int_t i=0; i<K-1; ++i) 
	  w[i] = w[i]*wp[i];
      }
   }
  
};

template<int_t M, typename T, int S, class W,
template<typename> class Complex>
class DFTk_x_Im_T<3,M,Complex<T>,S,W,false> {
   //typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   static const int_t N = 3*M;
   //static const int_t M2 = M*2;
   DFTk_inp<3,M,Complex<T>,S> spec_inp;
public:
   void apply(Complex<T>* data) 
   {
      spec_inp.apply(data);

      Complex<T> w[2];

      // W = (wpr1, wpi1)
//       LocalVType t = Sin<N,1,LocalVType>::value();
//       const LocalVType wpr1 = 1 - 2.0*t*t;
//       const LocalVType wpi1 = -S*Sin<N,2,LocalVType>::value();
      Complex<T> wp1(WR::value(), WI::value());
      
      // W^2 = (wpr2, wpi2)
      Complex<T> wp2(wp1*wp1);
      
      w[0] = wp1;
      w[1] = wp2;
      for (int_t i=1; i<M; i++) {
	spec_inp.apply(data+i, w);

        w[0] = w[0]*wp1;
        w[1] = w[1]*wp2;
      }
   }
};

template<int_t M, typename T, int S, class W,
template<typename> class Complex>
class DFTk_x_Im_T<2,M,Complex<T>,S,W,false> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Compute<typename W::Re,2> WR;
   typedef Compute<typename W::Im,2> WI;
   DFTk_inp<2,M,Complex<T>,S> spec_inp;
public:
   void apply(Complex<T>* data) 
   {
      spec_inp.apply(data);

//    LocalVType  t = Sin<N,1,LocalVType>::value();
//       const LocalVType wpr = 1-2.0*t*t;
//       const LocalVType wpi = -S*Sin<N,2,LocalVType>::value();
      Complex<T> wp(WR::value(), WI::value());

      Complex<T> w(wp);
      for (int_t i=1; i<M; i++) {
	spec_inp.apply(data+i, &w);

        w = w*wp;
      }
   }
};

// template<int_t N, typename Head, typename Tail, typename T, int S, class W1, int_t LastK>
// class InTime<N, Loki::Typelist<Head,Tail>, T, S, W1, LastK>
// {
//   // Not implemented, because not allowed
//   // Transforms in-place are allowed for powers of primes only!!!
// };

template<int_t N, typename Head, typename T, int S, class W1, int_t LastK,
template<typename> class Complex>
class InTime<N, Loki::Typelist<Head,Loki::NullType>, Complex<T>, S, W1, LastK>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Loki::NullType> NFactNext;
   InTime<M,NFactNext,Complex<T>,S,WK,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,Complex<T>,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T<K,M,Complex<T>,S,W1,false> dft_scaled;
public:
   void apply(Complex<T>* data) 
   {
      for (int_t m=0; m < N; m+=M) 
	dft_str.apply(data + m);

      dft_scaled.apply(data);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename T, int S, class W1, int_t LastK,
template<typename> class Complex>
class InTime<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, Complex<T>, S, W1, LastK>
: public InTime<N, Tail, Complex<T>, S, W1, LastK> {};


// Specialization for prime N
template<int_t N, typename T, int S, class W1, int_t LastK,
template<typename> class Complex>
class InTime<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,Complex<T>,S,W1,LastK> {
  DFTk_inp<N, 1, Complex<T>, S> spec_inp;
public:
  void apply(Complex<T>* data) 
  { 
    spec_inp.apply(data);
  }
};



template<int_t N, typename Head, typename Tail, typename T, int S, class W1, int_t LastK,
template<typename> class Complex>
class InTimeOOP<N, Loki::Typelist<Head,Tail>, Complex<T>, S, W1, LastK>
{
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int_t K = Head::first::value;
   static const int_t M = N/K;
//    static const int_t M2 = M*2;
//    static const int_t N2 = N*2;
//    static const int_t LastK2 = LastK*2;
   
   typedef typename IPowBig<W1,K>::Result WK;
   typedef Loki::Typelist<Pair<typename Head::first, SInt<Head::second::value-1> >, Tail> NFactNext;
   InTimeOOP<M,NFactNext,Complex<T>,S,WK,K*LastK> dft_str;
//   DFTk_x_Im_T<K,M,Complex<T>,S,W1,(N<=StaticLoopLimit)> dft_scaled;
   DFTk_x_Im_T<K,M,Complex<T>,S,W1,false> dft_scaled;
public:

   void apply(const Complex<T>* src, Complex<T>* dst) 
   {
      int_t lk = 0;
      for (int_t m = 0; m < N; m+=M, lk+=LastK)
        dft_str.apply(src + lk, dst + m);

      dft_scaled.apply(dst);
   }
};

// Take the next factor from the list
template<int_t N, int_t K, typename Tail, typename T, int S, class W1, int_t LastK,
template<typename> class Complex>
class InTimeOOP<N, Loki::Typelist<Pair<SInt<K>, SInt<0> >,Tail>, Complex<T>, S, W1, LastK>
: public InTimeOOP<N, Tail, Complex<T>, S, W1, LastK> {};


// Specialization for prime N
template<int_t N, typename T, int S, class W1, int_t LastK,
template<typename> class Complex>
class InTimeOOP<N,Loki::Typelist<Pair<SInt<N>, SInt<1> >, Loki::NullType>,Complex<T>,S,W1,LastK> {
  DFTk<N, LastK, 1, Complex<T>, S> spec;
public:
  void apply(const Complex<T>* src, Complex<T>* dst) 
  { 
    spec.apply(src, dst);
  }
};

  
}  //namespace DFT

#endif /*__gfftstdalg_h*/

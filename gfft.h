/***************************************************************************
 *   Copyright (C) 2007 by Volodymyr Myrnyy                                *
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

/** \file
    \brief Definition of main template classes, object factory and singleton
*/
/**
\mainpage
    Generic simple and efficient Fast Fourier Transforms (FFT) implementation
    using policy-based design and template metaprogramming. \n
    Features:
    - Transforms in-place
    - Complex FFT of power-of-two length. \n
      See Numerical recipes in C, Chap. 12.2 for theory \n
      http://www.nrbook.com/a/bookcpdf/c12-2.pdf
    - Single real function FFT of power-of-two length \n
      See Numerical recipes in C, Chap. 12.3 for theory and data format \n
      http://www.nrbook.com/a/bookcpdf/c12-3.pdf
    - High and cache-independent performance
    - No additional data is stored. \n
      You can use all available RAM for your transformed data
    - One-step transform. \n
      Many known FFT implementation perform to steps: initialization and transform.
      Initialization for a given length is usually computationally expensive,
      while transform is fast. GFFT needs only to create an object instance that
      includes FFT-algorithm, but no additional data or computation.

\section refs References

- Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T.
  Numerical recipes in C, 2th edition, 1992

*/

#ifndef __gfft_h
#define __gfft_h

#include "loki/Typelist.h"
#include "loki/Factory.h"
#include "loki/Singleton.h"
//#include "loki/Functor.h"

#include "finit.h"
#include "gfftalg.h"
#include "gfftpolicy.h"

namespace DFT {

enum FFTDirection { DIRECT, INVERSE };
enum FFTDecimation { INTIME, INFREQ };
enum FFTType { COMPLEX, REAL };

/// Generic Fast Fourier transform in-place
template<unsigned P, class T,
FFTType Type,                          // COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction=DIRECT,         // DIRECT, INVERSE
class FactoryPolicy=Empty>
class GFFT:public FactoryPolicy {
   enum { N = 1<<P };
   typedef GFFTswap<N,T> Swap;
   typedef typename Loki::Select<Direction==DIRECT,
           Direct<N,T>,Inverse<N,T> >::Result Dir;

   typedef InTime<N,T,Dir::Sign> InT;
   typedef InFreq<N,T,Dir::Sign> InF;
   typedef Separate<N,T,Dir::Sign> Sep;

   typedef typename Loki::Select<Decimation==INTIME,
           TYPELIST_3(Swap,InT,Dir),TYPELIST_3(InF,Swap,Dir)>::Result TList;
   typedef typename Loki::Select<Type==REAL,
      typename Loki::Select<Direction==DIRECT,
      typename Loki::TL::Append<TList,Sep>::Result,
               Loki::Typelist<Sep,TList> >::Result,TList>::Result Alg;

   PoliciesHandler<Alg> run;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new GFFT<P,T,Type,Decimation,Direction,FactoryPolicy>();
   }
   void fft(T* data) {
      run.apply(data);
   }
};

/*
/// Real-valued direct fast Fourier transform in-place
template<unsigned P, class T,
template<unsigned,class,class> class Decimation,  // InTime, InFreq
template<unsigned,class> class Direction,         // Direct, Inverse
class FactoryPolicy=Empty>
class RVGFFT : public FactoryPolicy {
   typedef Decimation<P,T,Direct<P,T> > Decim;
   enum { N = Decim::N };
   Decim gfft;
   Separate<N,T,1> sep;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new RVGFFT<P,T,Decimation,Direct,FactoryPolicy>();
   }
   void fft(T* data) {
      gfft.apply(data);
      sep.apply(data);
   }
};


/// Real-valued inverse fast Fourier transform in-place
template<unsigned P, class T,
template<unsigned,class,class> class Decimation,  // InTime, InFreq
class FactoryPolicy>
class RVGFFT<P,T,Decimation,Inverse,FactoryPolicy>:public FactoryPolicy {
   typedef Decimation<P,T,Inverse<P,T> > Decim;
   enum { N = Decim::N };
   Decim gfft;
   Separate<N,T,-1> sep;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new RVGFFT<P,T,Decimation,Inverse,FactoryPolicy>();
   }
   void fft(T* data) {
      sep.apply(data);
      gfft.apply(data);
   }
};
*/

/// Generator of a Typelist of GFFTs classes
template<unsigned Begin, unsigned End, typename T,
FFTType Type,                          // COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction=DIRECT,         // DIRECT, INVERSE
class FactoryPolicy=AbstractFFT<T> >
struct GFFTList {
   typedef Loki::Typelist<GFFT<Begin,T,Type,Decimation,Direction,FactoryPolicy>,
           typename GFFTList<Begin+1,End,T,Type,Decimation,Direction,
                    FactoryPolicy>::Result> Result;
};

template<unsigned End, typename T,
FFTType Type,                          // COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction,                // DIRECT, INVERSE
class FactoryPolicy>
struct GFFTList<End,End,T,Type,Decimation,Direction,FactoryPolicy> {
   typedef Loki::NullType Result;
};


/// Singleton for GFFT object factory
template<unsigned Min, unsigned Max, typename T,
FFTType Type,                          // COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction=DIRECT>
class GFFT_Singleton
: public Loki::SingletonHolder<Loki::Factory<AbstractFFT<T>,unsigned int> > {
   typedef Loki::Factory<AbstractFFT<T>,unsigned int> GFFT_Factory;
   typedef Loki::SingletonHolder<GFFT_Factory> Base;
   typedef FactoryInit<typename GFFTList<Min,Max,T,Type,Decimation,Direction>::Result> FInit;

   //Protection
   //GFFT_Singleton();
public:
   static GFFT_Factory& Instance() {
      FInit::apply(Base::Instance());
      return Base::Instance();
   }
};


}  // namespace DFT

#endif /*__gfft_h*/

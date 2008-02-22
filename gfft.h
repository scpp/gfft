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

\section start Getting started

GFFT consists now from only 5 headers, one sample programm gfft.cpp.
and some files from Loki-lib library by Andrei Alexandrescu.
It doesn't need any installation or additional packages.

\section license Licensing

  GFFT is open source software distributed under terms of GPL.

\section download Download

  Download the source code from project web site at SourceForge \n
  http://sourceforge.net/projects/gfft/download

\section refs References

- Myrnyy, V. A Simple and Efficient FFT Implementation in C++\n
  http://www.ddj.com/cpp/199500857
- Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T.
  Numerical recipes in C, 2th edition, 1992
- Alexandrescu, A. Modern C++ Design. Addison-Wesley, 2001
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

/// Direction of transform
enum FFTDirection { FORWARD, BACKWARD };

/// Type of decimation used in transform: in-time or in-frequency
enum FFTDecimation { INTIME, INFREQ };

/// Type of transform
enum FFTType { COMPLEX, REAL };


/// Generic Fast Fourier transform in-place
/**
\param P defines transform length, which is 2^P
\param T type of data element
\param Type type of transform: COMPLEX, REAL
\param Decimation in-time or in-frequency: INTIME, INFREQ
\param Direction transform direction: FORWARD, BACKWARD
\param FactoryPolicy policy used to create an object factory. Don't define it explicitely, if unsure
*/
template<unsigned P, class T,
FFTType Type,                          ///< COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction=FORWARD,         // FORWARD, BACKWARD
class FactoryPolicy=Empty>
class GFFT:public FactoryPolicy {
   enum { N = 1<<P };
   typedef GFFTswap<N,T> Swap;
   typedef typename Loki::Select<Direction==FORWARD,
           Forward<N,T>,Backward<N,T> >::Result Dir;

   typedef InTime<N,T,Dir::Sign> InT;
   typedef InFreq<N,T,Dir::Sign> InF;
   typedef Separate<N,T,Dir::Sign> Sep;

   typedef typename Loki::Select<Decimation==INTIME,
           TYPELIST_3(Swap,InT,Dir),TYPELIST_3(InF,Swap,Dir)>::Result TList;
   typedef typename Loki::Select<Type==REAL,
      typename Loki::Select<Direction==FORWARD,
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

/// Generator of a Typelist of GFFTs classes for data lengths from 2^Begin to 2^End
template<unsigned Begin, unsigned End, typename T,
FFTType Type,                          // COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction=FORWARD,         // FORWARD, BACKWARD
class FactoryPolicy=AbstractFFT<T> >
struct GFFTList {
   typedef Loki::Typelist<GFFT<Begin,T,Type,Decimation,Direction,FactoryPolicy>,
           typename GFFTList<Begin+1,End,T,Type,Decimation,Direction,
                    FactoryPolicy>::Result> Result;
};

template<unsigned End, typename T,
FFTType Type,                          // COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction,                // FORWARD, BACKWARD
class FactoryPolicy>
struct GFFTList<End,End,T,Type,Decimation,Direction,FactoryPolicy> {
   typedef Loki::NullType Result;
};


/// Singleton for GFFT object factory
template<unsigned Min, unsigned Max, typename T,
FFTType Type,                          // COMPLEX, REAL
FFTDecimation Decimation,              // INTIME, INFREQ
FFTDirection Direction=FORWARD>
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

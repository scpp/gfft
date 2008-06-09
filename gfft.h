/***************************************************************************
 *   Copyright (C) 2008 by Volodymyr Myrnyy                                *
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


#ifndef __gfft_h
#define __gfft_h

// General Doxygen documentation
#include "gfftdoc.h"

#include "loki/Typelist.h"
#include "loki/Factory.h"
#include "loki/Singleton.h"
//#include "loki/Functor.h"

#include "finit.h"
#include "gfftalg.h"
#include "gfftstdalg.h"
#include "gfftpolicy.h"

namespace DFT {

/// Direction of transform
enum FFTDirection { FORWARD, BACKWARD };

/// Type of decimation used in transform: in-time or in-frequency
enum FFTDecimation { INTIME, INFREQ };

/// Type of transform
enum FFTType { COMPLEX, REAL };


template<unsigned N, typename T,
FFTType Type, FFTDecimation Decimation, FFTDirection Dir>
class DecimationTrait;

template<unsigned N, typename T>
class DecimationTrait<N,T,COMPLEX,INTIME,FORWARD> {
   typedef Forward<N,T> Dir;
   typedef GFFTswap<N,T> Swap;
   typedef InTime<N,T,Dir::Sign> InT;
public:
   typename TYPELIST_3(Swap,InT,Dir) Result;
};

template<unsigned N, typename T>
class DecimationTrait<N,T,COMPLEX,INTIME,BACKWARD> {
   typedef Backward<N,T> Dir;
   typedef GFFTswap<N,T> Swap;
   typedef InTime<N,T,Dir::Sign> InT;
public:
   typename TYPELIST_3(Swap,InT,Dir) Result;
};

template<unsigned N, typename T>
class DecimationTrait<N,T,COMPLEX,INFREQ,FORWARD> {
   typedef Forward<N,T> Dir;
   typedef GFFTswap<N,T> Swap;
   typedef InFreq<N,T,Dir::Sign> InF;
public:
   typename TYPELIST_3(InF,Swap,Dir) Result;
};

template<unsigned N, typename T>
class DecimationTrait<N,T,COMPLEX,INFREQ,BACKWARD> {
   typedef Backward<N,T> Dir;
   typedef GFFTswap<N,T> Swap;
   typedef InFreq<N,T,Dir::Sign> InF;
public:
   typename TYPELIST_3(InF,Swap,Dir) Result;
};

template<unsigned N, typename T,
FFTDecimation Decimation>
class DecimationTrait<N,T,REAL,Decimation,FORWARD> {
   typedef Forward<N,T> Dir;
   typedef Separate<N,T,Dir::Sign> Sep;
   typedef typename DecimationTrait<N,T,COMPLEX,Decimation,FORWARD>::Result TList;
public:
   typedef typename Loki::TL::Append<TList,Sep>::Result Result;
};

template<unsigned N, typename T,
FFTDecimation Decimation>
class DecimationTrait<N,T,REAL,Decimation,BACKWARD> {
   typedef Forward<N,T> Dir;
   typedef Separate<N,T,Dir::Sign> Sep;
   typedef typename DecimationTrait<N,T,COMPLEX,Decimation,BACKWARD>::Result TList;
public:
   typedef Loki::Typelist<Sep,TList> Result;
};

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
FFTDecimation Decimation=INFREQ,              // INTIME, INFREQ
FFTDirection Direction=FORWARD,         // FORWARD, BACKWARD
class FactoryPolicy=Empty>
class GFFT:public FactoryPolicy {
   enum { N = 1<<P };
   typedef GFFTswap<N,T> Swap;
   typedef typename Loki::Select<Direction==BACKWARD,
           Backward<N,T>,Forward<N,T> >::Result Dir;

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
FFTDecimation Decimation=INFREQ,              // INTIME, INFREQ
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

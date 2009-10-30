/***************************************************************************
 *   Copyright (C) 2009 by Volodymyr Myrnyy                                *
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
    \brief Definition of main template classes
*/


#ifndef __gfft_h
#define __gfft_h

// General Doxygen documentation
#ifdef GFFTDOC
#include "gfftdoc.h"
#endif

#include "loki/Typelist.h"
#include "loki/Factory.h"
#include "loki/Singleton.h"

#include "finit.h"
#include "gfftomp.h"
#include "gfftstdalg.h"
#include "gfftpolicy.h"

#include "gfftcaller.h"

/** \brief Main namespace of GFFT

*/
namespace GFFT {


/** \class GFFT
\brief Generic Fast Fourier transform in-place (deprecated)
\tparam P defines transform length, which is 2^P
\tparam VType type of data element
\tparam Type type of transform: COMPLEX, REAL
\tparam Direction transform direction: FORWARD, BACKWARD
\tparam Parall parallelization
\tparam Decimation in-time or in-frequency: INTIME, INFREQ
\tparam FactoryPolicy policy used to create an object factory. Don't define it explicitely, if unsure
\deprecated
*/
template<unsigned P,
class VType,
class Type,                          ///< COMPLEX, REAL
class Direction,         // FORWARD, BACKWARD
class Parall,
class Decimation,              // INTIME, INFREQ
class FactoryPolicy=Empty,
unsigned IDN = P>
class GFFT:public FactoryPolicy {
   enum { N = 1<<P };
   typedef typename VType::ValueType T;
   //typedef GFFTswap<N,T> Swap;
   //typedef GFFTswap2<P,T> Swap;
   //typedef GFFTswapOMP<2,N,T> Swap;
   typedef typename Parall::template Swap<P,T>::Result Swap;
   typedef typename Direction::template Type<N,T> Dir;
   typedef Separate<N,T,Dir::Sign> Sep;
   typedef typename Decimation::template List<N,T,Swap,Dir,Parall::NParProc>::Result TList;
   typedef typename Type::template Algorithm<Direction,TList,Sep>::Result Alg;

   Caller<Loki::Typelist<Parall,Alg> > run;
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Decimation DecimationType;
   typedef Direction DirectionType;

   enum { ID = IDN, Len = N, PLen = P };

   static FactoryPolicy* Create() {
      return new GFFT<P,VType,Type,Direction,Parall,Decimation,FactoryPolicy>();
   }

   void fft(T* data) {
      run.apply(data);
   }
};


/** 
\brief Generic Fast Fourier transform in-place
\tparam Power2 defines transform length, which is 2^Power2
\tparam VType type of data element
\tparam Type type of transform: DFT, IDFT, RDFT, IRDFT
\tparam Dim dimension of transform, defined as SIntID<N>, N=1,2,...
\tparam Parall parallelization
\tparam Decimation in-time or in-frequency: INTIME, INFREQ
\tparam FactoryPolicy policy used to create an object factory. Don't define it explicitely, if unsure
*/
template<class Power2,
class VType,
class Type,                    ///< DFT, IDFT, RDFT, IRDFT
class Dim,
class Parall,
class Decimation,              // INTIME, INFREQ
class FactoryPolicy=Empty,
uint IDN = Power2::ID>
class Transform:public FactoryPolicy {
   enum { P = Type::template Length<Power2>::Value };
   enum { N = 1<<P };
   typedef typename VType::ValueType T;
   typedef typename Parall::template Swap<P,T>::Result Swap;
   typedef typename Type::template Direction<N,T> Dir;
   typedef Separate<N,T,Dir::Sign> Sep;
   typedef typename Decimation::template List<N,T,Swap,Dir,Parall::NParProc>::Result TList;
   typedef typename Type::template Algorithm<TList,Sep>::Result Alg;

   Caller<Loki::Typelist<Parall,Alg> > run;
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Decimation DecimationType;
//   typedef Direction DirectionType;

   enum { ID = IDN, Len = N, PLen = P };

   static FactoryPolicy* Create() {
      return new Transform<Power2,VType,Type,Dim,Parall,Decimation,FactoryPolicy>();
   }

   void fft(T* data) {
      run.apply(data);
   }
};


// Generator of a Typelist of GFFTs classes for data lengths from 2^Begin to 2^End
template<unsigned Begin, unsigned End, 
typename ValueType,
class TransformType,                 // COMPLEX, REAL
class Direction=FORWARD,             // FORWARD, BACKWARD
class Parall=Serial,
class Decimation=INFREQ,             // INTIME, INFREQ
class FactoryPolicy=AbstractFFT<typename ValueType::ValueType> >
struct GFFTList {
   typedef Loki::Typelist<GFFT<Begin,ValueType,TransformType,Direction,Parall,Decimation,FactoryPolicy>,
           typename GFFTList<Begin+1,End,ValueType,TransformType,Direction,Parall,Decimation,
                    FactoryPolicy>::Result> Result;
};

template<unsigned End, 
typename ValueType,
class TransformType,           // COMPLEX, REAL
class Direction,               // FORWARD, BACKWARD
class Parall,
class Decimation,              // INTIME, INFREQ
class FactoryPolicy>
struct GFFTList<End,End,ValueType,TransformType,Direction,Parall,Decimation,FactoryPolicy> {
   typedef Loki::NullType Result;
};

template<
class Type,                          // COMPLEX, REAL
class Direction=FORWARD,
class Parall=Serial,
class Decimation=INFREQ>
class GFFT_Options { };

/** \class GFFT_Singleton
\brief Singleton for GFFT object factory
\param Min minimum power of two 
\param Max maximum power of two
*/
template<unsigned Min, unsigned Max, 
typename ValueType,
class TransformType,                          // COMPLEX, REAL
class Direction=FORWARD,
class Parall=Serial,
class Decimation=INFREQ>              // INTIME, INFREQ
class GFFT_Singleton
: public Loki::SingletonHolder<Loki::Factory<AbstractFFT<typename ValueType::ValueType>,unsigned int>,
                               GFFT_Options<TransformType,Direction,Parall,Decimation> > {
   typedef Loki::Factory<AbstractFFT<typename ValueType::ValueType>,unsigned int> GFFT_Factory;
   typedef GFFT_Options<TransformType,Direction,Parall,Decimation> Options;
   typedef Loki::SingletonHolder<GFFT_Factory,Options> Base;
   typedef FactoryInit<typename GFFTList<Min,Max,ValueType,TransformType,Direction,Parall,Decimation>::Result> FInit;

   //Protection
   GFFT_Singleton();
public:
   static GFFT_Factory& Instance() {
      FInit::apply(Base::Instance());
      return Base::Instance();
   }
};

//
template<class VType>
struct GFFTFactory {
   typedef Loki::Factory<AbstractFFT<typename VType::ValueType>, unsigned int> _Factory;

};


}  // namespace GFFT

#endif /*__gfft_h*/

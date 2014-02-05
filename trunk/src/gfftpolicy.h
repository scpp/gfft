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

#ifndef __gfftpolicy_h
#define __gfftpolicy_h

/** \file
    \brief Policy classes
*/

#include "Typelist.h"

#include <complex>

#include <omp.h>

#include "sint.h"

namespace GFFT {

typedef unsigned int id_t;

/** \class AbstractFFT_inp
\brief Abstract interface class to build %GFFT object factory

This class represents basic interface for %GFFT classes.
In other words, it shares the function fft(T*) between
classes that represent FFT of different lengths and types.
*/
template<typename T>
class AbstractFFT_inp {
public:
   virtual void fft(T*) = 0;
   virtual ~AbstractFFT_inp() {}
};

template<typename T>
class AbstractFFT_oop {
public:
   virtual void fft(const T*, T*) = 0;
   virtual ~AbstractFFT_oop() {}
};

/** \class Empty
\brief Abstract empty base class

This class is passed instead of AbstractFFT,
if object factory is not needed.
The virtual function call is then avoided.
*/
class Empty { };


/*!
\defgroup gr_params GFFT parameters
\brief Classes substituted as template parameters to define transform 
*/

/*! \brief Double precision type representation
\ingroup gr_params
*/
struct DOUBLE {
   static const id_t ID = 0;
   typedef double ValueType;
#ifdef __x86_64
   static const int Accuracy = 2;
#else  
   static const int Accuracy = 4;
#endif
};

/*! \brief Single precision type representation
\ingroup gr_params
*/
struct FLOAT {
   static const id_t ID = 1;
   typedef float ValueType;
#ifdef __x86_64
   static const int Accuracy = 1;
#else
   static const int Accuracy = 2;
#endif
};

/*! \brief Complex number of double precision type representation
\ingroup gr_params
*/
struct COMPLEX_DOUBLE {
   static const id_t ID = 2;
   typedef std::complex<double> ValueType;
#ifdef __x86_64
   static const int Accuracy = 2;
#else  
   static const int Accuracy = 4;
#endif
};

/*! \brief Complex number of single precision type representation
\ingroup gr_params
*/
struct COMPLEX_FLOAT {
   static const id_t ID = 3;
   typedef std::complex<float> ValueType;
#ifdef __x86_64
   static const int Accuracy = 1;
#else
   static const int Accuracy = 2;
#endif
};

/*! \brief Decimation in-time
\ingroup gr_params
*/
struct INTIME {
   static const id_t ID = 0;
   template<int_t N, typename NFact, typename T,
            class Swap, class Direction, short_t NT, class W1>
   class List {
//      typedef InTime<N,NFact,T,Direction::Sign,W1> InT;
      typedef InTimeOMP<NT,N,NFact,T,Direction::Sign,W1> InT;
   public:
      typedef TYPELIST_3(Swap,InT,Direction) Result;
   };
};

struct INTIME_OOP {
   static const id_t ID = 0;
   template<int_t N, typename NFact, typename T,
            class Swap, class Direction, short_t NT, class W1>
   class List {
//      typedef InTimeOOP<N,NFact,T,Direction::Sign,W1> InT;
      typedef InTimeOOP_OMP<NT,N,NFact,T,Direction::Sign,W1> InT;
   public:
       typedef TYPELIST_3(Swap,InT,Direction) Result;
   };
};

/*! \brief Decimation in-frequency
\ingroup gr_params
*/
struct INFREQ {
   static const id_t ID = 1;
   template<int_t N, typename NFact, typename T,
            class Swap, class Direction, short_t NT, class W1>
   class List {
//      typedef InFreq<N,NFact,T,Direction::Sign,W1> InF;
      typedef InFreqOMP<NT,N,NFact,T,Direction::Sign,W1> InF;
   public:
      typedef TYPELIST_3(InF,Swap,Direction) Result;
   };
};

/*! \brief In-place algorithm 
\ingroup gr_params
*/
struct IN_PLACE {
   static const id_t ID = 0;
   template<class T>
   struct Interface {
     typedef AbstractFFT_inp<T> Result;  
   };
   template<int_t N, typename NFact, typename T,
            class Swap, class Direction, short_t NT, class W1>
   struct List {
      typedef typename INTIME::template List<N,NFact,T,Swap,Direction,NT,W1>::Result Result;
   };
};

/*! \brief Out-of-place algorithm
\ingroup gr_params
*/
struct OUT_OF_PLACE {
   static const id_t ID = 1;
   template<class T>
   struct Interface {
     typedef AbstractFFT_oop<T> Result;  
   };
   template<int_t N, typename NFact, typename T,
            class Swap, class Direction, short_t NT, class W1>
   struct List {
      typedef typename INTIME_OOP::template List<N,NFact,T,Swap,Direction,NT,W1>::Result Result;
   };
};

struct IDFT;
struct IRDFT;

/*! \brief Forward compex-valued discrete Fourier transform
\ingroup gr_params
*/
struct DFT {
   static const id_t ID = 0;
   typedef IDFT Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Forward<N,T> {};

   template<class TList, class Separator>
   struct Algorithm {
      typedef TList Result;
   };
};

/*! \brief Inverse compex-valued discrete Fourier transform
\ingroup gr_params
*/
struct IDFT {
   static const id_t ID = 1;
   typedef DFT Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Backward<N,T> {};

   template<class TList, class Separator>
   struct Algorithm {
      typedef TList Result;
   };
};

/*! \brief Forward real-valued discrete Fourier transform
\ingroup gr_params
*/
struct RDFT {
   static const id_t ID = 2;
   typedef IRDFT Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Forward<N,T> {};

   template<class TList, class Separator>
   struct Algorithm {
      typedef typename Loki::TL::Append<TList,Separator>::Result Result;
   };
};

/*! \brief Inverse real-valued discrete Fourier transform
\ingroup gr_params
*/
struct IRDFT {
   static const id_t ID = 3;
   typedef RDFT Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Backward<N,T> {};

   template<class TList, class Separator>
   struct Algorithm {
      typedef Loki::Typelist<Separator,TList> Result;
   };
};


/*! \brief %Serial (single-core) implementation of transform
\sa OpenMP
\ingroup gr_params
*/
struct Serial {
   static const id_t ID = 0;
   static const uint NParProc = 1;

   template<uint M, uint P, class T>
   struct Swap {
      typedef GFFTswap2<M,P,T> Result;
   };

   template<typename T>
   void apply(T*) { }

   template<typename T>
   void apply(const T*, T*) { }
};

/*! \brief %Transform is parallelized using %OpenMP standard
\tparam NT number of parallel threads
\sa Serial
\ingroup gr_params
*/
template<unsigned int NT>
struct OpenMP {
   static const id_t ID = NT-1;
   static const uint NParProc = NT;

   template<uint_t M, uint_t P, class T>
   struct Swap {
      typedef GFFTswap2<M,P,T> Result;
//       typedef GFFTswap2OMP<NT,M,P,T> Result;
   };

   template<typename T>
   void apply(T*) {
      //omp_set_dynamic(0);
      //omp_set_num_threads(NT);
      omp_set_nested(true);
   }

   template<typename T>
   void apply(const T*, T* d) { apply(d); }
};

template<>
struct OpenMP<0>:public Serial { };

template<>
struct OpenMP<1>:public Serial { };

  
}  //namespace GFFT

#endif /*__gfftpolicy_h*/

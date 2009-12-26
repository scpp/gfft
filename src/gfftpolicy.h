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

#ifndef __gfftpolicy_h
#define __gfftpolicy_h

/** \file
    \brief Policy classes
*/

#include "loki/Typelist.h"
#include <complex>


namespace GFFT {

typedef unsigned int id_t;

/** \class AbstractFFT
\brief Abstract interface class to build %GFFT object factory

This class represents basic interface for %GFFT classes.
In other words, it shares the function fft(T*) between
classes that represent FFT of different lengths and types.
*/
template<typename T>
class AbstractFFT {
public:
   virtual void fft(T*) = 0;
   virtual ~AbstractFFT() {}
};


/** \class Empty
\brief Abstract empty base class

This class is passed instead of AbstractFFT,
if object factory is not needed
to avoid a virtual function call penalty.
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
};

/*! \brief Single precision type representation
\ingroup gr_params
*/
struct FLOAT {
   static const id_t ID = 1;
   typedef float ValueType;
};

/*! \brief Complex number of double precision type representation
\ingroup gr_params
*/
struct COMPLEX_DOUBLE {
   static const id_t ID = 2;
   typedef std::complex<double> ValueType;
};

/*! \brief Complex number of single precision type representation
\ingroup gr_params
*/
struct COMPLEX_FLOAT {
   static const id_t ID = 3;
   typedef std::complex<float> ValueType;
};

/*! \brief Decimation in-time
\ingroup gr_params
*/
struct INTIME {
   static const id_t ID = 0;
   template<unsigned long N, typename T,
            class Swap,
            class Direction, unsigned int NT>
   class List {
//      typedef InTime<N,T,Direction::Sign> InT;
      typedef InTimeOMP<NT,N,T,Direction::Sign> InT;
   public:
      typedef TYPELIST_3(Swap,InT,Direction) Result;
   };
};

/*! \brief Decimation in-frequency
\ingroup gr_params
*/
struct INFREQ {
   static const id_t ID = 1;
   template<unsigned long N, typename T,
            class Swap,
            class Direction, unsigned int NT>
   class List {
//      typedef InFreq<N,T,Direction::Sign> InF;
      typedef InFreqOMP<NT,N,T,Direction::Sign> InF;
   public:
      typedef TYPELIST_3(InF,Swap,Direction) Result;
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

   template<typename TPower>
   struct Length {
      enum { Value = TPower::value };
   };

   template<unsigned long N, typename T>
   struct Direction : public Forward<N,T> {};

   template<class List, class Separator>
   struct Algorithm {
      typedef List Result;
   };
};

/*! \brief Inverse compex-valued discrete Fourier transform
\ingroup gr_params
*/
struct IDFT {
   static const id_t ID = 1;
   typedef DFT Inverse;

   template<typename TPower>
   struct Length {
      enum { Value = TPower::value };
   };

   template<unsigned long N, typename T>
   struct Direction : public Backward<N,T> {};

   template<class List, class Separator>
   struct Algorithm {
      typedef List Result;
   };
};

/*! \brief Forward real-valued discrete Fourier transform
\ingroup gr_params
*/
struct RDFT {
   static const id_t ID = 2;
   typedef IRDFT Inverse;

   template<typename TPower>
   struct Length {
      static const unsigned int Value = TPower::value-1;
   };

   template<unsigned long N, typename T>
   struct Direction : public Forward<N,T> {};

   template<class List, class Separator>
   struct Algorithm {
      typedef typename Loki::TL::Append<List,Separator>::Result Result;
   };
};

/*! \brief Inverse real-valued discrete Fourier transform
\ingroup gr_params
*/
struct IRDFT {
   static const id_t ID = 3;
   typedef RDFT Inverse;

   template<typename TPower>
   struct Length {
      enum { Value = TPower::value-1 };
   };

   template<unsigned long N, typename T>
   struct Direction : public Backward<N,T> {};

   template<class List, class Separator>
   struct Algorithm {
      typedef Loki::Typelist<Separator,List> Result;
   };
};


/*! \brief %Serial (single-core) implementation of transform
\sa OpenMP
\ingroup gr_params
*/
struct Serial {
   static const id_t ID = 0;
   static const unsigned int NParProc = 1;

   template<unsigned int P, class T>
   struct Swap {
      enum { N = 1<<P };
      typedef GFFTswap<N,T> Result;
   };

   template<typename T>
   void apply(T*) { }
};

/*! \brief %Transform is parallelized using %OpenMP standard
\tparam NT number of parallel threads
\sa Serial
\ingroup gr_params
*/
template<unsigned int NT>
struct OpenMP {
   static const id_t ID = NT-1;
   static const unsigned int NParProc = NT;

   template<unsigned int P, class T>
   struct Swap {
      //typedef GFFTswap<(1<<P),T> Result;
      typedef GFFTswap2OMP<NT,P,T> Result;
   };

   template<typename T>
   void apply(T*) {
      omp_set_num_threads(NT);
      omp_set_nested(true);
   }
};

template<>
struct OpenMP<0>:public Serial { };

}  //namespace GFFT

#endif /*__gfftpolicy_h*/

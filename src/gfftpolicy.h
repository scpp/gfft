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


/** \class AbstractFFT
\brief Abstract interface class to build GFFT object factory

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
\defgroup gr_params %GFFT parameters
\brief Classes substituted as template parameters to define transform 
*/

/*! \brief Double precision type representation
\ingroup gr_params
*/
struct DOUBLE {
   enum { ID = 0 };
   typedef double ValueType;
};

/*! \brief Single precision type representation
\ingroup gr_params
*/
struct FLOAT {
   enum { ID = 1 };
   typedef float ValueType;
};

/*! \brief Complex number of double precision type representation
\ingroup gr_params
*/
struct COMPLEX_DOUBLE {
   enum { ID = 2 };
   typedef std::complex<double> ValueType;
};

/*! \brief Complex number of single precision type representation
\ingroup gr_params
*/
struct COMPLEX_FLOAT {
   enum { ID = 3 };
   typedef std::complex<float> ValueType;
};

/*! \brief Forward direction of transform (deprecated)
\deprecated
*/
struct FORWARD {
   enum { ID = 0 };
   template<unsigned N, typename T>
   struct Type : public Forward<N,T> {};

   template<class List, class Separator>
   struct AddSeparator {
      typedef typename Loki::TL::Append<List,Separator>::Result Result;
   };
};

/*! \brief Backward direction of transform (deprecated)
\deprecated
*/
struct BACKWARD {
   enum { ID = 1 };
   template<unsigned N, typename T>
   struct Type : public Backward<N,T> {};

   template<class List, class Separator>
   struct AddSeparator {
      typedef Loki::Typelist<Separator,List> Result;
   };
};

/*! \brief Decimation in-time
\ingroup gr_params
*/
struct INTIME {
   enum { ID = 0 };
   template<unsigned N, typename T,
            class Swap,
            class Direction, unsigned NT>
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
   enum { ID = 1 };
   template<unsigned N, typename T,
            class Swap,
            class Direction, unsigned NT>
   class List {
//      typedef InFreq<N,T,Direction::Sign> InF;
      typedef InFreqOMP<NT,N,T,Direction::Sign> InF;
   public:
      typedef TYPELIST_3(InF,Swap,Direction) Result;
   };
};

/*! \brief Complex valued transform
\deprecated
*/
struct COMPLEX {
   enum { ID = 0 };
   template<class Direction, class List, class Separator>
   struct Algorithm {
      typedef List Result;
   };
};

/*! \brief Real valued transform
\deprecated
*/
struct REAL {
   enum { ID = 1 };
   template<class Direction, class List, class Separator>
   struct Algorithm {
      typedef typename Direction::
         template AddSeparator<List,Separator>::Result Result;
   };
};

struct REAL2 {
   enum { ID = 2 };
   template<class Direction, class List, class Separator>
   struct Algorithm {
      typedef typename Direction::
         template AddSeparator<List,Separator>::Result Result;
   };
};

/*! \brief Forward compex-valued discrete Fourier transform
\ingroup gr_params
*/
struct DFT {
   enum { ID = 0 };

   template<typename TPower>
   struct Length {
      enum { Value = TPower::value };
   };

   template<unsigned N, typename T>
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
   enum { ID = 1 };

   template<typename TPower>
   struct Length {
      enum { Value = TPower::value };
   };

   template<unsigned N, typename T>
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
   enum { ID = 2 };

   template<typename TPower>
   struct Length {
      static const unsigned int Value = TPower::value-1;
   };

   template<unsigned N, typename T>
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
   enum { ID = 3 };

   template<typename TPower>
   struct Length {
      enum { Value = TPower::value-1 };
   };

   template<unsigned N, typename T>
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
   enum { ID = 0 };
   static const unsigned NParProc = 1;

   template<unsigned P, class T>
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
template<unsigned NT>
struct OpenMP {
   enum { ID = NT-1 };
   static const unsigned NParProc = NT;

   template<unsigned P, class T>
   struct Swap {
      typedef GFFTswap2OMP<NT,P,T> Result;
   };

   template<typename T>
   void apply(T*) {
      omp_set_num_threads(NT);
      omp_set_nested(true);
   }
};

}  //namespace GFFT

#endif /*__gfftpolicy_h*/

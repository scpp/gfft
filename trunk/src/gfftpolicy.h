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
#include "twiddles.h"

static const int_t SwitchToOMP = (1<<8);

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
   typedef double base_type;
   typedef double ValueType;
   typedef long double TempType;
#ifdef __x86_64
   static const int Accuracy = 2;
#else  
   static const int Accuracy = 4;
#endif
   static const char* name() { return "double"; }
};

/*! \brief Single precision type representation
\ingroup gr_params
*/
struct FLOAT {
   static const id_t ID = 1;
   typedef float base_type;
   typedef float ValueType;
   typedef double TempType;
#ifdef __x86_64
   static const int Accuracy = 1;
#else
   static const int Accuracy = 2;
#endif
   static const char* name() { return "float"; }
};

/*! \brief Complex number of double precision type representation
\ingroup gr_params
*/
struct COMPLEX_DOUBLE {
   static const id_t ID = 2;
   typedef double base_type;
   typedef std::complex<double> ValueType;
   typedef std::complex<long double> TempType;
#ifdef __x86_64
   static const int Accuracy = 2;
#else  
   static const int Accuracy = 4;
#endif
   static const char* name() { return "std::complex<double>"; }
};

/*! \brief Complex number of single precision type representation
\ingroup gr_params
*/
struct COMPLEX_FLOAT {
   static const id_t ID = 3;
   typedef float base_type;
   typedef std::complex<float> ValueType;
   typedef std::complex<double> TempType;
#ifdef __x86_64
   static const int Accuracy = 1;
#else
   static const int Accuracy = 2;
#endif
   static const char* name() { return "std::complex<float>"; }
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

   template<int_t N, typename NFact, typename VType,
            typename Parall, typename Direction, typename Twiddles>
   class List {
      typedef typename VType::ValueType T;
      typedef typename Parall::template ActualParall<N>::Result NewParall;
      typedef typename NewParall::template Swap<NFact,T>::Result Swap;
      typedef typename GetFirstRoot<N,Direction::Sign,VType::Accuracy>::Result W1;
      typedef InTime_omp<NewParall::NParProc,N,NFact,VType,Direction::Sign,W1> InT;
   public:
      typedef TYPELIST_3(Swap,InT,Direction) Result;
   };
   
   template<typename FuncList, typename T>
   struct Function : public Interface<T>::Result
   {
      FuncList m_run;
    //   in-place transform
      void fft(T* data) 
      { 
	m_run.apply(data); 
      }
   };
   
   static const char* name() { return "in-place"; }
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

//    template<typename Parall, typename NFact, typename T>
//    struct Swap {
//       typedef Caller<Loki::NullType> Result;
//    };

   template<int_t N, typename NFact, typename VType,
            typename Parall, typename Direction, typename Twiddles>
   class List {
      typedef typename Parall::template ActualParall<N>::Result NewParall;
      typedef typename GetFirstRoot<N,Direction::Sign,VType::Accuracy>::Result W1;
      typedef InTimeOOP_omp<NewParall::NParProc,N,NFact,VType,Direction::Sign,W1> InT;
   public:
       typedef TYPELIST_2(InT,Direction) Result;
   };

   template<typename FuncList, typename T>
   struct Function : public Interface<T>::Result
   {
      FuncList m_run;
    // out-of-place transform
      void fft(const T* src, T* dst) 
      { 
	m_run.apply(src, dst); 
      }
   };

   static const char* name() { return "out-of-place"; }
};

struct IDFT;
struct IRDFT;
struct IDCT1;
struct IDCT2;

/*! \brief Forward compex-valued discrete Fourier transform
\ingroup gr_params
*/
struct DFT {
   static const id_t ID = 0;
   static const int Sign = 1;
   typedef IDFT Inverse;

   template<int_t N, typename NFact, typename VType,
            typename Parall, typename Place, typename Twiddles>
   class Algorithm {
      typedef typename VType::ValueType T;
      typedef Forward<N,T> Direction;
//      typedef typename GetFirstRoot<N,Direction::Sign,VType::Accuracy>::Result W1;
   public:
      typedef typename Place::template List<N,NFact,VType,Parall,Direction,Twiddles>::Result Result;
   };
};

/*! \brief Inverse compex-valued discrete Fourier transform
\ingroup gr_params
*/
struct IDFT {
   static const id_t ID = 1;
   static const int Sign = -1;
   typedef DFT Inverse;

   template<int_t N, typename NFact, typename VType,
            typename Parall, typename Place, typename Twiddles>
   class Algorithm {
      typedef typename VType::ValueType T;
      typedef Backward<N,T> Direction;
   public:
      typedef typename Place::template List<N,NFact,VType,Parall,Direction,Twiddles>::Result Result;
   };
};

/*! \brief Forward real-valued discrete Fourier transform
\ingroup gr_params
*/
struct RDFT {
   static const id_t ID = 2;
   static const int Sign = 1;
   typedef IRDFT Inverse;

   template<int_t N, typename NFact, typename VType,
            typename Parall, typename Place, typename Twiddles>
   class Algorithm {
      typedef typename VType::ValueType T;
      typedef Forward<N,T> Direction;
      typedef Separate<N,VType,Direction::Sign> Separator;
      typedef typename Place::template List<N,NFact,VType,Parall,Direction,Twiddles>::Result TList;
   public:
      typedef typename Loki::TL::Append<TList,Separator>::Result Result;
   };
};

/*! \brief Inverse real-valued discrete Fourier transform
\ingroup gr_params
*/
struct IRDFT {
   static const id_t ID = 3;
   static const int Sign = -1;
   typedef RDFT Inverse;

   template<int_t N, typename NFact, typename VType,
            typename Parall, typename Place, typename Twiddles>
   class Algorithm {
      typedef typename VType::ValueType T;
      typedef Backward<N,T> Direction;
      typedef Separate<N,VType,Direction::Sign> Separator;
      typedef typename Place::template List<N,NFact,VType,Parall,Direction,Twiddles>::Result TList;
   public:
      typedef Loki::Typelist<Separator,TList> Result;
   };
};

/*! \brief Forward discrete cosine transform, type 1
\ingroup gr_params
*/
struct DCT1 {
   static const id_t ID = 4;
   typedef IDCT1 Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Forward<N,T> {};

   template<int_t N, typename NFact, typename VType,
            class Parall, class Place>
   struct Algorithm {
//      typedef TList Result;
   };
};

/*! \brief Inverse discrete cosine transform, type 1
\ingroup gr_params
*/
struct IDCT1 {
   static const id_t ID = 5;
   typedef DCT1 Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Backward<N,T> {};

   template<int_t N, typename NFact, typename VType,
            class Parall, class Place>
   struct Algorithm {
//      typedef TList Result;
   };
};

/*! \brief Forward discrete cosine transform, type 2
\ingroup gr_params
*/
struct DCT2 {
   static const id_t ID = 6;
   typedef IDCT2 Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Forward<N,T> {};

   template<int_t N, typename NFact, typename VType,
            class Parall, class Place>
   struct Algorithm {
      //typedef TList Result;
   };
};

/*! \brief Inverse discrete cosine transform, type 2
\ingroup gr_params
*/
struct IDCT2 {
   static const id_t ID = 7;
   typedef DCT2 Inverse;

   template<unsigned long N, typename T>
   struct Direction : public Backward<N,T> {};

   template<int_t N, typename NFact, typename VType,
            class Parall, class Place>
   struct Algorithm {
      //typedef TList Result;
   };
};

/*! \brief %Serial (single-core) implementation of transform
\sa OpenMP
\ingroup gr_params
*/
struct Serial {
   static const id_t ID = 0;
   static const uint_t NParProc = 1;

   template<int_t N>
   struct ActualParall {
      typedef Serial Result;
   };

   // used for in-place transforms only
   template<typename NFact, typename T>
   struct Swap {
      static const uint_t M = NFact::Head::first::value;
      static const uint_t P = NFact::Head::second::value;
      typedef GFFTswap2<M,P,T> Result;
   };

   template<typename N>
   struct Factor : public Factorization<N, SInt> {};

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
   static const uint_t NParProc = NT;

   template<int_t N>
   struct ActualParall {
      static const bool C = ((N > NT*NT) && (N >= SwitchToOMP));
      typedef typename Loki::Select<C,OpenMP<NT>,Serial>::Result Result;
   };
   
   // used for in-place transforms only
   template<typename NFact, typename T>
   struct Swap {
      static const uint_t M = NFact::Tail::Head::first::value;
      static const uint_t P = IsMultipleOf<NFact::Head::first::value,M>::value 
                            + NFact::Tail::Head::second::value;
      typedef GFFTswap2<M,P,T> Result;
//       typedef GFFTswap2OMP<NT,M,P,T> Result;
   };

   template<typename N>
   struct Factor {
      static const int_t G = GCD<SInt<N::value>, SInt<NT> >::Result::value;
      typedef typename Factorization<SIntID<N::value/G>, SInt>::Result NFact1;
      typedef ExtractFactor<NT/G, NFact1> EF;
      typedef Pair<SInt<G*EF::value>,SInt<1> > NParall;
      typedef Loki::Typelist<NParall,typename EF::Result> Multithreaded;
      typedef typename Factorization<N, SInt>::Result Singlethreaded;
      static const bool C = ((N::value > NT*NT) && (N::value >= SwitchToOMP));
      typedef typename Loki::Select<C,   // Condition to turn on multithreaded mode
	  Multithreaded, Singlethreaded>::Result Result;
   };
   
   template<typename T>
   void apply(T*) {
      //omp_set_dynamic(0);
      //omp_set_num_threads(NT);
      //omp_set_nested(true);
      //omp_set_max_active_levels(static_cast<int>(log(NT)/log(2))+1);
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

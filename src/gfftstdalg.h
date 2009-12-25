/***************************************************************************
 *   Copyright (C) 2008-2009 by Volodymyr Myrnyy                           *
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
    \brief General algorithms and short-radix FFT specializations
           using STL class complex<>
*/

#include "metafunc.h"

namespace GFFT {

using namespace MF;

/// Specialization for complex-valued radix 2 FFT in-place
/// \tparam T is value type
/// \tparam Complex<T> is a generic type representing complex numbers (like std::complex)
/// \param data is the array containing two complex numbers of type Complex<T>.
template<typename T,
template<typename> class Complex>
inline void _spec2(Complex<T>* data) {
      Complex<T> t(data[1]);
      data[1] = data[0]-t;
      data[0] += t;
}

/// Danielson-Lanczos section of the decimation-in-time FFT version
/**
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam Complex<T> is a generic type representing complex numbers (like std::complex)

This template class implements resursive metaprogram, which
runs funciton apply() twice recursively at the beginning of the function apply()
with the half of the transform length N
until the simplest case N=2 has been reached. Then function \a _spec2 is called.
Therefore, it has two specializations for N=2 and N=1 (the trivial and empty case).
*/
template<unsigned long N, typename T, int S,
template<typename> class Complex>
class InTime<N,Complex<T>,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const unsigned long N2 = N/2;
   InTime<N2,Complex<T>,S> next;
public:
   void apply(Complex<T>* data) {
      next.apply(data);
      next.apply(data+N2);

      LocalVType wtemp = Sin<N,1,LocalVType>::value();
      Complex<T> temp;
      Complex<LocalVType> w(1.,0.);
      Complex<LocalVType> wp(-2.0*wtemp*wtemp,-S*Sin<N,2,LocalVType>::value());

      for (unsigned long i=0; i<N2; ++i) {
        // rewritten componentwise because of the different types of the components
        temp = Complex<T>(static_cast<T>(data[i+N2].real()*w.real() - data[i+N2].imag()*w.imag()),
                          static_cast<T>(data[i+N2].real()*w.imag() + data[i+N2].imag()*w.real()));
        //temp = data[i+N2]*Complex<T>(w);
        data[i+N2] = data[i]-temp;
        data[i] += temp;

        w += w*wp;
      }
   }
};

// Specialization for N=4, decimation-in-time
// template<typename T, int S,
// template<typename> class Complex>
// class InTime<4,Complex<T>,S> {
// public:
//    void apply(Complex<T>* data) {
//       Complex<T> t(data[1]);
//       data[1] = data[0]-t;
//       data[0] += t;
//       t = data[3];
//       data[3] = Complex<T>(S*(data[2].imag()-t.imag()),S*(t.real()-data[2].real()));
//       data[2] += t;
// 
//       t = data[2];
//       data[2] = data[0]-t;
//       data[0] += t;
//       t = data[3];
//       data[3] = data[1]-t;
//       data[1] += t;
//    }
// };

// Specialization for N=2, decimation-in-time
template<typename T, int S,
template<typename> class Complex>
class InTime<2,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { _spec2(data); }
};

// Specialization for N=1, decimation-in-time
template<typename T, int S,
template<typename> class Complex>
class InTime<1,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { }
};

/// Danielson-Lanczos section of the decimation-in-frequency FFT version
/**
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam Complex<T> is a generic type representing complex numbers (like std::complex)

This template class implements resursive metaprogram, which
runs funciton apply() twice recursively at the end of the function apply()
with the half of the transform length N
until the simplest case N=2 has been reached. Then function \a _spec2 is called.
Therefore, it has two specializations for N=2 and N=1 (the trivial and empty case).
*/
template<unsigned long N, typename T, int S,
template<typename> class Complex>
class InFreq<N,Complex<T>,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const unsigned long N2 = N/2;
   InFreq<N2,Complex<T>,S> next;
public:
   void apply(Complex<T>* data) {

      LocalVType wtemp = Sin<N,1,LocalVType>::value();
      Complex<T> temp;
      Complex<LocalVType> w(1.,0.);
      Complex<LocalVType> wp(-2.0*wtemp*wtemp,-S*Sin<N,2,LocalVType>::value());

      for (unsigned long i=0; i<N2; ++i) {
        temp = data[i]-data[i+N2];
        data[i] += data[i+N2];
        // rewritten componentwise because of the different types of the components
	data[i+N2] = Complex<T>(static_cast<T>(temp.real()*w.real() - temp.imag()*w.imag()),
	                        static_cast<T>(temp.imag()*w.real() + temp.real()*w.imag()));
        //data[i+N2] = temp*w;

        w += w*wp;
      }

      next.apply(data);
      next.apply(data+N2);
   }
};

// Specialization for N=4, decimation-in-frequency
template<typename T, int S,
template<typename> class Complex>
class InFreq<2,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { _spec2(data); }
};

// Specialization for N=1, decimation-in-frequency
template<typename T, int S,
template<typename> class Complex>
class InFreq<1,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { }
};

/// Binary reordering of array elements
/*!
\tparam N length of the data
\tparam T value type
*/
template<unsigned long N, typename T,
template<typename> class Complex>
class GFFTswap<N,Complex<T> > {
public:
   void apply(Complex<T>* data) {
     unsigned long m,j=0;
     for (unsigned long i=0; i<N; ++i) {
        if (j>i) {
            std::swap(data[j], data[i]);
        }
        m = N/2;
        while (m>=1 && j>=m) {
            j -= m;
            m >>= 1;
        }
        j += m;
     }
   }
};

template<unsigned int P, typename T,
template<typename> class Complex,
unsigned int I>
class GFFTswap2<P,Complex<T>,I> {
   static const unsigned long BN = 1<<I;
   static const unsigned long BR = 1<<(P-I-1);
   GFFTswap2<P,Complex<T>,I+1> next;
public:
   void apply(Complex<T>* data, unsigned long n=0, unsigned long r=0) {
      next.apply(data,n,r);
      next.apply(data,n|BN,r|BR);
   }
};

template<unsigned int P, typename T,
template<typename> class Complex>
class GFFTswap2<P,Complex<T>,P> {
public:
   void apply(Complex<T>* data, unsigned long n=0, unsigned long r=0) {
      if (n>r)
        swap(data[n],data[r]);
   }
};

template<unsigned int NThreads, unsigned int P, typename T,
template<typename> class Complex, unsigned int I>
class GFFTswap2OMP<NThreads,P,Complex<T>,I,true> {
   static const unsigned long BN = 1<<I;
   static const unsigned long BR = 1<<(P-I-1);
   GFFTswap2OMP<NThreads/2,P,Complex<T>,I+1> next;
public:
   void apply(Complex<T>* data, const unsigned long n=0, const unsigned long r=0) {
     #pragma omp parallel shared(data)
     {
       #pragma omp sections
       {
         #pragma omp section
         next.apply(data,n,r);

         #pragma omp section
         next.apply(data,n|BN,r|BR);
       }
     }
   }
};

template<unsigned int NThreads, unsigned int P, typename T,
template<typename> class Complex>
class GFFTswap2OMP<NThreads,P,Complex<T>,P,true> {
public:
   void apply(Complex<T>* data, const unsigned long n, const unsigned long r) {
      if (n>r)
        swap(data[n],data[r]);
   }
};

template<unsigned int P, typename T, unsigned int I,
template<typename> class Complex>
class GFFTswap2OMP<1,P,Complex<T>,I,true> : public GFFTswap2<P,Complex<T>,I> { };

template<unsigned int P, typename T,
template<typename> class Complex>
class GFFTswap2OMP<1,P,Complex<T>,P,true> : public GFFTswap2<P,Complex<T>,P> { };

template<unsigned int NThreads, unsigned int P, typename T, unsigned int I,
template<typename> class Complex>
class GFFTswap2OMP<NThreads,P,Complex<T>,I,false> : public GFFTswap2<P,Complex<T>,I> { };

/// Reordering of data for real-valued transforms
/*!
\tparam N length of the data
\tparam T value type
\tparam S sign of the transform: 1 - forward, -1 - backward
*/
template<unsigned long N, typename T, int S,
template<typename> class Complex>
class Separate<N,Complex<T>,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Complex<LocalVType> LocalComplex;
   static const int M = (S==1) ? 2 : 1;
public:
   void apply(Complex<T>* data) {
      unsigned long i,i1;
      LocalComplex h1,h2,h3;
      LocalVType wtemp = Sin<2*N,1,LocalVType>::value();
      LocalComplex wp(-2.*wtemp*wtemp,-S*Sin<N,1,LocalVType>::value());
      LocalComplex w(1.+wp.real(),wp.imag());

      for (i=1; i<N/2; ++i) {
        i1 = N-i;
        h1 = Complex<LocalVType>(static_cast<LocalVType>(0.5*(data[i].real()+data[i1].real())),
                                 static_cast<LocalVType>(0.5*(data[i].imag()-data[i1].imag())));
        h2 = Complex<LocalVType>(static_cast<LocalVType>( S*0.5*(data[i].imag()+data[i1].imag())),
                                 static_cast<LocalVType>(-S*0.5*(data[i].real()-data[i1].real())));
        h3 = w*h2;
        data[i] = h1 + h3;
        data[i1]= h1 - h3;
        data[i1] = Complex<T>(data[i1].real(), -data[i1].imag());

        w += w*wp;
      }
      wtemp = data[0].real();
      data[0] = Complex<T>(M*0.5*(wtemp + data[0].imag()), M*0.5*(wtemp - data[0].imag()));

      data[N/2] = Complex<T>(data[N/2].real(), -data[N/2].imag());
   }
};

template<unsigned int NThreads, unsigned long N, typename T, int S,
template<typename> class Complex>
class InTimeOMP<NThreads,N,Complex<T>,S,true> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const unsigned long N2 = N/2;
   InTimeOMP<NThreads/2,N2,Complex<T>,S> next;
public:
   void apply(Complex<T>* data) {

      #pragma omp parallel shared(data)
      {
        #pragma omp sections
        {
          #pragma omp section
          next.apply(data);

          #pragma omp section
          next.apply(data+N2);
        }
      }

      Complex<T> temp;
      //Complex<LocalVType> w, wp;
      LocalVType wtemp = Sin<N,1,LocalVType>::value();
      Complex<LocalVType> w(1.,0.);
      Complex<LocalVType> wp(-2.0*wtemp*wtemp, -S*Sin<N,2,LocalVType>::value());

      for (unsigned long i=0; i<N2; ++i) {
        // rewritten componentwise because of the different types of the components
        temp = Complex<T>(data[i+N2].real()*w.real() - data[i+N2].imag()*w.imag(),
                          data[i+N2].real()*w.imag() + data[i+N2].imag()*w.real());
        //temp = data[i+N2]*Complex<T>(w);
        data[i+N2] = data[i]-temp;
        data[i] += temp;

        w += w*wp;
      }
   }
};

template<unsigned long N, typename T, int S,
template<typename> class Complex>
class InTimeOMP<1,N,Complex<T>,S,true> : public InTime<N,Complex<T>,S> { };

template<unsigned int NThreads, unsigned long N, typename T, int S,
template<typename> class Complex>
class InTimeOMP<NThreads,N,Complex<T>,S,false> : public InTime<N,Complex<T>,S> { };



template<unsigned int NThreads, unsigned long N, typename T, int S,
template<typename> class Complex>
class InFreqOMP<NThreads,N,Complex<T>,S,true> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const unsigned long N2 = N/2;
   InFreqOMP<NThreads/2,N2,Complex<T>,S> next;
public:
   void apply(Complex<T>* data) {

      LocalVType wtemp;
      Complex<T> temp;

      wtemp = Sin<N,1,LocalVType>::value();
      Complex<LocalVType> w(1.,0.);
      Complex<LocalVType> wp(-2.0*wtemp*wtemp,-S*Sin<N,2,LocalVType>::value());
      for (unsigned long i=0; i<N2; ++i) {
        temp = data[i]-data[i+N2];
        data[i] += data[i+N2];
        // rewritten componentwise because of the different types of the components
	data[i+N2] = Complex<T>(temp.real()*w.real() - temp.imag()*w.imag(),
	                        temp.imag()*w.real() + temp.real()*w.imag());
        //data[i+N2] = temp*w;

        w += w*wp;
      }

      #pragma omp parallel shared(data)
      {
      #pragma omp sections
      {
        #pragma omp section
        next.apply(data);
        #pragma omp section
        next.apply(data+N2);
      }
      } // parallel
   }
};

template<unsigned long N, typename T, int S,
template<typename> class Complex>
class InFreqOMP<1,N,Complex<T>,S,true> : public InFreq<N,Complex<T>,S> { };

template<unsigned int NThreads, unsigned long N, typename T, int S,
template<typename> class Complex>
class InFreqOMP<NThreads,N,Complex<T>,S,false> : public InFreq<N,Complex<T>,S> { };


} //namespace

#endif

/***************************************************************************
 *   Copyright (C) 2014-2015 by Vladimir Mirnyy                            *
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

#ifndef __gff_accuracy_h
#define __gff_accuracy_h

/** \file
    \brief Check accuracy of GFFT results comparing to FFTW and DFT of double-double precision
*/

#include <iostream>

#include "gfft.h"
#include "direct.h"
#include "thirdparty.h"

namespace GFFT {

template<typename T>
T norm_inf(const T* data, const unsigned int n) {
   if (n<1) return 0.;
   T d=fabs(data[0]);
   for (unsigned int i=1; i<n; ++i) {
     if (fabs(data[i])>d) d=fabs(data[i]);
   }
   return d;
}

template<typename T>
T norm2(const T* data, const unsigned int n) {
   T s=0.;
   for (unsigned int i=0; i<n; ++i) {
     s += data[i]*data[i];
   }
   return sqrt(s);
}

// template<typename T, template<typename> class Complex>
// Complex<T> norm_inf(const Complex<T>* data, const unsigned int n) {
//    if (n<1) return 0.;
//    Complex<T> d(fabs(data[0].real()),fabs(data[0].imag()));
//    for (unsigned int i=1; i<n; ++i) {
//      if (fabs(data[i].real()) > d.real()) d.real() = fabs(data[i].real());
//      if (fabs(data[i].imag()) > d.imag()) d.imag() = fabs(data[i].imag());
//    }
//    return d;
// }
// 
// template<typename T, template<typename> class Complex>
// Complex<T> norm2(const Complex<T>* data, const unsigned int n) {
//    Complex<T> s(0.,0.);
//    for (unsigned int i=0; i<n; ++i) {
//      s.real() += data[i].real()*data[i].real();
//      s.imag() += data[i].imag()*data[i].imag();
//    }
//    return Complex<T>(sqrt(s.real()), sqrt(s.imag()));
// }

template<typename T, template<typename> class Complex>
T norm_inf(const Complex<T>* data, const unsigned int n) {
   if (n<1) return 0.;
   T d(std::max(fabs(data[0].real()),fabs(data[0].imag())));
   for (unsigned int i=1; i<n; ++i) {
     if (fabs(data[i].real()) > d) d = fabs(data[i].real());
     if (fabs(data[i].imag()) > d) d = fabs(data[i].imag());
   }
   return d;
}

template<typename T, template<typename> class Complex>
T norm2(const Complex<T>* data, const unsigned int n) {
   T s(0.);
   for (unsigned int i=0; i<n; ++i) {
     s += data[i].real()*data[i].real();
     s += data[i].imag()*data[i].imag();
   }
   return sqrt(s);
}

/////////////////////////////////////////////////////

template<class FFT1, class FFT2>
class FFTcompare 
{
  typedef typename FFT1::value_type T1;
  typedef typename FFT2::value_type T2;

public:
  void apply(long_t pmin, long_t pmax)
  {
    long_t i;

    double *data = new double [(1 << pmax)*2];

    srand(StartSeed);
    
    for (long_t p = pmin; p <= pmax; ++p)
    {
      long_t n = 1 << p;
      
      for (i=0; i < n; ++i) {
	data[2*i] = rand()/(double)RAND_MAX - 0.5;  // distribute in [-0.5;0.5] as in FFTW
	data[2*i+1] = rand()/(double)RAND_MAX - 0.5;
      }
//       for (i=0; i < n; ++i) {
//          data[2*i] = 2*i;
//          data[2*i+1] = 2*i+1; 
//       }
      
      FFT1 fft1(data, n);
      FFT2 fft2(data, n);

      fft1.apply();
      fft2.apply();

      T1* out = fft1.getdata();

      const int n2 = n+n;
      
      T1 d = norm_inf(out, n2);
      
      fft2.diff(out);
      
      T1 nr2 = norm2(out, n2);
      T1 nrinf = norm_inf(out, n2);
      std::cout << n << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << std::endl;
    }
    
    delete [] data;
  }
};

//============================================================

static double MaxRelError = 0;

template<class NList, class DFTClass, class Place>
class GFFTcheck;

template<class H, class Tail, class DFTClass>
class GFFTcheck<Loki::Typelist<H,Tail>, DFTClass, IN_PLACE> 
{
  typedef typename H::ValueType::ValueType T1;
  typedef typename H::ValueType::base_type BT;
  typedef typename DFTClass::value_type T2;
  GFFTcheck<Tail,DFTClass,IN_PLACE> next;

  static const int C = Loki::TypeTraits<T1>::isStdFundamental ? 2 : 1;
  static const long_t N = H::Len;
  static const long_t N2 = N*C;
  
  typename H::Instance gfft;
  
public:
  void apply() 
  {
    next.apply();
    
    long_t i;

    T1 *data = new T1 [N2];

    for (i=0; i < N; ++i) 
      GenInput<T1>::rand(data, i);  // distribute in [-0.5;0.5] as in FFTW
      
    // Another input dataset
//     for (i=0; i < N; ++i) {
//        GenInput<T1>::seq(data, i);
//     }
    
    DFTClass dft(data, N);

    // apply FFT in-place
    gfft.fft(data);
    
    // External DFT
    dft.apply();

    BT d = norm_inf(data, N2);

    // Subtract result of dft from data
    dft.diff(data);
    
    BT nr2 = norm2(data, N2);
    BT nrinf = norm_inf(data, N2);
#ifdef FOUT
    std::cout << N << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << std::endl;
#endif
    delete [] data;
    
    nrinf /= d;
    if (MaxRelError < nrinf) MaxRelError = nrinf;
  }
};

template<class H, class Tail, class DFTClass>
class GFFTcheck<Loki::Typelist<H,Tail>, DFTClass, OUT_OF_PLACE> {
  typedef typename H::ValueType::ValueType T1;
  typedef typename H::ValueType::base_type BT;
  typedef typename DFTClass::value_type T2;
  GFFTcheck<Tail,DFTClass,OUT_OF_PLACE> next;

  static const int C = Loki::TypeTraits<T1>::isStdFundamental ? 2 : 1;
  static const long_t N = H::Len;
  static const long_t N2 = N*C;
  
  typename H::Instance gfft;
  
public:
  void apply() 
  {
    next.apply();
    
    long_t i;

    T1 *data = new T1 [N2];
    T1 *dataout = new T1 [N2];

    for (i=0; i < N2; ++i) 
      dataout[i] = 0.;
    
    for (i=0; i < N; ++i) 
      GenInput<T1>::rand(data, i);  // distribute in [-0.5;0.5] as in FFTW
      
    // Another input dataset
//     for (i=0; i < N; ++i) {
//        GenInput<T1>::seq(data, i);
//     }
    
    DFTClass dft(data, N);

    // apply FFT out-of-place
    gfft.fft(data, dataout);
    
    // External DFT
    dft.apply();

    BT d = norm_inf(dataout, N2);

    // Subtract result of dft from dataout
    dft.diff(dataout);
    
    BT nr2 = norm2(dataout, N2);
    BT nrinf = norm_inf(dataout, N2);
#ifdef FOUT
    std::cout << N << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << std::endl;
#endif
    delete [] dataout;
    delete [] data;

    nrinf /= d;
    if (MaxRelError < nrinf) MaxRelError = nrinf;
  }
};

template<class DFTClass>
class GFFTcheck<Loki::NullType, DFTClass, IN_PLACE> {
public:
  void apply() { }
};

template<class DFTClass>
class GFFTcheck<Loki::NullType, DFTClass, OUT_OF_PLACE> {
public:
  void apply() { }
};

} // namespace GFFT

#endif

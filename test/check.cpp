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
    \brief Check GFFT results comparing to FFTW and NR
*/

#include <iostream>

#ifdef FFTW
#include <fftw3.h>
#include <boost/iterator/iterator_concepts.hpp>
#endif

#include "gfft.h"
#include "nrfft.h"

#include <arprec/mp_real.h>


using namespace std;

using namespace GFFT;

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


template<class T>
class DFT_wrapper 
{
  void dft1(T* output_data, const T* input_data, const int_t size, bool inverse)
  {
    T pi2 = (inverse) ? 2.0 * M_PI : -2.0 * M_PI;
    T a, ca, sa;
    T invs = 1.0 / size;
    for(unsigned int y = 0; y < size; y++) {
      output_data[2*y] = 0;
      output_data[2*y+1] = 0;
      for(unsigned int x = 0; x < size; x++) {
	a = pi2 * y * x * invs;
	ca = cos(a);
	sa = sin(a);
	output_data[2*y]   += input_data[2*x] * ca - input_data[2*x+1] * sa;
	output_data[2*y+1] += input_data[2*x] * sa + input_data[2*x+1] * ca;
      }
      if(inverse) {
	output_data[2*y]   *= invs;
	output_data[2*y+1] *= invs;
      }
    }
  }
  
  T* input_data;
  T* output_data;
  int_t size;
  
public:
  
  DFT_wrapper(const T* data, int_t n) : size(n)
  {
    init(data, n);
  }
  ~DFT_wrapper() 
  {
    delete [] input_data;
    delete [] output_data;
  }
  
  void init(const T* data, int_t n)
  {
    input_data = new T [n*2];
    output_data = new T [n*2];
   
    for (int_t i=0; i < n*2; ++i) {
       input_data[i] = data[i];
    }
  }
  
  void apply()
  {
    dft1(output_data, input_data, size, false);
  }
  
  void diff(T* data)
  {
    for (int_t i=0; i<size*2; ++i) 
      data[i] -= output_data[i];
  }
};


#ifdef FFTW
template<class T>
class FFTW_wrapper 
{
    fftw_complex* in;
    fftw_complex* out;
    fftw_plan plan;
    int_t size;
public:
  FFTW_wrapper(const T* data, int_t n) : size(n)
  {
    init(data, n);
  }
  ~FFTW_wrapper() 
  {
    fftw_free(in);
    fftw_free(out);
  }
  
  void init(const T* data, int_t n)
  {
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
   
    for (int_t i=0; i < n; ++i) {
       in[i][0] = data[2*i];
       in[i][1] = data[2*i+1];
    }
  }
  
  void apply()
  {
    plan = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
  }
  
  void diff(T* data)
  {
    for (int_t i=0; i<size; ++i) {
      data[2*i]   -= out[i][0];
      data[2*i+1] -= out[i][1];
    }
  }
};
#endif

template<class T>
class Arprec_wrapper 
{
    T* y;
    T* in;
    int_t size;
    int m, m1, m2;
public:
  Arprec_wrapper(const T* data, int_t n) : size(n)
  {
    init(data, n);
  }
  ~Arprec_wrapper() 
  {
    delete [] in;
    delete [] y;
  }
  
  void init(const T* data, int_t n)
  {
     in = new T [n * 2];

     m = log(n)/log(2);
     m1 = ((m + 1) / 2);
     m2 = (m - (m + 1) / 2);
     int n1 = 1 << m1;
//      y = new T [(n + n1 * mp::mpnsp1) * 2];
     y = new T [n * 4];
     mp_real::mpinix(n);
   
     for (int_t i=0; i < 2*n; ++i) {
       in[i] = data[i];
     }
  }
  
  void apply()
  {
    mp_real::mpfft1(1, m, m1, m2, in, y);
  }
  
  void diff(T* data)
  {
    for (int_t i=0; i<2*size; ++i) {
      data[i]   -= y[i];
    }
  }
};

//==============================================================

template<class NList, template<class> class FFTWrapper>
class GFFTcheck;

template<class H, class T, template<class> class FFTWrapper>
class GFFTcheck<Loki::Typelist<H,T>, FFTWrapper> {
  typedef typename H::ValueType::ValueType Tp;
  GFFTcheck<T,FFTWrapper> next;

  static const int_t N = H::Len;
  static const int_t N2 = 2*N;
  
  H gfft;
  
public:
  void apply() 
  {
    next.apply();
    
    int_t i;
    double d, nr2, nrinf;

    Tp *data = new Tp [N2];

    srand(17);
    
    for (i=0; i < N; ++i) {
	data[2*i] = rand()/(Tp)RAND_MAX - 0.5;  // distribute in [-0.5;0.5] as in FFTW
	data[2*i+1] = rand()/(Tp)RAND_MAX - 0.5;
    }

    FFTWrapper<Tp> dft(data, N);

// apply FFT in-place
    gfft.fft(data);
    
    dft.apply();

    d = norm_inf(data, N2);
    
    dft.diff(data);
    
    nr2 = norm2(data, N2);
    nrinf = norm_inf(data, N2);
    cout << N << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << endl;

    delete [] data;
  }

};

template<template<class> class FFTWrapper>
class GFFTcheck<Loki::NullType, FFTWrapper> {
public:
  void apply() { }
};

typedef DOUBLE VType;
typedef IN_PLACE Place;
//typedef OUT_OF_PLACE Place;

const unsigned Min = 2;
const unsigned Max = 10;

typedef TYPELIST_3(OpenMP<2>, OpenMP<3>, OpenMP<4>) ParallList;
//typedef GenNumList<2, 10, SIntID>::Result NList;
typedef GenPowerList<Min, Max, 2>::Result NList;
typedef GenerateTransform<NList, VType, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::Default, Place> Trans;



int main(int argc, char *argv[])
{
  //B<C,A> a;
  cout << "Forward transforms in place:" << endl;
  //GFFTcheck<typename Trans::Result, DFT_wrapper> check;
//   GFFTcheck<typename Trans::Result, FFTW_wrapper> check;
  GFFTcheck<typename Trans::Result, Arprec_wrapper> check;
  check.apply();

  
  return 0;
}


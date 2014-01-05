/***************************************************************************
 *   Copyright (C) 2014 by Vladimir Mirnyy                                 *
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
//#include <boost/iterator/iterator_concepts.hpp>
#endif

//#include "nrfft.h"

#include <arprec/mp_real.h>


template<class T>
class DFT_wrapper 
{
  void dft1(T* output_data, const T* input_data, const int_t size, bool inverse)
  {
    T pi2 = (inverse) ? 2.0 * M_PI : -2.0 * M_PI;
    T a, ca, sa;
    T invs = 1.0 / size;
    for(unsigned int y = 0; y < size; y++) {
      output_data[2*y] = 0.;
      output_data[2*y+1] = 0.;
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
  
  typedef T value_type;
  
  DFT_wrapper(const T* data, int_t n) : size(n)
  {
    init(data, n);
  }
  DFT_wrapper(const std::complex<double>* data, int_t n) : size(n)
  {
    init(data, n);
  }
  ~DFT_wrapper() 
  {
    delete [] input_data;
    delete [] output_data;
  }
  
  T* getdata() const { return output_data; }
  
  void init(const T* data, int_t n)
  {
    input_data = new T [n*2];
    output_data = new T [n*2];
   
    for (int_t i=0; i < n*2; ++i) {
       input_data[i] = data[i];
    }
  }
  void init(const std::complex<double>* data, int_t n)
  {
    input_data = new T [n*2];
    output_data = new T [n*2];
   
    for (int_t i=0; i < n; ++i) {
       input_data[2*i] = data[i].real();
       input_data[2*i+1] = data[i].imag();
    }
  }
  
  void apply()
  {
    dft1(output_data, input_data, size, false);
  }
  
  template<class Tp>
  void diff(Tp* data)
  {
    for (int_t i=0; i<size*2; ++i) 
      data[i] -= output_data[i];
  }
  template<class Tp>
  void diff(std::complex<Tp>* data)
  {
    for (int_t i=0; i<size; ++i) {
      data[i].real() -= output_data[2*i];
      data[i].imag() -= output_data[2*i+1];
    }
  }
};


#ifdef FFTW
template<class T>
class FFTW_wrapper 
{
    T* in;
    T* out;
    fftw_plan plan;
    int_t size;
public:

  typedef T value_type;

  FFTW_wrapper(const double* data, int_t n) : size(n)
  {
    init(data, n);
  }
  ~FFTW_wrapper() 
  {
    fftw_free(in);
    fftw_free(out);
  }
  
  void init(const double* data, int_t n)
  {
    in = (T*)fftw_malloc(sizeof(T)*n);
    out = (T*)fftw_malloc(sizeof(T)*n);
   
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
  
  template<typename OtherType>
  void diff(OtherType* data)
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
  
  T* getdata() const { return in; }
  
  void init(const T* data, int_t n)
  {
     in = new T [n * 2];

     m = log(n)/log(2);
     m1 = ((m + 1) / 2);
     m2 = (m - (m + 1) / 2);
     int n1 = 1 << m1;
     y = new T [(n + n1 * mp::mpnsp1) * 2];
//     y = new T [n * 4];
     mp_real::mpinix(n);
   
// print out sample data
//     cout<<"Input data:"<<endl;
//     for (int i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

     for (int_t i=0; i < 2*n; ++i) {
       in[i] = data[i];
     }
  }
  
  void apply()
  {
    mp_real::mpfft1(-1, m, m1, m2, in, y);

    // print out sample data
//     cout<<"Output data:"<<endl;
//     for (int i=0; i < size; ++i)
//       cout<<"("<<in[2*i]<<","<<in[2*i+1]<<")"<<endl;
  }
  
  void diff(T* data)
  {
    for (int_t i=0; i<2*size; ++i) {
      data[i]   -= in[i];
    }
  }
};


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
#endif

#ifdef QD
#include "qd/dd_real.h"
#include "qd/qd_real.h"
#endif

double to_double(const double d) { return d; }


template<class T>
class DFT_wrapper 
{
  void dft1(T* output_data, const T* input_data, const int_t size, bool inverse)
  {
    T pi2 = (inverse) ? 2.0 * M_PI : -2.0 * M_PI;
    T a, ca, sa;
    T invs = 1.0 / size;
    for(int_t y = 0; y < size; y++) {
      output_data[2*y] = 0.;
      output_data[2*y+1] = 0.;
      for(int_t x = 0; x < size; x++) {
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
  
  DFT_wrapper(const double* data, int_t n) : size(n)
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
  
  template<typename Tp>
  void init(const Tp* data, int_t n)
  {
    input_data = new T [n*2];
    output_data = new T [n*2];
   
    for (int_t i=0; i < n*2; ++i) {
       input_data[i] = data[i];
    }
  }
  template<typename Tp>
  void init(const std::complex<Tp>* data, int_t n)
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
      data[i] -= to_double(output_data[i]);
  }
  template<class Tp>
  void diff(std::complex<Tp>* data)
  {
    for (int_t i=0; i<size; ++i) {
      data[i] -= std::complex<Tp>(output_data[2*i], output_data[2*i+1]);
    }
  }
  
  T norm2() 
  { 
    T s = 0.;
    for (int_t i=0; i<size*2; ++i) {
      s += output_data[i]*output_data[i];
    }
    return sqrt(s);
  }
  
  T norm_inf() 
  {
    if (size<1) return 0.;
    T d=fabs(output_data[0]);
    for (int_t i=1; i<size*2; ++i) {
      if (fabs(output_data[i])>d) d=fabs(output_data[i]);
    }
    return d;
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
  
  T* getdata() const { return out; }

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
  
  template<typename AnotherType>
  void diff(AnotherType* data)
  {
    for (int_t i=0; i<size; ++i) {
      data[2*i]   -= out[i][0];
      data[2*i+1] -= out[i][1];
    }
  }
  
  double norm2() 
  {
    double s = 0.;
    for (int_t i=0; i<size; ++i) {
      s += out[i][0]*out[i][0];
      s += out[i][1]*out[i][1];
    }
    return sqrt(s);
  }

  double norm_inf() 
  {
    if (size<1) return 0.;
    double d=fabs(out[0][0]);
    if (fabs(out[0][1])>d) d=fabs(out[0][1]);
    for (int_t i=1; i<size; ++i) {
      if (fabs(out[i][0])>d) d=fabs(out[i][0]);
      if (fabs(out[i][1])>d) d=fabs(out[i][1]);
    }
    return d;
  }
};
#endif


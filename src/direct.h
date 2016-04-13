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

#ifndef __direct_h
#define __direct_h

/** \file
    \brief Check GFFT results comparing to simple DFT implementations
*/

#include <iostream>

#ifdef QD
// install libqd-dev to get these packages
#include "qd/dd_real.h"
#include "qd/qd_real.h"
#endif

double to_double(const double d) { return d; }

static const int StartSeed = 17;

template<typename T, bool C = Loki::TypeTraits<T>::isStdFundamental>
struct GenInput;

template<typename T>
struct GenInput<T,true>
{
  GenInput() { srand(17); }
  
  static void rand(T* data, const int i) 
  { // distribute in [-0.5;0.5] as in FFTW
    data[2*i]   = ::rand()/static_cast<T>(RAND_MAX) - 0.5; 
    data[2*i+1] = ::rand()/static_cast<T>(RAND_MAX) - 0.5; 
  }
  static void seq(T* data, const int i) 
  { 
    data[2*i]   = static_cast<T>(2*i); 
    data[2*i+1] = static_cast<T>(2*i+1); 
  }
};

template<typename T>
struct GenInput<T,false>
{
  GenInput() { srand(17); }
  
  static void rand(T* data, const int i) 
  { // distribute in [-0.5;0.5] as in FFTW
    T tmp(::rand()/static_cast<typename T::value_type>(RAND_MAX) - 0.5, 
	  ::rand()/static_cast<typename T::value_type>(RAND_MAX) - 0.5);
    data[i] = tmp;
  }
  static void seq(T* data, const int i) 
  { 
    T tmp(static_cast<typename T::value_type>(2*i), static_cast<typename T::value_type>(2*i+1)); 
    data[i] = tmp;
  }
};

/////////////////////////////////////////////////////

template<typename T, bool C = Loki::TypeTraits<T>::isStdFundamental>
struct GenOutput;

template<typename T>
struct GenOutput<T,true>
{
  const T re,im;
  GenOutput(const T* data, const int i) : re(data[2*i]), im(data[2*i+1]) {}
};

template<typename T>
struct GenOutput<T,false>
{
  const T c;
  GenOutput(const T* data, const int i) : c(data[i]) {}
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const GenOutput<T,true>& v)
{
  os << "(" << v.re << "," << v.im << ")";
  return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const GenOutput<T,false>& v)
{
  os << v.c;       // User-defined complex types should implement this operator themselves like std::complex
  return os;
}

/////////////////////////////////////////////////////

template<typename T, bool C = Loki::TypeTraits<T>::isStdFundamental>
struct ComplexWrapper;

template<typename T>
struct ComplexWrapper<T,true>
{
  T& re,im;
  ComplexWrapper(T* data, const int i) : re(data[2*i]), im(data[2*i+1]) {}

  T& real() { return re; }
  T& imag() { return im; }
  T real() const { return re; }
  T imag() const { return im; }
};

template<typename T>
struct ComplexWrapper<T,false>
{
  typedef typename T::value_type VT;
  T& c;
  ComplexWrapper(T* data, const int i) : c(data[i]) {}

  VT& real() { return c.real(); }
  VT& imag() { return c.imag(); }
  const VT& real() const { return c.real(); }
  const VT& imag() const { return c.imag(); }
};

/////////////////////////////////////////////////////

/** \class DFT_wrapper
\brief A wrapper class for a simple implementation of DFT 

This class contains a simple and straightforward implementation of DFT. 
It is very close to its definition using cos() and sin() trigonometric functions.
Therefore, it can take a while running this implementation for big a transform size.
It is used for a check of correctness and accuracy of GFFT.
*/
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
  DFT_wrapper(const float* data, int_t n) : size(n)
  {
    init(data, n);
  }
  DFT_wrapper(const std::complex<double>* data, int_t n) : size(n)
  {
    init(data, n);
  }
  DFT_wrapper(const std::complex<float>* data, int_t n) : size(n)
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
  
  template<typename Tp, template<typename> class Complex>
  void diff(Complex<Tp>* data)
  {
    for (int_t i=0; i<size; ++i) {
      data[i] -= Complex<Tp>(to_double(output_data[2*i]), to_double(output_data[2*i+1]));
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


struct DCT1 
{
  template<typename T>
  void operator()(T* output_data, const T* input_data, const int_t size, bool inverse)
  {
    T pi = (inverse) ? M_PI : -M_PI;
    T a, ca, sa;
    T invs = pi / (size - 1);
    for(int_t k = 0; k < size; k++) {
      output_data[k] = 0.;
      for(int_t l = 0; l < size; l++) {
	a = k * l * invs;
	ca = cos(a);
	output_data[k] += input_data[l] * ca;
      }
//       if(inverse) {
// 	output_data[k]   *= invs;
//       }
    }
  }
};

struct DCT2 
{
  template<typename T>
  void operator()(T* output_data, const T* input_data, const int_t size, bool inverse)
  {
    T pi = (inverse) ? M_PI : -M_PI;
    T a, ca, sa;
    T invs = pi / (T)size;
    for(int_t k = 0; k < size; k++) {
      output_data[k] = 0.;
      for(int_t l = 0; l < size; l++) {
	a = k * l * invs;
	ca = cos(a);
	output_data[k] += input_data[l] * ca;
      }
//       if(inverse) {
// 	output_data[k]   *= invs;
//       }
    }
  }
};

template<class T, class DCT>
class DCT_wrapper 
{
 
  T* input_data;
  T* output_data;
  int_t size;
  
public:
  
  typedef T value_type;
  
  DCT_wrapper(const double* data, int_t n) : size(n)
  {
    init(data, n);
  }
  ~DCT_wrapper() 
  {
    delete [] input_data;
    delete [] output_data;
  }
  
  T* getdata() const { return output_data; }
  
  template<typename Tp>
  void init(const Tp* data, int_t n)
  {
    input_data = new T [n];
    output_data = new T [n];
   
    for (int_t i=0; i < n; ++i) {
       input_data[i] = data[i];
    }
  }
  
  void apply()
  {
    dct1(output_data, input_data, size, false);
  }
  
  template<class Tp>
  void diff(Tp* data)
  {
    for (int_t i=0; i<size; ++i) 
      data[i] -= to_double(output_data[i]);
  }
  
  T norm2() 
  { 
    T s = 0.;
    for (int_t i=0; i<size; ++i) {
      s += output_data[i]*output_data[i];
    }
    return sqrt(s);
  }
  
  T norm_inf() 
  {
    if (size<1) return 0.;
    T d=fabs(output_data[0]);
    for (int_t i=1; i<size; ++i) {
      if (fabs(output_data[i])>d) d=fabs(output_data[i]);
    }
    return d;
  }
};

typedef DCT_wrapper<double,DCT1> DCT1d;

#endif

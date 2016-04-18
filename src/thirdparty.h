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

#ifndef __thirdparty_h
#define __thirdparty_h

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


#ifdef FFTW
/** \class FFTW_wrapper
\brief A wrapper class for FFTW library 

This class contains an interface for the 1D complex DFT from FFTW library.
*/
template<class T>
class FFTW_wrapper 
{
    T* in;
    T* out;
    fftw_plan plan;
    long_t size;
public:

  typedef T value_type;

  FFTW_wrapper(const double* data, long_t n) : size(n)
  {
    init(data, n);
  }
  ~FFTW_wrapper() 
  {
    fftw_free(in);
    fftw_free(out);
  }
  
  T* getdata() const { return out; }

  void init(const double* data, long_t n)
  {
    in = (T*)fftw_malloc(sizeof(T)*n);
    out = (T*)fftw_malloc(sizeof(T)*n);
   
    for (long_t i=0; i < n; ++i) {
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
    for (long_t i=0; i<size; ++i) {
      data[2*i]   -= out[i][0];
      data[2*i+1] -= out[i][1];
    }
  }
  
  double norm2() 
  {
    double s = 0.;
    for (long_t i=0; i<size; ++i) {
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
    for (long_t i=1; i<size; ++i) {
      if (fabs(out[i][0])>d) d=fabs(out[i][0]);
      if (fabs(out[i][1])>d) d=fabs(out[i][1]);
    }
    return d;
  }
};
#endif

#endif

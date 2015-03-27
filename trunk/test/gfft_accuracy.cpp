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
    \brief Check accuracy of GFFT results comparing to FFTW and DFT of double-double precision
*/

#include <iostream>

#include "gfft.h"
#include "direct.h"
#include "thirdparty.h"

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

//==============================================================

template<class FFT1, class FFT2>
class FFTcompare 
{
  typedef typename FFT1::value_type T1;
  typedef typename FFT2::value_type T2;

public:
  void apply(int_t pmin, int_t pmax) 
  {
    int_t i;

    double *data = new double [(1 << pmax)*2];

    srand(17);
    
    for (int_t p = pmin; p <= pmax; ++p) 
    {
      int_t n = 1 << p;
      
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
      cout << n << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << endl;
    }
    
    delete [] data;
  }
};

//============================================================

template<class NList, class DFTClass, class Place>
class GFFTcheck;

template<class H, class T, class DFTClass>
class GFFTcheck<Loki::Typelist<H,T>, DFTClass, IN_PLACE> 
{
  typedef typename H::ValueType::ValueType T1;
  typedef typename DFTClass::value_type T2;
  GFFTcheck<T,DFTClass,IN_PLACE> next;

  static const int_t N = H::Len;
  static const int_t N2 = 2*N;
  
  typename H::Instance gfft;
  
public:
  void apply() 
  {
    next.apply();
    
    int_t i;

    T1 *data = new T1 [N2];

    srand(17);
    
    for (i=0; i < N; ++i) {
      data[2*i] = rand()/(T1)RAND_MAX - 0.5;  // distribute in [-0.5;0.5] as in FFTW
      data[2*i+1] = rand()/(T1)RAND_MAX - 0.5;
    }
    // Another input dataset
//     for (i=0; i < N; ++i) {
//        data[2*i] = 2*i;
//        data[2*i+1] = 2*i+1; 
//     }
    
    DFTClass dft(data, N);

    // apply FFT in-place
    gfft.fft(data);
    
    // External DFT
    dft.apply();

    T1 d = norm_inf(data, N2);

    // Subtract result of dft from data
    dft.diff(data);
    
    T1 nr2 = norm2(data, N2);
    T1 nrinf = norm_inf(data, N2);
    cout << N << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << endl;

    delete [] data;
  }
};

template<class H, class T, class DFTClass>
class GFFTcheck<Loki::Typelist<H,T>, DFTClass, OUT_OF_PLACE> {
  typedef typename H::ValueType::ValueType T1;
  typedef typename DFTClass::value_type T2;
  GFFTcheck<T,DFTClass,OUT_OF_PLACE> next;

  static const int_t N = H::Len;
  static const int_t N2 = 2*N;
  
  typename H::Instance gfft;
  
public:
  void apply() 
  {
    next.apply();
    
    int_t i;

    T1 *data = new T1 [N2];
    T1 *dataout = new T1 [N2];

    srand(17);
    
    for (i=0; i < 2*N; ++i) 
      dataout[i] = 0.;
    
    for (i=0; i < N; ++i) {
      data[2*i] = rand()/(T1)RAND_MAX - 0.5;  // distribute in [-0.5;0.5] as in FFTW
      data[2*i+1] = rand()/(T1)RAND_MAX - 0.5;
    }
    // Another input dataset
//     for (i=0; i < N; ++i) {
//        data[2*i] = 2*i;
//        data[2*i+1] = 2*i+1; 
//     }
    
    DFTClass dft(data, N);

    // apply FFT out-of-place
    gfft.fft(data, dataout);
    
    // External DFT
    dft.apply();

    T1 d = norm_inf(dataout, N2);

    // Subtract result of dft from dataout
    dft.diff(dataout);
    
    T1 nr2 = norm2(dataout, N2);
    T1 nrinf = norm_inf(dataout, N2);
    cout << N << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << endl;

    delete [] dataout;
    delete [] data;
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


typedef DOUBLE VType;
typedef IN_PLACE Place;
//typedef OUT_OF_PLACE Place;

const unsigned Min = 2;
const unsigned Max = 10;

typedef TYPELIST_3(OpenMP<2>, OpenMP<4>, OpenMP<8>) ParallList;
//typedef GenNumList<2, 15, SIntID>::Result NList;
typedef GenPowerList<Min, Max, 2>::Result NList;
typedef GenerateTransform<NList, VType, TransformTypeGroup::Default, SIntID<1>, ParallList, Place> Trans;

ostream& operator<<(ostream& os, const dd_real& v)
{
  os << v.to_string(16);
  return os;
}

ostream& operator<<(ostream& os, const qd_real& v)
{
  os << v.to_string(16);
  return os;
}

void print_header() 
{
  cout<<"------------------------------------------------------------------------------"<<endl;
  cout<<" N              Norm2                   NormInf           Relative NormInf    "<<endl;
  cout<<"------------------------------------------------------------------------------"<<endl;
}


int main(int argc, char *argv[])
{
  unsigned int oldcw;
  fpu_fix_start(&oldcw);
  
  cout.precision(16);

//   cout << "double-double DFT vs. FFTW:" << endl;
//   FFTcompare<DFT_wrapper<dd_real>, FFTW_wrapper<fftw_complex> > comp;
//   print_header();
//   comp.apply(Min,Max);
//   cout<<"------------------------------------------------------------------------------"<<endl<<endl;
  
  cout << "double-double DFT vs. GFFT:" << endl;
  GFFTcheck<Trans::Result, DFT_wrapper<dd_real>, Place> check_dft;
  print_header();
  check_dft.apply();
  cout<<"------------------------------------------------------------------------------"<<endl<<endl;
  
//   cout << "GFFT vs. FFTW:" << endl;
//   GFFTcheck<Trans::Result, FFTW_wrapper<fftw_complex>, Place> check_fftw;
//   print_header();
//   check_fftw.apply();
//   cout<<"------------------------------------------------------------------------------"<<endl<<endl;

  fpu_fix_end(&oldcw);
  return 0;
}


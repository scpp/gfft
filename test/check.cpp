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
#endif

#include "gfft.h"
#include "nrfft.h"

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


void dft1(double* output_data, const double* input_data, const unsigned int size, bool inverse)
{
  double pi2 = (inverse) ? 2.0 * M_PI : -2.0 * M_PI;
  double a, ca, sa;
  double invs = 1.0 / size;
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

//==============================================================

template<class NList>
class GFFTcheck;

template<class H, class T>
class GFFTcheck<Loki::Typelist<H,T> > {
  typedef typename H::ValueType::ValueType Tp;
  GFFTcheck<T> next;

  static const int_t N = H::Len;
  static const int_t N2 = 2*N;
  
  H gfft;
  
public:
  void check_vs_dft() 
  {
    next.check_vs_dft();
    
    int_t i;
    double d, nr2, nrinf;

    srand(17);
    
    Tp *data;
    Tp *datain, *dataout;

// sample data
    data = new Tp [N2];
    datain  = new Tp [N2];
    dataout = new Tp [N2];
    
    for (i=0; i < N; ++i) {
       data[2*i] = rand()/(Tp)RAND_MAX - 0.5;  // distribute in [-0.5;0.5] as in FFTW
       data[2*i+1] = rand()/(Tp)RAND_MAX - 0.5;
       datain[2*i] = data[2*i];
       datain[2*i+1] = data[2*i+1];
    }

// apply FFT in-place
    gfft.fft(data);
    dft1(dataout, datain, N, false);
//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<nrdata[2*i]<<","<<data[2*i+1]<<")"<<endl;

    d = norm_inf(data, N2);
    for (i=0; i<N2; ++i) 
      data[i] -= dataout[i];
    nr2 = norm2(data, N2);
    nrinf = norm_inf(data, N2);
    cout << N << "\t" << nr2 << "\t" << nrinf << "\t" << nrinf/d << endl;

    delete [] data;
    delete [] datain;
    delete [] dataout;
  }

#ifdef FFTW
  void check_vs_fftw() 
  {
    unsigned int i,p;
    double d, nr2, nrinf, dftn2, dftninf;
    unsigned n;

    srand(17);
    
    Tp *data = new Tp [N2];

    fftw_complex* in;
    fftw_complex* out;
    fftw_plan plan;
    double fftwinf, fftw2;

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    for (i=0; i < N; ++i) {
       data[2*i] = rand()/(Tp)RAND_MAX - 0.5;  // distribute in [-0.5;0.5] as in FFTW
       data[2*i+1] = rand()/(Tp)RAND_MAX - 0.5;
       in[i][0] = data[2*i];
       in[i][1] = data[2*i+1];
    }

// apply FFT in-place
    gfft.fft(data);
//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<nrdata[2*i]<<","<<data[2*i+1]<<")"<<endl;

    plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//     plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    d = norm_inf(data, N2);
    for (i=0; i<N; ++i) {
      data[2*i]   -= out[i][0];
      data[2*i+1] -= out[i][1];
    }
    fftw2 = norm2(data, N2);
    fftwinf = norm_inf(data, N2);
    cout << fftw2 << "\t" << fftwinf << "\t" << fftwinf/d << endl;

    delete [] data;
    fftw_free(in);
    fftw_free(out);
  }
#endif
};

template<>
class GFFTcheck<Loki::NullType> {
public:
  void check_vs_dft() { }
  void check_vs_fftw() { }
};

typedef DOUBLE VType;
typedef IN_PLACE Place;
//typedef OUT_OF_PLACE Place;

const unsigned Min = 10;
const unsigned Max = 14;

typedef TYPELIST_3(OpenMP<2>, OpenMP<3>, OpenMP<4>) ParallList;
//typedef GenNumList<2, 10, SIntID>::Result NList;
typedef GenPowerList<2, 20, 2>::Result NList;
typedef GenerateTransform<NList, VType, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::Default, Place> Trans;


void check_complex_gfft(const int sign) 
{
    unsigned int i,p;
    double d, nr2, nrinf, dftn2, dftninf;
    unsigned n;

    srand(17);
    
    VType::ValueType *data, *nrdata;
    VType::ValueType *datain, *dataout;

#ifdef FFTW
    fftw_complex* in;
    fftw_complex* out;
    fftw_plan plan;
    double fftwinf, fftw2;
#endif

    for (p=Min; p<Max; ++p) {
//    for (p=1; p<2; ++p) {
    n = 1<<p;

    //AbstractFFT<VType>* fftobj = gfft->Instance().CreateObject(p);
    Trans::ObjectType* fftobj;
    if (sign == 1)
      fftobj = gfft.CreateTransformObject(p, VType::ID, DFT::ID, 1, OpenMP<4>::ID);
    else
      fftobj = gfft.CreateTransformObject(p, VType::ID, IDFT::ID, 1, OpenMP<4>::ID);

// sample data
    data = new VType::ValueType [2*n];
    nrdata = new VType::ValueType [2*n];
    datain  = new VType::ValueType [2*n];
    dataout = new VType::ValueType [2*n];
    
#ifdef FFTW
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
#endif
    for (i=0; i < n; ++i) {
       data[2*i] = rand()/(VType::ValueType)RAND_MAX - 0.5;  // distribute in [-0.5;0.5] as in FFTW
       nrdata[2*i] = data[2*i];
       datain[2*i] = data[2*i];
       data[2*i+1] = rand()/(VType::ValueType)RAND_MAX - 0.5;
       nrdata[2*i+1] = data[2*i+1];
       datain[2*i+1] = data[2*i+1];
#ifdef FFTW
       in[i][0] = data[2*i];
       in[i][1] = data[2*i+1];
#endif
    }

// apply FFT in-place
    fftobj->fft(data);
    four1(nrdata,n,sign);
    dft1(dataout, datain, n, (sign == -1));
    if (sign==-1) {
       for (i=0; i < 2*n; ++i) nrdata[i]/=(VType::ValueType)n;
    }
//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<nrdata[2*i]<<","<<data[2*i+1]<<")"<<endl;

#ifdef FFTW
    if (sign==1)
       plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    else
       plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
#endif

    d = norm_inf(data, 2*n);
    for (i=0; i<2*n; ++i) data[i]-=nrdata[i];
    nr2 = norm2(data, 2*n);
    nrinf = norm_inf(data, 2*n);
    cout<<"L2:"<<nr2<<"  Max:"<<nrinf<<"  Rel:"<<nrinf/d<<"  ";

#ifdef FFTW
    for (i=0; i<n; ++i) {
      data[2*i]   += nrdata[2*i] - out[i][0];
      data[2*i+1] += nrdata[2*i+1] - out[i][1];
    }
    fftw2=norm2(data,2*n);
    fftwinf=norm_inf(data,2*n);
    cout<<"L2:"<<fftw2<<"  Max:"<<fftwinf<<"  Rel:"<<fftwinf/d<<"  ";

    for (i=0; i<n; ++i) {
      data[2*i]   += out[i][0] - dataout[2*i];
      data[2*i+1] += out[i][1] - dataout[2*i+1];
    }
#else
    for (i=0; i<2*n; ++i) 
      data[i]   += nrdata[i] - dataout[i];
#endif
    dftn2 = norm2(data,2*n);
    dftninf = norm_inf(data,2*n);
    cout<<"L2:"<<dftn2<<"  Max:"<<dftninf<<"  Rel:"<<dftninf/d;

    cout<<endl;

    delete [] nrdata;
    delete [] data;
#ifdef FFTW
    fftw_free(in);
#endif

    }
}



int main(int argc, char *argv[])
{
  cout << "Forward transforms in place:" << endl;
  GFFTcheck<Trans::Result> check;
  check.check_vs_dft();
    
//     cout << "Backward transforms:" << endl;
//     check_complex_gfft(-1);

  return 0;
}


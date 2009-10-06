/***************************************************************************
 *   Copyright (C) 2007 by Volodymyr Myrnyy                                *
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

#include <cmath>
#include <iostream>

#include <fftw3.h>

#include "gfft.h"
#include "nrfft.h"

using namespace std;

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


typedef float VType;

int main(int argc, char *argv[])
{

    const unsigned Min = 1;
    const unsigned Max = 21;

    int p;
    int i;
    double d,d1;
    unsigned n;

    VType *data, *data1;
//     Loki::Factory<DFT::AbstractFFT<Tp>,unsigned int> fft;
//     FactoryInit<DFT::GFFTList<DFT::GFFTf,Min,Max>::Result>::apply(fft);
//     DFT::AbstractFFT<Tp>* mf = fft.CreateObject(k);

//    gfft_init();
//FactoryInit<DFT::GFFTList<DFT::GFFTf,Min,Max>::Result>::apply(DFT::gfft);
    DFT::GFFT_Singleton<Min,Max,VType,DFT::COMPLEX,DFT::INFREQ,DFT::FORWARD>* gfft;

   fftw_complex* in;
//   fftw_plan plan;

    for (p=1; p<Max; ++p) {
    n=1<<p;

    DFT::AbstractFFT<VType>* fftobj = gfft->Instance().CreateObject(p);

// sample data
    data = new VType [2*n];
    data1 = new VType [2*n];
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
    for (i=0; i < n; ++i) {
       data[2*i] = rand();
       data1[2*i] = data[2*i];
       data[2*i+1] = rand();
       data1[2*i+1] = data[2*i+1];
       in[i][0] = data[2*i];
       in[i][1] = data[2*i+1];
    }


// print out sample data
//     for (i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT in-place
    fftobj->fft(data);
    four1(data1,n,1);

//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
//     cout<<"Result of transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<data1[2*i]/(VType)n<<","<<data1[2*i+1]/(VType)n<<")"<<endl;

//     fftw_plan plan = fftw_plan_dft_1d(n, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
//     fftw_execute(plan);

    d1=norm_inf(data,2*n);
    for (i=0; i<2*n; ++i) data[i]-=data1[i];
//     for (i=0; i<n; ++i) {
//       data[2*i]-=in[i][0];
//       data[2*i+1]-=in[i][1];
//     }

    d=norm_inf(data,2*n);
    cout<<"L2:"<<norm2(data,2*n)<<"  Max:"<<d<<"  Rel:"<<d/d1<<endl;

    delete [] data1;
    delete [] data;
    fftw_free(in);

    }
// print out transformed data
//     for (i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;


}

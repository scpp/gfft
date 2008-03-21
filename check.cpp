/***************************************************************************
 *   Copyright (C) 2008 by Volodymyr Myrnyy                                *
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

typedef double VType;

const unsigned Min = 1;
const unsigned Max = 21;

template<class GFFTSingleton>
void check_complex_gfft(GFFTSingleton* gfft, const int sign) {

    unsigned int i,p;
    double d,nr2,nrinf;
    unsigned n;

    VType *data, *nrdata;

#ifdef FFTW
    fftw_complex* in;
    fftw_plan plan;
    double fftwinf,fftw2;
#endif

    for (p=1; p<Max; ++p) {
    n=1<<p;

    DFT::AbstractFFT<VType>* fftobj = gfft->Instance().CreateObject(p);

// sample data
    data = new VType [2*n];
    nrdata = new VType [2*n];
#ifdef FFTW
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
#endif
    for (i=0; i < n; ++i) {
       data[2*i] = rand();
       nrdata[2*i] = data[2*i];
       data[2*i+1] = rand();
       nrdata[2*i+1] = data[2*i+1];
#ifdef FFTW
       in[i][0] = data[2*i];
       in[i][1] = data[2*i+1];
#endif
    }

// apply FFT in-place
    fftobj->fft(data);
    four1(nrdata,n,sign);
    if (sign==-1) {
       for (i=0; i < 2*n; ++i) data[i]/=n;
    }

#ifdef FFTW
    if (sign==1)
       plan = fftw_plan_dft_1d(n, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
    else
       plan = fftw_plan_dft_1d(n, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
#endif

    d=norm_inf(data,2*n);
    for (i=0; i<2*n; ++i) data[i]-=nrdata[i];
    nr2=norm2(data,2*n);
    nrinf=norm_inf(data,2*n);
    cout<<"L2:"<<nr2<<"  Max:"<<nrinf<<"  Rel:"<<nrinf/d<<"  ";

#ifdef FFTW
    for (i=0; i<n; ++i) {
      data[2*i]   += nrdata[2*i]-in[i][0];
      data[2*i+1] += nrdata[2*i+1]-in[i][1];
    }
    fftw2=norm2(data,2*n);
    fftwinf=norm_inf(data,2*n);
    cout<<"L2:"<<fftw2<<"  Max:"<<fftwinf<<"  Rel:"<<fftwinf/d;
#endif
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

    static DFT::GFFT_Singleton<Min,Max,VType,DFT::COMPLEX,DFT::INFREQ,DFT::FORWARD>* gfft;
    static DFT::GFFT_Singleton<Min,Max,VType,DFT::COMPLEX,DFT::INFREQ,DFT::BACKWARD>* gfftback;

    check_complex_gfft(gfft,1);
    check_complex_gfft(gfftback,-1);

    return 0;
}


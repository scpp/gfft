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
/*
#include <iostream>

#include "gfft.h"
#include "nrfft.h"

using namespace std;

typedef double Tp;

int main(int argc, char *argv[])
{

    const unsigned Min = 1;
    const unsigned Max = 30;

    int k=2;
    int i,it;
    double t,mt;
    unsigned n= 1<<k;
    it = int(10000000./n)+1;

//     Loki::Factory<DFT::AbstractFFT<Tp>,unsigned int> fft;
//     FactoryInit<DFT::GFFTList<DFT::GFFTf,Min,Max>::Result>::apply(fft);
//     DFT::AbstractFFT<Tp>* mf = fft.CreateObject(k);

//    gfft_init();
//FactoryInit<DFT::GFFTList<DFT::GFFTf,Min,Max>::Result>::apply(DFT::gfft);

    DFT::GFFT_Singleton<double>* gfft;
    DFT::AbstractFFT<Tp>* mf = gfft->Instance().CreateObject(k);


// sample data
    Tp* data = new Tp [2*n*it];
    for (unsigned int j=0; j<it; ++j) {
      for (i=0; i < n; ++i) {
        data[2*(n*j+i)] = 2*i+j;
        data[2*(n*j+i)+1] = 2*i+j+1;
      }
    }

// print out sample data
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT in-place
    mf->fft(data);
//    four1(data,n,1);

// print out transformed data
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;


}
*/

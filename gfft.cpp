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

#include <iostream>

#include "gfft.h"

#include "nrfft.h"


using namespace std;

typedef double ValueType;

const unsigned Min = 1;
const unsigned Max = 30;

int main(int argc, char *argv[])
{

    unsigned int i,k=3;
    unsigned int n= 1<<k;

// Singleton holds the object factory for GFFT
    DFT::GFFT_Singleton<DFT::RVGFFTf,Min,Max,ValueType>* gfft;
    DFT::AbstractFFT<ValueType>* fftobj = gfft->Instance().CreateObject(k);

// or create the object factory directly
//    Loki::Factory<DFT::AbstractFFT<ValueType>,unsigned int> gfft;
//    FactoryInit<DFT::GFFTList<DFT::GFFTf,Min,Max>::Result>::apply(gfft);
//    DFT::AbstractFFT<ValueType>* fftobj = gfft.CreateObject(k);

// create FFT object of specific length, if the length is known at compile-time
//    DFT::GFFTf<3,ValueType>* fftobj = new DFT::GFFTf<3,ValueType>;

// create sample data
    ValueType* data = new ValueType [2*n];
    for (i=0; i < n; ++i) {
       data[2*i] = 2*i;
       data[2*i+1] = 2*i+1;
    }

// print out sample data
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT in-place
//    fftobj->fft(data);
    realft(data,n,1);

// print out transformed data
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

}


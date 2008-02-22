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

/** \file
    \brief Sample program to represent GFFT usage
*/

#include <iostream>

#include "gfft.h"

#include "nrfft.h"


using namespace std;

typedef double ValueType;

const unsigned Min = 1;
const unsigned Max = 24;

int main(int argc, char *argv[])
{

    unsigned int i,p=2;
    unsigned int n= 1<<p;

// There are three ways to create object to perform FFT of the length 2^p
// 1) Singleton holds the object factory for GFFT
    DFT::GFFT_Singleton<Min,Max,ValueType,DFT::RVGFFT,DFT::InFreq,DFT::Direct>* gfft;
    DFT::AbstractFFT<ValueType>* fftobj = gfft->Instance().CreateObject(p);

// 2) Create the object factory without singleton
//    Loki::Factory<DFT::AbstractFFT<ValueType>,unsigned int> gfft;
//    FactoryInit<DFT::GFFTList<Min,Max,ValueType,DFT::GFFT,DFT::InFreq,DFT::Direct>::Result>::apply(gfft);
//    DFT::AbstractFFT<ValueType>* fftobj = gfft.CreateObject(p);

// 3) create FFT object of specific length, if the length is known at compile-time
//     typedef DFT::GFFT<2,ValueType,DFT::InFreq,DFT::Direct> MyGFFT;
//     MyGFFT* fftobj = new MyGFFT;

// create sample data
    ValueType* data = new ValueType [2*n];
    for (i=0; i < n; ++i) {
       data[2*i] = 2*i;
       data[2*i+1] = 2*i+1;
    }

// print out sample data
    cout<<"Input data:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT in-place
    fftobj->fft(data);
//      four1(data,n,1);
//    realft(data,n,1);

// print out transformed data
    cout<<"Result of transform:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

}


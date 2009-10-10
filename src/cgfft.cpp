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
    \brief Sample program to represent GFFT usage with std::complex
*/

#include <iostream>

#include "gfftconf.h"

#include <omp.h>

using namespace std;

using namespace GFFT;

typedef COMPLEX_DOUBLE ValueType;
typedef AbstractFFT<ValueType::ValueType> AbstFFT;

const unsigned Min = 1;
const unsigned Max = 27;

int main(int argc, char *argv[])
{

    unsigned int i,p=2;
    unsigned int n= 1<<p;

    typedef Generate<1,3,ValueType> List;
    Loki::Factory<AbstFFT,unsigned int> gfft;
    FactoryInit<List::Result>::apply(gfft);
    unsigned int id1[5] = {p,ValueType::ID,RDFT::ID, Serial::ID,INFREQ::ID};
    unsigned int id2[5] = {p,ValueType::ID,IRDFT::ID,Serial::ID,INFREQ::ID};
    unsigned int p1 = List::trans_id(id1);
    unsigned int p2 = List::trans_id(id2);
    cout<<p1<<" "<<p2<<endl;
    AbstFFT* fftobj  = gfft.CreateObject(p1);
    AbstFFT* ifftobj = gfft.CreateObject(p2);

// There are three ways to create object to perform FFT of the length 2^p
// 1) Singleton holds the object factory for GFFT
//     DFT::GFFT_Singleton<Min,Max,ValueType,DFT::COMPLEX,DFT::INTIME,DFT::FORWARD>* gfft;
//     DFT::AbstractFFT<ValueType>* fftobj = gfft->Instance().CreateObject(p);
//
//     DFT::GFFT_Singleton<Min,Max,ValueType,DFT::COMPLEX,DFT::INTIME,DFT::BACKWARD>* igfft;
//     DFT::AbstractFFT<ValueType>* ifftobj = igfft->Instance().CreateObject(p);

// 2) Create the object factory without singleton
//    Loki::Factory<DFT::AbstractFFT<ValueType>,unsigned int> gfft;
//    FactoryInit<DFT::GFFTList<Min,Max,ValueType,DFT::COMPLEX,DFT::INTIME,DFT::FORWARD>::Result>::apply(gfft);
//    DFT::AbstractFFT<ValueType>* fftobj = gfft.CreateObject(p);
// 
//    Loki::Factory<DFT::AbstractFFT<ValueType>,unsigned int> igfft;
//    FactoryInit<DFT::GFFTList<Min,Max,ValueType,DFT::COMPLEX,DFT::INTIME,DFT::BACKWARD>::Result>::apply(igfft);
//    DFT::AbstractFFT<ValueType>* ifftobj = igfft.CreateObject(p);

// 3) create FFT object of specific length, if the length is known at compile-time
//     typedef DFT::GFFT<2,ValueType,DFT::REAL,DFT::INTIME,DFT::FORWARD> MyGFFT;
//     MyGFFT* fftobj = new MyGFFT;

// create sample data
    ValueType::ValueType* data = new ValueType::ValueType [2*n];
    for (i=0; i < n; ++i) {
       data[i].real() = 2*i;
       data[i].imag() = 2*i+1;
    }

// print out sample data
    cout<<"Input data:"<<endl;
    for (i=0; i < n; ++i)
      cout<<data[i]<<endl;

// apply FFT in-place
    fftobj->fft(data);

// print out transformed data
    cout<<"Result of transform:"<<endl;
    for (i=0; i < n; ++i)
      cout<<data[i]<<endl;

    ifftobj->fft(data);

// print out transformed data
    cout<<"Result of backward transform:"<<endl;
    for (i=0; i < n; ++i)
      cout<<data[i]<<endl;

}


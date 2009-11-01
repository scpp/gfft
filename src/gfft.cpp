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
    \brief Sample program to represent %GFFT usage
*/

#include <iostream>

#include "gfft.h"

#include <omp.h>

using namespace std;

using namespace GFFT;

typedef DOUBLE ValueType;
typedef AbstractFFT<ValueType::ValueType> AbstFFT;

const unsigned Min = 1;
const unsigned Max = 7;

int main(int argc, char *argv[])
{

    unsigned int i,p=2;
    unsigned int n= 1<<p;

// There are three ways to create object to perform FFT of the length 2^p
   //typedef GenList<1,7,DOUBLE,COMPLEX,FORWARD,OpenMP<2>,INFREQ> List;
//     typedef GenList<1,7,DOUBLE> List;
//     Loki::Factory<AbstFFT,unsigned int> gfft;
//     FactoryInit<List::Result>::apply(gfft);
//     unsigned int id1[6] = {p-1,DOUBLE::ID,COMPLEX::ID,FORWARD::ID,OpenMP<2>::ID,INFREQ::ID};
//     unsigned int id2[6] = {p-1,DOUBLE::ID,COMPLEX::ID,BACKWARD::ID,OpenMP<2>::ID,INFREQ::ID};
//     unsigned int p1 = List::trans_id(id1);
//     unsigned int p2 = List::trans_id(id2);
//     cout<<p1<<" "<<p2<<endl;
//     AbstFFT* fftobj  = gfft.CreateObject(p1);
//     AbstFFT* ifftobj = gfft.CreateObject(p2);

    typedef GenerateTransform<1,3,ValueType> Trans;

/*    Loki::Factory<AbstFFT,unsigned int> gfft;
    FactoryInit<Trans::Result>::apply(gfft);
    unsigned int id1[6] = {p,ValueType::ID,DFT::ID, SIntID<1>::ID,Serial::ID,INFREQ::ID};
    unsigned int id2[6] = {p,ValueType::ID,IDFT::ID,SIntID<1>::ID,Serial::ID,INFREQ::ID};
    unsigned int p1 = Trans::trans_id(id1);
    unsigned int p2 = Trans::trans_id(id2);
    cout<<p1<<" "<<p2<<endl;
    AbstFFT* fftobj  = gfft.CreateObject(p1);
    AbstFFT* ifftobj = gfft.CreateObject(p2);*/

    Trans gfft;
    Trans::Abstract* fftobj  = gfft.CreateTransformObject(p, ValueType::ID, DFT::ID);
    Trans::Abstract* ifftobj = gfft.CreateTransformObject(p, ValueType::ID, IDFT::ID);

 // 1) Singleton holds the object factory for GFFT
//    GFFT_Singleton<Min,Max,ValueType,COMPLEX,FORWARD,OpenMP<2>,INTIME>* gfft;
//    AbstFFT* fftobj = gfft->Instance().CreateObject(p);
//
//    GFFT_Singleton<Min,Max,ValueType,COMPLEX,BACKWARD,OpenMP<2>,INTIME>* igfft;
//    AbstFFT* ifftobj = igfft->Instance().CreateObject(p);

// 2) Create the object factory without singleton
//    Loki::Factory<AbstFFT,unsigned int> gfft;
//    FactoryInit<GFFTList<Min,Max,ValueType,COMPLEX,FORWARD,OpenMP<2>,INFREQ>::Result>::apply(gfft);
//    AbstFFT* fftobj = gfft.CreateObject(p);
//
//    Loki::Factory<AbstFFT,unsigned int> igfft;
//    FactoryInit<GFFTList<Min,Max,ValueType,COMPLEX,BACKWARD,OpenMP<2>,INFREQ>::Result>::apply(igfft);
//    AbstFFT* ifftobj = igfft.CreateObject(p);

// 3) create FFT object of specific length, if the length is known at compile-time
//     typedef GFFT<3,ValueType,COMPLEX,FORWARD,OpenMP<2>,INTIME> MyGFFT;
//     typedef GFFT<3,ValueType,COMPLEX,BACKWARD,OpenMP<2>,INTIME> MyIGFFT;
//     MyGFFT* fftobj = new MyGFFT;
//     MyIGFFT* ifftobj = new MyIGFFT;

// create sample data
    ValueType::ValueType* data = new ValueType::ValueType [2*n];
    for (i=0; i < n; ++i) {
       data[2*i] = 2*i;
       data[2*i+1] = 2*i+1; //2*i+1;
    }

// print out sample data
    cout<<"Input data:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT in-place
    fftobj->fft(data);

// print out transformed data
    cout<<"Result of transform:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

    ifftobj->fft(data);

// print out transformed data
    cout<<"Result of backward transform:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

}


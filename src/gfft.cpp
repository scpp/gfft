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

using namespace std;

using namespace GFFT;


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


typedef DOUBLE ValueType;
//typedef typename GenNumList<2, 3>::Result NList;
typedef TYPELIST_5(SIntID<2>, SIntID<3>, SIntID<4>, SIntID<5>, SIntID<7>) NList;
typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, ParallelizationGroup::Default, DecimationGroup::FullList > TransformSet;

int main(int argc, char *argv[])
{
//     unsigned int p = 2;
//     unsigned long i, n = (TransformType::ID == RDFT::ID) ? (1<<(p-1)) : (1<<p);
    unsigned int i, n = 7;
    
    typedef DFT TransformType;
/*
    TransformSet gfft;
    TransformSet::ObjectType* fftobj  = gfft.CreateTransformObject(n, ValueType::ID, TransformType::ID, 1, 
								   ParallelizationGroup::Default::ID, INFREQ::ID);
    TransformSet::ObjectType* ifftobj = gfft.CreateTransformObject(n, ValueType::ID, TransformType::Inverse::ID, 1, 
								   ParallelizationGroup::Default::ID, INFREQ::ID);
    
// create sample data
    ValueType::ValueType* data = new ValueType::ValueType [2*n];
    ValueType::ValueType* dataout = new ValueType::ValueType [2*n];
    ValueType::ValueType* data1 = new ValueType::ValueType [2*n];
    ValueType::ValueType* dataout1 = new ValueType::ValueType [2*n];
    for (i=0; i < n; ++i) {
       data[2*i] = 2*i;
       data[2*i+1] = 2*i+1; 
       data1[2*i] = 2*i;
       data1[2*i+1] = 2*i+1; 
       dataout[2*i] = 0;
       dataout[2*i+1] = 0; 
       dataout1[2*i] = 0;
       dataout1[2*i+1] = 0; 
    }

// print out sample data
    cout<<"Input data:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT in-place
//     fftobj->fft(data);
//     for (i=0; i < 2*n; ++i)
//       dataout[i] = data[i];

//      fftobj->fft(data);
    fftobj->fft(data, dataout);
    dft1(dataout1, data1, n, false);

// print out transformed data
    cout<<"Result of transform:"<<endl;
//      for (i=0; i < n; ++i)
//        cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
   for (i=0; i < n; ++i)
     cout<<"("<<dataout[2*i]<<","<<dataout[2*i+1]<<")"<<endl;

//    ifftobj->fft(data);
   ifftobj->fft(dataout, data);
   dft1(data1, dataout1, n, true);

// print out transformed data
    cout<<"Result of backward transform:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
    
    cout<<"Check against DFT:"<<endl;
    cout.precision(2);
    for (i=0; i < n; ++i)
      cout<<"("<<fabs(data[2*i]-data1[2*i])<<","<<fabs(data[2*i+1]-data1[2*i+1])<<")   \t("
      <<fabs(dataout[2*i]-dataout1[2*i])<<","<<fabs(dataout[2*i+1]-dataout1[2*i+1])<<")"<<endl;
*/

//    typedef Print<Factorization<SInt<111234> >::Result>::Result TTT;  // 2*3*18539
    typedef Print<Factorization<SInt<1024*169*1999> >::Result>::Result TTT;
//    typedef Print<Factorization<SInt<1024> >::Result>::Result TTT;
//    typedef Print<FactorizationLoop<13,2>::Result>::Result TTT;

}


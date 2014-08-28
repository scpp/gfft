/***************************************************************************
 *   Copyright (C) 2009-2014 by Vladimir Mirnyy                            *
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

#include "static_check.h"

#include "gfft.h"

using namespace std;

using namespace GFFT;


void dct1(double* output_data, const double* input_data, const unsigned int size, bool inverse)
{
  double pi = (inverse) ? M_PI : -M_PI;
  double a, ca;
  double invs = 1.0 / (size-1);
  for(unsigned int y = 0; y < size; y++) {
    output_data[y] = 0;
    for(unsigned int x = 1; x < size-1; x++) {
      a = pi * y * x * invs;
      ca = cos(a);
      output_data[y] += input_data[x] * ca;
    }
    output_data[y] += 0.5 * (input_data[0] + (y%2 > 0 ? -input_data[size-1] : input_data[size-1]));
    if(inverse) 
      output_data[y] *= invs;
  }
}

void dct2(double* output_data, const double* input_data, const unsigned int size, bool inverse)
{
  double pi = (inverse) ? M_PI : -M_PI;
  double a, ca;
  double invs = 1.0 / size;
  for(unsigned int y = 0; y < size; y++) {
    output_data[y] = 0;
    for(unsigned int x = 0; x < size; x++) {
      a = pi * y * (x+0.5) * invs;
      ca = cos(a);
      output_data[y] += input_data[x] * ca;
    }
    output_data[y] *= sqrt(2*invs);
    if(inverse) 
      output_data[y] *= invs;
  }
  output_data[0] /= sqrt(2);
}



typedef DOUBLE ValueType;
//typedef IN_PLACE Place;
typedef OUT_OF_PLACE Place;

static const int_t N = 4;
//typedef typename GenNumList<2, 3>::Result NList;
//typedef TYPELIST_4(SIntID<2>, SIntID<3>, SIntID<4>, SIntID<5>) NList;
typedef TYPELIST_1(SIntID<N>) NList;
typedef GenerateTransform<NList, ValueType, DFT, SIntID<1>, ParallelizationGroup::Default, Place> TransformSet;
//typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, OpenMP<2>, Place> TransformSet;

void print(ValueType::ValueType* d, const int_t n)
{
  cout << d[0];
  for (int_t i=1; i < n; ++i)
    cout << ", " << d[i];
  cout << endl;
}


int main(int argc, char *argv[])
{
    cout.precision(16);
//     unsigned int p = 2;
//     unsigned long i, n = (TransformType::ID == RDFT::ID) ? (1<<(p-1)) : (1<<p);
    int_t i, n = N;
    //cin >> n;
   
    typedef DFT TransformType;

    TransformSet gfft;
    TransformSet::ObjectType* fftobj  = gfft.CreateTransformObject(n, ValueType::ID, TransformType::ID, 1, 
								   ParallelizationGroup::Default::ID, Place::ID);
//     TransformSet::ObjectType* ifftobj = gfft.CreateTransformObject(n, ValueType::ID, TransformType::Inverse::ID, 1, 
// 								   ParallelizationGroup::Default::ID, Place::ID);
    
// create sample data
    ValueType::ValueType* data = new ValueType::ValueType [2*n];
    ValueType::ValueType* dataout = new ValueType::ValueType [2*n];
    ValueType::ValueType* data1 = new ValueType::ValueType [4*n];
    ValueType::ValueType* dataout1 = new ValueType::ValueType [4*n];
    for (i=0; i < n; ++i) {
//        data[2*i]   = rand()/(double)RAND_MAX - 0.5;   // distribute in [-0.5;0.5] as in FFTW
//        data[2*i+1] = rand()/(double)RAND_MAX - 0.5;
       data[i] = i;
       data1[i]   = data[i];
       dataout[i] = 0;
       dataout1[i] = 0;
    }

// print out sample data
    cout<<"Input : ";
    print(data, n);

// apply FFT in-place
//     fftobj->fft(data);
//     for (i=0; i < 2*n; ++i)
//       dataout[i] = data[i];

//       GFFTswap2<5,2,double> swp;
//       swp.apply(data);
//    fftobj->fft(data);
  fftobj->fft(data, dataout);
    dct2(dataout1, data1, n, false);
    
    //realft(data1, n, 1);

// print out transformed data
    cout.precision(3);
    cout<<"Result: ";
    print(dataout1, n);

   //ifftobj->fft(data);
//   ifftobj->fft(dataout, data);
   //dft1(data1, dataout1, n, true);

// print out transformed data
//     cout<<"Result of backward transform:"<<endl;
//     for (i=0; i < n; ++i)
//       cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
    
//     cout<<"Check against DFT:"<<endl;
//     double mx1(-1);
//     for (i=0; i < n; ++i) {
//       cout<<"("<<fabs(data[2*i]-data1[2*i])<<","<<fabs(data[2*i+1]-data1[2*i+1])<<")"<<endl;
//       mx1 = max(mx1, fabs(data[2*i]-data1[2*i]));
//       mx1 = max(mx1, fabs(data[2*i+1]-data1[2*i+1]));
//     }
//     cout<<"---------------------------------------------"<<endl;
//     cout << mx1 << endl;


//  cout << DOUBLE::Accuracy << endl;
}


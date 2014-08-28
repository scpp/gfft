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

#include "gfft.h"
#include "direct.h"

using namespace std;

using namespace GFFT;


typedef DOUBLE ValueType;
typedef OUT_OF_PLACE Place;

static const int_t N = 30;
//typedef typename GenNumList<2, 3>::Result NList;
//typedef TYPELIST_4(SIntID<2>, SIntID<3>, SIntID<4>, SIntID<5>) NList;
typedef TYPELIST_1(SIntID<N>) NList;
typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, ParallelizationGroup::Default, Place> TransformSet;
//typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, OpenMP<2>, Place> TransformSet;

int main(int argc, char *argv[])
{
    cout.precision(16);
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
    for (i=0; i < n; ++i) {
//        data[2*i]   = rand()/(double)RAND_MAX - 0.5;   // distribute in [-0.5;0.5] as in FFTW
//        data[2*i+1] = rand()/(double)RAND_MAX - 0.5;
       data[2*i] = 2*i;
       data[2*i+1] = 2*i+1; 
       dataout[2*i] = 0;
       dataout[2*i+1] = 0; 
    }

    DFT_wrapper<ValueType::ValueType> dft(data, n);

 // print out sample data
    cout<<"Input data:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT out-of-place
   fftobj->fft(data, dataout);

// do simple dft
   dft.apply();
    
   ValueType::ValueType* dataout1 = dft.getdata();

// print out transformed data
   cout.precision(3);
   cout<<"Result of transform:"<<endl;
   for (i=0; i < n; ++i)
     cout<<"("<<dataout[2*i]<<","<<dataout[2*i+1]<<")   \t("<<dataout1[2*i]<<","<<dataout1[2*i+1]<<") \t"<<endl;

   dft.diff(dataout);

   cout<<"Check against DFT:"<<endl;
   double mx(-1);
   for (i=0; i < n; ++i) {
      cout<<"("<<fabs(dataout[2*i])<<","<<fabs(dataout[2*i+1])<<")"<<endl;
      mx = max(mx, fabs(dataout[2*i]));
      mx = max(mx, fabs(dataout[2*i+1]));
   }
   cout<<"---------------------------------------------"<<endl;
   cout << mx << endl;

   delete [] data;
   delete [] dataout;
}


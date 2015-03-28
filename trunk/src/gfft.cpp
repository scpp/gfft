/***************************************************************************
 *   Copyright (C) 2009-2015 by Vladimir Mirnyy                            *
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
    \brief Sample program to show %GFFT usage
    
    The variables NUM, TYPE, PLACE and NUMTHREADS are defined in CMakeLists.txt file 
    in the same folder. They can also be supplied to cmake command
    Example: cmake -DNUM=8 -DTYPE=FLOAT -DPLACE=IN_PLACE -DNUMTHREADS=4 .
*/

#include <iostream>

#include "gfft.h"
#include "direct.h"

using namespace std;

using namespace GFFT;

typedef TYPE ValueType;
typedef PLACE Place;

static const int_t N = NUM;
static const int_t NThreads = NUMTHREADS;
//typedef TYPELIST_4(SIntID<2>, SIntID<3>, SIntID<4>, SIntID<5>) NList;
typedef TYPELIST_1(SIntID<N>) NList;
//typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, ParallelizationGroup::Default, Place> TransformSet;
typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, OpenMP<NThreads>, Place> TransformSet;


int main(int argc, char *argv[])
{
    cout.precision(16);
    int_t i, n = N;
    //cin >> n;
   
    typedef DFT TransformType;
    typedef ValueType::ValueType T;

    TransformSet gfft;
    TransformSet::ObjectType* fftobj  = gfft.CreateTransformObject(n, ValueType::ID, TransformType::ID, 1, 
								   OpenMP<NThreads>::ID, Place::ID);
//     TransformSet::ObjectType* ifftobj = gfft.CreateTransformObject(n, ValueType::ID, TransformType::Inverse::ID, 1, 
// 								   ParallelizationGroup::Default::ID, Place::ID);
    
// create sample data
    T* data = new T [2*n];
#if PLACE == OUT_OF_PLACE
    T* dataout = new T [2*n];
#endif
    for (i=0; i < n; ++i) {
//        data[2*i]   = rand()/(double)RAND_MAX - 0.5;   // distribute in [-0.5;0.5] as in FFTW
//        data[2*i+1] = rand()/(double)RAND_MAX - 0.5;
       data[2*i] = 2*i;
       data[2*i+1] = 2*i+1; 
#if PLACE == OUT_OF_PLACE
       dataout[2*i] = 0;
       dataout[2*i+1] = 0; 
#endif
    }

    DFT_wrapper<T> dft(data, n);

 // print out sample data
#ifdef FOUT
    cout<<"Input data:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
#endif

#if PLACE == OUT_OF_PLACE
// apply FFT out-of-place
   fftobj->fft(data, dataout);
#else
// apply FFT in-place
   fftobj->fft(data);
   T* dataout = data;
#endif
   
// do simple dft
   dft.apply();

   T* dataout1 = dft.getdata();

// print out transformed data
   cout.precision(3);
#ifdef FOUT
   cout<<"Result of transform:"<<endl;
   for (i=0; i < n; ++i)
     cout<<"("<<dataout[2*i]<<","<<dataout[2*i+1]<<")   \t("<<dataout1[2*i]<<","<<dataout1[2*i+1]<<") \t"<<endl;
#endif

   dft.diff(dataout);

   cout<<"Check against DFT:"<<endl;
   double mx(-1);
   double s = 0.;
   for (i=0; i < n; ++i) {
#ifdef FOUT
      cout<<"("<<fabs(dataout[2*i])<<","<<fabs(dataout[2*i+1])<<")"<<endl;
#endif
      mx = max(mx, fabs(dataout[2*i]));
      mx = max(mx, fabs(dataout[2*i+1]));
      s += dataout[2*i]*dataout[2*i];
      s += dataout[2*i+1]*dataout[2*i+1];
   }
   cout<<"---------------------------------------------"<<endl;
   cout << N << ": " << ValueType::name() << ", " << Place::name() << endl;
   cout << mx << "  " << sqrt(s) << endl;

   delete [] data;
#if PLACE == OUT_OF_PLACE
   delete [] dataout;
#endif
}


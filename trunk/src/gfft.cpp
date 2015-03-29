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
    
    This sample program is generalized to handle both simple arrays containing consequent pairs 
    of real and imaginary parts and std::complex. It would also handle any user-defined complex 
    types, if they define functions real() and imag() and the output operator to std::ostream
    
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
    int_t i;
   
    typedef DFT TransformType;
    typedef ValueType::ValueType T;
    typedef ValueType::base_type BT;
    static const int C = Loki::TypeTraits<T>::isStdFundamental ? 2 : 1;

    TransformSet gfft;
    TransformSet::ObjectType* fftobj  = gfft.CreateTransformObject(N, ValueType::ID, TransformType::ID, 1, 
								   OpenMP<NThreads>::ID, Place::ID);
//     TransformSet::ObjectType* ifftobj = gfft.CreateTransformObject(n, ValueType::ID, TransformType::Inverse::ID, 1, 
// 								   ParallelizationGroup::Default::ID, Place::ID);
    
// create sample data
    T* data = new T [C*N];
#if PLACE == OUT_OF_PLACE
    T* dataout = new T [C*N];
#endif
    for (i=0; i < N; ++i) {
//      GenInput<T>::rand(data, i);  // distribute in [-0.5;0.5] as in FFTW
       GenInput<T>::seq(data, i);
    }

    DFT_wrapper<BT> dft(data, N);

 // print out sample data
#ifdef FOUT
    cout<<"Input data:"<<endl;
    for (i=0; i < N; ++i)
      cout << GenOutput<T>(data,i) << endl;
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

   BT* dataout1 = dft.getdata();

// print out transformed data
   cout.precision(3);
#ifdef FOUT
   cout<<"Result of transform:"<<endl;
   for (i=0; i < N; ++i)
     cout << GenOutput<T>(dataout,i) << "   \t"<< GenOutput<BT>(dataout1,i) <<" \t"<<endl;
#endif

   dft.diff(dataout);

   cout<<"Check against DFT:"<<endl;
   double mx(-1);
   double s = 0.;
   for (i=0; i < N; ++i) {
#ifdef FOUT
      cout << GenOutput<T>(dataout,i) << endl;
#endif
      double re = ComplexWrapper<T>(dataout,i).real();
      double im = ComplexWrapper<T>(dataout,i).imag();
      mx = max(mx, fabs(re));
      mx = max(mx, fabs(im));
      s += re*re;
      s += im*im;
   }
   cout<<"---------------------------------------------"<<endl;
   cout << N << ": " << ValueType::name() << ", " << Place::name() << endl;
   cout << mx << "  " << sqrt(s) << endl;

   delete [] data;
#if PLACE == OUT_OF_PLACE
   delete [] dataout;
#endif
}


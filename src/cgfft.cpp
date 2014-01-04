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
    \brief Sample program to represent %GFFT usage with std::complex
*/

#include <iostream>

#include "gfft.h"

#include <omp.h>

using namespace std;

using namespace GFFT;

typedef COMPLEX_DOUBLE ValueType;

int main(int argc, char *argv[])
{

    unsigned int i,p=2;
    unsigned int n= 1<<p;

    typedef TYPELIST_1(OpenMP<2>) ParallList;
    typedef GenerateTransform<1,3,ValueType, TransformTypeGroup::FullList, SIntID<1>, OpenMP<2>, INTIME> TransformSet;
    TransformSet gfft;
    TransformSet::ObjectType* fftobj  = gfft.CreateTransformObject(p, ValueType::ID, RDFT::ID, 1, OpenMP<2>::ID, INTIME::ID);
    TransformSet::ObjectType* ifftobj = gfft.CreateTransformObject(p, ValueType::ID, IRDFT::ID, 1, OpenMP<2>::ID, INTIME::ID);

// create sample data
    ValueType::ValueType* data = new ValueType::ValueType [2*n];
    for (i=0; i < n; ++i) {
       data[i] = ValueType::ValueType(2*i, 2*i+1);
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


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

typedef DOUBLE ValueType;
typedef GenerateTransform<1, 4, ValueType, TransformTypeGroup::FullList, SIntID<1> > TransformSet;


int main(int argc, char *argv[])
{
    unsigned int p = 2;
    typedef DFT TransformType;

    TransformSet gfft;
    TransformSet::ObjectType* fftobj  = gfft.CreateTransformObject(p, ValueType::ID, TransformType::ID, 1);
    TransformSet::ObjectType* ifftobj = gfft.CreateTransformObject(p, ValueType::ID, TransformType::Inverse::ID, 1);

    unsigned long i, n = (TransformType::ID == RDFT::ID) ? (1<<(p-1)) : (1<<p);

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


/***************************************************************************
 *   Copyright (C) 2014 by Vladimir Mirnyy                                 *
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
    \brief Check accuracy of GFFT results comparing to FFTW and DFT of double-double precision
*/

#include "gfft_accuracy.h"

using namespace std;

using namespace GFFT;

typedef COMPLEX_DOUBLE VType;
//typedef IN_PLACE Place;
typedef OUT_OF_PLACE Place;

const unsigned Min = 2;
const unsigned Max = 11;

typedef TYPELIST_3(OpenMP<2>, OpenMP<4>, OpenMP<8>) ParallList;
//typedef GenNumList<2, 15, SIntID>::Result NList;
typedef GenPowerList<Min, Max, 2>::Result NList;
typedef GenerateTransform<NList, VType, TransformTypeGroup::Default, SIntID<1>, ParallList, Place> Trans;

ostream& operator<<(ostream& os, const dd_real& v)
{
  os << v.to_string(16);
  return os;
}

ostream& operator<<(ostream& os, const qd_real& v)
{
  os << v.to_string(16);
  return os;
}

void print_header() 
{
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<" N         Norm2           NormInf       Relative NormInf    "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
}


int main(int argc, char *argv[])
{
  unsigned int oldcw;
  fpu_fix_start(&oldcw);
  
  cout.precision(4);

//   cout << "double-double DFT vs. FFTW:" << endl;
//   FFTcompare<DFT_wrapper<dd_real>, FFTW_wrapper<fftw_complex> > comp;
//   print_header();
//   comp.apply(Min,Max);
//   cout<<"------------------------------------------------------------------------------"<<endl<<endl;
  
  cout << "double-double DFT vs. GFFT:" << endl;
  GFFTcheck<Trans::Result, DFT_wrapper<dd_real>, Place> check_dft;
  print_header();
  check_dft.apply();
  cout<<"------------------------------------------------------------------------------"<<endl<<endl;
  
//   cout << "GFFT vs. FFTW:" << endl;
//   GFFTcheck<Trans::Result, FFTW_wrapper<fftw_complex>, Place> check_fftw;
//   print_header();
//   check_fftw.apply();
//   cout<<"------------------------------------------------------------------------------"<<endl<<endl;

  fpu_fix_end(&oldcw);
  return 0;
}


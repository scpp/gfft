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
    \brief Program to demonstrate compile-time computation of PI with arbitrary accuracy
    
    Two last digits of the resulted PI may not be correct due to roundings!
*/

#include <iostream>

#include "gfft.h"

using namespace GFFT;

using namespace std;

static const int Accuracy = 2; // *9 - 64-bit; *4 - 32 bit

template<int N>
struct MetaSqrtTest 
{
  typedef typename SqrtDecAcc<SInt<N>,Accuracy>::Result TSqrtDec;
  
  static void apply() 
  {
    MetaSqrtTest<N-1>::apply();
    cout << ">>> Sqrt(" << N << ") = " <<endl;
    cout << "Compile-time: ";
    Cout<TSqrtDec>::apply(cout);
    cout << endl;
    cout << "Run time    : " << sqrt(static_cast<double>(N)) << endl;
  }
};

template<>
struct MetaSqrtTest<0>
{
  static void apply() {}
};

int main(int argc, char *argv[])
{
  cout.precision(16);
//   typedef SqrtInitGuessDec<SInt<N>,Accuracy> TSqrtGuess;
//   typedef SqrtDecimal<0,SInt<N>,Loki::NullType,Accuracy> Step0;
//   typedef SqrtDecimal<0,SInt<N>,Step0::ResultAux,Accuracy> Step1;
//   typedef SqrtDecimal<0,SInt<N>,Step1::ResultAux,Accuracy> Step2;
//   typedef SqrtDecimal<0,SInt<N>,Step2::ResultAux,Accuracy> Step3;
//   typedef SqrtDecimal<0,SInt<N>,Step3::ResultAux,Accuracy> Step4;

//   Cout<TSqrtGuess::Result>::apply(cout);
//   cout << endl;
//   Cout<Step1::Result>::apply(cout);
//   cout << endl;
//   Cout<Step2::Result>::apply(cout);
//   cout << endl;
//   Cout<Step3::Result>::apply(cout);
//   cout << endl;
//   Cout<Step4::Result>::apply(cout);
//   cout << endl;

  static const int_t N = 125348; // example from wikipedia
  typedef typename SqrtDecAcc<SInt<N>,Accuracy>::Result TSqrtDec;
  
  cout << ">>> Sqrt(" << N << ") = " <<endl;
  cout << "Compile-time: ";
  Cout<TSqrtDec>::apply(cout);
  cout << endl;
  cout << "Run time    : " << sqrt(static_cast<double>(N)) << endl;

  //MetaSqrtTest<10>::apply();
}


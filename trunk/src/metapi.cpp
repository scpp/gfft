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
    \brief Program to demonstrate compile-time computation of PI with arbitrary accuracy
    
    Two last digits of the resulted PI may not be correct due to roundings!
*/

#include <iostream>

#include "gfft.h"

using namespace GFFT;

using namespace std;

static const int Accuracy = 7; // *9 - 64-bit; *4 - 32 bit

int main(int argc, char *argv[])
{
  cout.precision(16);
  typedef PiDecAcc<Accuracy>::Result TPiDec;
  cout << "Compile-time PI: ";
  Cout<TPiDec>::apply(cout);
  cout << endl;
  cout << "           M_PI: " << M_PI << endl;
}


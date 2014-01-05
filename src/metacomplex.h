/***************************************************************************
 *   Copyright (C) 2007-2014 by Vladimir Mirnyy                            *
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

#ifndef __metacomplex_h
#define __metacomplex_h

/** \file
    \brief Compile-time computing of complex numbers
*/

#include <cmath>

#include "srational.h"



template<class Real, class Imag> 
struct MComplex {
  typedef Real Re;
  typedef Imag Im;
};


template<class Re1, class Im1, class Re2, class Im2>
class Add<MComplex<Re1,Im1>, MComplex<Re2,Im2> > {
  typedef typename Add<Re1,Re2>::Result Re;
  typedef typename Add<Im1,Im2>::Result Im;
public:
  typedef MComplex<Re,Im> Result;
};
  
template<class Re1, class Im1, class Re2, class Im2>
class Sub<MComplex<Re1,Im1>, MComplex<Re2,Im2> > {
  typedef typename Sub<Re1,Re2>::Result Re;
  typedef typename Sub<Im1,Im2>::Result Im;
public:
  typedef MComplex<Re,Im> Result;
};

template<class Re1, class Im1, class Re2, class Im2>
class Mult<MComplex<Re1,Im1>, MComplex<Re2,Im2> > {
  typedef typename Mult<Re1,Re2>::Result PRR;
  typedef typename Mult<Im1,Im2>::Result PII;
  typedef typename Mult<Im1,Re2>::Result PIR;
  typedef typename Mult<Re1,Im2>::Result PRI;
  typedef typename Sub<PRR,PII>::Result Re;
  typedef typename Add<PIR,PRI>::Result Im;
public:
  typedef MComplex<Re,Im> Result;
};


#endif /*__metacomplex_h*/

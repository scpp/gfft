/***************************************************************************
 *   Copyright (C) 2007 by Volodymyr Myrnyy                                *
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

#ifndef __gfftpolicy_h
#define __gfftpolicy_h

#include "loki/Typelist.h"

/** \file
    \brief Policy classes
*/


namespace DFT {


template<class TList>
struct PoliciesHandler;

/// Instantiates each type of the typelist and calls function apply
/**
    All the classes in TList must include function apply(T*) with one parameter
    as a pointer on type T.
*/

template<class Head, class Tail>
struct PoliciesHandler<Loki::Typelist<Head,Tail> > {
   template<typename T>
   void apply(T* data) {
      obj.apply(data);
      next.apply(data);
   }
private:
   Head obj;
   PoliciesHandler<Tail> next;
};

template<>
struct PoliciesHandler<Loki::NullType> {
   template<typename T>
   void apply(T*) { }
};


/// Abstract base class to build an object factory
/** It shares the function fft(T*) between
    classes that represent FFT of different length.
*/
template<typename T>
class AbstractFFT {
public:
   virtual void fft(T*) = 0;
};

/// Abstract empty base class
/** This class is passed instead of AbstractFFT,
    if object factory is not needed
    to avoid a virtual function call penalty
*/
class Empty { };


/// Policy for a definition of forward FFT
template<unsigned N, typename T>
struct Forward {
   enum { Sign = 1 };
   void apply(T*) { }
};

/// Policy for a definition of backward FFT
template<unsigned N, typename T>
struct Backward {
   enum { Sign = -1 };
   void apply(T* data) {
      for (T* i=data; i<data+2*N; ++i) *i/=N;
   }
};


}  //namespace DFT

#endif /*__gfftpolicy_h*/

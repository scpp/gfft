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

/** \file
    \brief Policy classes
*/


namespace DFT {

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


/// Policy for decimation-in-time FFT
template<unsigned P, class T, class Direction>
class InTime {
public:
   enum { Power = P, N = 1<<P };
   typedef T ValueType;
   typedef Direction Dir;

   void apply(T* data) {
      scramble.apply(data);
      recursion.apply(data);
      direction.apply(data);
   }
private:
   GFFTswap<N,T> scramble;
   DLTime<N,T,Direction::Sign> recursion;
   Direction direction;
};


/// Policy for decimation-in-frequency FFT
template<unsigned P, class T, class Direction>
class InFreq {
public:
   enum { Power = P, N = 1<<P };
   typedef T ValueType;
   typedef Direction Dir;

   void apply(T* data) {
      recursion.apply(data);
      scramble.apply(data);
      direction.apply(data);
   }
private:
   GFFTswap<N,T> scramble;
   DLFreq<N,T,Direction::Sign> recursion;
   Direction direction;
};


/// Policy for a definition of direct FFT
template<unsigned N, typename T>
struct Direct {
   enum { Sign = 1 };
   void apply(T*) { }
};

/// Policy for a definition of inverse FFT
template<unsigned N, typename T>
struct Inverse {
   enum { Sign = -1 };
   void apply(T* data) {
      for (T* i=data; i<data+2*N; ++i) *i/=N;
   }
};


}  //namespace DFT

#endif /*__gfftpolicy_h*/

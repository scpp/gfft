/***************************************************************************
 *   Copyright (C) 2007-2015 by Vladimir Mirnyy                            *
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

#ifndef __metapow_h
#define __metapow_h

/** \file
    \brief Compile-time computing of power function
*/


namespace MF {


template<int_t A, int_t P, class RetType = long double>
struct DPow {
  static RetType value() {
    return static_cast<RetType>(A) * DPow<A,P-1,RetType>::value();
  }
};

template<int_t A, class RetType>
struct DPow<A,1,RetType> {
  static RetType value() { return A; }
};

template<int_t A, class RetType>
struct DPow<A,0,RetType> {
  static RetType value() { return 1; }
};


template<unsigned N, unsigned P>
struct IPow {
  static const unsigned long value = IPow<N,P-1>::value * N;
};

template<unsigned N>
struct IPow<N,1> {
  static const unsigned long value = N;
};

template<unsigned N>
struct IPow<N,0> {
  static const unsigned long value = 1;
};


template<class N, int_t P>
struct IPowBig {
  typedef typename Mult<N, typename IPowBig<N,P-1>::Result>::Result Result;
};

template<class N>
struct IPowBig<N,1> {
  typedef N Result;
};

template<class N>
struct IPowBig<N,0> {
  typedef SInt<1> Result;
};



} // namespace MF

#endif /*__metapow_h*/

/***************************************************************************
 *   Copyright (C) 2008 by Volodymyr Myrnyy                                *
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

#ifndef __sfraction_h
#define __sfraction_h

#include "sbigint.h"

template<class Num, class Den>
struct SFraction {
   typedef Num Numer;
   typedef Den Denom;
};

template<class C>
class Simplify;

///Greatest common divisor of A and B
template<class A, class B>
struct GCD;

template<int N1, int N2>
struct GCD<SInt<N1>, SInt<N2> > {
   typedef typename GCD<SInt<N2>,SInt<N1%N2> >::Result Result;
};

template<int N>
struct GCD<SInt<N>, SInt<0> > {
   typedef SInt<N> Result;
};

template<class N, class D>
class Simplify<SFraction<N,D> > {
   typedef typename GCD<N,D>::Result T;
   typedef typename Div<N,T>::Result Num;
   typedef typename Div<D,T>::Result Den;
public:
   typedef SFraction<Num,Den> Result;
};

/// Multiplication of two compile-time Fractions
template<class N1, class D1, class N2, class D2>
class Mult<SFraction<N1,D1>, SFraction<N2,D2> > {
   typedef typename Mult<N1,N2>::Result Num;
   typedef typename Mult<D1,D2>::Result Den;
public:
   typedef SFraction<Num,Den> Result;
};


//////////////////////////////////////////////

template<class N1, class D1, class N2, class D2>
class Add<SFraction<N1,D1>, SFraction<N2,D2> > {
   typedef typename Mult<N1,D2>::Result T1;
   typedef typename Mult<N2,D1>::Result T2;
   typedef typename Add<T1,T2>::Result Num;
   typedef typename Mult<D1,D2>::Result Den;
public:
   typedef SFraction<Num,Den> Result;
};



#endif

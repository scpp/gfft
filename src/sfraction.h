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

template<class BigInt, int_t NDecPlaces, base_t DecBase>
struct SDecimalFraction {
   typedef BigInt Num;
   static const int_t NDec = NDecPlaces;
   static const base_t Base = DecBase;
};



///Greatest common divisor of A and B
/*
template<class A, class B>
struct GCD;

template<bool S1, class N1, bool S2, class N2, base_t B>
struct GCD<SBigInt<S1,N1,B>, SBigInt<S2,N2,B> > {
  typedef SBigInt<S1,N1,B> TN1;
  typedef SBigInt<S2,N2,B> TN2;
  static const int C = NL::Compare<TN1,TN2>::value;
  typedef typename Loki::Select<(C==0), TN1,
          typename Loki::Select<(C>0), typename GCD<TN2, typename Simplify<typename Sub<TN1,TN2>::Result>::Result>::Result,
                                       typename GCD<TN1, typename Simplify<typename Sub<TN2,TN1>::Result>::Result>::Result>::Result>::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct GCD<SBigInt<S,N1,Base>, SInt<N2> > {
  typedef SBigInt<S,N1,Base> TN1;
  typedef SInt<N2> TN2;
  static const int C = NL::Compare<TN1,TN2>::value;
  typedef typename Loki::Select<(C==0), TN1,
          typename Loki::Select<(C>0), typename GCD<TN2, typename Simplify<typename Sub<TN1,TN2>::Result>::Result>::Result,
                                       typename GCD<TN1, typename Simplify<typename Sub<TN2,TN1>::Result>::Result>::Result>::Result>::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct GCD<SInt<N2>, SBigInt<S,N1,Base> > 
: public GCD<SBigInt<S,N1,Base>, SInt<N2> > {};

template<int_t N1, int_t N2>
struct GCD<SInt<N1>, SInt<N2> > {
   typedef typename Loki::Select<(N1==N2), SInt<N1>, 
           typename Loki::Select<(N1>N2), typename GCD<SInt<N2>, SInt<N1-N2> >::Result,
                                          typename GCD<SInt<N1>, SInt<N2-N1> >::Result>::Result>::Result Result;
};

template<bool S, class N, base_t Base>
struct GCD<SBigInt<S,N,Base>, SInt<0> > {
   typedef SBigInt<S,N,Base> Result;
};

template<int_t N>
struct GCD<SInt<N>, SInt<0> > {
   typedef SInt<N> Result;
};
*/

template<class A, class B,
int C = (int)NL::Compare<A,B>::value>
struct GCD;

template<class A, class B>
struct GCD<A,B,0> {
  typedef A Result;
};

template<bool S1, class N1, bool S2, class N2, base_t B>
struct GCD<SBigInt<S1,N1,B>, SBigInt<S2,N2,B>, 1> {
  typedef SBigInt<S1,N1,B> TN1;
  typedef SBigInt<S2,N2,B> TN2;
  typedef typename Simplify<typename Div<TN1,TN2>::ModResult>::Result Mod;
//  typedef Loki::TL::Print<Mod> T;
  typedef typename GCD<TN2,Mod>::Result Result;
};

template<bool S1, class N1, base_t B1, bool S2, class N2, base_t B2>
struct GCD<SBigInt<S1,N1,B1>, SBigInt<S2,N2,B2>, -1>
: public GCD<SBigInt<S2,N2,B2>, SBigInt<S1,N1,B1> > {};

template<bool S, class N1, base_t Base, int_t N2>
struct GCD<SBigInt<S,N1,Base>, SInt<N2>, 1> {
  typedef typename Simplify<typename Div<SBigInt<S,N1,Base>,SInt<N2> >::ModResult>::Result Mod;
  typedef typename GCD<SInt<N2>,Mod>::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct GCD<SBigInt<S,N1,Base>, SInt<N2>, -1> 
: public GCD<SInt<N2>, SBigInt<S,N1,Base> > {};

template<bool S, class N1, base_t Base, int_t N2>
struct GCD<SInt<N2>, SBigInt<S,N1,Base>, 1> 
{
   static const int_t N = Evaluate2Int<SBigInt<S,N1,Base>, int_t>::Value;
   typedef typename GCD<SInt<N2>, SInt<N> >::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct GCD<SInt<N2>, SBigInt<S,N1,Base>, -1> 
: public GCD<SBigInt<S,N1,Base>, SInt<N2> > {};

template<int_t N1, int_t N2>
struct GCD<SInt<N1>, SInt<N2>, 1> {
   typedef typename GCD<SInt<N2>,SInt<N1%N2> >::Result Result;
};

template<int_t N1, int_t N2>
struct GCD<SInt<N1>, SInt<N2>, -1> 
: public GCD<SInt<N2>, SInt<N1> > {};

template<bool S1, bool S2, class N, base_t Base>
struct GCD<SBigInt<S1,N,Base>, SBigInt<S2,Loki::NullType,Base>, 1> {
   typedef SBigInt<S1,N,Base> Result;
};

template<bool S, class N, base_t Base>
struct GCD<SBigInt<S,N,Base>, SInt<0>, 1> {
   typedef SBigInt<S,N,Base> Result;
};

template<int_t N>
struct GCD<SInt<N>, SInt<0>, 1> {
   typedef SInt<N> Result;
};


template<class N, class D>
class Simplify<SFraction<N,D> > {
   typedef typename GCD<N,D>::Result T;
   typedef typename Div<N,T>::DivResult Num;
   typedef typename Div<D,T>::DivResult Den;
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

template<class N, class D, int Num>
class Mult<SFraction<N,D>,SInt<Num> > {
   typedef typename Mult<N,SInt<Num> >::Result Numer;
public:
   typedef SFraction<Numer,D> Result;
};

template<class N, class D, bool S, class NList, unsigned int Base>
class Mult<SFraction<N,D>,SBigInt<S,NList,Base> > {
   typedef typename Mult<N,SBigInt<S,NList,Base> >::Result Numer;
public:
   typedef SFraction<Numer,D> Result;
};

/////////////////////////////////////////////////////////////

template<class Numer, class Denom>
struct Negate<SFraction<Numer,Denom> > {
   typedef typename Negate<Numer>::Result NewNumer;
   typedef SFraction<NewNumer,Denom> Result;
};

/////////////////////////////////////////////////////////////

template<class N1, class D1, class N2, class D2>
class Add<SFraction<N1,D1>, SFraction<N2,D2> > {
   typedef typename Mult<N1,D2>::Result T1;
   typedef typename Mult<N2,D1>::Result T2;
   typedef typename Add<T1,T2>::Result Num;
   typedef typename Mult<D1,D2>::Result Den;
public:
   typedef SFraction<Num,Den> Result;
};

template<class N, class D, int Num>
class Add<SFraction<N,D>,SInt<Num> > {
   typedef SInt<Num> T;
   typedef typename Add<
           typename Mult<D,T>::Result,N>::Result Numer;
public:
   typedef SFraction<Numer,D> Result;
};

template<class N, class D, bool S, class NList, unsigned int Base>
class Add<SFraction<N,D>,SBigInt<S,NList,Base> > {
   typedef SBigInt<S,NList,Base> T;
   typedef typename Add<
           typename Mult<D,T>::Result,N>::Result Numer;
public:
   typedef SFraction<Numer,D> Result;
};



#endif

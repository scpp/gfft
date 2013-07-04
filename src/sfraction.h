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

typedef SFraction<SInt<1>,SInt<1> > UnitFraction;

template<class BigInt, int_t NDecPlaces, base_t DecBase>
struct SDecimalFraction {
   typedef BigInt Num;
   static const int_t NDec = NDecPlaces;
   static const base_t Base = DecBase;
};



///Greatest common divisor of A and B

template<class A, class B, char C>
struct __GCD;

template<class A, class B>
struct __GCD<A,B,0> {
  typedef A Result;
};

template<bool S1, class N1, bool S2, class N2, base_t B>
struct __GCD<SBigInt<S1,N1,B>, SBigInt<S2,N2,B>, 1> {
  typedef SBigInt<S1,N1,B> TN1;
  typedef SBigInt<S2,N2,B> TN2;
  typedef typename Div<TN1,TN2>::ModResult Mod;
  static const char C = NL::Compare<TN2,Mod>::value;
  typedef typename __GCD<TN2,Mod,C>::Result Result;
};

template<bool S1, class N1, base_t B1, bool S2, class N2, base_t B2>
struct __GCD<SBigInt<S1,N1,B1>, SBigInt<S2,N2,B2>, -1>
: public __GCD<SBigInt<S2,N2,B2>, SBigInt<S1,N1,B1>, 1> {};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SBigInt<S,N1,Base>, SInt<N2>, 1> {
  typedef typename Div<SBigInt<S,N1,Base>,SInt<N2> >::ModResult Mod;
  static const char C = NL::Compare<SInt<N2>,Mod>::value;
  typedef typename __GCD<SInt<N2>,Mod,C>::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SBigInt<S,N1,Base>, SInt<N2>, -1> 
: public __GCD<SInt<N2>, SBigInt<S,N1,Base>, 1> {};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SInt<N2>, SBigInt<S,N1,Base>, 1> 
{
   static const int_t N = Evaluate2Int<SBigInt<S,N1,Base>, int_t>::Value;
   static const char C = (N2<N) ? -1 : (N2>N) ? 1 : 0;
   typedef typename __GCD<SInt<N2>,SInt<N>,C>::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SInt<N2>, SBigInt<S,N1,Base>, -1> 
: public __GCD<SBigInt<S,N1,Base>, SInt<N2>, 1> {};

template<int_t N1, int_t N2>
struct __GCD<SInt<N1>, SInt<N2>, 1> {
   static const int_t N = N1%N2;
   static const char C = (N2<N) ? -1 : (N2>N) ? 1 : 0;
   typedef typename __GCD<SInt<N2>,SInt<N>,C>::Result Result;
};

template<int_t N1, int_t N2>
struct __GCD<SInt<N1>, SInt<N2>, -1> 
: public __GCD<SInt<N2>, SInt<N1>, 1> {};

template<bool S1, bool S2, class N, base_t Base>
struct __GCD<SBigInt<S1,N,Base>, SBigInt<S2,Loki::NullType,Base>, 1> {
   typedef SBigInt<S1,N,Base> Result;
};

template<bool S, class N, base_t Base>
struct __GCD<SBigInt<S,N,Base>, SInt<0>, 1> {
   typedef SBigInt<S,N,Base> Result;
};

template<int_t N>
struct __GCD<SInt<N>, SInt<0>, 1> {
   typedef SInt<N> Result;
};


template<class N1, class N2>
class GCD {
  typedef typename Abs<N1>::Result AN1;
  typedef typename Abs<N2>::Result AN2;
  static const char C = NL::Compare<AN1,AN2>::value;
public:
  typedef typename __GCD<AN1,AN2,C>::Result Result;
};


template<class N, class D>
class Simplify<SFraction<N,D> > {
   typedef typename Abs<N>::Result AN;
   typedef typename GCD<AN,D>::Result T;
   typedef typename Div<AN,T>::DivResult AbsNum;
   typedef typename Div<D,T>::DivResult Den;

   static const int SN = Sign<N>::value;
   typedef typename Loki::Select<(SN >= 0), AbsNum,
           typename Negate<AbsNum>::Result>::Result Num;
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

template<class N, class D, int_t Num>
class Mult<SFraction<N,D>,SInt<Num> > {
   typedef typename Mult<N,SInt<Num> >::Result Numer;
public:
   typedef SFraction<Numer,D> Result;
};

template<class N, class D, int_t Num>
class Mult<SInt<Num>, SFraction<N,D> > : public Mult<SFraction<N,D>,SInt<Num> > {};

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

template<class N, class D, int_t Num>
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

///////////////////////////////////////////////

template<class N, class D>
struct Check<SFraction<N,D> >
{
  typedef typename Check<N>::Result CheckN;
  typedef typename Check<D>::Result CheckD;
  typedef SFraction<CheckN,CheckD> Result;
};

#endif

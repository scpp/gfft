/***************************************************************************
 *   Copyright (C) 2008-2014 by Vladimir Mirnyy                            *
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

#ifndef __srational_h
#define __srational_h

#include "sbigint.h"

template<class Num, class Den>
struct SRational {
   typedef Num Numer;
   typedef Den Denom;
};

typedef SRational<SInt<1>,SInt<1> > UnitRational;


///Greatest common divisor of A and B

template<class A, class B, int C>
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
  static const int C = NL::Compare<TN2,Mod>::value;
  typedef typename __GCD<TN2,Mod,C>::Result Result;
};

template<bool S1, class N1, base_t B1, bool S2, class N2, base_t B2>
struct __GCD<SBigInt<S1,N1,B1>, SBigInt<S2,N2,B2>, -1>
: public __GCD<SBigInt<S2,N2,B2>, SBigInt<S1,N1,B1>, 1> {};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SBigInt<S,N1,Base>, SInt<N2>, 1> {
  typedef typename Div<SBigInt<S,N1,Base>,SInt<N2> >::ModResult Mod;
  static const int C = NL::Compare<SInt<N2>,Mod>::value;
  typedef typename __GCD<SInt<N2>,Mod,C>::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SBigInt<S,N1,Base>, SInt<N2>, -1> 
: public __GCD<SInt<N2>, SBigInt<S,N1,Base>, 1> {};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SInt<N2>, SBigInt<S,N1,Base>, 1> 
{
   static const int_t N = Evaluate2Int<SBigInt<S,N1,Base>, int_t>::Value;
   static const int C = (N2<N) ? -1 : (N2>N) ? 1 : 0;
   typedef typename __GCD<SInt<N2>,SInt<N>,C>::Result Result;
};

template<bool S, class N1, base_t Base, int_t N2>
struct __GCD<SInt<N2>, SBigInt<S,N1,Base>, -1> 
: public __GCD<SBigInt<S,N1,Base>, SInt<N2>, 1> {};

template<int_t N1, int_t N2>
struct __GCD<SInt<N1>, SInt<N2>, 1> {
   static const int_t N = N1%N2;
   static const int C = (N2<N) ? -1 : (N2>N) ? 1 : 0;
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
struct __GCD<SBigInt<S,N,Base>, SInt<1>, 1> {
   typedef SInt<1> Result;
};

template<int_t N>
struct __GCD<SInt<N>, SInt<1>, 1> {
   typedef SInt<1> Result;
};

template<bool S, class N, base_t Base>
struct __GCD<SBigInt<S,N,Base>, SInt<0>, 1> {
   typedef SBigInt<S,N,Base> Result;
};

template<int_t N>
struct __GCD<SInt<N>, SInt<0>, 1> {
   typedef typename __SwitchToBigInt<N>::Result Result;
};


template<class N1, class N2>
class GCD {
  typedef typename Abs<N1>::Result AN1;
  typedef typename Abs<N2>::Result AN2;
  static const int C = NL::Compare<AN1,AN2>::value;
public:
  typedef typename __GCD<AN1,AN2,C>::Result Result;
};

template<bool S1, class T1, bool S2, class T2, base_t B>
class GCD<SBigInt<S1,Loki::Typelist<SInt<0>,T1>,B>,
          SBigInt<S2,Loki::Typelist<SInt<0>,T2>,B> > {
  typedef typename GCD<SBigInt<S1,T1,B>,SBigInt<S2,T2,B> >::Result Next;
  typedef typename CreateBigInt<Next,B>::Result NextBig;
  typedef Loki::Typelist<SInt<0>,typename NextBig::Num> NList;
public:
  typedef SBigInt<true,NList,B> Result;
};



template<class N, class D>
class Simplify<SRational<N,D> > {
   typedef typename Abs<N>::Result AN;
   typedef typename GCD<AN,D>::Result T;
   typedef typename Div<N,T>::DivResult Num;
   typedef typename Div<D,T>::DivResult Den;
public:
   typedef SRational<Num,Den> Result;
};

/// Multiplication of two compile-time Rationals
template<class N1, class D1, class N2, class D2>
class Mult<SRational<N1,D1>, SRational<N2,D2> > {
   typedef typename GCD<N1,D2>::Result G1;
   typedef typename GCD<N2,D1>::Result G2;
   typedef typename Div<N1,G1>::DivResult NN1;
   typedef typename Div<N2,G2>::DivResult NN2;
   typedef typename Div<D2,G1>::DivResult DD2;
   typedef typename Div<D1,G2>::DivResult DD1;
   typedef typename Mult<NN1,NN2>::Result Num;
   typedef typename Mult<DD1,DD2>::Result Den;
public:
   typedef SRational<Num,Den> Result;
};

template<class N, class D>
class Mult<SRational<N,D>, SRational<N,D> > {
   typedef typename Mult<N,N>::Result Num;
   typedef typename Mult<D,D>::Result Den;
public:
   typedef SRational<Num,Den> Result;
};

template<int_t N1, int_t D1, int_t N2, int_t D2>
class Mult<SRational<SInt<N1>,SInt<D1> >, SRational<SInt<N2>,SInt<D2> > > {
   static const int_t Num = N1*N2;
   static const int_t Den = D1*D2;
   typedef typename GCD<SInt<Num>, SInt<Den> >::Result T;  // T must remain SInt
   typedef typename __SwitchToBigInt<Num/T::value>::Result NumT;
   typedef typename __SwitchToBigInt<Den/T::value>::Result DenT;
public:
   typedef SRational<NumT,DenT> Result;
};

template<int_t N, int_t D>
class Mult<SRational<SInt<N>,SInt<D> >, SRational<SInt<N>,SInt<D> > > {
   static const int_t Num = N*N;
   static const int_t Den = D*D;
   typedef typename __SwitchToBigInt<Num>::Result NumT;
   typedef typename __SwitchToBigInt<Den>::Result DenT;
public:
   typedef SRational<NumT,DenT> Result;
};

template<class N, class D, int_t Num>
class Mult<SRational<N,D>,SInt<Num> > {
   typedef typename GCD<SInt<Num>,D>::Result G;
   typedef typename Div<SInt<Num>,G>::DivResult Num1;
   typedef typename Div<D,G>::DivResult D1;
   typedef typename Mult<N,Num1>::Result Numer;
public:
   typedef SRational<Numer,D1> Result;
};

template<class N, class D, int_t Num>
class Mult<SInt<Num>, SRational<N,D> > : public Mult<SRational<N,D>,SInt<Num> > {};

template<class N, class D, bool S, class NList, unsigned int Base>
class Mult<SRational<N,D>,SBigInt<S,NList,Base> > {
   typedef SBigInt<S,NList,Base> Num;
   typedef typename GCD<Num,D>::Result G;
   typedef typename Div<Num,G>::DivResult Num1;
   typedef typename Div<D,G>::DivResult D1;
   typedef typename Mult<N,Num1>::Result Numer;
public:
   typedef SRational<Numer,D1> Result;
};

template<class N, class D, bool S, class NList, unsigned int Base>
class Mult<SBigInt<S,NList,Base>,SRational<N,D> > : public Mult<SRational<N,D>,SBigInt<S,NList,Base> > {};

/////////////////////////////////////////////////////////////

template<class Numer, class Denom>
class Negate<SRational<Numer,Denom> > {
   typedef typename Negate<Numer>::Result NewNumer;
public:
   typedef SRational<NewNumer,Denom> Result;
};

/////////////////////////////////////////////////////////////

template<class N1, class D1, class N2, class D2>
class Add<SRational<N1,D1>, SRational<N2,D2> > {
   typedef typename GCD<D1,D2>::Result G;
   typedef typename Div<D1,G>::DivResult DD1;
   typedef typename Div<D2,G>::DivResult DD2;
   typedef typename Mult<N1,DD2>::Result T1;
   typedef typename Mult<N2,DD1>::Result T2;
   typedef typename Add<T1,T2>::Result Num;
   typedef typename Mult<DD1,D2>::Result Den;
public:
   typedef SRational<Num,Den> Result;
};

template<class N, class D, int_t Num>
class Add<SRational<N,D>,SInt<Num> > {
   typedef SInt<Num> T;
   typedef typename Add<
           typename Mult<D,T>::Result,N>::Result Numer;
public:
   typedef SRational<Numer,D> Result;
};

template<class N, class D, int_t Num>
class Add<SInt<Num>, SRational<N,D> > : public Add<SRational<N,D>,SInt<Num> > {};


template<class N, class D, bool S, class NList, base_t Base>
class Add<SRational<N,D>,SBigInt<S,NList,Base> > {
   typedef SBigInt<S,NList,Base> T;
   typedef typename Add<
           typename Mult<D,T>::Result,N>::Result Numer;
public:
   typedef SRational<Numer,D> Result;
};

template<class N, class D, bool S, class NList, base_t Base>
class Add<SBigInt<S,NList,Base>,SRational<N,D> > : public Add<SRational<N,D>,SBigInt<S,NList,Base> > {};

///////////////////////////////////////////////

template<class N1, class D1, class N2, class D2>
class Sub<SRational<N1,D1>, SRational<N2,D2> > 
: public Add<SRational<N1,D1>, typename Negate<SRational<N2,D2> >::Result> { };

template<class N, class D, int_t Num>
class Sub<SRational<N,D>,SInt<Num> > : public Add<SRational<N,D>,SInt<-Num> > { };

template<class N, class D, int_t Num>
class Sub<SInt<Num>, SRational<N,D> > : public Add<SInt<Num>, typename Negate<SRational<N,D> >::Result> { };

///////////////////////////////////////////////

template<class N, class D>
struct Check<SRational<N,D> >
{
  typedef typename Check<N>::Result CheckN;
  typedef typename Check<D>::Result CheckD;
  typedef SRational<CheckN,CheckD> Result;
};

#endif

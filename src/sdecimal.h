/***************************************************************************
 *   Copyright (C) 2008-2015 by Vladimir Mirnyy                            *
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

#ifndef __sdecimal_h
#define __sdecimal_h

#include "sbigint.h"


template<class BigInt, long_t NDecPlaces, base_t DecBase = DefaultDecimalBase>
struct SDecimal {
   typedef BigInt Num;
   static const long_t NDec = NDecPlaces;
   static const base_t Base = DecBase;
};


template<class BI, unsigned int N, unsigned int I=N>
struct ShiftLeftRound;

template<bool S, class H1, class H2, class T, base_t Base, unsigned int N, unsigned int I>
struct ShiftLeftRound<SBigInt<S,Loki::Typelist<H1,Loki::Typelist<H2,T> >,Base>,N,I>
{
  static const base_t HalfBase = Base >> 1;
  typedef typename Loki::Select<(H1::value >= HalfBase),
     Loki::Typelist<long_<H2::value+1>,T>,Loki::Typelist<H2,T> >::Result TList;
  typedef typename ShiftLeftRound<SBigInt<S,TList,Base>,N,I-1>::Result Result;
};

template<bool S, class H1, class H2, class T, base_t Base, unsigned int N>
struct ShiftLeftRound<SBigInt<S,Loki::Typelist<H1,Loki::Typelist<H2,T> >,Base>,N,0>
{
  typedef SBigInt<S,Loki::Typelist<H1,Loki::Typelist<H2,T> >,Base> Result;
};

template<bool S, class H1, base_t Base, unsigned int N, unsigned int I>
struct ShiftLeftRound<SBigInt<S,Loki::Typelist<H1,Loki::NullType>,Base>,N,I>
{
  static const base_t HalfBase = Base >> 1;
  typedef typename Loki::Select<(H1::value >= HalfBase),
     Loki::Typelist<long_<1>,Loki::NullType>,Loki::NullType>::Result TList;
  typedef typename ShiftLeftRound<SBigInt<S,TList,Base>,N,I-1>::Result Result;
};

template<bool S, class H1, base_t Base, unsigned int N>
struct ShiftLeftRound<SBigInt<S,Loki::Typelist<H1,Loki::NullType>,Base>,N,0>
{
  typedef SBigInt<S,Loki::Typelist<H1,Loki::NullType>,Base> Result;
};

template<bool S, class TList, base_t Base, unsigned int N>
struct ShiftLeftRound<SBigInt<S,TList,Base>,N,0>
{
   typedef SBigInt<S,TList,Base> Result;
};

template<bool S, base_t Base, unsigned int N, unsigned int I>
struct ShiftLeftRound<SBigInt<S,Loki::NullType,Base>,N,I>
{
   typedef SBigInt<S,Loki::NullType,Base> Result;
};

template<bool S, base_t Base, unsigned int N>
struct ShiftLeftRound<SBigInt<S,Loki::NullType,Base>,N,0>
{
   typedef SBigInt<S,Loki::NullType,Base> Result;
};


template<class BI, long_t ND, int Accuracy, base_t Base>
struct Reduce<SDecimal<BI,ND,Base>,Accuracy,Base> {
//   typedef typename BI::Num NList;
//   typedef typename Loki::Select<(ND>Accuracy),
//           typename Loki::TL::ShiftLeft<NList,ND-Accuracy>::Result,NList>::Result NewList;
//   typedef SDecimal<SBigInt<BI::isPositive,NewList,BI::Base>,Accuracy,Base> Result;

  typedef typename Loki::Select<(ND>Accuracy),
          typename ShiftLeftRound<BI,ND-Accuracy>::Result,BI>::Result NewBI;
  typedef SDecimal<NewBI,Accuracy,Base> Result;
};

template<long_t N, long_t ND, int Accuracy, base_t Base>
struct Reduce<SDecimal<long_<N>,ND,Base>,Accuracy,Base> {
  typedef SDecimal<long_<N>,Accuracy,Base> Result;
};

template<long_t ND, int Accuracy, base_t Base>
struct Reduce<SDecimal<long_<0>,ND,Base>,Accuracy,Base> {
  typedef long_<0> Result;
};


template<class BI1, long_t ND1, class BI2, long_t ND2, base_t DecBase>
class Mult<SDecimal<BI1,ND1,DecBase>,SDecimal<BI2,ND2,DecBase> > {
  static const long_t MaxND = (ND1 > ND2) ? ND1 : ND2;
  typedef typename Mult<BI1,BI2>::Result Prod;
  typedef SDecimal<Prod,ND1+ND2,DecBase> NewDec;
public:
  typedef typename Reduce<NewDec,MaxND,DecBase>::Result Result;
};

template<long_t N, class BI, long_t ND, base_t DecBase>
class Mult<long_<N>,SDecimal<BI,ND,DecBase> > {
  typedef typename Mult<long_<N>,BI>::Result Prod;
public:
  typedef SDecimal<Prod,ND,DecBase> Result;
};

template<long_t N, class BI, long_t ND, base_t DecBase>
class Mult<SDecimal<BI,ND,DecBase>,long_<N> > 
: public Mult<long_<N>,SDecimal<BI,ND,DecBase> > {};


/////////////////////////////////////////////////////////////

template<class BI1, class BI2, long_t ND, base_t DecBase>
class Add<SDecimal<BI1,ND,DecBase>,SDecimal<BI2,ND,DecBase> > {
  typedef typename Add<BI1,BI2>::Result Sum;
public:
  typedef SDecimal<Sum,ND,DecBase> Result;
};

template<class BI, long_t ND, base_t DecBase, long_t N>
class Add<SDecimal<BI,ND,DecBase>,long_<N> > {
  typedef typename CreateBigInt<long_<N>,DecBase>::Result BI1;
  typedef typename Loki::TL::ShiftRight<typename BI1::Num,ND,long_<0> >::Result NList;
  typedef SBigInt<BI1::isPositive,NList,DecBase> NewBI;
  typedef typename Add<BI,NewBI>::Result Sum;
public:
  typedef SDecimal<Sum,ND,DecBase> Result;
};

template<class BI, long_t ND, base_t DecBase, long_t N>
class Add<long_<N>, SDecimal<BI,ND,DecBase> > 
: public Add<SDecimal<BI,ND,DecBase>,long_<N> > {};

template<class BI, long_t ND, base_t DecBase>
class Add<long_<0>, SDecimal<BI,ND,DecBase> > 
{
public:
    typedef SDecimal<BI,ND,DecBase> Result;
};

template<class BI, long_t ND, base_t DecBase>
class Add<SDecimal<BI,ND,DecBase>,long_<0> >
: public Add<long_<0>, SDecimal<BI,ND,DecBase> > {};

///////////////////////////////////////////////

template<class BI1, class BI2, long_t ND, base_t DecBase>
class Sub<SDecimal<BI1,ND,DecBase>,SDecimal<BI2,ND,DecBase> > {
  typedef typename Sub<BI1,BI2>::Result Dif;
public:
  typedef SDecimal<Dif,ND,DecBase> Result;
};

template<class BI, long_t ND, base_t DecBase, long_t N>
class Sub<long_<N>, SDecimal<BI,ND,DecBase> > 
: public Add<long_<N>, typename Negate<SDecimal<BI,ND,DecBase> >::Result> {};

template<class BI, long_t ND, base_t DecBase, long_t N>
class Sub<SDecimal<BI,ND,DecBase>, long_<N> > 
: public Add<SDecimal<BI,ND,DecBase>, long_<-N> > {};

///////////////////////////////////////////////

template<class BI, long_t ND, base_t DecBase>
class Negate<SDecimal<BI,ND,DecBase> > {
  typedef typename Negate<BI>::Result NewBI;
public:
  typedef SDecimal<NewBI,ND,DecBase> Result;
};

///////////////////////////////////////////////

template<class BI, long_t ND, base_t DecBase>
struct Check<SDecimal<BI,ND,DecBase> > : public Check<BI> {};

///////////////////////////////////////////////

template<class SDec>
struct DoubleAccuracy;

template<class BI, long_t ND, base_t DecBase>
struct DoubleAccuracy<SDecimal<BI,ND,DecBase> > {
  typedef typename Loki::TL::ShiftRight<typename BI::Num,ND,long_<0> >::Result NList;
  typedef SBigInt<BI::isPositive,NList,DecBase> NewBI;
  typedef SDecimal<NewBI,ND+ND,DecBase> Result;
};

#endif

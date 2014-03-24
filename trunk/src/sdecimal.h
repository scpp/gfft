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

#ifndef __sdecimal_h
#define __sdecimal_h

#include "sbigint.h"


template<class BigInt, int_t NDecPlaces, base_t DecBase>
struct SDecimal {
   typedef BigInt Num;
   static const int_t NDec = NDecPlaces;
   static const base_t Base = DecBase;
};


template<class BI, unsigned int N, unsigned int I=N>
struct ShiftLeftRound;

template<bool S, class H1, class H2, class T, base_t Base, unsigned int N, unsigned int I>
struct ShiftLeftRound<SBigInt<S,Loki::Typelist<H1,Loki::Typelist<H2,T> >,Base>,N,I>
{
  static const base_t HalfBase = Base >> 1;
  typedef typename Loki::Select<(H1::value >= HalfBase),
     Loki::Typelist<SInt<H2::value+1>,T>,Loki::Typelist<H2,T> >::Result TList;
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
     Loki::Typelist<SInt<1>,Loki::NullType>,Loki::NullType>::Result TList;
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


template<class BI, int_t ND, int Accuracy, base_t Base>
struct Reduce<SDecimal<BI,ND,Base>,Accuracy,Base> {
//   typedef typename BI::Num NList;
//   typedef typename Loki::Select<(ND>Accuracy),
//           typename Loki::TL::ShiftLeft<NList,ND-Accuracy>::Result,NList>::Result NewList;
//   typedef SDecimal<SBigInt<BI::isPositive,NewList,BI::Base>,Accuracy,Base> Result;

  typedef typename Loki::Select<(ND>Accuracy),
          typename ShiftLeftRound<BI,ND-Accuracy>::Result,BI>::Result NewBI;
  typedef SDecimal<NewBI,Accuracy,Base> Result;
};

template<int_t N, int_t ND, int Accuracy, base_t Base>
struct Reduce<SDecimal<SInt<N>,ND,Base>,Accuracy,Base> {
  typedef SDecimal<SInt<N>,Accuracy,Base> Result;
};

template<class BI1, int_t ND1, class BI2, int_t ND2, base_t DecBase>
class Mult<SDecimal<BI1,ND1,DecBase>,SDecimal<BI2,ND2,DecBase> > {
  static const int_t MaxND = (ND1 > ND2) ? ND1 : ND2;
  typedef typename Mult<BI1,BI2>::Result Prod;
  typedef SDecimal<Prod,ND1+ND2,DecBase> NewDec;
public:
  typedef typename Reduce<NewDec,MaxND,DecBase>::Result Result;
};

template<int_t N, class BI, int_t ND, base_t DecBase>
class Mult<SInt<N>,SDecimal<BI,ND,DecBase> > {
  typedef typename Mult<SInt<N>,BI>::Result Prod;
public:
  typedef SDecimal<Prod,ND,DecBase> Result;
};

template<int_t N, class BI, int_t ND, base_t DecBase>
class Mult<SDecimal<BI,ND,DecBase>,SInt<N> > 
: public Mult<SInt<N>,SDecimal<BI,ND,DecBase> > {};

/////////////////////////////////////////////////////////////

template<class BI1, class BI2, int_t ND, base_t DecBase>
class Add<SDecimal<BI1,ND,DecBase>,SDecimal<BI2,ND,DecBase> > {
  typedef typename Add<BI1,BI2>::Result Sum;
public:
  typedef SDecimal<Sum,ND,DecBase> Result;
};

template<class BI, int_t ND, base_t DecBase, int_t N>
class Add<SDecimal<BI,ND,DecBase>,SInt<N> > {
  typedef typename CreateBigInt<SInt<N>,DecBase>::Result BI1;
  typedef typename Loki::TL::ShiftRight<typename BI1::Num,ND,SInt<0> >::Result NList;
  typedef SBigInt<BI1::isPositive,NList,DecBase> NewBI;
  typedef typename Add<BI,NewBI>::Result Sum;
public:
  typedef SDecimal<Sum,ND,DecBase> Result;
};

///////////////////////////////////////////////

template<class BI1, class BI2, int_t ND, base_t DecBase>
class Sub<SDecimal<BI1,ND,DecBase>,SDecimal<BI2,ND,DecBase> > {
  typedef typename Sub<BI1,BI2>::Result Dif;
public:
  typedef SDecimal<Dif,ND,DecBase> Result;
};

///////////////////////////////////////////////

template<class BI, int_t ND, base_t DecBase>
class Negate<SDecimal<BI,ND,DecBase> > {
  typedef typename Negate<BI>::Result NewBI;
public:
  typedef SDecimal<NewBI,ND,DecBase> Result;
};

///////////////////////////////////////////////

template<class BI, int_t ND, base_t DecBase>
struct Check<SDecimal<BI,ND,DecBase> > : public Check<BI> {};

#endif

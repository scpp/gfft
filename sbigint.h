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

#ifndef __sbigint_h
#define __sbigint_h

typedef unsigned int uint;
#define IntT uint
#include "Numlist.h"

#include "loki/TypeManip.h"

template<int N>
struct SInt {
   enum { Value = N };
};

template<bool S, class NList,
         unsigned int B/*=(1<<(sizeof(IntT)*4))*/>
struct SBigInt {
   enum { isPositive = S };
   enum { Base = B };
   typedef NList Num;
};

template<class F1, class F2>
class Mult;

template<class F1, class F2>
class Add;

template<class F1, class F2>
class Sub;

template<class F>
class Negate;

template<class F1, class F2>
class Div;

template<class F1, class F2>
class Mod;

/////////////////////////////////////////////////////////////

template<class NList, unsigned int Base, IntT Rest=0> struct Align;

template<IntT H, class T, unsigned int Base, IntT Rest>
struct Align<NL::Numlist<H,T>,Base,Rest> {
   enum { A = H+Rest };
   typedef NL::Numlist<A%Base,
      typename Align<T,Base,A/Base>::Result> Result;
};

template<unsigned int Base, IntT Rest>
struct Align<NL::NullType,Base,Rest> {
   typedef NL::Numlist<Rest%Base,
      typename Align<NL::NullType,Base,Rest/Base>::Result>  Result;
};

template<unsigned int Base>
struct Align<NL::NullType,Base,0> {
   typedef NL::NullType Result;
};

/////////////////////////////////////////////////////////////

template<int N, bool isSmallEnough=(N<(1<<(sizeof(IntT)*4)))>
struct __SwitchToBigInt;

template<int N>
struct __SwitchToBigInt<N,true> {
   typedef SInt<N> Result;
};

template<int N>
struct __SwitchToBigInt<N,false> {
   enum { Base = (1<<(sizeof(IntT)*4)) };
   typedef typename Align<
      NL::Numlist<N,NL::NullType>,Base>::Result Result;
};

template<int N1, int N2>
class Add<SInt<N1>, SInt<N2> > {
public:
   typedef typename __SwitchToBigInt<N1+N2>::Result Result;
};

template<int N1, int N2>
class Sub<SInt<N1>, SInt<N2> > {
public:
   typedef SInt<N1-N2> Result;
};

template<int N>
class Negate<SInt<N> > {
public:
   typedef SInt<-N> Result;
};

template<int N1, int N2>
class Mult<SInt<N1>, SInt<N2> > {
public:
   typedef typename __SwitchToBigInt<N1*N2>::Result Result;
};

template<int N1, int N2>
class Div<SInt<N1>, SInt<N2> > {
public:
   typedef SInt<N1/N2> Result;
};

template<int N1, int N2>
class Mod<SInt<N1>, SInt<N2> > {
public:
   typedef SInt<N1%N2> Result;
};

/////////////////////////////////////////////////////////////
/*
template<bool S, class NList>
class Negate<SBigInt<S,NList> > {
   enum { Base = 1<<(sizeof(IntT)*4) };
   static const unsigned L = NL::Length<NList>::value;
   static const bool B = NL::NumAt<NList,L-1>::value!=(Base-1);
   typedef typename Loki::Select<(B && S),typename NL::Append<NList,0>::Result,NList>::Result NList1;
   static const unsigned L1 = NL::Length<NList1>::value;
   typedef typename NL::InitNumlist<L1,Base-1>::Result MaxList;
   typedef typename NL::AddAt<MaxList,0,1>::Result MaxList1;
   typedef typename NL::Sub<MaxList1,NList1>::Result NListModif;
   typedef typename Align<NListModif>::Result ANList;
public:
   typedef SBigInt<!S,ANList> Result;
};
*/

/////////////////////////////////////////////////////////////
/*
template<class B1, class B2, bool C>
class __Add;

template<bool S1, class NList1, bool S2, class NList2>
class __Add<SBigInt<S1,NList1>,SBigInt<S2,NList2>,true> {
   typedef typename NL::Add<NList1,NList2>::Result Sum;
   typedef typename Align<Sum>::Result ASum;
public:
   typedef SBigInt<S1,ASum> Result;
};

template<bool S1, class NList1, bool S2, class NList2>
class __Add<SBigInt<S1,NList1>,SBigInt<S2,NList2>,false> {
public:
   enum { Base = 1<<(sizeof(IntT)*4) };
   static const unsigned L1 = NL::Length<NList1>::value;
   static const unsigned L2 = NL::Length<NList2>::value;
   static const unsigned ND = L1-L2;
   typedef typename NL::FillNumlist<NList1,L2,0>::Result NList11;
   static const unsigned ML = (L1>L2) ? L1 : L2;
   typedef typename NL::AddAt<
           typename NL::AddConst<NList11,Base-1>::Result,0,1>::Result NList1modif;
   typedef typename NL::Sub<NList1modif,NList2>::Result Dif;
   typedef typename Align<Dif>::Result ADif;
   static const bool S = (NL::Length<ADif>::value>ML);
   typedef typename NL::EraseAt<ADif,ML>::Result ADifMod;

   typedef SBigInt<S,ADifMod> Result;
//   typedef SBigInt<S,ADif> Result;
};
*/
template<bool S, class NList, unsigned int Base>
class Negate<SBigInt<S,NList,Base> > {
public:
   typedef SBigInt<!S,NList,Base> Result;
};

template<class B1, class B2,
         bool C=NL::Greater<typename B1::Num,typename B2::Num>::value>
class __Add;

template<class NList1, class NList2, unsigned int Base>
class __Add<SBigInt<true,NList1,Base>,SBigInt<true,NList2,Base>,true> {
   typedef typename NL::Add<NList1,NList2>::Result Sum;
   typedef typename Align<Sum,Base>::Result ASum;
public:
   typedef SBigInt<true,ASum,Base> Result;
};

template<class NList1, class NList2, unsigned int Base>
class __Add<SBigInt<true,NList1,Base>,SBigInt<false,NList2,Base>,true> {
//   enum { Base = 1<<(sizeof(IntT)*4) };
   enum { L = NL::Length<NList1>::value };
   typedef typename NL::AddAt<
           typename NL::AddConst<NList1,Base-1>::Result,0,1>::Result NList12;
   typedef typename NL::Sub<NList12,NList2>::Result Dif;
   typedef typename Align<Dif,Base>::Result ADif;
   typedef typename NL::EraseAt<ADif,L>::Result ADif1;
public:
   typedef SBigInt<true,ADif1,Base> Result;
};

template<class NList1, class NList2, unsigned int Base>
class __Add<SBigInt<false,NList1,Base>,SBigInt<true,NList2,Base>,true> {
   typedef SBigInt<true,NList1,Base> BI1;
   typedef SBigInt<false,NList2,Base> BI2;
public:
   typedef typename Negate<
           typename __Add<BI1,BI2>::Result>::Result Result;
};

template<class NList1, class NList2, unsigned int Base>
class __Add<SBigInt<false,NList1,Base>,SBigInt<false,NList2,Base>,true> {
   typedef SBigInt<true,NList1,Base> BI1;
   typedef SBigInt<true,NList2,Base> BI2;
public:
   typedef typename Negate<
           typename __Add<BI1,BI2>::Result>::Result Result;
};

template<bool S1, bool S2, class NList1, class NList2, unsigned int Base>
class __Add<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base>,false> {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<S2,NList2,Base> BI2;
public:
   typedef typename __Add<BI2,BI1,true>::Result Result;
};

//////////////////////////////////////////////////////////

template<bool S1, class NList1, bool S2, class NList2, unsigned int Base>
class Add<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<S2,NList2,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, unsigned int Base, int N>
class Add<SBigInt<S,NList,Base>,SInt<N> > {
   static const IntT A = (N<0) ? -N : N;
   static const bool S1 = (N>=0);
   typedef NL::Numlist<A,NL::NullType> NList1;
   typedef SBigInt<S,NList,Base> BI1;
   typedef SBigInt<S1,NList1,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, unsigned int Base, int N>
class Add<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Add<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

template<bool S1, class NList1, bool S2, class NList2, unsigned int Base>
class Sub<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<!S2,NList2,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, unsigned int Base, int N>
class Sub<SBigInt<S,NList,Base>,SInt<N> > {
   static const IntT A = (N<0) ? -N : N;
   static const bool S1 = (N>=0);
   typedef NL::Numlist<A,NL::NullType> NList1;
   typedef SBigInt<S,NList,Base> BI1;
   typedef SBigInt<!S1,NList1,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, unsigned int Base, int N>
class Sub<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Negate<
           typename Sub<SBigInt<S,NList,Base>,SInt<N> >
              ::Result>::Result Result;
};

///////////////////////////////////////////////////////////////

template<class NList1, class NList2, int I=0>
struct __MultLoop;

template<IntT Num, class Tail, class NList2, int I>
struct __MultLoop<NL::Numlist<Num,Tail>,NList2,I> {
private:
   typedef typename NL::MultConst<NList2,Num>::Result Prod;
   typedef typename NL::InitNumlist<I,0>::Result Shift;
   typedef typename NL::AppendList<Shift,Prod>::Result ShiftedProd;
public:
   typedef typename NL::Add<ShiftedProd,
           typename __MultLoop<Tail,NList2,I+1>::Result>::Result Result;
};

template<class NList2, int I>
struct __MultLoop<NL::NullType,NList2,I> {
   typedef NL::NullType Result;
};

template<bool S1, class NList1, bool S2, class NList2, unsigned int Base>
class Mult<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef typename __MultLoop<NList1,NList2>::Result NListProd;
   typedef typename Align<NListProd,Base>::Result ANListProd;
public:
   typedef SBigInt<(S1==S2),ANListProd,Base> Result;
};

template<bool S, class NList, unsigned int Base, int N>
class Mult<SBigInt<S,NList,Base>,SInt<N> > {
   enum { A = (N<0) ? -N : N };
   static const bool S1 = (N>=0);
   typedef typename NL::MultConst<NList,A>::Result Prod;
   typedef typename Align<Prod,Base>::Result AProd;
public:
   typedef SBigInt<(S1==S),AProd,Base> Result;
};

template<bool S, class NList, unsigned int Base, int N>
class Mult<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Mult<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

////////////////////////////////////////////////////////

template<class B, unsigned int NewBase>
class Translate;

template<bool S, IntT H, class T, unsigned int Base, unsigned int NewBase>
class Translate<SBigInt<S,NL::Numlist<H,T>,Base>,NewBase> {
   typedef SBigInt<S,T,NewBase> BI;
   typedef typename Translate<BI,NewBase>::Result Next;
public:
   typedef typename Add<typename Mult<Next,SInt<Base> >::Result,SInt<H> >::Result Result;
};

template<bool S, unsigned int Base, unsigned int NewBase>
class Translate<SBigInt<S,NL::NullType,Base>,NewBase> {
public:
   typedef SBigInt<S,NL::NullType,NewBase> Result;
};

#endif

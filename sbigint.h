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

/*! \file
    \brief Compile-time big integer implementation based on Numlist compile-time arrays
*/

typedef unsigned int uint;
#define IntT uint
#include "Numlist.h"

#include "sint.h"

#include "loki/TypeManip.h"


/// \brief Big integer number metacontainer.
/// \param S sign of the big integer
/// \param NList Numlist containing series of digits in the Base numerical system,
///              first element is least significant
/// \param Base Numerical system base for the big integer representation
/////////////////////////////////////////////////////////////////////////
template<bool S, class NList,
         unsigned int B=(1<<(sizeof(IntT)*4))>
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

template<class F>
class Abs;

template<class C>
class Simplify;

/// \class Align
/// \brief Division of every element in NList by Base.
/// Reminder is assigned to a digit and quotient is added to a higher digit.
/// This common operation is needed after addition and multiplication of
/// big integers.
/// \param NList Numlist containing series of digits
/// \param Base Numerical system base
/// \param Rest Quotient to carry to a higher digit
/// \return Numlist with the digits of a big integer in system Base
/////////////////////////////////////////////////////////////////////////
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

/// \brief Compile-time addition of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing sum of N1 and N2.
/// This can be big integer, if N1+N2 exceeds the predefined type IntT
////////////////////////////////////////////////////////////
template<int N1, int N2>
class Add<SInt<N1>, SInt<N2> > {
public:
   typedef typename __SwitchToBigInt<N1+N2>::Result Result;
};

/// \brief Compile-time subtraction of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing difference of N1 and N2.
////////////////////////////////////////////////////////////
template<int N1, int N2>
class Sub<SInt<N1>, SInt<N2> > {
public:
   typedef SInt<N1-N2> Result;
};

/// \brief Compile-time negation of an integer.
/// \param N an integer
/// \return container class representing (-N)
////////////////////////////////////////////////////////////
template<int N>
class Negate<SInt<N> > {
public:
   typedef SInt<-N> Result;
};

/// \brief Compile-time multiplication of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing product of N1 and N2.
/// This can be big integer, if N1*N2 exceeds the predefined type IntT
////////////////////////////////////////////////////////////
template<int N1, int N2>
class Mult<SInt<N1>, SInt<N2> > {
public:
   typedef typename __SwitchToBigInt<N1*N2>::Result Result;
};

/// \brief Quotient from division of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing quotient from division of N1 by N2.
////////////////////////////////////////////////////////////
template<int N1, int N2>
class Div<SInt<N1>, SInt<N2> > {
public:
   typedef SInt<N1/N2> Result;
};

/// \brief Reminder from division of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing reminder from division of N1 by N2.
////////////////////////////////////////////////////////////
template<int N1, int N2>
class Mod<SInt<N1>, SInt<N2> > {
public:
   typedef SInt<N1%N2> Result;
};

/// \brief Absolute value of an integer.
/// \param N an integer
/// \return container class representing absolute value of N.
////////////////////////////////////////////////////////////
template<int N>
class Abs<SInt<N> > {
   enum { AN = (N>0) ? N : -N };
public:
   typedef SInt<AN> Result;
};

/// \brief Absolute value of a big integer.
/// \param S sign of the big integer
/// \param NList Numlist containing series of digits
/// \param Base Numerical system base
/// \return container class representing absolute value of the big integer
////////////////////////////////////////////////////////////
template<bool S, class NList, unsigned int Base>
class Abs<SBigInt<S,NList,Base> > {
public:
   typedef SBigInt<true,NList,Base>  Result;
};

/// \brief Compile-time negation of a big integer.
/// \param S sign of the big integer
/// \param NList Numlist containing series of digits
/// \param Base Numerical system base
/// \return container class representing the
/// negative big integer SBigInt<-S,NList,Base>
////////////////////////////////////////////////////////////
template<bool S, class NList, unsigned int Base>
class Negate<SBigInt<S,NList,Base> > {
public:
   typedef SBigInt<!S,NList,Base> Result;
};

template<bool S, IntT H, class T, unsigned int Base>
class Simplify<SBigInt<S,NL::Numlist<H,T>,Base> > {
public:
   typedef typename Simplify<SBigInt<S,T,Base> >::Result Result;
};

template<bool S, IntT H, unsigned int Base>
class Simplify<SBigInt<S,NL::Numlist<H,NL::NullType>,Base> > {
public:
   typedef SBigInt<S,NL::Numlist<H,NL::NullType>,Base> Result;
};

template<bool S, unsigned int Base>
class Simplify<SBigInt<S,NL::Numlist<0,NL::NullType>,Base> > {
public:
   typedef SBigInt<S,NL::NullType,Base> Result;
};

template<bool S, unsigned int Base>
class Simplify<SBigInt<S,NL::NullType,Base> > {
public:
   typedef SBigInt<S,NL::NullType,Base> Result;
};

///////////////////////////////////////////////////////

template<class B1, class B2,
         bool C=(NL::Compare<typename B1::Num,typename B2::Num>::value>0)>
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

template<class NList1, class NList2, unsigned int Base, int I>
class __Div;

template<IntT H, class T, class NList2, unsigned int Base, int I>
class __Div<NL::Numlist<H,T>,NList2,Base,I> {
   typedef __Div<T,NList2,Base,I-1> Next;
   typedef NL::Numlist<H,typename Next::UList> NList1;
   enum { L1 = NL::Length<NList1>::value };
   enum { L2 = NL::Length<NList2>::value };
   enum { V1 = NL::NumAt<NList2,L2-1>::value };
   enum { V2 = NL::NumAt<NList2,L2-2>::value };
   enum { U0 = NL::NumAtNonStrict<NList1,L1-I>::value };
   enum { U1 = NL::NumAtNonStrict<NList1,L1-I-1>::value };
   enum { U2 = NL::NumAtNonStrict<NList1,L1-I-2>::value };
   enum { U = U0*Base+U1 };
   static const unsigned Q = U/V1;  // trial quotient
   static const unsigned R = U%V1;  // trial reminder
   // Test trial Q
   static const unsigned Q1 = ((Q==Base) || (Q*V2>R*Base+U2)) ? Q-1 : Q;
   static const unsigned R1 = (Q1<Q) ? R+V1 : R;
   // Repeat, if (R1<Base)
   static const unsigned Q2 = ((R1<Base) && ((Q1==Base) || (Q1*V2>R1*Base+U2))) ? Q1-1 : Q1;
   typedef SBigInt<true,NList1,Base> BI1;
   typedef SBigInt<true,NList2,Base> BI2;
   typedef typename Sub<BI1,typename Mult<BI2,SInt<Q2> >::Result>::Result Dif;
public:
//typedef typename NL::Print<T1>::Result DBG;
//typedef typename NL::Print<SInt<Q> >::Result DBG;
   typedef NL::Numlist<Q2,typename Next::DivList> DivList;
   typedef typename Dif::Num UList;
};

template<IntT H, class T, class NList2, unsigned int Base>
class __Div<NL::Numlist<H,T>,NList2,Base,0> {
public:
   typedef NL::NullType DivList;
   typedef NL::Numlist<H,T> UList;
};

template<bool S1, class NList1, bool S2, class NList2, unsigned int Base>
class Div<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
public:
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<S2,NList2,Base> BI2;
   static const unsigned L1 = NL::Length<NList1>::value;
   static const unsigned L2 = NL::Length<NList2>::value;
   // Normalization
   enum { D = Base/(NL::NumAt<NList2,L2-1>::value+1) };
   typedef typename Mult<BI1,SInt<D> >::Result U;
   static const unsigned ULen = NL::Length<typename U::Num>::value;
   typedef typename Loki::Select<(ULen==L1),
           typename NL::Append<typename U::Num,0>::Result,
           typename U::Num>::Result UList;
   typedef typename Mult<BI2,SInt<D> >::Result V;
   // Loop
   typedef __Div<UList,typename V::Num,Base,L1-L2+1> Loop;

/*  typedef typename Div<SBigInt<true,
           typename NL::Range<NList1,0,L2-1>::Result,Base>,SInt<D> >::DivResult ModResult;*/
   typedef typename Div<SBigInt<true,typename Loop::UList,Base>,SInt<D> >::DivResult ModResult;
   typedef SBigInt<(S1==S2),typename Loop::DivList,Base> DivResult;
//   typedef SBigInt<true,typename Loop::ModList,Base> ModResult;
};

template<bool S1, bool S2, class NList2, unsigned int Base>
class Div<SBigInt<S1,NL::NullType,Base>,SBigInt<S2,NList2,Base> > {
public:
   typedef SBigInt<true,NL::NullType,Base> Result;
};

//////////////

template<bool S, IntT H, class T, unsigned int Base, int N>
class Div<SBigInt<S,NL::Numlist<H,T>,Base>,SInt<N> > {
   typedef Div<SBigInt<S,T,Base>,SInt<N> > Next;
   enum { AN = (N>0) ? N : -N };
   enum { B = Next::Q*Base + H };
   enum { R = B/AN };
public:
   enum { Q = B%AN };
   typedef SBigInt<(S==(N>0)),NL::Numlist<R,
     typename Next::DivResult::Num>,Base> DivResult;
   typedef SInt<Q> ModResult;
};

template<bool S, unsigned int Base, int N>
class Div<SBigInt<S,NL::NullType,Base>,SInt<N> > {
public:
   enum { Q = 0 };
   typedef SBigInt<true,NL::NullType,Base> DivResult;
};

template<bool S, class NList, unsigned int Base, int N>
class Div<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Div<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

////////////////////////////////////////////////////////

template<bool S1, class NList1, bool S2, class NList2, unsigned int Base>
class Mod<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {

};

template<bool S, class NList, unsigned int Base, int N>
class Mod<SBigInt<S,NList,Base>,SInt<N> > {

};

template<bool S, class NList, unsigned int Base, int N>
class Mod<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Mod<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
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

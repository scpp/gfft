/***************************************************************************
 *   Copyright (C) 2008-2013 by Volodymyr Myrnyy                                *
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

#include "sint.h"

#include "numtypelist.h"
#include "typelistext.h"

#include "loki/Typelist.h"
#include "loki/TypeManip.h"
#include "loki/TypeTraits.h"

#include <iostream>

typedef unsigned int base_t;
typedef long IntT;

//static const base_t DefaultBase = (1<<(sizeof(IntT)*4));
// for 32bit
//static const base_t DefaultBase = 65536;
// for 64bit
static const base_t DefaultBase = 2147483648;

/// \brief Big integer number metacontainer.
/// \param S sign of the big integer
/// \param NList Numlist containing series of digits in the Base numerical system,
///              first element is least significant
/// \param Base Numerical system base for the big integer representation
/////////////////////////////////////////////////////////////////////////
template<bool S, class NList,
         base_t B = DefaultBase>
struct SBigInt {
   static const bool isPositive = S;
   static const base_t Base = B;
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



/////////////////////////////////////////////////////////////////////////
// bool isInt = Loki::TypeTraits<RetType>::isIntegral,
// bool isFloat = Loki::TypeTraits<RetType>::isFloat
template<class BInt, class RetType, unsigned long AccumBase = 1>
struct Evaluate2Int;

template<bool S, class H, class T, base_t Base, class RetType, unsigned long AccumBase>
struct Evaluate2Int<SBigInt<S,Loki::Typelist<H,T>,Base>,RetType,AccumBase>
{
  static const RetType Next = Evaluate2Int<SBigInt<S,T,Base>,RetType,AccumBase*Base>::Value;
  static const RetType Value = H::Value*AccumBase + Next;
};

template<bool S, base_t Base, class RetType, unsigned long AccumBase>
struct Evaluate2Int<SBigInt<S,Loki::NullType,Base>,RetType,AccumBase>
{
  static const RetType Value = 0;
};

template<int_t N, class RetType>
struct Evaluate2Int<SInt<N>, RetType>
{
  static const RetType Value = N;
};


template<class BInt, class RetType>
struct Evaluate2Float;

template<bool S, class H, class T, base_t Base, class RetType>
struct Evaluate2Float<SBigInt<S,Loki::Typelist<H,T>,Base>,RetType>
{
  static RetType value(const RetType AccumBase = 1) 
  {
    return (RetType)H::Value*AccumBase + Evaluate2Float<SBigInt<S,T,Base>,RetType>::value(AccumBase*Base);
  }
};

template<bool S, base_t Base, class RetType>
struct Evaluate2Float<SBigInt<S,Loki::NullType,Base>,RetType>
{
  static RetType value(const RetType AccumBase = 1) { return 0; }
};

template<int_t N, class RetType>
struct Evaluate2Float<SInt<N>,RetType>
{
  static RetType value(const RetType AccumBase = 1) { return N; }
};


namespace NL 
{

template<bool S1, class N1, bool S2, class N2, base_t B>
struct Compare<SBigInt<S1,N1,B>, SBigInt<S2,N2,B> > 
{
   static const char c = Compare<N1,N2>::value;
   static const char value = (S1 && S2) ? c : (!S1 && !S2) ? -c 
                        : (S1 && !S2) ? 1 : -1;
};
  
template<int_t N1, int_t N2>
struct Compare<SInt<N1>, SInt<N2> > 
{
   static const char value = (N1 > N2) ? 1 : (N1 < N2) ? -1 : 0; 
};
  
template<bool S1, class N1, base_t B1, int_t N2>
struct Compare<SBigInt<S1,N1,B1>, SInt<N2> > 
{
   static const int_t N = Evaluate2Int<SBigInt<S1,N1,B1>, int_t>::Value;
   static const char value = Compare<SInt<N>, SInt<N2> >::value;
};

template<bool S1, class N1, base_t B1, int_t N2>
struct Compare<SInt<N2>, SBigInt<S1,N1,B1> > 
{
   static const int_t N = Evaluate2Int<SBigInt<S1,N1,B1>, int_t>::Value;
   static const char value = Compare<SInt<N2>, SInt<N> >::value;
};
  
}



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
template<class NList, base_t Base, IntT Rest=0> 
struct Align;

template<class H, class T, base_t Base, IntT Rest>
struct Align<Loki::Typelist<H,T>,Base,Rest> {
   static const IntT A = H::Value+Rest;
   typedef Loki::Typelist<SInt<A%Base>,
      typename Align<T,Base,A/Base>::Result> Result;
};

template<base_t Base, IntT Rest>
struct Align<Loki::NullType,Base,Rest> {
   typedef Loki::Typelist<SInt<Rest%Base>,
      typename Align<Loki::NullType,Base,Rest/Base>::Result>  Result;
};

template<base_t Base>
struct Align<Loki::NullType,Base,0> {
   typedef Loki::NullType Result;
};

/////////////////////////////////////////////////////////////

template<IntT N, bool isSmallEnough=(N<DefaultBase)>
struct __SwitchToBigInt;

template<IntT N>
struct __SwitchToBigInt<N,true> {
   typedef SInt<N> Result;
};

template<IntT N>
struct __SwitchToBigInt<N,false> {
   static const bool S = (N>=0);
   static const base_t Base = DefaultBase;
   typedef typename Align<
      Loki::Typelist<SInt<N>,Loki::NullType>,Base>::Result NList;
   typedef SBigInt<S,NList,Base> Result;
};

/// \brief Compile-time addition of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing sum of N1 and N2.
/// This can be big integer, if N1+N2 exceeds the predefined type IntT
////////////////////////////////////////////////////////////
template<IntT N1, IntT N2>
struct Add<SInt<N1>, SInt<N2> > {
   typedef typename __SwitchToBigInt<N1+N2>::Result Result;
};
template<int_t N1, int_t N2>
struct Add<SInt<N1>, SInt<N2> > {
   typedef typename __SwitchToBigInt<N1+N2>::Result Result;
};

/// \brief Compile-time subtraction of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing difference of N1 and N2.
////////////////////////////////////////////////////////////
template<IntT N1, IntT N2>
class Sub<SInt<N1>, SInt<N2> > {
public:
   typedef SInt<N1-N2> Result;
};

/// \brief Compile-time negation of an integer.
/// \param N an integer
/// \return container class representing (-N)
////////////////////////////////////////////////////////////
template<IntT N>
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
template<IntT N1, IntT N2>
struct Mult<SInt<N1>, SInt<N2> > {
   typedef typename __SwitchToBigInt<N1*N2>::Result Result;
};

template<int_t N1, int_t N2>
struct Mult<SInt<N1>, SInt<N2> > {
   typedef typename __SwitchToBigInt<N1*N2>::Result Result;
};

/// \brief Quotient from division of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing quotient from division of N1 by N2.
////////////////////////////////////////////////////////////
template<IntT N1, IntT N2>
struct Div<SInt<N1>, SInt<N2> > {
   typedef SInt<N1/N2> DivResult;
   typedef SInt<N1%N2> ModResult;
};

template<int_t N1, int_t N2>
struct Div<SInt<N1>, SInt<N2> > {
   typedef SInt<N1/N2> DivResult;
   typedef SInt<N1%N2> ModResult;
};

/// \brief Reminder from division of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing reminder from division of N1 by N2.
////////////////////////////////////////////////////////////
template<IntT N1, IntT N2>
struct Mod<SInt<N1>, SInt<N2> > {
   typedef SInt<N1%N2> Result;
};

/// \brief Absolute value of an integer.
/// \param N an integer
/// \return container class representing absolute value of N.
////////////////////////////////////////////////////////////
template<IntT N>
class Abs<SInt<N> > {
   static const IntT AN = (N>0) ? N : -N ;
public:
   typedef SInt<AN> Result;
};

/// \brief Absolute value of a big integer.
/// \param S sign of the big integer
/// \param NList Numlist containing series of digits
/// \param Base Numerical system base
/// \return container class representing absolute value of the big integer
////////////////////////////////////////////////////////////
template<bool S, class NList, base_t Base>
struct Abs<SBigInt<S,NList,Base> > {
   typedef SBigInt<true,NList,Base>  Result;
};

/// \brief Compile-time negation of a big integer.
/// \param S sign of the big integer
/// \param NList Numlist containing series of digits
/// \param Base Numerical system base
/// \return container class representing the
/// negative big integer SBigInt<-S,NList,Base>
////////////////////////////////////////////////////////////
template<bool S, class NList, base_t Base>
struct Negate<SBigInt<S,NList,Base> > {
   typedef SBigInt<!S,NList,Base> Result;
};

//////////////////////////////////////////////////////
template<typename NList>
struct NZeros;

template<int_t N, typename Tail>
struct NZeros<Loki::Typelist<SInt<N>,Tail> > {
  static const int_t Value = NZeros<Tail>::Value;
};

template<typename Tail>
struct NZeros<Loki::Typelist<SInt<0>,Tail> > {
  static const int_t Value = NZeros<Tail>::Value + 1;
};

template<>
struct NZeros<Loki::NullType> {
  static const int_t Value = 0;
};


template<class NList,
int L = Loki::TL::Length<NList>::value - NZeros<NList>::Value>
struct SimplifyLoop;

template<class H, class Tail, int L>
struct SimplifyLoop<Loki::Typelist<H,Tail>,L> {
  typedef Loki::Typelist<H, typename SimplifyLoop<Tail>::Result> Result;
};

template<class H, class Tail>
struct SimplifyLoop<Loki::Typelist<H,Tail>,0> {
  typedef Loki::NullType Result;
};

template<>
struct SimplifyLoop<Loki::NullType,0> {
  typedef Loki::NullType Result;
};

template<class C>
struct Simplify;

template<bool S, typename NList, base_t Base>
struct Simplify<SBigInt<S,NList,Base> > {
   typedef SBigInt<S,typename SimplifyLoop<NList>::Result,Base> Result;
};

template<IntT N>
struct Simplify<SInt<N> > {
   typedef SInt<N> Result;
};

template<int_t N>
struct Simplify<SInt<N> > {
   typedef SInt<N> Result;
};

///////////////////////////////////////////////////////

template<class B1, class B2,
         bool C=(NL::Compare<typename B1::Num,typename B2::Num>::value>=0)>
class __Add;

template<class NList1, class NList2, base_t Base>
class __Add<SBigInt<true,NList1,Base>,SBigInt<true,NList2,Base>,true> {
   typedef typename NL::Add<NList1,NList2>::Result Sum;
   typedef typename Align<Sum,Base>::Result ASum;
public:
   typedef SBigInt<true,ASum,Base> Result;
};

template<class NList1, class NList2, base_t Base>
class __Add<SBigInt<true,NList1,Base>,SBigInt<false,NList2,Base>,true> {
   static const IntT L = Loki::TL::Length<NList1>::value;
   typedef typename NL::AddAt<
           typename NL::AddConst<NList1,SInt<Base-1> >::Result,0,SInt<1> >::Result NList12;
   typedef typename NL::Sub<NList12,NList2>::Result Dif;
   typedef typename Align<Dif,Base>::Result ADif;
   typedef typename Loki::TL::EraseAt<ADif,L>::Result ADif1;
public:
   typedef SBigInt<true,ADif1,Base> Result;
};

template<class NList1, class NList2, base_t Base>
class __Add<SBigInt<false,NList1,Base>,SBigInt<true,NList2,Base>,true> {
   typedef SBigInt<true,NList1,Base> BI1;
   typedef SBigInt<false,NList2,Base> BI2;
public:
   typedef typename Negate<
           typename __Add<BI1,BI2>::Result>::Result Result;
};

template<class NList1, class NList2, base_t Base>
class __Add<SBigInt<false,NList1,Base>,SBigInt<false,NList2,Base>,true> {
   typedef SBigInt<true,NList1,Base> BI1;
   typedef SBigInt<true,NList2,Base> BI2;
public:
   typedef typename Negate<
           typename __Add<BI1,BI2>::Result>::Result Result;
};

template<bool S1, bool S2, class NList1, class NList2, base_t Base>
class __Add<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base>,false> {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<S2,NList2,Base> BI2;
public:
   typedef typename __Add<BI2,BI1>::Result Result;
};

//////////////////////////////////////////////////////////

template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Add<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<S2,NList2,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, base_t Base, int_t N>
class Add<SBigInt<S,NList,Base>,SInt<N> > {
   static const IntT A = (N<0) ? -N : N;
   static const bool S1 = (N>=0);
   typedef Loki::Typelist<SInt<A>,Loki::NullType> NList1;
   typedef SBigInt<S,NList,Base> BI1;
   typedef SBigInt<S1,NList1,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, base_t Base, IntT N>
class Add<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Add<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};
template<bool S, class NList, base_t Base, int_t N>
class Add<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Add<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Sub<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<!S2,NList2,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, base_t Base, IntT N>
class Sub<SBigInt<S,NList,Base>,SInt<N> > {
   static const IntT A = (N<0) ? -N : N;
   static const bool S1 = (N>=0);
   typedef Loki::Typelist<SInt<A>,Loki::NullType> NList1;
   typedef SBigInt<S,NList,Base> BI1;
   typedef SBigInt<!S1,NList1,Base> BI2;
public:
   typedef typename __Add<BI1,BI2>::Result Result;
};

template<bool S, class NList, base_t Base, IntT N>
class Sub<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Negate<
           typename Sub<SBigInt<S,NList,Base>,SInt<N> >
              ::Result>::Result Result;
};

///////////////////////////////////////////////////////////////

template<class NList1, class NList2, IntT I=0>
struct __MultLoop;

template<class Num, class Tail, class NList2, IntT I>
struct __MultLoop<Loki::Typelist<Num,Tail>,NList2,I> {
private:
   typedef typename NL::MultConst<NList2,Num>::Result Prod;
   typedef typename Loki::TL::Repeat<SInt<0>,I>::Result Shift;
   typedef typename Loki::TL::Append<Shift,Prod>::Result ShiftedProd;
public:
   typedef typename NL::Add<ShiftedProd,
           typename __MultLoop<Tail,NList2,I+1>::Result>::Result Result;
};

template<class NList2, IntT I>
struct __MultLoop<Loki::NullType,NList2,I> {
   typedef Loki::NullType Result;
};

template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Mult<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef typename __MultLoop<NList1,NList2>::Result NListProd;
   typedef typename Align<NListProd,Base>::Result ANListProd;
public:
   typedef SBigInt<(S1==S2),ANListProd,Base> Result;
};

template<bool S, class NList, base_t Base, IntT N>
class Mult<SBigInt<S,NList,Base>,SInt<N> > {
   static const IntT A = (N<0) ? -N : N ;
   static const bool S1 = (N>=0);
   typedef typename NL::MultConst<NList,SInt<A> >::Result Prod;
   typedef typename Align<Prod,Base>::Result AProd;
public:
   typedef SBigInt<(S1==S),AProd,Base> Result;
};

template<bool S, class NList, base_t Base, int_t N>
class Mult<SBigInt<S,NList,Base>,SInt<N> > {
   static const IntT A = (N<0) ? -N : N ;
   static const bool S1 = (N>=0);
   typedef typename NL::MultConst<NList,SInt<A> >::Result Prod;
   typedef typename Align<Prod,Base>::Result AProd;
public:
   typedef SBigInt<(S1==S),AProd,Base> Result;
};

template<bool S, class NList, base_t Base, IntT N>
class Mult<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Mult<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};
template<bool S, class NList, base_t Base, int_t N>
class Mult<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Mult<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

////////////////////////////////////////////////////////

template<class NList1, class NList2, base_t Base, IntT I>
class __Div;

template<class H, class T, class NList2, base_t Base, IntT I>
class __Div<Loki::Typelist<H,T>,NList2,Base,I> {
   typedef __Div<T,NList2,Base,I-1> Next;
   typedef Loki::Typelist<H,typename Next::UList> NList1;
   static const unsigned int L1 = Loki::TL::Length<NList1>::value;
   static const unsigned int L2 = Loki::TL::Length<NList2>::value;
   static const int_t V1 = Loki::TL::TypeAt<NList2,L2-1>::Result::value;
   static const int_t V2 = Loki::TL::TypeAt<NList2,L2-2>::Result::value;
   static const int_t U0 = Loki::TL::TypeAtNonStrict<NList1,L1-I,SInt<0> >::Result::value;
   static const int_t U1 = Loki::TL::TypeAtNonStrict<NList1,L1-I-1,SInt<0> >::Result::value;
   static const int_t U2 = Loki::TL::TypeAtNonStrict<NList1,L1-I-2,SInt<0> >::Result::value;
   static const int_t U = U0*Base+U1;
   static const int_t Q = U/V1;  // trial quotient
   static const int_t R = U%V1;  // trial reminder
   // Test trial Q
   static const int_t Q1 = ((Q==Base) || (Q*V2>R*Base+U2)) ? Q-1 : Q;
   static const int_t R1 = (Q1<Q) ? R+V1 : R;
   // Repeat, if (R1<Base)
   static const int_t Q2 = ((R1<Base) && ((Q1==Base) || (Q1*V2>R1*Base+U2))) ? Q1-1 : Q1;
   typedef SBigInt<true,NList1,Base> BI1;
   typedef SBigInt<true,NList2,Base> BI2;
   typedef typename Sub<BI1,typename Mult<BI2,SInt<Q2> >::Result>::Result Dif;
public:
//typedef typename Loki::Print<T1>::Result DBG;
//typedef typename Loki::Print<SInt<Q> >::Result DBG;
   typedef Loki::Typelist<SInt<Q2>,typename Next::DivList> DivList;
   typedef typename Dif::Num UList;
};

template<class H, class T, class NList2, base_t Base>
class __Div<Loki::Typelist<H,T>,NList2,Base,0> {
public:
   typedef Loki::NullType DivList;
   typedef Loki::Typelist<H,T> UList;
};



template<class B1, class B2>
class BigDiv;


template<class B1, class B2,
int_t L1 = Loki::TL::Length<typename B1::Num>::value,
int_t L2 = Loki::TL::Length<typename B2::Num>::value>
struct DivSelect {
  typedef BigDiv<B1,B2> T;
  typedef typename T::DivResult DivResult;
  typedef typename T::ModResult ModResult;
};

template<class B1, class B2, int_t L>
struct DivSelect<B1,B2,1,L> {
  typedef SInt<Evaluate2Int<B1,int_t>::Value> I1;
  typedef Div<I1,B2> T;
  typedef typename T::DivResult DivResult;
  typedef typename T::ModResult ModResult;
};

template<class B1, class B2, int_t L>
struct DivSelect<B1,B2,L,1> {
  typedef SInt<Evaluate2Int<B2,int_t>::Value> I2;
  typedef Div<B1,I2> T;
  typedef typename T::DivResult DivResult;
  typedef typename T::ModResult ModResult;
};

template<class B1, class B2>
struct DivSelect<B1,B2,1,1> {
  typedef SInt<Evaluate2Int<B1,int_t>::Value> I1;
  typedef SInt<Evaluate2Int<B2,int_t>::Value> I2;
  typedef Div<I1,I2> T;
  typedef typename T::DivResult DivResult;
  typedef typename T::ModResult ModResult;
};


template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class BigDiv<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
public:
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<S2,NList2,Base> BI2;
   static const unsigned L1 = Loki::TL::Length<NList1>::value;
   static const unsigned L2 = Loki::TL::Length<NList2>::value;
   // Normalization
   static const IntT D = Base/(Loki::TL::TypeAt<NList2,L2-1>::Result::value+1);
   typedef typename Mult<BI1,SInt<D> >::Result U;
   static const unsigned ULen = Loki::TL::Length<typename U::Num>::value;
   typedef typename Loki::Select<(ULen==L1),
           typename Loki::TL::Append<typename U::Num,SInt<0> >::Result,
           typename U::Num>::Result UList;
   typedef typename Mult<BI2,SInt<D> >::Result V;
   // Loop
   typedef __Div<UList,typename V::Num,Base,L1-L2+1> Loop;
   //typedef SInt<Evaluate2Int<BI2,int_t>::Value> I2;

/*  typedef typename Div<SBigInt<true,
           typename Loki::Range<NList1,0,L2-1>::Result,Base>,SInt<D> >::DivResult ModResult;*/
   typedef typename Div<SBigInt<true,typename Loop::UList,Base>,SInt<D> >::DivResult ModResult;
	     
   typedef SBigInt<(S1==S2),typename Loop::DivList,Base> DivResult;
//   typedef SBigInt<true,typename Loop::ModList,Base> ModResult;
};



template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Div<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > 
{
public:
  typedef DivSelect<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > T;
  typedef typename Simplify<typename T::DivResult>::Result DivResult;
  typedef typename Simplify<typename T::ModResult>::Result ModResult;
};


template<bool S1, bool S2, class NList2, base_t Base>
class Div<SBigInt<S1,Loki::NullType,Base>,SBigInt<S2,NList2,Base> > {
public:
//    typedef SBigInt<true,Loki::NullType,Base> Result;
   typedef SInt<0> ModResult;
   typedef SInt<0> DivResult;
};

template<bool S1, class NList1, bool S2, base_t Base>
class Div<SBigInt<S1,NList1,Base>,SBigInt<S2,Loki::NullType,Base> > {}; // error, division by zero

//////////////

template<bool S, class H, class T, base_t Base, IntT N>
class Div<SBigInt<S,Loki::Typelist<H,T>,Base>,SInt<N> > {
   typedef Div<SBigInt<S,T,Base>,SInt<N> > Next;
   static const IntT AN = (N>0) ? N : -N;
   static const IntT B = Next::Q*Base + H::Value;
   static const IntT R = B/AN;
public:
   static const IntT Q = B%AN;
   typedef SBigInt<(S==(N>0)),Loki::Typelist<SInt<R>,
     typename Next::DivResult::Num>,Base> DivResult;
   typedef SInt<Q> ModResult;
};

template<bool S, class H, class T, base_t Base, int_t N>
class Div<SBigInt<S,Loki::Typelist<H,T>,Base>,SInt<N> > {
   typedef Div<SBigInt<S,T,Base>,SInt<N> > Next;
   static const IntT AN = (N>0) ? N : -N;
   static const IntT B = Next::Q*Base + H::Value;
   static const IntT R = B/AN;
public:
   static const IntT Q = B%AN;
   typedef SBigInt<(S==(N>0)),Loki::Typelist<SInt<R>,
     typename Next::DivResult::Num>,Base> DivResult;
   typedef SInt<Q> ModResult;
};
template<bool S, class H, class T, base_t Base>
class Div<SBigInt<S,Loki::Typelist<H,T>,Base>,SInt<0> > {};

template<bool S, base_t Base, IntT N>
class Div<SBigInt<S,Loki::NullType,Base>,SInt<N> > {
public:
   static const IntT Q = 0;
//    typedef SBigInt<true,Loki::NullType,Base> DivResult;
   typedef SInt<0> ModResult;
   typedef SInt<0> DivResult;
};
template<bool S, base_t Base, int_t N>
class Div<SBigInt<S,Loki::NullType,Base>,SInt<N> > {
public:
   static const IntT Q = 0;
   typedef SBigInt<true,Loki::NullType,Base> DivResult;
//    typedef SInt<0> ModResult;
//    typedef SInt<0> DivResult;
};

template<bool S, base_t Base, IntT N>
class Div<SInt<N>, SBigInt<S,Loki::NullType,Base> > {}; // error


template<bool S, class NList, base_t Base, IntT N>
class Div<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Div<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};
template<bool S, class NList, base_t Base, int_t N>
class Div<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Div<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

////////////////////////////////////////////////////////

template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Mod<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {

};

template<bool S, class NList, base_t Base, IntT N>
class Mod<SBigInt<S,NList,Base>,SInt<N> > {

};

template<bool S, class NList, base_t Base, IntT N>
class Mod<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Mod<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

////////////////////////////////////////////////////////

template<class B, base_t NewBase>
class Translate;

template<bool S, class H, class T, base_t Base, base_t NewBase>
class Translate<SBigInt<S,Loki::Typelist<H,T>,Base>,NewBase> {
   typedef SBigInt<S,T,Base> BI;
   typedef typename Translate<BI,NewBase>::Result Next;
public:
   typedef typename Add<typename Mult<Next,SInt<Base> >::Result,H>::Result Result;
};

template<bool S, base_t Base, base_t NewBase>
class Translate<SBigInt<S,Loki::NullType,Base>,NewBase> {
public:
   typedef SBigInt<S,Loki::NullType,NewBase> Result;
};


////////////////////////////////////////////////////////

template<class BigInt>
struct Cout;

template<bool S, class H, class T, base_t Base>
struct Cout<SBigInt<S,Loki::Typelist<H,T>,Base> > 
{
  typedef Cout<SBigInt<S,T,Base> > Next;
  
  static void apply(std::ostream& os, const int_t width = 0) { 
    int_t w = width;
    if (w == 0) {
      base_t b = Base-1;
      while (b > 0) { b/=10; ++w; }
    }
    Next::apply(os,w);
    os.fill('0');
    os.width(w);
    os << std::right << H::Value;
  }
};

template<bool S, class H, base_t Base>
struct Cout<SBigInt<S,Loki::Typelist<H,Loki::NullType>,Base> > {
  static void apply(std::ostream& os, const int_t width = 0) { 
    os << H::Value;
  }
};

#endif

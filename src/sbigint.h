/***************************************************************************
 *   Copyright (C) 2008-2013 by Volodymyr Myrnyy                           *
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

//typedef unsigned int base_t;
typedef int_t base_t;

//static const base_t DefaultBase = (1<<(sizeof(int_t)*4));

// for 32bit
//static const base_t DefaultBase        = 65536;
//static const base_t DefaultDecimalBase = 10000;

// for 64bit
//static const base_t DefaultBase        = 2147483648;
static const base_t DefaultBase        = 1000000000;
static const base_t DefaultDecimalBase = 1000000000;


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

template<typename T1, typename T2>
struct Pair {
  typedef T1 first;
  typedef T2 second;
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

template<class F>
class Sign;

template<class N, base_t Base>
struct CreateBigInt;


template<class T, int Accuracy, base_t Base = DefaultBase>
struct Reduce;

template<int_t N, int Accuracy, base_t Base>
struct Reduce<SInt<N>,Accuracy,Base> {
  typedef SInt<N> Result;
};


/////////////////////////////////////////////////////////////////////////
// bool isInt = Loki::TypeTraits<RetType>::isIntegral,
// bool isFloat = Loki::TypeTraits<RetType>::isFloat
template<class BInt, class RetType, unsigned long AccumBase = 1>
struct Evaluate2IntLoop;

template<bool S, class H, class T, base_t Base, class RetType, unsigned long AccumBase>
struct Evaluate2IntLoop<SBigInt<S,Loki::Typelist<H,T>,Base>,RetType,AccumBase>
{
  static const RetType Next = Evaluate2IntLoop<SBigInt<S,T,Base>,RetType,AccumBase*Base>::value;
  static const RetType Value = H::value*AccumBase + Next;
};

template<bool S, base_t Base, class RetType, unsigned long AccumBase>
struct Evaluate2IntLoop<SBigInt<S,Loki::NullType,Base>,RetType,AccumBase>
{
  static const RetType Value = 0;
};

template<class BInt, class RetType>
struct Evaluate2Int;

template<bool S, class NList, base_t Base, class RetType>
struct Evaluate2Int<SBigInt<S,NList,Base>,RetType>
{
  static const RetType v = Evaluate2IntLoop<SBigInt<S,NList,Base>,RetType>::value;
  static const RetType value = S ? v : -v;
};

template<int_t N, class RetType>
struct Evaluate2Int<SInt<N>, RetType>
{
  static const RetType value = N;
};


template<class BInt, class RetType>
struct EvaluateToFloatLoop;

template<bool S, class H, class T, base_t Base, class RetType>
struct EvaluateToFloatLoop<SBigInt<S,Loki::Typelist<H,T>,Base>,RetType>
{
  static RetType value() 
  {
    return (RetType)H::value + Base * EvaluateToFloatLoop<SBigInt<S,T,Base>,RetType>::value();
  }
};

template<bool S, class H, base_t Base, class RetType>
struct EvaluateToFloatLoop<SBigInt<S,Loki::Typelist<H,Loki::NullType>,Base>,RetType>
{
  static RetType value() { return H::value; }
};

template<bool S, base_t Base, class RetType>
struct EvaluateToFloatLoop<SBigInt<S,Loki::NullType,Base>,RetType>
{
  static RetType value() { return 0; }
};


template<class BInt, class RetType>
struct EvaluateToFloat;

template<bool S, class NList, base_t Base, class RetType>
struct EvaluateToFloat<SBigInt<S,NList,Base>,RetType> 
{
  static RetType value() 
  { 
    RetType v = EvaluateToFloatLoop<SBigInt<S,NList,Base>,RetType>::value();
    return S ? v : -v; 
  }
};

template<int_t N, class RetType>
struct EvaluateToFloat<SInt<N>,RetType>
{
  static RetType value() { return static_cast<RetType>(N); }
};


namespace NL 
{

template<bool S1, class N1, bool S2, class N2, base_t B>
struct Compare<SBigInt<S1,N1,B>, SBigInt<S2,N2,B> > 
{
   static const int c = Compare<N1,N2>::value;
   static const int value = (S1 && S2) ? c : (!S1 && !S2) ? -c 
                        : (S1 && !S2) ? 1 : -1;
};
  
template<int_t N1, int_t N2>
struct Compare<SInt<N1>, SInt<N2> > 
{
   static const int value = (N1 > N2) ? 1 : (N1 < N2) ? -1 : 0; 
};
  
template<bool S1, class N1, base_t Base, int_t N2>
struct Compare<SBigInt<S1,N1,Base>, SInt<N2> > 
{
  typedef typename CreateBigInt<SInt<N2>,Base>::Result BI;
  static const int value = Compare<SBigInt<S1,N1,Base>,BI>::value;
};

template<bool S1, class N1, base_t B1, int_t N2>
struct Compare<SInt<N2>, SBigInt<S1,N1,B1> > 
{
   static const int value = -Compare<SBigInt<S1,N1,B1>, SInt<N2> >::value;
};
  
template<class N1, class N2>
struct CompareAbs : public Compare<typename Abs<N1>::Result, typename Abs<N2>::Result> {};


template<class Num>
struct Length;

template<bool S, class NList, base_t B>
struct Length<SBigInt<S,NList,B> > {
  static const int value = Loki::TL::Length<NList>::value;
};

template<int_t N>
struct Length<SInt<N> > {
  static const int value = 1;
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
template<class NList, base_t Base, int_t Rest=0> 
struct Align;

template<class H, class T, base_t Base, int_t Rest>
struct Align<Loki::Typelist<H,T>,Base,Rest> {
   static const int_t A = H::value+Rest;
   typedef Loki::Typelist<SInt<A%Base>,
      typename Align<T,Base,A/Base>::Result> Result;
};

template<base_t Base, int_t Rest>
struct Align<Loki::NullType,Base,Rest> {
   typedef Loki::Typelist<SInt<Rest%Base>,
      typename Align<Loki::NullType,Base,Rest/Base>::Result>  Result;
};

template<base_t Base>
struct Align<Loki::Typelist<SInt<0>,Loki::NullType>,Base,0> {
   typedef Loki::NullType Result;
};

template<base_t Base>
struct Align<Loki::NullType,Base,0> {
   typedef Loki::NullType Result;
};

/////////////////////////////////////////////////////////////

template<class N, base_t Base>
struct CreateBigInt;

template<bool S, class N, base_t Base>
struct CreateBigInt<SBigInt<S,N,Base>,Base> {
  typedef SBigInt<S,N,Base> Result;
};

template<int_t N, base_t Base>
struct CreateBigInt<SInt<N>,Base> {
   static const bool S = (N>=0);
   typedef typename Abs<SInt<N> >::Result AN;
   typedef typename Align<
      Loki::Typelist<AN,Loki::NullType>,Base>::Result NList;
   typedef SBigInt<S,NList,Base> Result;
};

template<base_t Base>
struct CreateBigInt<SInt<0>,Base> {
  typedef SBigInt<true,Loki::NullType,DefaultBase> Result;
};

/////////////////////////////////////////////////////////////

template<int_t N, base_t Base = DefaultBase,
bool isSmallEnough=((N>=0 && N<Base) || (N<0 && -N<Base))>
struct __SwitchToBigInt;

template<int_t N, base_t Base>
struct __SwitchToBigInt<N,Base,true> {
   typedef SInt<N> Result;
};

template<int_t N, base_t Base>
struct __SwitchToBigInt<N,Base,false> {
   typedef typename CreateBigInt<SInt<N>,Base>::Result Result;
};

/////////////////////////////////////////////////////////////

template<class N, 
int_t L = NL::Length<N>::value>
struct __SwitchToInt {
  typedef N Result;
};

template<class N>
struct __SwitchToInt<N,1> {
   typedef SInt<Evaluate2Int<N,int_t>::value> Result;
};

template<class N>
struct __SwitchToInt<N,0> {
   typedef SInt<0> Result;
};

/// \brief Compile-time addition of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing sum of N1 and N2.
/// This can be big integer, if N1+N2 exceeds the predefined type int_t
////////////////////////////////////////////////////////////
template<int_t N1, int_t N2>
struct Add<SInt<N1>, SInt<N2> > {
   typedef typename __SwitchToBigInt<N1+N2>::Result Result;
};

/// \brief Compile-time subtraction of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing difference of N1 and N2.
////////////////////////////////////////////////////////////
template<int_t N1, int_t N2>
struct Sub<SInt<N1>, SInt<N2> > : public Add<SInt<N1>, SInt<-N2> > {};

/// \brief Compile-time negation of an integer.
/// \param N an integer
/// \return container class representing (-N)
////////////////////////////////////////////////////////////
template<int_t N>
struct Negate<SInt<N> > {
   typedef SInt<-N> Result;
};

/// \brief Compile-time multiplication of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing product of N1 and N2.
/// This can be big integer, if N1*N2 exceeds the predefined type int_t
////////////////////////////////////////////////////////////
template<int_t N1, int_t N2>
struct Mult<SInt<N1>, SInt<N2> > {
   typedef typename __SwitchToBigInt<N1*N2>::Result Result;
};

/// \brief Quotient from division of two integers.
/// \param N1 an integer
/// \param N2 an integer
/// \return container class representing quotient from division of N1 by N2.
////////////////////////////////////////////////////////////
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
template<int_t N1, int_t N2>
struct Mod<SInt<N1>, SInt<N2> > {
   typedef SInt<N1%N2> Result;
};

/// \brief Absolute value of an integer.
/// \param N an integer
/// \return container class representing absolute value of N.
////////////////////////////////////////////////////////////
template<int_t N>
struct Abs<SInt<N> > {
   static const int_t AN = (N>0) ? N : -N ;
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
template<int_t N>
struct Sign<SInt<N> > {
   static const char value = (N>0) ? 1 : ((N<0) ? -1 : 0);
};

template<bool S, class NList, base_t Base>
struct Sign<SBigInt<S,NList,Base> > {
   static const char value = S ? 1 : -1;
};

template<bool S, base_t Base>
struct Sign<SBigInt<S,Loki::NullType,Base> > {
   static const char value = 0;
};

///////////////////////////////////////////////////////

template<class C>
struct Simplify;

template<bool S, typename NList, base_t Base>
struct Simplify<SBigInt<S,NList,Base> > {
  typedef typename NL::CutTrailingZeros<NList>::Result NewList;
  typedef typename __SwitchToInt<SBigInt<S,NewList,Base> >::Result Result;
};

template<int_t N>
struct Simplify<SInt<N> > {
  typedef SInt<N> Result;
};

///////////////////////////////////////////////////////

template<class B1, class B2, int Comparison>
struct __Add;

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<true,NList1,Base>,SBigInt<true,NList2,Base>,1> {
private:
   typedef typename NL::Add<NList1,NList2>::Result Sum;
   typedef typename Align<Sum,Base>::Result ASum;
public:
   typedef SBigInt<true,ASum,Base> Result;
};

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<true,NList1,Base>,SBigInt<false,NList2,Base>,1> {
private:
   static const int_t L = Loki::TL::Length<NList1>::value;
   typedef typename NL::AddAt<
           typename NL::AddConst<NList1,SInt<Base-1> >::Result,0,SInt<1> >::Result NList12;
   typedef typename NL::Sub<NList12,NList2>::Result Dif;
   typedef typename Align<Dif,Base>::Result ADif;
   typedef typename Loki::TL::EraseAt<ADif,L>::Result ADif1;
public:
//   typedef typename Simplify<SBigInt<true,ADif1,Base> >::Result Result;
   typedef SBigInt<true,ADif1,Base> Result;
};

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<false,NList1,Base>,SBigInt<true,NList2,Base>,1> {
   typedef typename Negate<
           typename __Add<SBigInt<true,NList1,Base>,
	                  SBigInt<false,NList2,Base>,1>::Result>::Result Result;
};

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<false,NList1,Base>,SBigInt<false,NList2,Base>,1> {
   typedef typename Negate<
           typename __Add<SBigInt<true,NList1,Base>,
	                  SBigInt<true,NList2,Base>,1>::Result>::Result Result;
};

template<bool S1, bool S2, class NList1, class NList2, base_t Base>
struct __Add<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base>,-1> {
   typedef typename __Add<SBigInt<S2,NList2,Base>,SBigInt<S1,NList1,Base>,1>::Result Result;
};

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<true,NList1,Base>,SBigInt<true,NList2,Base>,0> 
: public __Add<SBigInt<true,NList1,Base>,SBigInt<true,NList2,Base>,1> {};

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<true,NList1,Base>,SBigInt<false,NList2,Base>,0> {
   typedef SInt<0> Result;
};

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<false,NList1,Base>,SBigInt<true,NList2,Base>,0> {
   typedef SInt<0> Result;
};

template<class NList1, class NList2, base_t Base>
struct __Add<SBigInt<false,NList1,Base>,SBigInt<false,NList2,Base>,0> 
: public __Add<SBigInt<false,NList1,Base>,SBigInt<false,NList2,Base>,1> {};

//////////////////////////////////////////////////////////

template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Add<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<S2,NList2,Base> BI2;
   static const char C = NL::Compare<NList1,NList2>::value;
public:
   typedef typename __Add<BI1,BI2,C>::Result Result;
};

template<bool S, class NList, base_t Base, int_t N>
class Add<SBigInt<S,NList,Base>,SInt<N> > {
   static const int C = NL::Compare<SBigInt<true,NList,Base>,typename Abs<SInt<N> >::Result>::value;
   typedef SBigInt<S,NList,Base> BI1;
   typedef typename CreateBigInt<SInt<N>,Base>::Result BI2;
public:
   typedef typename __Add<BI1,BI2,C>::Result Result;
};

template<bool S, class NList, base_t Base, int_t N>
class Add<SInt<N>,SBigInt<S,NList,Base> > : public Add<SBigInt<S,NList,Base>,SInt<N> > {};

template<bool S, class NList, base_t Base>
class Add<SBigInt<S,NList,Base>,SInt<0> > {
public:
  typedef SBigInt<S,NList,Base> Result;
};

template<bool S, base_t Base, int_t N>
class Add<SBigInt<S,Loki::NullType,Base>,SInt<N> > {
public:
  typedef SInt<N> Result;
};

  
template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Sub<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef SBigInt<S1,NList1,Base> BI1;
   typedef SBigInt<!S2,NList2,Base> BI2;
   static const char C = NL::Compare<NList1,NList2>::value;
public:
   typedef typename __Add<BI1,BI2,C>::Result Result;
};

template<bool S, class NList, base_t Base, int_t N>
class Sub<SBigInt<S,NList,Base>,SInt<N> > : public Add<SBigInt<S,NList,Base>,SInt<-N> > {};

template<bool S, class NList, base_t Base, int_t N>
class Sub<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Negate<
           typename Sub<SBigInt<S,NList,Base>,SInt<N> >::Result>::Result Result;
};

///////////////////////////////////////////////////////////////

template<class NList1, class NList2, base_t Base, int_t I=0>
struct __MultLoop;

template<class Num, class Tail, class NList2, base_t Base, int_t I>
struct __MultLoop<Loki::Typelist<Num,Tail>,NList2,Base,I> {
private:
   typedef typename NL::MultConst<NList2,Num>::Result Prod;
   typedef typename Loki::TL::Repeat<SInt<0>,I>::Result Shift;
   typedef typename Loki::TL::Append<Shift,Prod>::Result ShiftedProd;
public:
   typedef typename Align<typename NL::Add<ShiftedProd,
           typename __MultLoop<Tail,NList2,Base,I+1>::Result>::Result,Base>::Result Result;
};

template<class NList2, base_t Base, int_t I>
struct __MultLoop<Loki::NullType,NList2,Base,I> {
   typedef Loki::NullType Result;
};



template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Mult<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {
   typedef typename __MultLoop<NList1,NList2,Base>::Result NListProd;
   //typedef typename Align<NListProd,Base>::Result ANListProd;
public:
   typedef SBigInt<(S1==S2),NListProd,Base> Result;
};

template<bool S, class NList, base_t Base, int_t N>
class Mult<SBigInt<S,NList,Base>,SInt<N> > {
   static const int_t A = (N<0) ? -N : N ;
   static const bool S1 = (N>=0);
   typedef typename NL::MultConst<NList,SInt<A> >::Result Prod;
   typedef typename Align<Prod,Base>::Result AProd;
public:
   typedef SBigInt<(S1==S),AProd,Base> Result;
};

template<bool S, class NList, base_t Base>
class Mult<SBigInt<S,NList,Base>,SInt<Base> > {
public:
   typedef SBigInt<S,Loki::Typelist<SInt<0>,NList>,Base> Result;
};

template<bool S, class NList, base_t Base>
class Mult<SBigInt<S,NList,Base>,SInt<-Base> > {
public:
   typedef SBigInt<!S,Loki::Typelist<SInt<0>,NList>,Base> Result;
};

template<bool S, class NList, base_t Base>
class Mult<SBigInt<S,NList,Base>,SInt<1> > {
public:
   typedef SBigInt<S,NList,Base> Result;
};

template<bool S, class NList, base_t Base>
class Mult<SBigInt<S,NList,Base>,SInt<-1> > {
public:
   typedef SBigInt<!S,NList,Base> Result;
};

template<bool S, class NList, base_t Base>
class Mult<SBigInt<S,NList,Base>,SInt<0> > {
public:
   typedef SInt<0> Result;
};

template<bool S, class NList, base_t Base, int_t N>
class Mult<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Mult<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

////////////////////////////////////////////////////////

template<class C>
struct SelectUList;

template<>
struct SelectUList<SInt<0> > {
  typedef Loki::NullType Result;
};

template<bool S, class NList, base_t Base>
struct SelectUList<SBigInt<S,NList,Base> > {
  typedef NList Result;
};


template<class NList1, class NList2, base_t Base, int_t I,
unsigned int NList1Len, int_t V1, int_t V2>
class __Div;

template<class H, class T, class NList2, base_t Base, int_t I,
unsigned int NList1Len, int_t V1, int_t V2>
class __Div<Loki::Typelist<H,T>,NList2,Base,I,NList1Len,V1,V2> {
   typedef __Div<T,NList2,Base,I-1,NList1Len-1,V1,V2> Next;
   typedef Loki::Typelist<H,typename Next::UList> NList1;
   typedef typename Loki::TL::Next<NList1,NList1Len-I-2>::Result UT;
   static const int_t U0 = Loki::TL::TypeAtNonStrict<UT,2,SInt<0> >::Result::value;
   static const int_t U1 = Loki::TL::TypeAtNonStrict<UT,1,SInt<0> >::Result::value;
   static const int_t U2 = Loki::TL::TypeAtNonStrict<UT,0,SInt<0> >::Result::value;
//    static const int_t U0 = Loki::TL::TypeAtNonStrict<NList1,NList1Len-I,SInt<0> >::Result::value;
//    static const int_t U1 = Loki::TL::TypeAtNonStrict<NList1,NList1Len-I-1,SInt<0> >::Result::value;
//    static const int_t U2 = Loki::TL::TypeAtNonStrict<NList1,NList1Len-I-2,SInt<0> >::Result::value;
   static const int_t U = U0*Base+U1;
   static const uint_t Q = U/V1;  // trial quotient
   static const uint_t R = U%V1;  // trial reminder
   // Test trial Q
   static const uint_t Q1 = ((Q==Base) || (Q*V2>R*Base+U2)) ? Q-1 : Q;
   static const uint_t R1 = (Q1<Q) ? R+V1 : R;
   // Repeat, if (R1<Base)
   static const int_t Q2 = ((R1<Base) && ((Q1==Base) || (Q1*V2>R1*Base+U2))) ? Q1-1 : Q1;
   typedef SBigInt<true,NList1,Base> BI1;
   typedef SBigInt<true,NList2,Base> BI2;
   //typedef typename __SwitchToBigInt<Q2>::Result Q2T;
   //typedef typename Sub<BI1,typename Mult<BI2,Q2T>::Result>::Result Dif;
   typedef typename Sub<BI1,typename Mult<BI2,SInt<Q2> >::Result>::Result Dif;
public:
//typedef typename Loki::Print<T1>::Result DBG;
//typedef typename Loki::Print<SInt<Q> >::Result DBG;
   typedef Loki::Typelist<SInt<Q2>,typename Next::DivList> DivList;
   typedef typename SelectUList<Dif>::Result UList;
};

template<class H, class T, class NList2, base_t Base,
unsigned int NList1Len, int_t V1, int_t V2>
class __Div<Loki::Typelist<H,T>,NList2,Base,0,NList1Len,V1,V2> {
public:
   typedef Loki::NullType DivList;
   typedef Loki::Typelist<H,T> UList;
};



template<class B1, class B2>
class BigDiv;


template<class B1, class B2, int_t L1, int_t L2, char C>
struct DivSelect;

template<class B1, class B2, int_t L1, int_t L2>
struct DivSelect<B1,B2,L1,L2,1> {
  typedef typename Simplify<B1>::Result BS1;
  typedef typename Simplify<B2>::Result BS2;
  typedef BigDiv<BS1,BS2> T;
//   typedef BigDiv<B1,B2> T;
  typedef typename T::DivResult DivResult;
  typedef typename T::ModResult ModResult;
};

template<class B1, class B2, int_t L1, int_t L2>
struct DivSelect<B1,B2,L1,L2,-1> {
  typedef SInt<0> DivResult;
  typedef B1 ModResult;         //<<< TODO: modify for negative B1
};

template<class B1, class B2, int_t L1, int_t L2>
struct DivSelect<B1,B2,L1,L2,0> {
  typedef SInt<1> DivResult;
  typedef SInt<0> ModResult;               
};

template<class B1, class B2, int_t L>
struct DivSelect<B1,B2,2,L,1> {
  typedef SInt<Evaluate2Int<B1,int_t>::value> I1;
  typedef SInt<Evaluate2Int<B2,int_t>::value> I2;
  typedef Div<I1,I2> T;
  typedef typename T::DivResult DivResult;
  typedef typename T::ModResult ModResult;
};

template<class B1, class B2>
struct DivSelect<B1,B2,2,1,1> {
  typedef SInt<Evaluate2Int<B1,int_t>::value> I1;
  typedef SInt<Evaluate2Int<B2,int_t>::value> I2;
  typedef Div<I1,I2> T;
  typedef typename T::DivResult DivResult;
  typedef typename T::ModResult ModResult;
};

template<class B1, class B2, int_t L>
struct DivSelect<B1,B2,L,1,1> {
  typedef SInt<Evaluate2Int<B2,int_t>::value> I2;
  typedef Div<B1,I2> T;
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
//typedef typename Loki::TL::Print<SInt<L2> >::Result A;

   // Normalization
   static const int_t D = Base/(Loki::TL::TypeAt<NList2,L2-1>::Result::value+1);
   typedef typename Mult<BI1,SInt<D> >::Result U;
   static const unsigned int ULen = Loki::TL::Length<typename U::Num>::value;
   typedef typename Loki::Select<(ULen==L1),
           typename Loki::TL::Append<typename U::Num,SInt<0> >::Result,
           typename U::Num>::Result UList;
   typedef typename Mult<BI2,SInt<D> >::Result V;
   
   // Loop
   static const unsigned int UListLen = (ULen==L1) ? ULen+1 : ULen;
   typedef typename V::Num VList;
   static const unsigned int VLen = Loki::TL::Length<VList>::value;
   static const int_t V1 = Loki::TL::TypeAt<VList,VLen-1>::Result::value;
   static const int_t V2 = Loki::TL::TypeAt<VList,VLen-2>::Result::value;
   typedef __Div<UList,VList,Base,L1-L2+1,UListLen,V1,V2> Loop;
   //typedef SInt<Evaluate2Int<BI2,int_t>::value> I2;

/*  typedef typename Div<SBigInt<true,
           typename Loki::Range<NList1,0,L2-1>::Result,Base>,SInt<D> >::DivResult ModResult;*/
   typedef SBigInt<true,typename Loop::UList,Base> ModT;
   typedef typename Div<ModT,SInt<D> >::DivResult ModResult;
	     
   typedef SBigInt<(S1==S2),typename Loop::DivList,Base> DivResult;
//   typedef SBigInt<true,typename Loop::ModList,Base> ModResult;
};

template<class List, class N, base_t Base>
struct __DivN;

template<class H, class T, int_t N, base_t Base>
struct __DivN<Loki::Typelist<H,T>,SInt<N>,Base> {
   typedef __DivN<T,SInt<N>,Base> Next;
   static const int_t AN = (N>0) ? N : -N;
   static const int_t B = Next::Q*Base + H::value;
   static const int_t R = B/AN;
   static const int_t Q = B%AN;
   typedef Loki::Typelist<SInt<R>,typename Next::Result> Result;
};

template<int_t N, base_t Base>
struct __DivN<Loki::NullType,SInt<N>,Base> {
   static const int_t Q = 0;
   typedef Loki::NullType Result;
};

//////////////////

template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Div<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > 
{
  static const int_t L1 = Loki::TL::Length<NList1>::value;
  static const int_t L2 = Loki::TL::Length<NList2>::value;
  static const char C = NL::Compare<NList1,NList2>::value;
  
  typedef DivSelect<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base>,L1,L2,C> T;
public:
  typedef typename Simplify<typename T::DivResult>::Result DivResult;
  typedef typename Simplify<typename T::ModResult>::Result ModResult;
//   typedef typename T::DivResult DivResult;
//   typedef typename T::ModResult ModResult;
};


template<bool S1, bool S2, class NList2, base_t Base>
class Div<SBigInt<S1,Loki::NullType,Base>,SBigInt<S2,NList2,Base> > {
public:
   typedef SInt<0> ModResult;
   typedef SInt<0> DivResult;
};

template<bool S1, class NList1, bool S2, base_t Base>
class Div<SBigInt<S1,NList1,Base>,SBigInt<S2,Loki::NullType,Base> > {}; // error, division by zero

//////////////

template<bool S, class H, class T, base_t Base, int_t N>
class Div<SBigInt<S,Loki::Typelist<H,T>,Base>,SInt<N> > {
  typedef __DivN<Loki::Typelist<H,T>,typename Abs<SInt<N> >::Result,Base> DivLoop;
  typedef SBigInt<(S==(N>0)),typename DivLoop::Result,Base> BI;
public:
   typedef typename Simplify<BI>::Result DivResult;
   typedef SInt<DivLoop::Q> ModResult;
};

template<bool S, class H, class T, base_t Base>
class Div<SBigInt<S,Loki::Typelist<H,T>,Base>,SInt<1> > {
public:  
   typedef SBigInt<S,Loki::Typelist<H,T>,Base> DivResult;
   typedef SInt<0> ModResult;
};

template<bool S, class H, class T, base_t Base>
class Div<SBigInt<S,Loki::Typelist<H,T>,Base>,SInt<0> > {};  // error, division by zero

template<bool S, class H, class T, base_t Base>
class Div<SInt<0>, SBigInt<S,Loki::Typelist<H,T>,Base> > {
public:
   typedef SInt<0> DivResult;
   typedef SInt<0> ModResult;
};  

template<bool S, base_t Base, int_t N>
class Div<SBigInt<S,Loki::NullType,Base>,SInt<N> > {
public:
   typedef SInt<0> DivResult;
   typedef SInt<0> ModResult;
};

template<bool S, base_t Base, int_t N>
class Div<SInt<N>, SBigInt<S,Loki::NullType,Base> > {}; // error


template<bool S, class NList, base_t Base, int_t N>
class Div<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Div<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};

////////////////////////////////////////////////////////
/*
template<bool S1, class NList1, bool S2, class NList2, base_t Base>
class Mod<SBigInt<S1,NList1,Base>,SBigInt<S2,NList2,Base> > {

};

template<bool S, class NList, base_t Base, int_t N>
class Mod<SBigInt<S,NList,Base>,SInt<N> > {

};

template<bool S, class NList, base_t Base, int_t N>
class Mod<SInt<N>,SBigInt<S,NList,Base> > {
public:
   typedef typename Mod<SBigInt<S,NList,Base>,SInt<N> >::Result Result;
};
*/

// Returns number of digits in N in the Base-system (Base=2 for binary)
template<unsigned N, unsigned Base>
struct NDigits {
  static const unsigned value = NDigits<N/Base, Base>::value + 1;
};

template<unsigned Base>
struct NDigits<0, Base> {
  static const unsigned value = 0;
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

template<class NList, base_t Base>
struct DoubleBaseLoop;

template<class H1, class H2, class T, base_t Base>
struct DoubleBaseLoop<Loki::Typelist<H1,Loki::Typelist<H2,T> >,Base> {
  static const int_t H = H2::value*Base + H1::value;
  typedef Loki::Typelist<SInt<H>,typename DoubleBaseLoop<T,Base>::Result> Result;
};

template<class H1, base_t Base>
struct DoubleBaseLoop<Loki::Typelist<H1,Loki::NullType>,Base> {
  typedef Loki::Typelist<H1,Loki::NullType> Result;
};

template<base_t Base>
struct DoubleBaseLoop<Loki::NullType,Base> {
  typedef Loki::NullType Result;
};


template<class BI>
struct DoubleBase;

template<bool S, class N, base_t Base>
struct DoubleBase<SBigInt<S,N,Base> > {
  typedef typename DoubleBaseLoop<N,Base>::Result List;
  typedef SBigInt<S,List,Base*Base> Result;
};

template<int_t N>
struct DoubleBase<SInt<N> > {
  typedef SInt<N> Result;
};

////////////////////////////////////////////////////////

template<class H, base_t Base, 
bool C = (Abs<H>::Result::value < Base)>
struct CheckStep;

template<class H, base_t Base>
struct CheckStep<H,Base,true> {
  typedef H Result;
};

template<class H, base_t Base>
struct CheckStep<H,Base,false> { }; // error


template<class T>
struct Check;

template<bool S, class H, class T, base_t Base>
struct Check<SBigInt<S,Loki::Typelist<H,T>,Base> >
{
  typedef Loki::Typelist<typename CheckStep<H,Base>::Result, 
                         typename Check<SBigInt<S,T,Base> >::Result> Result;
};

template<bool S, base_t Base>
struct Check<SBigInt<S,Loki::NullType,Base> >
{
  typedef Loki::NullType Result;
};

template<int_t N>
struct Check<SInt<N> > : public CheckStep<SInt<N>,DefaultBase> {};

#endif

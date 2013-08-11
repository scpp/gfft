/***************************************************************************
 *   Copyright (C) 2007-2013 by Volodymyr Myrnyy                           *
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

#ifndef __metafunc_h
#define __metafunc_h

/** \file
    \brief Compile-time computing of trigonometric functions
*/

#include <cmath>

#include "sfraction.h"
#include "sdecimal.h"

#include "pseudometafunc.h"

template<typename T>
struct TempTypeTrait;

template<>
struct TempTypeTrait<float> {
   typedef double Result;
};

template<>
struct TempTypeTrait<double> {
   typedef long double Result;
};

template<typename T,
template<typename> class Complex>
struct TempTypeTrait<Complex<T> > {
   typedef typename TempTypeTrait<T>::Result Result;
};

template<typename T, typename A,
template<typename,typename> class Complex>
struct TempTypeTrait<Complex<T,A> > {
   typedef typename TempTypeTrait<T>::Result Result;
};

// template<typename T, typename A,
// template<typename,typename> class Complex>
// struct TempTypeTrait<Complex<T,A> > {
//    typedef T Result;
// };

//// To save compile time
// typedef TYPELIST_1(SInt<314159265>) NL1;
// typedef TYPELIST_1(SInt<100000000>) DL1;
// typedef SFraction<SBigInt<true,NL1,DefaultBase>,SBigInt<true,DL1,DefaultBase> > TPi1;
// 
// typedef TYPELIST_2(SInt<358979324>,SInt<314159265>) NL2;
// typedef TYPELIST_2(SInt<0>,SInt<100000000>) DL2;
// typedef SFraction<SBigInt<true,NL2,DefaultBase>,SBigInt<true,DL2,DefaultBase> > TPi2;

typedef TYPELIST_3(SInt<589793238>,SInt<141592653>,SInt<3>) NL21;
typedef SDecimalFraction<SBigInt<true,NL21,DefaultDecimalBase>,2,DefaultDecimalBase> TPi2Dec;

/// Metafunctions template classes
/*!
Template classes under this namespace are dedicated to the calculation of
different mathematical functions at compile time. Return values can be static constants, 
if they can be represented as integers, or static functions returning floating point values.
*/
namespace MF {


// Works with SFraction of decimal bases (10^n) only
// TODO: change that
template<class Fraction, int_t NDigits, base_t DecBase=DefaultDecimalBase>
struct FractionToDecimal;

template<class Numer, class Denom, int_t NDigits, base_t DecBase>
struct FractionToDecimal<SFraction<Numer,Denom>,NDigits,DecBase> {
  typedef typename IPowBig<SInt<DecBase>,NDigits>::Result D;
  typedef typename Mult<Numer,D>::Result NewNumer;
  typedef typename Div<NewNumer,Denom>::DivResult AllDecimals;
  typedef SDecimalFraction<AllDecimals,NDigits,DecBase> Result;
};

template<bool S1, class N1, class Denom, int_t NDigits, base_t Base>
struct FractionToDecimal<SFraction<SBigInt<S1,N1,Base>,Denom>,NDigits,Base> {
  typedef typename Loki::TL::ShiftRight<N1,NDigits,SInt<0> >::Result NList;
  typedef SBigInt<S1,NList,Base> NewNumer;
//typedef typename NL::Print<NewNumer>::Result TT2;
  typedef typename Div<NewNumer,Denom>::DivResult AllDecimals;
  typedef SDecimalFraction<AllDecimals,NDigits,Base> Result;
};

template<int_t N, int_t NDigits, base_t DecBase>
struct FractionToDecimal<SInt<N>,NDigits,DecBase> {
  typedef SDecimalFraction<SInt<N>,0,DecBase> Result;
};

template<class BI, int_t ND, base_t Base, int_t NDigits, base_t DecBase>
struct FractionToDecimal<SDecimalFraction<BI,ND,Base>,NDigits,DecBase> {
  typedef BI AllDecimals;
  typedef SDecimalFraction<BI,ND,Base> Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int_t Count> 
struct FuncSeries
{
  typedef FuncSeries<X,FuncStep,Accum,Count-1> NextIter;
  typedef FuncStep<Count-1,X,typename NextIter::ResultAux> FStep;
  typedef typename FStep::Result Step;
  typedef typename FStep::ResultAux ResultAux;

  typedef typename Accum<Step,typename NextIter::Result>::Result Result;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum> 
struct FuncSeries<X,FuncStep,Accum,1>
{
  typedef FuncStep<0,X,Loki::NullType> FStep;
  typedef typename FStep::Result Result;
  typedef typename FStep::ResultAux ResultAux;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum> 
struct FuncSeries<X,FuncStep,Accum,0> {};  // Error

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Accuracy, int I, 
class Value, class Dec1, class Dec2, class Aux,
bool C = (NL::Compare<Dec1,Dec2>::value == 0)>
class FuncAccuracyLoop;

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Accuracy, int I, class Value, class Dec1, class Dec2, class Aux>
struct FuncAccuracyLoop<X,FuncStep,Accum,Accuracy,I,Value,Dec1,Dec2,Aux,true>
{
  typedef Dec2 NextDecimal;
  typedef Value Result;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Accuracy, int I, class Value, class Dec1, class Dec2, class Aux>
struct FuncAccuracyLoop<X,FuncStep,Accum,Accuracy,I,Value,Dec1,Dec2,Aux,false>
{
  typedef FuncStep<I,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accum<Value,NextStep>::Result NextValue;
  
  typedef typename FractionToDecimal<NextValue,Accuracy,DefaultDecimalBase>::AllDecimals NextDecimal;
  typedef typename FuncAccuracyLoop<X,FuncStep,Accum,Accuracy,I+1,NextValue,Dec2,NextDecimal,NextAux>::Result Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Len, int I, class Value1, class Value2, class Aux,
bool C = (NL::Length<typename Value2::Numer>::value > Len 
       || NL::Length<typename Value2::Denom>::value > Len)>
class FuncLengthLoop;

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Len, int I, class Value1, class Value2, class Aux>
struct FuncLengthLoop<X,FuncStep,Accum,Len,I,Value1,Value2,Aux,true>
{
  typedef Value1 Result;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Len, int I, class Value1, class Value2, class Aux>
struct FuncLengthLoop<X,FuncStep,Accum,Len,I,Value1,Value2,Aux,false>
{
  typedef FuncStep<I,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accum<Value2,NextStep>::Result NextValue;
  typedef typename FuncLengthLoop<X,FuncStep,Accum,Len,I+1,Value2,NextValue,NextAux>::Result Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,  // One step of the series to compute the function
template<class,class> class Accumulator,      // How the steps are accumulated, normally Add or Mult
int Accuracy,                 // in powers of DefaultBase
int NStartingSteps>
struct GenericAccuracyBasedFunc;

template<class N, class D,
template<int,class,class> class FuncStep,  // One step of the series to compute the function
template<class,class> class Accumulator,      // How the steps are accumulated, normally Add or Mult
int Accuracy,                 // in powers of DefaultBase
int NStartingSteps>
struct GenericAccuracyBasedFunc<SFraction<N,D>,FuncStep,Accumulator,Accuracy,NStartingSteps>
{
  typedef SFraction<N,D> X;
  typedef FuncSeries<X,FuncStep,Accumulator,NStartingSteps> Sum;
  typedef typename Sum::Result StartValue;
  typedef typename Sum::ResultAux Aux;
  typedef typename FractionToDecimal<StartValue,Accuracy,DefaultBase>::AllDecimals StartDecimal;

  typedef FuncStep<NStartingSteps,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accumulator<NextStep,StartValue>::Result NextValue;
  typedef typename FractionToDecimal<NextValue,Accuracy,DefaultBase>::AllDecimals NextDecimal;
  
  typedef FuncAccuracyLoop<X,FuncStep,Accumulator,Accuracy,NStartingSteps+1,
                           NextValue,StartDecimal,NextDecimal,NextAux> Loop;
  //typedef SDecimalFraction<typename Loop::NextDecimal,Accuracy,DefaultDecimalBase> ResultDecimal;
  typedef typename Loop::Result Result;
};

template<class BI, int_t ND, base_t Base,
template<int,class,class> class FuncStep,  // One step of the series to compute the function
template<class,class> class Accumulator,      // How the steps are accumulated, normally Add or Mult
int Accuracy,                 // in powers of DefaultBase
int NStartingSteps>
struct GenericAccuracyBasedFunc<SDecimalFraction<BI,ND,Base>,FuncStep,Accumulator,Accuracy,NStartingSteps>
{
  typedef SDecimalFraction<BI,ND,Base> X;
  typedef FuncSeries<X,FuncStep,Accumulator,NStartingSteps> Sum;
  typedef typename Sum::Result StartValue;
  typedef typename Sum::ResultAux Aux;
  //typedef typename FractionToDecimal<StartValue,Accuracy,DefaultBase>::AllDecimals StartDecimal;

  typedef FuncStep<NStartingSteps,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accumulator<NextStep,StartValue>::Result NextValue;
  //typedef typename FractionToDecimal<NextValue,Accuracy,DefaultBase>::AllDecimals NextDecimal;
  
  typedef FuncAccuracyLoop<X,FuncStep,Accumulator,Accuracy,NStartingSteps+1,
                           NextValue,typename StartValue::Num,typename NextValue::Num,NextAux> Loop;
  //typedef SDecimalFraction<typename Loop::NextDecimal,Accuracy,DefaultDecimalBase> ResultDecimal;
  typedef typename Loop::Result Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,  // One step of the series to compute the function
template<class,class> class Accumulator,      // How the steps are accumulated, normally Add or Mult
int Length,    // in digits of DefaultBase
int NStartingSteps>
struct GenericLengthBasedFunc
{
  typedef FuncSeries<X,FuncStep,Accumulator,NStartingSteps> Sum;
  typedef typename Sum::Result StartValue;
  typedef typename Sum::ResultAux Aux;

  typedef FuncStep<NStartingSteps,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accumulator<NextStep,StartValue>::Result NextValue;

  typedef typename FuncLengthLoop<X,FuncStep,Accumulator,Length,NStartingSteps+1,StartValue,NextValue,NextAux>::Result Result;
};

/////////////////////////////////////////////////////

template<class T, int Accuracy, base_t Base = DefaultBase>
struct Reduce;

template<class N, class D, int Accuracy, base_t Base>
struct Reduce<SFraction<N,D>,Accuracy,Base> {
  typedef typename IPowBig<SInt<Base>,Accuracy>::Result Denom;
  typedef typename FractionToDecimal<SFraction<N,D>,Accuracy,Base>::AllDecimals Decimals;
  typedef typename Simplify<SFraction<Decimals,Denom> >::Result Result;
};

template<int_t N, int Accuracy, base_t Base>
struct Reduce<SInt<N>,Accuracy,Base> {
  typedef SInt<N> Result;
};

template<class BI, int_t ND, int Accuracy, base_t Base>
struct Reduce<SDecimalFraction<BI,ND,Base>,Accuracy,Base> {
  typedef typename BI::Num NList;
  typedef typename Loki::Select<(ND>Accuracy),
          typename Loki::TL::ShiftLeft<NList,ND-Accuracy>::Result,NList>::Result NewList;
  typedef SDecimalFraction<SBigInt<BI::isPositive,NewList,BI::Base>,Accuracy,Base> Result;
};


template<class SFrac, int Accuracy, class RetType = long double>
struct Compute;

template<class Numer, class Denom, int Accuracy, class RetType>
struct Compute<SFraction<Numer,Denom>,Accuracy,RetType> {
  typedef SFraction<Numer,Denom> Value;
  typedef typename FractionToDecimal<Value,Accuracy,DefaultDecimalBase>::Result TDec;
  typedef typename DoubleBase<typename TDec::Num>::Result BigInt;
  
  static RetType value() {
    return EvaluateToFloat<BigInt,RetType>::value()
         / DPow<DefaultDecimalBase,Accuracy,RetType>::value();
  }
};

template<int_t N, int Accuracy, class RetType>
struct Compute<SInt<N>,Accuracy,RetType> {
  typedef SInt<N> BigInt;
  static RetType value() { return static_cast<RetType>(N); }
};

template<class BI, int_t ND, base_t Base, int Accuracy, class RetType>
struct Compute<SDecimalFraction<BI,ND,Base>,Accuracy,RetType> {
  typedef SDecimalFraction<BI,ND,Base> Value;
  typedef typename Reduce<Value,Accuracy,Base>::Result TDec;
  typedef typename DoubleBase<typename TDec::Num>::Result BigInt;
  
  static RetType value() {
    return EvaluateToFloat<BigInt,RetType>::value() * 1e-18;
      //   / DPow<Base,Accuracy,RetType>::value();
  }
};

template<int_t N, int_t ND, base_t Base, int Accuracy, class RetType>
struct Compute<SDecimalFraction<SInt<N>,ND,Base>,Accuracy,RetType> {
  typedef SDecimalFraction<SInt<N>,ND,Base> Value;
  
  static RetType value() {
    return static_cast<RetType>(N) * 1e-18;
      //   / DPow<Base,ND,RetType>::value();
  }
};

////////////////////////////////////////////////////////

template<int K,class C,class Aux>
struct PiFraction
{
  typedef typename IPowBig<SInt<16>,K>::Result PBig;
  typedef SInt<8*K+1> TK1;
  typedef SInt<4*K+2> TK2;
  typedef SInt<8*K+5> TK3;
  typedef SInt<8*K+6> TK4;
  typedef typename Mult<typename Mult<typename Mult<TK1,TK2>::Result, 
                        typename Mult<TK3,TK4>::Result>::Result,PBig>::Result Denom;
  typedef typename Add<SInt<188>, typename Mult<SInt<4*K>,SInt<120*K+151> >::Result>::Result Numer;
  
  typedef typename Simplify<SFraction<Numer,Denom> >::Result Result;
//  typedef SFraction<Numer,Denom> Result;
  typedef Loki::NullType ResultAux;
};

template<class C,class Aux>
struct PiFraction<0,C,Aux>
{
  typedef SFraction<SInt<47>, SInt<15> > Result;
  typedef Loki::NullType ResultAux;
};



template<int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct PiAcc : public GenericAccuracyBasedFunc<Loki::NullType,PiFraction,Add,Accuracy,NStartingSteps> 
{};

template<int NStartingSteps>  
struct PiAcc<1,NStartingSteps> {
  static const base_t Base = DefaultDecimalBase;
  typedef TYPELIST_1(SInt<314159265>) NL1;
  typedef TYPELIST_1(SInt<Base/10>) DL1;
  typedef SFraction<SBigInt<true,NL1,Base>,SBigInt<true,DL1,Base> > Result;
};

template<int NStartingSteps>  
struct PiAcc<2,NStartingSteps> {
  static const base_t Base = DefaultDecimalBase;
  typedef TYPELIST_2(SInt<358979323>,SInt<314159265>) NL2;
  typedef TYPELIST_2(SInt<0>,SInt<Base/10>) DL2;
  typedef SFraction<SBigInt<true,NL2,Base>,SBigInt<true,DL2,Base> > Result;
};

template<int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct PiLen : public GenericLengthBasedFunc<Loki::NullType,PiFraction,Add,Len,NStartingSteps> 
{};

////////////////////////////////////////////////////////

template<class X, class Step, class Aux = Loki::NullType>
struct SinCosAux 
{
  typedef Pair<typename Aux::first, Step> Result;
};

template<class X, class Step>
struct SinCosAux<X,Step,Loki::NullType> 
{
  typedef typename Mult<X,X>::Result XX;
  typedef Pair<XX,Step> Result;
};

// Aux is Pair, where the first type is X^2 and the second is the previous series member
template<int K, class X, class Aux, int_t D>   // D=1 (for cos);   D=2 (for sin)
struct SinCosFraction 
{
  static const int_t M = 2*(K-1)+D;
  typedef SFraction<SInt<1>,SInt<M*(M+1)> > Divider;
  typedef typename Mult<typename Aux::first,Divider>::Result XX;
  typedef typename Mult<XX,typename Aux::second>::Result XP;
  typedef typename Negate<XP>::Result Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
//  typedef typename NL::Print<Result>::Result TT2;
};

template<class X, class Aux, int_t D>   // D=1 (for cos);   D=2 (for sin)
struct SinCosFraction<0,X,Aux,D>
{
  typedef typename Aux::second Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
};


template<int K, class X, class Aux, int_t D>   // D=1 (for cos);   D=2 (for sin)
struct SinCosDecimal 
{
  static const int Accuracy = 2;
  static const int_t M = 2*(K-1)+D;
  typedef SFraction<SInt<1>,SInt<M*(M+1)> > Divider;
  typedef typename FractionToDecimal<Divider,Accuracy>::Result DividerDec;
  typedef typename Mult<typename Aux::first,DividerDec>::Result XX;
  typedef typename Mult<XX,typename Aux::second>::Result XP;
  typedef typename Reduce<XP,Accuracy>::Result XPR;
  typedef typename Negate<XPR>::Result Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
};

template<class X, class Aux, int_t D>   // D=1 (for cos);   D=2 (for sin)
struct SinCosDecimal<0,X,Aux,D>
{
  typedef typename Aux::second Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
};


template<int K, class X, class Aux>
struct CosFraction : public SinCosFraction<K,X,Aux,1> {};

template<int K, class X>
struct CosFraction<K,X,Loki::NullType>
: public SinCosFraction<K,X,typename SinCosAux<X,SInt<1> >::Result,1> {};


template<int K, class X, class Aux>
struct CosDecimal : public SinCosDecimal<K,X,Aux,1> {};

template<int K, class X>
struct CosDecimal<K,X,Loki::NullType>
: public SinCosDecimal<K,X,typename SinCosAux<X,SInt<1> >::Result,1> {};



template<int K, class X, class Aux>
struct SinFraction : public SinCosFraction<K,X,Aux,2> {};

template<int K, class X>
struct SinFraction<K,X,Loki::NullType> 
: public SinCosFraction<K,X,typename SinCosAux<X,X>::Result,2> {};


template<int K, class X, class Aux>
struct SinDecimal : public SinCosDecimal<K,X,Aux,2> {};

template<int K, class X>
struct SinDecimal<K,X,Loki::NullType> 
: public SinCosDecimal<K,X,typename SinCosAux<X,X>::Result,2> {};


template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct CosAcc : public GenericAccuracyBasedFunc<X,CosFraction,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct CosDecAcc : public GenericAccuracyBasedFunc<X,CosDecimal,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct CosLen : public GenericLengthBasedFunc<X,CosFraction,Add,Len,NStartingSteps> {};


template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinAcc : public GenericAccuracyBasedFunc<X,SinFraction,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinDecAcc : public GenericAccuracyBasedFunc<X,SinDecimal,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct SinLen : public GenericLengthBasedFunc<X,SinFraction,Add,Len,NStartingSteps> {};


template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __SinPiFrac {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename PiAcc<Accuracy,NStartingSteps>::Result TPi;
   typedef typename Mult<TPi,F>::Result X;
   typedef typename SinAcc<X,Accuracy>::Result Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiFrac<A,1,Accuracy,NStartingSteps> {
  typedef SInt<0> Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiFrac<A,2,Accuracy,NStartingSteps> {
  typedef typename Loki::Select<(A%4 == 1),SInt<1>,SInt<-1> >::Result Result;
};

template<int_t A, int NStartingSteps>  
struct __SinPiFrac<A,3,2,NStartingSteps> {
  static const int_t R = A%6;
  typedef TYPELIST_2(SInt<784438645>,SInt<866025403>) NList;
  typedef TYPELIST_3(SInt<0>,SInt<0>,SInt<1>) DList;
  typedef SBigInt<(R==1 || R==2),NList,DefaultDecimalBase> Numer;
  typedef SBigInt<true,DList,DefaultDecimalBase> Denom;
  typedef SFraction<Numer,Denom> Result;
};

template<int_t A, int NStartingSteps>  
struct __SinPiFrac<A,4,2,NStartingSteps> {
  static const int_t R = A%8;
  typedef TYPELIST_2(SInt<186547523>,SInt<707106781>) NList;
  typedef TYPELIST_3(SInt<0>,SInt<0>,SInt<1>) DList;
  typedef SBigInt<(R==1 || R==3),NList,DefaultDecimalBase> Numer;
  typedef SBigInt<true,DList,DefaultDecimalBase> Denom;
  typedef SFraction<Numer,Denom> Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiFrac<A,6,Accuracy,NStartingSteps> {
  static const int_t R = A%12;
  typedef SFraction<SInt<1>,SInt<2> >  V1;
  typedef SFraction<SInt<-1>,SInt<2> > V2;
  typedef typename Loki::Select<(R==1 || R==5),V1,V2>::Result Result;
};

template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinPiFrac {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __SinPiFrac<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};


template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __SinPiDecimal {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename FractionToDecimal<F,Accuracy>::Result FDec;
   //typedef typename PiAcc<Accuracy,NStartingSteps>::Result TPi;
   typedef TPi2Dec TPi;   // <<< TODO: make general accuracy
   typedef typename Mult<TPi,FDec>::Result X;
   typedef typename Reduce<X,Accuracy>::Result XR;
   typedef typename SinDecAcc<XR,Accuracy>::Result Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiDecimal<A,1,Accuracy,NStartingSteps> {
  typedef SInt<0> Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiDecimal<A,2,Accuracy,NStartingSteps> {
  typedef typename Loki::Select<(A%4 == 1),SInt<1>,SInt<-1> >::Result Result;
};

template<int_t A, int NStartingSteps>  
struct __SinPiDecimal<A,3,2,NStartingSteps> {
  static const int_t R = A%6;
  typedef TYPELIST_2(SInt<784438645>,SInt<866025403>) NList;
  typedef SBigInt<(R==1 || R==2),NList,DefaultDecimalBase> Numer;
  typedef SDecimalFraction<Numer,2,DefaultDecimalBase> Result;
};

template<int_t A, int NStartingSteps>  
struct __SinPiDecimal<A,4,2,NStartingSteps> {
  static const int_t R = A%8;
  typedef TYPELIST_2(SInt<186547523>,SInt<707106781>) NList;
  typedef SBigInt<(R==1 || R==3),NList,DefaultDecimalBase> Numer;
  typedef SDecimalFraction<Numer,2,DefaultDecimalBase> Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiDecimal<A,6,Accuracy,NStartingSteps> {
  static const int_t R = A%12;
  typedef Loki::Typelist<SInt<500000000>,Loki::NullType> NList1;
  typedef typename Loki::TL::ShiftRight<NList1,Accuracy-1,SInt<0> >::Result NList;
  typedef SBigInt<(R==1 || R==5),NList,DefaultDecimalBase> Numer;
  typedef SDecimalFraction<Numer,2,DefaultDecimalBase> Result;
};

template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinPiDecimal {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __SinPiDecimal<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};



template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __CosPiFrac {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename PiAcc<Accuracy,NStartingSteps>::Result TPi;
   typedef typename Mult<TPi,F>::Result X;
   typedef typename CosAcc<X,Accuracy>::Result Result;
};
  
template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __CosPiFrac<A,1,Accuracy,NStartingSteps> {
  typedef typename Loki::Select<(A%2==0),SInt<1>,SInt<-1> >::Result Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __CosPiFrac<A,2,Accuracy,NStartingSteps> {
  typedef SInt<0> Result;
};

template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct CosPiFrac {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __CosPiFrac<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};



template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __CosPiDecimal {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename FractionToDecimal<F,Accuracy>::Result FDec;
   //typedef typename PiAcc<Accuracy,NStartingSteps>::Result TPi;
   typedef TPi2Dec TPi;   // <<< TODO: make general accuracy
   typedef typename Mult<TPi,FDec>::Result X;
   typedef typename Reduce<X,Accuracy>::Result XR;
   typedef typename CosDecAcc<XR,Accuracy>::Result Result;
};
  
template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __CosPiDecimal<A,1,Accuracy,NStartingSteps> {
  typedef typename Loki::Select<(A%2==0),SInt<1>,SInt<-1> >::Result Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __CosPiDecimal<A,2,Accuracy,NStartingSteps> {
  typedef SInt<0> Result;
};

template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct CosPiDecimal {
   typedef SFraction<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __CosPiDecimal<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};
  
////////////////////////////////////////////////////////

template<class T>
struct Cout;

template<int_t N>
struct Cout<SInt<N> > 
{
  static void apply(std::ostream& os) { 
    os << N;
  }
};

template<bool S, class H, class T, base_t Base>
struct Cout<SBigInt<S,Loki::Typelist<H,T>,Base> > 
{
  static const int_t W = NDigits<Base-1,10>::value;
  typedef Cout<SBigInt<S,T,Base> > Next;
  
  static void apply(std::ostream& os) { 
    Next::apply(os);
    os.fill('0');
    os.width(W);
    os << std::right << H::value << " ";
  }
};

template<bool S, class H, base_t Base>
struct Cout<SBigInt<S,Loki::Typelist<H,Loki::NullType>,Base> > {
  static void apply(std::ostream& os) { 
//     os << H::Value << " ";
    if (!S)
      os << "-";
    os << H::value << " ";
  }
};

template<class N, class D>
struct Cout<SFraction<N,D> > 
{
  typedef Cout<N> CN;
  typedef Cout<D> CD;
  
  static void apply(std::ostream& os) { 
    CN::apply(os);
    os << " / ";
    CD::apply(os);
  }
};

//////////////////////////////////////////

template<bool S, class H, class T, base_t Base, int_t NDecPlaces, base_t DecBase>
struct Cout<SDecimalFraction<SBigInt<S,Loki::Typelist<H,T>,Base>,NDecPlaces,DecBase> >
{
  static const int_t W = NDigits<Base-1,10>::value;
  static const int_t DW = NDigits<DecBase-1,10>::value;
  static const int_t Len = NL::Length<SBigInt<S,Loki::Typelist<H,T>,Base> >::value;
  typedef Cout<SDecimalFraction<SBigInt<S,T,Base>,NDecPlaces,DecBase> > Next;
  
  static void apply(std::ostream& os, const int_t len = 0) { 
    Next::apply(os,len+W);
    os.fill('0');
    if (NDecPlaces < len+W && NDecPlaces > len) {
      int_t d = 1;
      for (int i = 0; i < NDecPlaces-len; ++i) d *= 10;
      os.width(W-NDecPlaces+len);
      os << std::right << H::value/d << "." << H::value%d;
    }
    else {
      os.width(W);
      os << std::right << H::value;
    }
    if (NDecPlaces == len)
      os << ".";
  }
};

template<bool S, class H, base_t Base, int_t NDecPlaces, base_t DecBase>
struct Cout<SDecimalFraction<SBigInt<S,Loki::Typelist<H,Loki::NullType>,Base>,NDecPlaces,DecBase> > 
{
  static const int_t HW = NDigits<H::value,10>::value;
  static void apply(std::ostream& os, const int_t len = 0) { 
    if (!S)
      os << "-";
    if (NDecPlaces >= len+HW) {
      os << "0.";
      os.fill('0');
      os.width(NDecPlaces-len);
      os << std::right << H::value;
    }
    else if (NDecPlaces < len+HW && NDecPlaces > len) {
      int_t d = 1;
      for (int i = 0; i < NDecPlaces-len; ++i) d *= 10;
      os << H::value/d << "." << H::value%d;
    }
    else
      os << H::value;
    if (NDecPlaces == len)
      os << ".";
  }
};

template<int_t N, int_t NDecPlaces, base_t DecBase>
struct Cout<SDecimalFraction<SInt<N>,NDecPlaces,DecBase> > 
{
  static const bool S = (N>=0);
  static const int_t AN = S ? N : -N;
  static const int_t HW = NDigits<AN,10>::value;
  static void apply(std::ostream& os, const int_t len = 0) { 
    if (!S)
      os << "-";
    if (NDecPlaces >= len+HW) {
      os << "0.";
      os.fill('0');
      os.width(NDecPlaces-len);
      os << std::right << AN;
    }
    else if (NDecPlaces < len+HW && NDecPlaces > len) {
      int_t d = 1;
      for (int i = 0; i < NDecPlaces-len; ++i) d *= 10;
      os << AN/d << "." << AN%d;
    }
    else
      os << AN;
  }
};

} // namespace MF

#endif /*__metafunc_h*/

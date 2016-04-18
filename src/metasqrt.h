/***************************************************************************
 *   Copyright (C) 2014-2015 by Vladimir Mirnyy                            *
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

#ifndef __metasqrt_h
#define __metasqrt_h

/** \file
    \brief Compile-time computing of square root
*/

#include <cmath>

#include "metafunc.h"

namespace MF {

template<class X, int Accuracy>
class SqrtInitGuessDec;

template<long_t N, int Accuracy>
class SqrtInitGuessDec<long_<N>,Accuracy> {
  static const long_t ND2 = NDigits<N/100,10>::value/2;
  static const long_t NP = IPow<10,ND2>::value;
  static const long_t X0 = (N<10) ? 2*NP : 6*NP;
  //static const long_t D = N - X0*X0;
  static const base_t DecBase = DefaultDecimalBase;
  typedef typename CreateBigInt<long_<X0>,DecBase>::Result BI1;
  typedef typename Loki::TL::ShiftRight<typename BI1::Num,Accuracy,long_<0> >::Result NList;
  typedef SBigInt<BI1::isPositive,NList,DecBase> NewBI;
public:
  typedef SDecimal<NewBI,Accuracy,DecBase> Result;
};

// TODO: extend for any X (now works with long_)
template<int K, class X, class Aux, int Accuracy>
struct SqrtDecimal {
  typedef Aux XN;
  typedef typename DoubleAccuracy<XN>::Result XND;
  typedef typename Add<typename Mult<XND,XND>::Result,X>::Result Numerator;
  typedef typename Add<XN,XN>::Result Denominator;
  typedef typename Div<typename Numerator::Num,typename Denominator::Num>::DivResult Res;
  typedef SDecimal<Res,XN::NDec,XN::Base> ResultAux;
  typedef ResultAux Result;
};

template<int K, class X, int Accuracy>
struct SqrtDecimal<K,X,Loki::NullType,Accuracy> {
  typedef typename SqrtInitGuessDec<X,Accuracy>::Result ResultAux;
  typedef ResultAux Result;
};

struct SqrtDecimalFunc {
  template<int K, class X, class Aux, int Accuracy>
  struct Value : public SqrtDecimal<K,X,Aux,Accuracy> {};
};

template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct SqrtDecAcc : public GenericAccuracyBasedFunc<X,SqrtDecimalFunc,EmptyOperation,Accuracy,NStartingSteps> {};

//////////// Rational ///////////////////////

template<class X>
class SqrtInitGuess;

template<long_t N>
class SqrtInitGuess<long_<N> > {
  static const long_t ND2 = NDigits<N,10>::value/2;
  static const long_t NP = IPow<10,ND2>::value;
  static const long_t X0 = (N<10) ? 2*NP : 6*NP;
public:
  typedef SRational<long_<X0>,long_<1> > Result;
};

template<int K, class X, class Aux>
struct SqrtRational {
  typedef Aux XN;
  typedef typename DoubleAccuracy<XN>::Result XND;
  typedef typename Add<typename Mult<XND,XND>::Result,X>::Result Numerator;
  typedef typename Add<XN,XN>::Result Denominator;
  typedef typename Div<typename Numerator::Num,typename Denominator::Num>::DivResult Res;
  typedef SDecimal<Res,XN::NDec,XN::Base> Result;
  typedef Result ResultAux;
};

template<int K, class X>
struct SqrtRational<K,X,Loki::NullType> {
  typedef typename SqrtInitGuess<X>::Result Result;
  typedef Result ResultAux;
};

struct SqrtRationalFunc {
  template<int K, class X, class Aux>
  struct Value : public SqrtRational<K,X,Aux> {};
};

template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SqrtAcc : public GenericAccuracyBasedFunc<X,SqrtRationalFunc,EmptyOperation,Accuracy,NStartingSteps> {};



} // namespace MF

#endif /*__metasqrt_h*/

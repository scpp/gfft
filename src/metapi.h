/***************************************************************************
 *   Copyright (C) 2007-2015 by Vladimir Mirnyy                            *
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

#ifndef __metapi_h
#define __metapi_h

/** \file
    \brief Compile-time computing of Pi
*/

#include "metafunc.h"


//// To save compile time
// typedef TYPELIST_1(SInt<314159265>) NL1;
// typedef TYPELIST_1(SInt<100000000>) DL1;
// typedef SRational<SBigInt<true,NL1,DefaultBase>,SBigInt<true,DL1,DefaultBase> > TPi1;
// 
// typedef TYPELIST_2(SInt<358979324>,SInt<314159265>) NL2;
// typedef TYPELIST_2(SInt<0>,SInt<100000000>) DL2;
// typedef SRational<SBigInt<true,NL2,DefaultBase>,SBigInt<true,DL2,DefaultBase> > TPi2;

typedef TYPELIST_2(SInt<141592653>,SInt<3>) NL11;
typedef SDecimal<SBigInt<true,NL11,DefaultDecimalBase>,1,DefaultDecimalBase> TPi1Dec;

typedef TYPELIST_3(SInt<589793238>,SInt<141592653>,SInt<3>) NL21;
typedef SDecimal<SBigInt<true,NL21,DefaultDecimalBase>,2,DefaultDecimalBase> TPi2Dec;

typedef TYPELIST_4(SInt<462643383>,SInt<589793238>,SInt<141592653>,SInt<3>) NL31;
typedef SDecimal<SBigInt<true,NL31,DefaultDecimalBase>,3,DefaultDecimalBase> TPi3Dec;


namespace MF {


template<int K,class C,class Aux>
struct PiRational
{
  typedef typename IPowBig<SInt<16>,K>::Result PBig;
  typedef SInt<8*K+1> TK1;
  typedef SInt<4*K+2> TK2;
  typedef SInt<8*K+5> TK3;
  typedef SInt<8*K+6> TK4;
  typedef typename Mult<typename Mult<typename Mult<TK1,TK2>::Result, 
                        typename Mult<TK3,TK4>::Result>::Result,PBig>::Result Denom;
  typedef typename Add<SInt<188>, typename Mult<SInt<4*K>,SInt<120*K+151> >::Result>::Result Numer;
  
  typedef typename Simplify<SRational<Numer,Denom> >::Result Result;
//  typedef SRational<Numer,Denom> Result;
  typedef Loki::NullType ResultAux;
};

template<class C,class Aux>
struct PiRational<0,C,Aux>
{
  typedef SRational<SInt<47>, SInt<15> > Result;
  typedef Loki::NullType ResultAux;
};

struct PiRationalFunc 
{
  template<int K,class C,class Aux,int Accuracy>
  struct Value : public PiRational<K,C,Aux> {};
};


template<int K,class C,class Aux,int Accuracy>
struct PiDecimal
{
  typedef typename IPowBig<SInt<16>,K>::Result PBig;
  typedef SInt<8*K+1> TK1;
  typedef SInt<4*K+2> TK2;
  typedef SInt<8*K+5> TK3;
  typedef SInt<8*K+6> TK4;
  typedef typename Mult<typename Mult<typename Mult<TK1,TK2>::Result, 
                        typename Mult<TK3,TK4>::Result>::Result,PBig>::Result Denom;
  typedef typename Add<SInt<188>, typename Mult<SInt<4*K>,SInt<120*K+151> >::Result>::Result Numer;
  
  typedef SRational<Numer,Denom> F;
  typedef typename RationalToDecimal<F,Accuracy,DefaultDecimalBase>::Result Result;
  typedef Loki::NullType ResultAux;
};

template<class C,class Aux,int Accuracy>
struct PiDecimal<0,C,Aux,Accuracy>
{
  typedef SRational<SInt<47>, SInt<15> > F;
  typedef typename RationalToDecimal<F,Accuracy,DefaultDecimalBase>::Result Result;
  typedef Loki::NullType ResultAux;
};

struct PiDecimalFunc 
{
  template<int K,class C,class Aux,int Accuracy>
  struct Value : public PiDecimal<K,C,Aux,Accuracy> {};
};


template<int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct PiAcc : public GenericAccuracyBasedFunc<Loki::NullType,PiRationalFunc,Add,Accuracy,NStartingSteps> 
{};

template<int NStartingSteps>  
struct PiAcc<1,NStartingSteps> {
  static const base_t Base = DefaultDecimalBase;
  typedef TYPELIST_1(SInt<314159265>) NL1;
  typedef TYPELIST_1(SInt<Base/10>) DL1;
  typedef SRational<SBigInt<true,NL1,Base>,SBigInt<true,DL1,Base> > Result;
};

template<int NStartingSteps>  
struct PiAcc<2,NStartingSteps> {
  static const base_t Base = DefaultDecimalBase;
  typedef TYPELIST_2(SInt<358979323>,SInt<314159265>) NL2;
  typedef TYPELIST_2(SInt<0>,SInt<Base/10>) DL2;
  typedef SRational<SBigInt<true,NL2,Base>,SBigInt<true,DL2,Base> > Result;
};


template<int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct PiDecAcc : public GenericAccuracyBasedFunc<Loki::NullType,PiDecimalFunc,Add,Accuracy,NStartingSteps> {};


template<int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct PiLen : public GenericLengthBasedFunc<Loki::NullType,PiRationalFunc,Add,Len,NStartingSteps> 
{};

} // namespace MF

#endif /*__metapi_h*/

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

#ifndef __metasincos_h
#define __metasincos_h

/** \file
    \brief Compile-time computing of sine and cosine functions
*/

#include "metapi.h"

namespace MF {

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
struct SinCosRational 
{
  static const int_t M = 2*(K-1)+D;
  typedef SRational<SInt<1>,SInt<M*(M+1)> > Divider;
  typedef typename Mult<typename Aux::first,Divider>::Result XX;
  typedef typename Mult<XX,typename Aux::second>::Result XP;
  typedef typename Negate<XP>::Result Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
//  typedef typename NL::Print<Result>::Result TT2;
};

template<class X, class Aux, int_t D>   // D=1 (for cos);   D=2 (for sin)
struct SinCosRational<0,X,Aux,D>
{
  typedef typename Aux::second Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
};


template<int K, class X, class Aux, int_t D, // D=1 (for cos);   D=2 (for sin)
int Accuracy>   
struct SinCosDecimal 
{
  static const int_t M = 2*(K-1)+D;
  typedef SRational<SInt<1>,SInt<M*(M+1)> > Divider;
  typedef typename RationalToDecimal<Divider,Accuracy>::Result DividerDec;
  typedef typename Mult<typename Aux::first,DividerDec>::Result XX;
  typedef typename Mult<XX,typename Aux::second>::Result XP;
//  typedef typename Reduce<XP,Accuracy>::Result XPR;
  typedef typename Negate<XP>::Result Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
};

template<class X, class Aux, int_t D,   // D=1 (for cos);   D=2 (for sin)
int Accuracy>  
struct SinCosDecimal<0,X,Aux,D,Accuracy>
{
  typedef typename Aux::second Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
};


template<int K, class X, class Aux>
struct CosRational : public SinCosRational<K,X,Aux,1> {};

template<int K, class X>
struct CosRational<K,X,Loki::NullType>
: public SinCosRational<K,X,typename SinCosAux<X,SInt<1> >::Result,1> {};

struct CosRationalFunc {
  template<int K, class X, class Aux, int Accuracy>
  struct Value : public CosRational<K,X,Aux> {};
};


template<int K, class X, class Aux, int Accuracy>
struct CosDecimal : public SinCosDecimal<K,X,Aux,1,Accuracy> {};

template<int K, class X, int Accuracy>
struct CosDecimal<K,X,Loki::NullType,Accuracy>
: public SinCosDecimal<K,X,typename SinCosAux<X,SInt<1> >::Result,1,Accuracy> {};

struct CosDecimalFunc {
  template<int K, class X, class Aux, int Accuracy>
  struct Value : public CosDecimal<K,X,Aux,Accuracy> {};
};


template<int K, class X, class Aux>
struct SinRational : public SinCosRational<K,X,Aux,2> {};

template<int K, class X>
struct SinRational<K,X,Loki::NullType> 
: public SinCosRational<K,X,typename SinCosAux<X,X>::Result,2> {};

struct SinRationalFunc {
  template<int K, class X, class Aux, int Accuracy>
  struct Value : public SinRational<K,X,Aux> {};
};


template<int K, class X, class Aux, int Accuracy>
struct SinDecimal : public SinCosDecimal<K,X,Aux,2,Accuracy> {};

template<int K, class X, int Accuracy>
struct SinDecimal<K,X,Loki::NullType,Accuracy> 
: public SinCosDecimal<K,X,typename SinCosAux<X,X>::Result,2,Accuracy> {};


struct SinDecimalFunc {
  template<int K, class X, class Aux, int Accuracy>
  struct Value : public SinDecimal<K,X,Aux,Accuracy> {};
};


template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct CosAcc : public GenericAccuracyBasedFunc<X,CosRationalFunc,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct CosDecAcc : public GenericAccuracyBasedFunc<X,CosDecimalFunc,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct CosLen : public GenericLengthBasedFunc<X,CosRationalFunc,Add,Len,NStartingSteps> {};


template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinAcc : public GenericAccuracyBasedFunc<X,SinRationalFunc,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinDecAcc : public GenericAccuracyBasedFunc<X,SinDecimalFunc,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct SinLen : public GenericLengthBasedFunc<X,SinRationalFunc,Add,Len,NStartingSteps> {};


template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __SinPiFrac {
   typedef SRational<SInt<A>,SInt<B> > F;
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
  typedef SRational<Numer,Denom> Result;
};

template<int_t A, int NStartingSteps>  
struct __SinPiFrac<A,4,2,NStartingSteps> {
  static const int_t R = A%8;
  typedef TYPELIST_2(SInt<186547523>,SInt<707106781>) NList;
  typedef TYPELIST_3(SInt<0>,SInt<0>,SInt<1>) DList;
  typedef SBigInt<(R==1 || R==3),NList,DefaultDecimalBase> Numer;
  typedef SBigInt<true,DList,DefaultDecimalBase> Denom;
  typedef SRational<Numer,Denom> Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiFrac<A,6,Accuracy,NStartingSteps> {
  static const int_t R = A%12;
  typedef SRational<SInt<1>,SInt<2> >  V1;
  typedef SRational<SInt<-1>,SInt<2> > V2;
  typedef typename Loki::Select<(R==1 || R==5),V1,V2>::Result Result;
};

template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinPiFrac {
   typedef SRational<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __SinPiFrac<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};


template<int_t A, int_t B, 
int Accuracy,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __SinPiDecimal {
   typedef SRational<SInt<A>,SInt<B> > F;
   typedef typename RationalToDecimal<F,Accuracy>::Result FDec;
   typedef typename PiDecAcc<Accuracy,NStartingSteps>::Result TPi;
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
  typedef SDecimal<Numer,2,DefaultDecimalBase> Result;
};

template<int_t A, int NStartingSteps>  
struct __SinPiDecimal<A,4,2,NStartingSteps> {
  static const int_t R = A%8;
  typedef TYPELIST_2(SInt<186547523>,SInt<707106781>) NList;
  typedef SBigInt<(R==1 || R==3),NList,DefaultDecimalBase> Numer;
  typedef SDecimal<Numer,2,DefaultDecimalBase> Result;
};

template<int_t A, 
int Accuracy, int NStartingSteps>  
struct __SinPiDecimal<A,6,Accuracy,NStartingSteps> {
  static const int_t R = A%12;
  typedef Loki::Typelist<SInt<500000000>,Loki::NullType> NList1;
  typedef typename Loki::TL::ShiftRight<NList1,Accuracy-1,SInt<0> >::Result NList;
  typedef SBigInt<(R==1 || R==5),NList,DefaultDecimalBase> Numer;
  typedef SDecimal<Numer,2,DefaultDecimalBase> Result;
};

template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinPiDecimal {
   typedef SRational<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __SinPiDecimal<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};



template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __CosPiFrac {
   typedef SRational<SInt<A>,SInt<B> > F;
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
   typedef SRational<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __CosPiFrac<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};


template<int_t A, int_t B, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct __CosPiDecimal {
   typedef SRational<SInt<A>,SInt<B> > F;
   typedef typename RationalToDecimal<F,Accuracy>::Result FDec;
   typedef typename PiDecAcc<Accuracy,NStartingSteps>::Result TPi;
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
   typedef SRational<SInt<A>,SInt<B> > F;
   typedef typename Simplify<F>::Result SF;
   typedef typename __CosPiDecimal<SF::Numer::value,SF::Denom::value,Accuracy,NStartingSteps>::Result Result;
};


} // namespace MF

#endif /*__metasincos_h*/

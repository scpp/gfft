/***************************************************************************
 *   Copyright (C) 2014 by Vladimir Mirnyy                                 *
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

#ifndef __metaroot_h
#define __metaroot_h

/** \file
    \brief Compile-time computing of complex roots of unity
*/

#include <cmath>

#include "metasincos.h"
#include "metasqrt.h"

namespace MF {
  
// template<class W, class W1, int_t N, int_t I, int S>
// struct GetNextRoot;
// 
// template<class H, class T, class H1, class T1, int_t N, int_t I, int S>
// struct GetNextRoot<Loki::Typelist<H,T>,Loki::Typelist<H1,T1>,N,I,S>
// {
//    typedef typename Mult<H,H1>::Result W;
//    typedef typename GetNextRoot<T,T1,N,I,S>::Result Next;
//    typedef Loki::Typelist<W,Next> Result;
// };
// 
// template<int_t N, int_t I, int S>
// struct GetNextRoot<Loki::NullType,Loki::NullType,N,I,S>
// {
//    typedef Loki::NullType Result;
// };

template<int_t I, int_t M, class W1, class W, int Accuracy>
struct GetNextRoot {
  typedef typename Mult<W1,W>::Result Result;
};
/*
template<class W1, class W, int_t I, int Accuracy>
struct GetNextRoot<I,2,W1,W,Accuracy> {
  typedef typename CosPiFrac<I,2,Accuracy>::Result Re;
  typedef typename SinPiFrac<I,2,Accuracy>::Result Im;
  typedef typename RationalToDecimal<Re,Accuracy>::Result ReDec;
  typedef typename RationalToDecimal<Im,Accuracy>::Result ImDec;
  typedef MComplex<ReDec,ImDec> Result;
};

template<class W1, class W, int_t I, int Accuracy>
struct GetNextRoot<I,1,W1,W,Accuracy> {
  typedef typename CosPiFrac<I,1,Accuracy>::Result Re;
  typedef typename SinPiFrac<I,1,Accuracy>::Result Im;
  typedef typename RationalToDecimal<Re,Accuracy>::Result ReDec;
  typedef typename RationalToDecimal<Im,Accuracy>::Result ImDec;
  typedef MComplex<ReDec,ImDec> Result;
};
*/

template<class W1, int_t N, int Accuracy, class W, int_t Count, int_t I = 2>
struct __RootListLoop {
  typedef typename Simplify<SRational<SInt<2*I>,SInt<N> > >::Result SF;
//typedef typename NL::Print<SF>::Result TTT;
  typedef typename GetNextRoot<SF::Numer::value,SF::Denom::value,W1,W,Accuracy>::Result WW;
  typedef Compute<typename WW::Re,Accuracy> CRe;
  typedef Compute<typename WW::Im,Accuracy> CIm;
  typedef typename __RootListLoop<W1,N,Accuracy,WW,Count,I+1>::Result Next;
  typedef Loki::Typelist<Pair<CRe,CIm>,Next> Result;
};

template<class W1, int_t N, int Accuracy, class W, int_t Count>
struct __RootListLoop<W1,N,Accuracy,W,Count,Count> {
  typedef Loki::NullType Result;
};

template<class RList>
struct GenerateSymmetricPart;

template<class H, class T>
struct GenerateSymmetricPart<Loki::Typelist<H,T> > {
  typedef typename H::first T1;
  typedef typename Negate<typename H::second::BigInt>::Result T2;
  typedef typename GenerateSymmetricPart<T>::Result Next;
  typedef Loki::Typelist<Pair<T1,T2>,Next> Result;
};

template<>
struct GenerateSymmetricPart<Loki::NullType> {
  typedef Loki::NullType Result;
};


template<int_t N, int S, int Accuracy>
class GetFirstRoot {
  //typedef typename SinPiDecimal<1,N,Accuracy>::Result Sin1;
  typedef typename SinPiDecimal<2,N,Accuracy>::Result Sin2;
  
  typedef typename Loki::Select<(S<0),Sin2,
          typename Negate<Sin2>::Result>::Result WI;
  //typedef typename RationalToDecimal<WI,Accuracy>::Result WIDec;
  
//   typedef typename Sub<SInt<1>,typename Mult<SInt<2>,
//           typename Mult<Sin1,Sin1>::Result>::Result>::Result WR;
  typedef typename CosPiDecimal<2,N,Accuracy>::Result WR;
  //typedef typename RationalToDecimal<WR,Accuracy>::Result WRDec;
  
public:
  typedef MComplex<WR,WI> Result;
};



template<int_t N, int S, int Accuracy>
class GenerateRootList 
{
  typedef typename GetFirstRoot<N,-S,Accuracy>::Result W1;
  
public:
  typedef Compute<typename W1::Re,Accuracy> CRe;
  typedef Compute<typename W1::Im,Accuracy> CIm;
  typedef Loki::Typelist<Pair<CRe,CIm>,typename __RootListLoop<W1,N,Accuracy,W1,(N%2==0) ? N/2 : N/2+1>::Result> FirstHalf;
  typedef typename Loki::Select<(N%2==0),
          typename Loki::TL::Append<FirstHalf,typename GenerateRootList<2,S,Accuracy>::Result>::Result,FirstHalf>::Result FirstHalfMod;
//  typedef typename Loki::TL::Reverse<typename GenerateSymmetricPart<FirstHalf>::Result>::Result SecondHalf;

//public:
//  typedef typename Loki::TL::Append<FirstHalfMod,SecondHalf>::Result Result;
  typedef FirstHalf Result;
};

template<int S, int Accuracy>
class GenerateRootList<2,S,Accuracy> {
  typedef Compute<SInt<-1>,Accuracy> CRe;
  typedef Compute<SInt<0>, Accuracy> CIm;
public:
  typedef Loki::Typelist<Pair<CRe,CIm>,Loki::NullType> Result;
};

template<int S, int Accuracy>
class GenerateRootList<1,S,Accuracy> {
public:
  typedef Loki::NullType Result;
};

} // namespace

#endif /*__metaroot_h*/

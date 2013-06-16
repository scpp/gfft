/***************************************************************************
 *   Copyright (C) 2006-2013 by Volodymyr Myrnyy                           *
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

#ifndef __gfftfactor_h
#define __gfftfactor_h

/** \file
    \brief Factorization meta-algorithms 
*/

#include "gfftspec.h"


namespace GFFT {


template<class TList> struct Print;

template<> struct Print<Loki::NullType> { };

template<class Head, class Tail>
struct Print<Loki::Typelist<Head,Tail> > {
   typedef typename Print<Tail>::Result Result;
};


template<typename T1, typename T2>
struct Pair {
  typedef T1 first;
  typedef T2 second;
};

template<int_t N, int_t Factor, 
bool C = (N % Factor == 0)>
struct IsMultipleOf;
  
template<int_t N, int_t Factor>
struct IsMultipleOf<N, Factor, true> {
  static const int_t value = IsMultipleOf<N/Factor, Factor>::value + 1;
};

template<int_t N, int_t Factor> 
struct IsMultipleOf<N, Factor, false> {
  static const int_t value = 0;
};

template<int_t Factor> 
struct IsMultipleOf<0, Factor, true> {
  static const int_t value = 0;
};

template<int_t N, int_t K,  
template<int_t> class IntHolder = SInt,
bool C1 = ((6*K+1)*(6*K+1) <= N),
bool C2 = ((N % (6*K+1) == 0) || (N % (6*K+5) == 0))>
struct FactorizationLoop;

template<int_t N, int_t K,
template<int_t> class IntHolder>
struct FactorizationLoop<N, K, IntHolder, true, true>
{
  static const int_t Candidate1 = 6*K + 1;
  static const int_t Candidate2 = 6*K + 5;
  static const int_t P1 = IsMultipleOf<N, Candidate1>::value;
  static const int_t P2 = IsMultipleOf<N, Candidate2>::value;
  static const int_t F1 = IPow<Candidate1, P1>::value;
  static const int_t F2 = IPow<Candidate2, P2>::value;
  typedef Pair<IntHolder<Candidate1>, IntHolder<P1> > T1;
  typedef Pair<IntHolder<Candidate2>, IntHolder<P2> > T2;
  
  static const int_t NextN = N/F1/F2;
  typedef typename FactorizationLoop<NextN, K+1>::Result NextIter;
  
  typedef typename Loki::Select<(P1>0) && (P2>0), 
    Loki::Typelist<T1, Loki::Typelist<T2, NextIter> >, 
    typename Loki::Select<(P1>0), Loki::Typelist<T1, NextIter>, 
    typename Loki::Select<(P2>0), Loki::Typelist<T2, NextIter>, NextIter>::Result>::Result>::Result Result;
};

template<int_t N, int_t K,
template<int_t> class IntHolder>
struct FactorizationLoop<N, K, IntHolder, true, false> : public FactorizationLoop<N, K+1> {};

template<int_t N, int_t K, 
template<int_t> class IntHolder, bool C>
struct FactorizationLoop<N, K, IntHolder, false, C>
{
  typedef Pair<IntHolder<N>, IntHolder<1> > T;
  typedef Loki::Typelist<T, Loki::NullType> Result;
};


typedef TYPELIST_5(SInt<2>, SInt<3>, SInt<5>, SInt<7>, SInt<11>) InitialPrimesList;

template<typename Num,
template<int_t> class IntHolder = SInt,
typename StartList = InitialPrimesList>
struct Factorization;

template<typename Num, template<int_t> class IntHolder, typename H, typename Tail>
struct Factorization<Num, IntHolder, Loki::Typelist<H,Tail> >
{
  static const int_t P = IsMultipleOf<Num::value, H::value>::value;
  typedef IntHolder<Num::value / IPow<H::value,P>::value> NextNum;
  typedef typename Factorization<NextNum,IntHolder,Tail>::Result Next;
  typedef typename Loki::Select<(P > 0), 
     Loki::Typelist<Pair<IntHolder<H::value>, IntHolder<P> >, Next>, Next>::Result Result;
};

template<typename Num, template<int_t> class IntHolder>
struct Factorization<Num, IntHolder, Loki::NullType> : public FactorizationLoop<Num::value, 2> {};

template<template<int_t> class IntHolder>
struct Factorization<IntHolder<1>, IntHolder, Loki::NullType> {
  typedef Loki::NullType Result;
};

}  //namespace DFT

#endif /*__gfftfactor_h*/
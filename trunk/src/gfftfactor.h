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

#ifndef __gfftfactor_h
#define __gfftfactor_h

/** \file
    \brief Factorization meta-algorithms 
*/

#include "sint.h"

namespace GFFT {


template<class TList> struct Print;

template<> struct Print<Loki::NullType> { };

template<class Head, class Tail>
struct Print<Loki::Typelist<Head,Tail> > {
   typedef typename Print<Tail>::Result Result;
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
int_t AddPower = 0,
bool C1 = ((6*K+1)*(6*K+1) <= N),
bool C2 = ((N % (6*K+1) == 0) || (N % (6*K+5) == 0))>
struct FactorizationLoop;

template<int_t N, int_t K,
template<int_t> class IntHolder, int_t AddPower>
struct FactorizationLoop<N, K, IntHolder, AddPower, true, true>
{
  static const int_t Candidate1 = 6*K + 1;
  static const int_t Candidate2 = 6*K + 5;
  static const int_t P1 = IsMultipleOf<N, Candidate1>::value + AddPower;
  static const int_t P2 = IsMultipleOf<N, Candidate2>::value + AddPower;
  static const int_t F1 = IPow<Candidate1, P1>::value;
  static const int_t F2 = IPow<Candidate2, P2>::value;
  typedef Pair<IntHolder<Candidate1>, IntHolder<P1> > T1;
  typedef Pair<IntHolder<Candidate2>, IntHolder<P2> > T2;
  
  static const int_t NextN = N/F1/F2;
  typedef typename FactorizationLoop<NextN, K+1, IntHolder, AddPower>::Result NextIter;
  
  typedef typename Loki::Select<(P1>0) && (P2>0), 
    Loki::Typelist<T1, Loki::Typelist<T2, NextIter> >, 
    typename Loki::Select<(P1>0), Loki::Typelist<T1, NextIter>, 
    typename Loki::Select<(P2>0), Loki::Typelist<T2, NextIter>, NextIter>::Result>::Result>::Result Result;
};

template<int_t N, int_t K,
template<int_t> class IntHolder, int_t AddPower>
struct FactorizationLoop<N, K, IntHolder, AddPower, true, false> 
: public FactorizationLoop<N, K+1, IntHolder, AddPower> {};

template<int_t N, int_t K, 
template<int_t> class IntHolder, int_t AddPower, bool C>
struct FactorizationLoop<N, K, IntHolder, AddPower, false, C>
{
  typedef Pair<IntHolder<N>, IntHolder<1+AddPower> > T;
  typedef Loki::Typelist<T, Loki::NullType> Result;
};

template<int_t K, 
template<int_t> class IntHolder, int_t AddPower, bool C>
struct FactorizationLoop<1, K, IntHolder, AddPower, false, C>
{
  typedef Loki::NullType Result;
};


typedef TYPELIST_5(SInt<2>, SInt<3>, SInt<5>, SInt<7>, SInt<11>) InitialPrimesList;

template<typename Num,
template<int_t> class IntHolder = SInt,
typename StartList = InitialPrimesList, 
int_t AddPower = 0>
struct Factorization;

// Factorization using trial deletion from InitialPrimesList
template<int_t N, template<int_t> class IntHolder, typename H, typename Tail, int_t AddPower>
struct Factorization<SIntID<N>, IntHolder, Loki::Typelist<H,Tail>, AddPower>
{
  //static const int_t N = Num::value;
  static const int_t P = IsMultipleOf<N, H::value>::value;
  typedef SIntID<N / IPow<H::value,P>::value> NextNum;
  typedef typename Factorization<NextNum,IntHolder,Tail,AddPower>::Result Next;
  typedef typename Loki::Select<(P > 0), 
     Loki::Typelist<Pair<IntHolder<H::value>, IntHolder<P+AddPower> >, Next>, Next>::Result Result;
};

// Further factorization 
template<int_t N, template<int_t> class IntHolder, int_t AddPower>
struct Factorization<SIntID<N>, IntHolder, Loki::NullType, AddPower> 
: public FactorizationLoop<N, 2, IntHolder, AddPower> {};

// End of factorization
template<template<int_t> class IntHolder, typename H, typename Tail, int_t AddPower>
struct Factorization<SIntID<1>, IntHolder, Loki::Typelist<H,Tail>, AddPower> {
  typedef Loki::NullType Result;
};


template<int_t M, int_t P>
struct PowerHolder;

// The power is predefined, no factorization needed
template<int_t M, int_t P,
template<int_t> class IntHolder, typename StartList, int_t AddPower>
struct Factorization<PowerHolder<M,P>, IntHolder, StartList, AddPower> 
: public Factorization<SIntID<M>, IntHolder, StartList, P-1> { };



template<int_t N, typename FactorList, int_t Accum = 1, bool C = (N > Accum)>
struct ExtractFactor;

template<int_t N, typename NT, int_t P, typename Tail, int_t Accum>
struct ExtractFactor<N, Loki::Typelist<Pair<NT,SInt<P> >,Tail>, Accum, true>
{
  typedef ExtractFactor<N, Tail, Accum> T1;
  typedef ExtractFactor<N, Loki::Typelist<Pair<NT,SInt<P-1> >,Tail>, Accum*NT::value> Next;
  typedef typename Next::Result Result;
  static const int_t value = Next::value;
};

template<int_t N, typename NT, typename Tail, int_t Accum>
struct ExtractFactor<N, Loki::Typelist<Pair<NT,SInt<0> >,Tail>, Accum, true>
{
  typedef ExtractFactor<N, Tail, Accum> Next;
  typedef typename Next::Result Result;
  static const int_t value = Next::value;
};

template<int_t N, int_t Accum>
struct ExtractFactor<N, Loki::NullType, Accum, true>
{
  typedef ExtractFactor<N, Loki::NullType, Accum, false> Next;
  typedef typename Next::Result Result;
  static const int_t value = Next::value;
};
  
template<int_t N, typename NT, int_t P, typename Tail, int_t Accum>
struct ExtractFactor<N, Loki::Typelist<Pair<NT,SInt<P> >,Tail>, Accum, false>
{
  typedef Loki::Typelist<Pair<NT,SInt<P> >,Tail> Result;
  static const int_t value = Accum;
};

template<int_t N, typename NT, typename Tail, int_t Accum>
struct ExtractFactor<N, Loki::Typelist<Pair<NT,SInt<0> >,Tail>, Accum, false>
{
  typedef Tail Result;
  static const int_t value = Accum;
};

template<int_t N, int_t Accum>
struct ExtractFactor<N, Loki::NullType, Accum, false>
{
  typedef Loki::NullType Result;
  static const int_t value = Accum;
};

  
}  //namespace DFT

#endif /*__gfftfactor_h*/

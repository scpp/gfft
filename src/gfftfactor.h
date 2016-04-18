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

#include "metafactor.h"

namespace GFFT {

template<class TList> struct Print;

template<> struct Print<Loki::NullType> { };

template<class Head, class Tail>
struct Print<Loki::Typelist<Head,Tail> > {
   typedef typename Print<Tail>::Result Result;
};


template<long_t N, long_t Factor,
bool C = (N % Factor == 0)>
struct IsMultipleOf;
  
template<long_t N, long_t Factor>
struct IsMultipleOf<N, Factor, true> {
  static const long_t value = IsMultipleOf<N/Factor, Factor>::value + 1;
};

template<long_t N, long_t Factor>
struct IsMultipleOf<N, Factor, false> {
  static const long_t value = 0;
};

template<long_t Factor>
struct IsMultipleOf<0, Factor, true> {
  static const long_t value = 0;
};

/*
template<long_t N, long_t K,
template<long_t> class IntHolder = ulong_,
long_t AddPower = 0,
bool C1 = ((6*K+1)*(6*K+1) <= N),
bool C2 = ((N % (6*K+1) == 0) || (N % (6*K+5) == 0))>
struct FactorizationLoop;

template<long_t N, long_t K,
template<long_t> class IntHolder, long_t AddPower>
struct FactorizationLoop<N, K, IntHolder, AddPower, true, true>
{
  static const long_t Candidate1 = 6*K + 1;
  static const long_t Candidate2 = 6*K + 5;
  static const long_t P1 = IsMultipleOf<N, Candidate1>::value + AddPower;
  static const long_t P2 = IsMultipleOf<N, Candidate2>::value + AddPower;
  static const long_t F1 = IPow<Candidate1, P1>::value;
  static const long_t F2 = IPow<Candidate2, P2>::value;
  typedef pair_<IntHolder<Candidate1>, IntHolder<P1> > T1;
  typedef pair_<IntHolder<Candidate2>, IntHolder<P2> > T2;
  
  static const long_t NextN = N/F1/F2;
  typedef typename FactorizationLoop<NextN, K+1, IntHolder, AddPower>::Result NextIter;
  
  typedef typename Loki::Select<(P1>0) && (P2>0), 
    Loki::Typelist<T1, Loki::Typelist<T2, NextIter> >, 
    typename Loki::Select<(P1>0), Loki::Typelist<T1, NextIter>, 
    typename Loki::Select<(P2>0), Loki::Typelist<T2, NextIter>, NextIter>::Result>::Result>::Result Result;
};

template<long_t N, long_t K,
template<long_t> class IntHolder, long_t AddPower>
struct FactorizationLoop<N, K, IntHolder, AddPower, true, false> 
: public FactorizationLoop<N, K+1, IntHolder, AddPower> {};

template<long_t N, long_t K,
template<long_t> class IntHolder, long_t AddPower, bool C>
struct FactorizationLoop<N, K, IntHolder, AddPower, false, C>
{
  typedef pair_<IntHolder<N>, IntHolder<1+AddPower> > T;
  typedef Loki::Typelist<T, Loki::NullType> Result;
};

template<long_t K,
template<long_t> class IntHolder, long_t AddPower, bool C>
struct FactorizationLoop<1, K, IntHolder, AddPower, false, C>
{
  typedef Loki::NullType Result;
};


typedef TYPELIST_5(ulong_<2>, ulong_<3>, ulong_<5>, ulong_<7>, ulong_<11>) InitialPrimesList;

template<typename Num,
template<long_t> class IntHolder = ulong_,
typename StartList = InitialPrimesList, 
long_t AddPower = 0>
struct Factorization;

// Factorization using trial deletion from InitialPrimesList
template<long_t N, template<long_t> class IntHolder, typename H, typename Tail, long_t AddPower>
struct Factorization<ulong_ID<N>, IntHolder, Loki::Typelist<H,Tail>, AddPower>
{
  //static const long_t N = Num::value;
  static const long_t P = IsMultipleOf<N, H::value>::value;
  typedef ulong_ID<N / IPow<H::value,P>::value> NextNum;
  typedef typename Factorization<NextNum,IntHolder,Tail,AddPower>::Result Next;
  typedef typename Loki::Select<(P > 0), 
     Loki::Typelist<pair_<IntHolder<H::value>, IntHolder<P+AddPower> >, Next>, Next>::Result Result;
};

// Further factorization 
template<long_t N, template<long_t> class IntHolder, long_t AddPower>
struct Factorization<ulong_ID<N>, IntHolder, Loki::NullType, AddPower>
: public FactorizationLoop<N, 2, IntHolder, AddPower> {};

// End of factorization
template<template<long_t> class IntHolder, typename H, typename Tail, long_t AddPower>
struct Factorization<ulong_ID<1>, IntHolder, Loki::Typelist<H,Tail>, AddPower> {
  typedef Loki::NullType Result;
};
*/


template<typename TList, long_t AddP>
struct AddPower;

template<typename NT, ulong_t P, typename Tail, long_t AddP>
struct AddPower<Loki::Typelist<pair_<NT,ulong_<P> >,Tail>, AddP>
{
    typedef typename AddPower<Tail,AddP>::Result Next;
    typedef Loki::Typelist<pair_<NT, ulong_<P + AddP> >, Next> Result;
};

template<long_t AddP>
struct AddPower<Loki::NullType, AddP>
{
    typedef Loki::NullType Result;
};


template<long_t M, long_t P>
struct PowerHolder;


template<typename N>
struct Factorize : public Factorization<N> { };

// The power is predefined, no factorization needed
template<long_t M, long_t P>
struct Factorize<PowerHolder<M,P> >  
{ 
    typedef typename Factorization<ulong_<M> >::Result F;
    typedef typename AddPower<F, P-1>::Result Result;
};


template<long_t N, typename FactorList, long_t Accum = 1, bool C = (N > Accum)>
struct ExtractFactor;

template<long_t N, typename NT, long_t P, typename Tail, long_t Accum>
struct ExtractFactor<N, Loki::Typelist<pair_<NT,ulong_<P> >,Tail>, Accum, true>
{
  typedef ExtractFactor<N, Tail, Accum> T1;
  typedef ExtractFactor<N, Loki::Typelist<pair_<NT,ulong_<P-1> >,Tail>, Accum*NT::value> Next;
  typedef typename Next::Result Result;
  static const long_t value = Next::value;
};

template<long_t N, typename NT, typename Tail, long_t Accum>
struct ExtractFactor<N, Loki::Typelist<pair_<NT,ulong_<0> >,Tail>, Accum, true>
{
  typedef ExtractFactor<N, Tail, Accum> Next;
  typedef typename Next::Result Result;
  static const long_t value = Next::value;
};

template<long_t N, long_t Accum>
struct ExtractFactor<N, Loki::NullType, Accum, true>
{
  typedef ExtractFactor<N, Loki::NullType, Accum, false> Next;
  typedef typename Next::Result Result;
  static const long_t value = Next::value;
};
  
template<long_t N, typename NT, typename P, typename Tail, long_t Accum>
struct ExtractFactor<N, Loki::Typelist<pair_<NT,P>,Tail>, Accum, false>
{
  typedef Loki::Typelist<pair_<NT,P>,Tail> Result;
  static const long_t value = Accum;
};

template<long_t N, typename NT, typename Tail, long_t Accum>
struct ExtractFactor<N, Loki::Typelist<pair_<NT,ulong_<0> >,Tail>, Accum, false>
{
  typedef Tail Result;
  static const long_t value = Accum;
};

template<long_t N, long_t Accum>
struct ExtractFactor<N, Loki::NullType, Accum, false>
{
  typedef Loki::NullType Result;
  static const long_t value = Accum;
};

  
}  //namespace DFT

#endif /*__gfftfactor_h*/

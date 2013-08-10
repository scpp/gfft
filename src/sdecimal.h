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

#ifndef __sdecimal_h
#define __sdecimal_h

#include "sbigint.h"


template<class BigInt, int_t NDecPlaces, base_t DecBase>
struct SDecimalFraction {
   typedef BigInt Num;
   static const int_t NDec = NDecPlaces;
   static const base_t Base = DecBase;
};


template<class BI1, int_t ND1, class BI2, int_t ND2, base_t DecBase>
class Mult<SDecimalFraction<BI1,ND1,DecBase>,SDecimalFraction<BI2,ND2,DecBase> > {
  typedef typename Mult<BI1,BI2>::Result Prod;
public:
  typedef SDecimalFraction<Prod,ND1+ND2,DecBase> Result;
};

template<int_t N, class BI, int_t ND, base_t DecBase>
class Mult<SInt<N>,SDecimalFraction<BI,ND,DecBase> > {
  typedef typename Mult<SInt<N>,BI>::Result Prod;
public:
  typedef SDecimalFraction<Prod,ND,DecBase> Result;
};

template<int_t N, class BI, int_t ND, base_t DecBase>
class Mult<SDecimalFraction<BI,ND,DecBase>,SInt<N> > 
: public Mult<SInt<N>,SDecimalFraction<BI,ND,DecBase> > {};

/////////////////////////////////////////////////////////////

template<class BI1, class BI2, int_t ND, base_t DecBase>
class Add<SDecimalFraction<BI1,ND,DecBase>,SDecimalFraction<BI2,ND,DecBase> > {
  typedef typename Add<BI1,BI2>::Result Sum;
public:
  typedef SDecimalFraction<Sum,ND,DecBase> Result;
};

///////////////////////////////////////////////

template<class BI1, class BI2, int_t ND, base_t DecBase>
class Sub<SDecimalFraction<BI1,ND,DecBase>,SDecimalFraction<BI2,ND,DecBase> > {
  typedef typename Sub<BI1,BI2>::Result Dif;
public:
  typedef SDecimalFraction<Dif,ND,DecBase> Result;
};

///////////////////////////////////////////////

template<class BI, int_t ND, base_t DecBase>
struct Negate<SDecimalFraction<BI,ND,DecBase> > {
  typedef typename Negate<BI>::Result NewBI;
  typedef SDecimalFraction<NewBI,ND,DecBase> Result;
};

///////////////////////////////////////////////

template<class BI, int_t ND, base_t DecBase>
struct Check<SDecimalFraction<BI,ND,DecBase> > : public Check<BI> {};

#endif

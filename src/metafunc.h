/***************************************************************************
 *   Copyright (C) 2007-2009 by Volodymyr Myrnyy                           *
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

/// Metafunctions template classes
/*!
Template classes under this namespace are dedicated to the calculation of
different mathematical functions at compile time. Return values can be static constants, 
if they can be represented as integers, or static functions returning floating point values.
*/
namespace MF {

/// Common series to compile-time calculation of sine and cosine functions
/*!
\tparam M is the starting counter of members in the series (2 for Sin function and 1 for Cos function)
\tparam N is the number of last member in the series
\tparam A numerator
\tparam B denominator

Using theory of Taylor series sine and cosine functions can be defined as infinite series, 
which are valid for all real numbers \e x:
\f[
\sin(x) = \sum_{n=0}^{\infty}\frac{(-1)^n x^{2n+1}}{(2n+1)!} =
          x - \frac{x^3}{3!} + \frac{x^5}{5!} -  \frac{x^7}{7!} + ... = 
          x(1 - \frac{x^2}{2\cdot 3}(1 - \frac{x^2}{4\cdot 5}(1 - \frac{x^2}{6\cdot 7}(1 - ... )))) \approx
          x S(x,M,N), \quad M=2 ,
\f]
\f[
\cos(x) = \sum_{n=0}^{\infty}\frac{(-1)^n x^{2n}}{(2n)!} =
          1 - \frac{x^2}{2!} + \frac{x^4}{4!} -  \frac{x^6}{6!} + ... = 
          1 - \frac{x^2}{1\cdot 2}(1 - \frac{x^2}{3\cdot 4}(1 - \frac{x^2}{5\cdot 6}(1 - ... ))) \approx
          S(x,M,N), \quad M=1 .
\f]
Both series contain common series \e S :
\f[
S(x,M,N) = 1 - \frac{x^2}{M(M+1)}(1 - \frac{x^2}{(M+2)(M+3)}(1 - ... \frac{x^2}{N(N+1)}))
\f]
which can be parametrized by the starting denominator coefficient \e M and parameter \e N as the stopping criterium \e M = \e N.
This template class implements the common series \e S for the argument \f$ x = \frac{A\pi}{B} \f$.
*/
template<unsigned M, unsigned N, unsigned B, unsigned A>
struct SinCosSeries {
   static long double value() {
      return 1-(A*M_PI/B)*(A*M_PI/B)/M/(M+1)
               *SinCosSeries<M+2,N,B,A>::value();
   }
};

template<unsigned N, unsigned B, unsigned A>
struct SinCosSeries<N,N,B,A> {
   static long double value() { return 1.; }
};



/** \class {MF::Sin}
\brief Sine function

Compile-time calculation of \f$ \sin(\frac{A\pi}{B})\f$ function.
The function is computed as convergent series. Number of used series entries
is dependent on necessary accuracy. Therefore, this template class
is specialized for float, double and long double types.
\sa SinCosSeries
*/
template<unsigned B, unsigned A, typename T=double>
struct Sin;

template<unsigned B, unsigned A>
struct Sin<B,A,float> {
   static float value() {
      return (A*M_PI/B)*SinCosSeries<2,24,B,A>::value();
   }
};

template<unsigned B, unsigned A>
struct Sin<B,A,double> {
   static double value() {
      return (A*M_PI/B)*SinCosSeries<2,34,B,A>::value();
   }
};

template<unsigned B, unsigned A>
struct Sin<B,A,long double> {
   static long double value() {
      return (A*M_PI/B)*SinCosSeries<2,60,B,A>::value();
   }
};

/** \class {MF::Cos}
\brief Cosine function

Compile-time calculation of \f$ \cos(\frac{A\pi}{B})\f$ function.
The function is computed as convergent series. Number of used series entries
is dependent on necessary accuracy. Therefore, this template class
is specialized for float, double and long double types.
\sa SinCosSeries
*/
template<unsigned B, unsigned A, typename T=double>
struct Cos;

template<unsigned B, unsigned A>
struct Cos<B,A,float> {
   static float value() {
      return SinCosSeries<1,23,B,A>::value();
   }
};

template<unsigned B, unsigned A>
struct Cos<B,A,double> {
   static double value() {
      return SinCosSeries<1,33,B,A>::value();
   }
};

template<unsigned B, unsigned A>
struct Cos<B,A,long double> {
   static long double value() {
      return SinCosSeries<1,59,B,A>::value();
   }
};

// Returns number of digits in N in the Base-system (Base=2 for binary)
template<unsigned N, unsigned Base>
struct NDigits {
  static const unsigned value = NDigits<N/Base, Base>::value + 1;
};

template<unsigned Base>
struct NDigits<0, Base> {
  static const unsigned value = 0;
};

template<unsigned N, unsigned P>
struct IPow {
  static const unsigned long value = IPow<N,P-1>::value * N;
};

template<unsigned N>
struct IPow<N,1> {
  static const unsigned long value = N;
};

template<unsigned N>
struct IPow<N,0> {
  static const unsigned long value = 1;
};

template<int_t N, int_t P>
struct IPowBig {
  typedef typename Mult<typename IPowBig<N,P-1>::Result, SInt<N> >::Result Result;
};

template<int_t N>
struct IPowBig<N,1> {
  typedef SInt<N> Result;
};

template<int_t N>
struct IPowBig<N,0> {
  typedef SInt<1> Result;
};

template<unsigned N, unsigned I>
class SqrtSeries {
public:
   static long double value() {
      static const long double XI = SqrtSeries<N,I-1>::value();
      return 0.5*(XI + N/XI);
   }
};

template<unsigned N>
class SqrtSeries<N,0> {
  static const unsigned ND = NDigits<N, 2>::value;
  static const unsigned X0 = IPow<2, ND/2>::value;
public:
   static long double value() {
     return 0.5*(X0 + N/static_cast<long double>(X0));
   }
};

template<unsigned N, typename T=double>
struct Sqrt;

template<unsigned N>
struct Sqrt<N, float> {
   static float value() {
      return SqrtSeries<N,5>::value();
   }
};

template<unsigned N>
struct Sqrt<N, double> {
   static double value() {
      return SqrtSeries<N,6>::value();
   }
};

template<unsigned N>
struct Sqrt<N, long double> {
   static long double value() {
      return SqrtSeries<N,6>::value();
   }
};


template<int K = 10>
struct Pi
{
  static const unsigned long P16 = IPow<16,K>::value;
  static const int_t K1 = 8*K+1;
  static const int_t K2 = 4*K+2;
  static const int_t K3 = 8*K+5;
  static const int_t K4 = 8*K+6;
  
  typedef double T;
  static T value() 
  {
    return (4./static_cast<T>(K1) - 2./static_cast<T>(K2) - 1./static_cast<T>(K3) - 1./static_cast<T>(K4))
            /static_cast<T>(P16) + Pi<K-1>::value();
  }
};

template<>
struct Pi<0>
{
  static double value() { return 47./15.; }
};




namespace EX {
  
template<class Fraction, int_t NDigits, base_t DecBase>
struct FractionToDecimal;

template<class Numer, class Denom, int_t NDigits, base_t DecBase>
struct FractionToDecimal<SFraction<Numer,Denom>,NDigits,DecBase> {
  typedef typename IPowBig<DecBase,NDigits>::Result M;
  typedef typename Mult<Numer,M>::Result NewNumer;
  typedef typename Div<NewNumer,Denom>::DivResult AllDecimals;
  typedef SDecimalFraction<AllDecimals,NDigits,DecBase> Result;
};

/////////////////////////////////////////////////////

template<int K = 10>
struct Pi
{
  typedef typename IPowBig<16,K>::Result PBig;
  typedef SInt<8*K+1> TK1;
  typedef SInt<4*K+2> TK2;
  typedef SInt<8*K+5> TK3;
  typedef SInt<8*K+6> TK4;
  typedef typename Mult<typename Mult<typename Mult<TK1,TK2>::Result, 
                        typename Mult<TK3,TK4>::Result>::Result,PBig>::Result Denom;
  typedef typename Add<SInt<188>, typename Mult<SInt<4*K>,SInt<120*K+151> >::Result>::Result Numer;
  
//   typedef typename Simplify<SFraction<Numer,Denom> >::Result Fraction;
  typedef SFraction<Numer,Denom> Fraction;
  typedef typename Add<typename Pi<K-1>::Result, Fraction>::Result Result;
  //typedef typename Simplify<typename Add<typename Pi<K-1>::Result, Fraction>::Result>::Result Result;
};

template<>
struct Pi<0>
{
  typedef SFraction<SInt<47>, SInt<15> > Result;
};


} // namespcae EX

} // namespace MF

#endif /*__metafunc_h*/

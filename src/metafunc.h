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
   static double value() {
      return 1-(A*M_PI/B)*(A*M_PI/B)/M/(M+1)
               *SinCosSeries<M+2,N,B,A>::value();
   }
};

template<unsigned N, unsigned B, unsigned A>
struct SinCosSeries<N,N,B,A> {
   static double value() { return 1.; }
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
   static double value() {
      return SinCosSeries<1,60,B,A>::value();
   }
};

} // namespace MF

#endif /*__metafunc_h*/

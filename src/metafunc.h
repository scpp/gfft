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

//// To save compile time
typedef TYPELIST_1(SInt<314159265>) NL1;
typedef TYPELIST_1(SInt<100000000>) DL1;
typedef SFraction<SBigInt<true,NL1,DefaultBase>,SBigInt<true,DL1,DefaultBase> > TPi1;

typedef TYPELIST_2(SInt<358979323>,SInt<314159265>) NL2;
typedef TYPELIST_2(SInt<0>,SInt<100000000>) DL2;
typedef SFraction<SBigInt<true,NL2,DefaultBase>,SBigInt<true,DL2,DefaultBase> > TPi2;


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
  
// Works with SFraction of decimal bases (10^n) only
// TODO: change that
template<class Fraction, int_t NDigits, base_t DecBase>
struct FractionToDecimal;

template<class Numer, class Denom, int_t NDigits, base_t DecBase>
struct FractionToDecimal<SFraction<Numer,Denom>,NDigits,DecBase> {
  typedef typename IPowBig<DecBase,NDigits>::Result M;
  typedef typename Mult<Numer,M>::Result NewNumer;
  typedef typename Div<NewNumer,Denom>::DivResult AllDecimals;
  typedef SDecimalFraction<AllDecimals,NDigits,DecBase> Result;
};

template<bool S1, class N1, class Denom, int_t NDigits, base_t Base>
struct FractionToDecimal<SFraction<SBigInt<S1,N1,Base>,Denom>,NDigits,Base> {
  typedef typename Loki::TL::ShiftRight<N1,NDigits,SInt<0> >::Result NList;
  typedef SBigInt<S1,NList,Base> NewNumer;
//typedef typename NL::Print<NewNumer>::Result TT2;
  typedef typename Div<NewNumer,Denom>::DivResult AllDecimals;
  typedef SDecimalFraction<AllDecimals,NDigits,Base> Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int_t Count> 
struct FuncSeries
{
  typedef FuncSeries<X,FuncStep,Accum,Count-1> NextIter;
  typedef FuncStep<Count-1,X,typename NextIter::ResultAux> FStep;
  typedef typename FStep::Result Step;
  typedef typename FStep::ResultAux ResultAux;

//   typedef FuncStep<Start,X,Aux> FStep;
//   typedef typename FStep::Result Step;
//   typedef typename FStep::ResultAux StepAux;
  //typedef typename NextIter::ResultAux::second LastStep;

  //typedef typename FuncAux<X,LastStep,StepAux>::Result ResultAux;
  typedef typename Accum<Step,typename NextIter::Result>::Result Result;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum> 
struct FuncSeries<X,FuncStep,Accum,1>
{
  typedef FuncStep<0,X,Loki::NullType> FStep;
  typedef typename FStep::Result Result;
  typedef typename FStep::ResultAux ResultAux;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum> 
struct FuncSeries<X,FuncStep,Accum,0> {};  // Error

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Accuracy, int I, 
class Value, class Dec1, class Dec2, class Aux,
bool C = (NL::Compare<Dec1,Dec2>::value == 0)>
class FuncAccuracyLoop;

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Accuracy, int I, class Value, class Dec1, class Dec2, class Aux>
struct FuncAccuracyLoop<X,FuncStep,Accum,Accuracy,I,Value,Dec1,Dec2,Aux,true>
{
  typedef Dec2 NextDecimal;
  typedef Value Result;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Accuracy, int I, class Value, class Dec1, class Dec2, class Aux>
struct FuncAccuracyLoop<X,FuncStep,Accum,Accuracy,I,Value,Dec1,Dec2,Aux,false>
{
  typedef FuncStep<I,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accum<Value,NextStep>::Result NextValue;
  
  typedef typename FractionToDecimal<NextValue,Accuracy,DefaultDecimalBase>::AllDecimals NextDecimal;
  typedef typename FuncAccuracyLoop<X,FuncStep,Accum,Accuracy,I+1,NextValue,Dec2,NextDecimal,NextAux>::Result Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Len, int I, class Value1, class Value2, class Aux,
bool C = (NL::Length<typename Value2::Numer>::value > Len 
       || NL::Length<typename Value2::Denom>::value > Len)>
class FuncLengthLoop;

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Len, int I, class Value1, class Value2, class Aux>
struct FuncLengthLoop<X,FuncStep,Accum,Len,I,Value1,Value2,Aux,true>
{
  typedef Value1 Result;
};

template<class X,
template<int,class,class> class FuncStep,
template<class,class> class Accum,
int Len, int I, class Value1, class Value2, class Aux>
struct FuncLengthLoop<X,FuncStep,Accum,Len,I,Value1,Value2,Aux,false>
{
  typedef FuncStep<I,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accum<Value2,NextStep>::Result NextValue;
  typedef typename FuncLengthLoop<X,FuncStep,Accum,Len,I+1,Value2,NextValue,NextAux>::Result Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,  // One step of the series to compute the function
template<class,class> class Accumulator,      // How the steps are accumulated, normally Add or Mult
int Accuracy,                 // in powers of DefaultBase
int NStartingSteps>
struct GenericAccuracyBasedFunc
{
  typedef FuncSeries<X,FuncStep,Accumulator,NStartingSteps> Sum;
  typedef typename Sum::Result StartValue;
  typedef typename Sum::ResultAux Aux;
  typedef typename FractionToDecimal<StartValue,Accuracy,DefaultBase>::AllDecimals StartDecimal;

  typedef FuncStep<NStartingSteps,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
  typedef typename Accumulator<NextStep,StartValue>::Result NextValue;
  typedef typename FractionToDecimal<NextValue,Accuracy,DefaultBase>::AllDecimals NextDecimal;
  
  typedef FuncAccuracyLoop<X,FuncStep,Accumulator,Accuracy,NStartingSteps+1,
                           NextValue,StartDecimal,NextDecimal,NextAux> Loop;
  //typedef SDecimalFraction<typename Loop::NextDecimal,Accuracy,DefaultDecimalBase> ResultDecimal;
  typedef typename Loop::Result Result;
};

/////////////////////////////////////////////////////

template<class X,
template<int,class,class> class FuncStep,  // One step of the series to compute the function
template<class,class> class Accumulator,      // How the steps are accumulated, normally Add or Mult
int Length,    // in digits of DefaultBase
int NStartingSteps>
struct GenericLengthBasedFunc
{
  typedef FuncSeries<X,FuncStep,Accumulator,NStartingSteps> Sum;
//   typedef typename Simplify<typename Sum::Result>::Result StartValue;
  typedef typename Sum::Result StartValue;
  typedef typename Sum::ResultAux Aux;

  typedef FuncStep<NStartingSteps,X,Aux> FStep;
  typedef typename FStep::Result NextStep;
  typedef typename FStep::ResultAux NextAux;
//   typedef typename Simplify<typename Accumulator<NextStep,StartValue>::Result>::Result NextValue;
  typedef typename Accumulator<NextStep,StartValue>::Result NextValue;

  typedef typename FuncLengthLoop<X,FuncStep,Accumulator,Length,NStartingSteps+1,StartValue,NextValue,NextAux>::Result Result;
};

////////////////////////////////////////////////////////

template<int K,class C,class Aux>
struct PiFraction
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
  typedef SFraction<Numer,Denom> Result;
  typedef Loki::NullType ResultAux;
};

template<class C,class Aux>
struct PiFraction<0,C,Aux>
{
  typedef SFraction<SInt<47>, SInt<15> > Result;
  typedef Loki::NullType ResultAux;
};



template<int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 6>  
struct PiAcc : public GenericAccuracyBasedFunc<Loki::NullType,PiFraction,Add,Accuracy,NStartingSteps> 
{};

template<int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct PiLen : public GenericLengthBasedFunc<Loki::NullType,PiFraction,Add,Len,NStartingSteps> 
{};

////////////////////////////////////////////////////////

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
struct SinCosFraction 
{
  static const int_t M = 2*(K-1)+D;
  typedef SFraction<SInt<1>,SInt<M*(M+1)> > Divider;
//   typedef typename Mult<typename Mult<X,X>::Result,Divider>::Result XX;
//   typedef typename Mult<XX,Aux>::Result XP;
  typedef typename Mult<typename Aux::first,Divider>::Result XX;
  typedef typename Mult<XX,typename Aux::second>::Result XP;
  typedef typename Negate<XP>::Result Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
//  typedef typename NL::Print<Result>::Result TT2;
};

template<class X, class Aux, int_t D>   // D=1 (for cos);   D=2 (for sin)
struct SinCosFraction<0,X,Aux,D>
{
  typedef typename Aux::second Result;
  typedef typename SinCosAux<X,Result,Aux>::Result ResultAux;
};


template<int K, class X, class Aux>
struct CosFraction : public SinCosFraction<K,X,Aux,1> {};

template<int K, class X>
struct CosFraction<K,X,Loki::NullType>
: public SinCosFraction<K,X,typename SinCosAux<X,UnitFraction>::Result,1> {};


template<int K, class X, class Aux>
struct SinFraction : public SinCosFraction<K,X,Aux,2> {};

template<int K, class X>
struct SinFraction<K,X,Loki::NullType> 
: public SinCosFraction<K,X,typename SinCosAux<X,X>::Result,2> {};


template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct CosAcc : public GenericAccuracyBasedFunc<X,CosFraction,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct CosLen : public GenericLengthBasedFunc<X,CosFraction,Add,Len,NStartingSteps> {};


template<class X, 
int Accuracy = 2,    // in powers of DefaultBase
int NStartingSteps = 5>  
struct SinAcc : public GenericAccuracyBasedFunc<X,SinFraction,Add,Accuracy,NStartingSteps> {};

template<class X, 
int Len = 2,    // in powers of DefaultBase
int NStartingSteps = 3>  
struct SinLen : public GenericLengthBasedFunc<X,SinFraction,Add,Len,NStartingSteps> {};

  
} // namespcae EX

////////////////////////////////////////////////////////

template<class T>
struct Cout;

template<int_t N>
struct Cout<SInt<N> > 
{
  static void apply(std::ostream& os) { 
    os << N;
  }
};

template<bool S, class H, class T, base_t Base>
struct Cout<SBigInt<S,Loki::Typelist<H,T>,Base> > 
{
  static const int_t W = NDigits<Base-1,10>::value;
  typedef Cout<SBigInt<S,T,Base> > Next;
  
  static void apply(std::ostream& os) { 
    Next::apply(os);
    os.fill('0');
    os.width(W);
    os << std::right << H::value << " ";
  }
};

template<bool S, class H, base_t Base>
struct Cout<SBigInt<S,Loki::Typelist<H,Loki::NullType>,Base> > {
  static void apply(std::ostream& os) { 
//     os << H::Value << " ";
    if (!S)
      os << "-";
    os << H::value << " ";
  }
};

template<class N, class D>
struct Cout<SFraction<N,D> > 
{
  typedef Cout<N> CN;
  typedef Cout<D> CD;
  
  static void apply(std::ostream& os) { 
    CN::apply(os);
    os << " / ";
    CD::apply(os);
  }
};

//////////////////////////////////////////

template<bool S, class H, class T, base_t Base, int_t NDecPlaces, base_t DecBase>
struct Cout<SDecimalFraction<SBigInt<S,Loki::Typelist<H,T>,Base>,NDecPlaces,DecBase> >
{
  static const int_t W = NDigits<Base-1,10>::value;
  static const int_t DW = NDigits<DecBase-1,10>::value;
  static const int_t Len = NL::Length<SBigInt<S,Loki::Typelist<H,T>,Base> >::value;
  typedef Cout<SDecimalFraction<SBigInt<S,T,Base>,NDecPlaces,DecBase> > Next;
  
  static void apply(std::ostream& os, const int_t len = 0) { 
    Next::apply(os,len+W);
    os.fill('0');
    if (NDecPlaces < len+W && NDecPlaces > len) {
      int_t d = 1;
      for (int i = 0; i < NDecPlaces-len; ++i) d *= 10;
      os.width(W-NDecPlaces+len);
      os << std::right << H::value/d << "." << H::value%d;
    }
    else {
      os.width(W);
      os << std::right << H::value;
    }
    if (NDecPlaces == len)
      os << ".";
  }
};

template<bool S, class H, base_t Base, int_t NDecPlaces, base_t DecBase>
struct Cout<SDecimalFraction<SBigInt<S,Loki::Typelist<H,Loki::NullType>,Base>,NDecPlaces,DecBase> > 
{
  static const int_t HW = NDigits<H::value,10>::value;
  static void apply(std::ostream& os, const int_t len = 0) { 
    if (!S)
      os << "-";
    if (NDecPlaces >= len+HW) {
      os << "0.";
      os.fill('0');
      os.width(NDecPlaces-len);
      os << std::right << H::value;
    }
    else if (NDecPlaces < len+HW && NDecPlaces > len) {
      int_t d = 1;
      for (int i = 0; i < NDecPlaces-len; ++i) d *= 10;
      os << H::value/d << "." << H::value%d;
    }
    else
      os << H::value;
    if (NDecPlaces == len)
      os << ".";
  }
};

template<int_t N, int_t NDecPlaces, base_t DecBase>
struct Cout<SDecimalFraction<SInt<N>,NDecPlaces,DecBase> > 
{
  static const bool S = (N>=0);
  static const int_t AN = S ? N : -N;
  static const int_t HW = NDigits<AN,10>::value;
  static void apply(std::ostream& os, const int_t len = 0) { 
    if (!S)
      os << "-";
    if (NDecPlaces >= len+HW) {
      os << "0.";
      os.fill('0');
      os.width(NDecPlaces-len);
      os << std::right << AN;
    }
    else if (NDecPlaces < len+HW && NDecPlaces > len) {
      int_t d = 1;
      for (int i = 0; i < NDecPlaces-len; ++i) d *= 10;
      os << AN/d << "." << AN%d;
    }
    else
      os << AN;
  }
};

} // namespace MF

#endif /*__metafunc_h*/

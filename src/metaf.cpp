/***************************************************************************
 *   Copyright (C) 2009-2013 by Volodymyr Myrnyy                                *
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

/** \file
    \brief Sample program to test meta functions
*/

#include <iostream>

#include "gfft.h"

using namespace std;

using namespace GFFT;


template<class List, int_t N, int_t I=2>
struct PrintRootList;

template<class H, class T, int_t N, int_t I>
struct PrintRootList<Loki::Typelist<H,T>,N,I> {
  static void apply() {
    const long double a = H::first::value();
    const long double b = H::second::value();
    cout.precision(18);
    cout << "(" << a << ", " << b << ") --- ";
    cout.precision(2);
    cout << "(" << fabs(a-cos(I*M_PI/(long double)N)) << ", " << fabs(b-sin(I*M_PI/(long double)N)) << ")" << endl;
    PrintRootList<T,N,I+2>::apply();
  }
};

template<int_t N, int_t I>
struct PrintRootList<Loki::NullType,N,I> {
  static void apply() {}
};
//typedef SFraction<SInt<157>,SInt<100> > X;
//typedef SInt<1> X;
// typedef TPi2 X;
// typedef Mult<TPi2,SFraction<SInt<1>,SInt<4> > >::Result X4;


typedef DOUBLE ValueType;
//typedef typename GenNumList<2, 3>::Result NList;
//typedef TYPELIST_4(SIntID<2>, SIntID<3>, SIntID<4>, SIntID<5>) NList;
// typedef TYPELIST_1(SIntID<8>) NList;
// typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, ParallelizationGroup::Default, INFREQ> TransformSet;

int main(int argc, char *argv[])
{
    cout.precision(18);
//     unsigned int p = 2;
//     unsigned long i, n = (TransformType::ID == RDFT::ID) ? (1<<(p-1)) : (1<<p);
//    int_t i, n = 8;
//    cin >> n;
   
//    typedef Print<Factorization<SInt<2>, SIntID>::Result>::Result TTT;  // 2*3*18539
//    typedef Print<Factorization<SInt<111234> >::Result>::Result TTT;  // 2*3*18539
//    typedef Print<Factorization<SInt<1024*169*1999> >::Result>::Result TTT;
//    typedef Print<Factorization<SInt<1024> >::Result>::Result TTT;
//    typedef Print<FactorizationLoop<13,2>::Result>::Result TTT;

// typedef TYPELIST_3(SInt<1>,SInt<2>,SInt<3>) LL1;
// typedef TYPELIST_2(SInt<999999999>,SInt<999999999>) LL2;
// typedef TYPELIST_2(SInt<999999999>,SInt<999999999>) LL3;
// typedef SBigInt<true,LL1,10> BI1;
// typedef SBigInt<true,LL2,DefaultBase> BI2;
// typedef Check<BI2>::Result a;
// typedef SBigInt<true,LL3,DefaultBase> BI3;
// typedef Mult<BI2,BI3>::Result M;
// Cout<M>::apply(cout);

// typedef Div<BI2,SInt<3> >::DivResult DD;
// typedef Simplify<Div<BI2,SInt<30> >::ModResult>::Result DM;
// typedef Add<BI1,BI2>::Result A;
// typedef Sub<BI1,BI2>::Result S;
// typedef Sub<S,BI2>::Result SS;
//typedef Simplify<Div<M,BI2>::ModResult>::Result D;
//typedef Print<DM>::Result TTT;

//typedef unsigned long int RetType;
typedef long double RetType;
//typedef Simplify<SF>::Result SI;
//typedef Print<G1>::Result TTT;
//cout<< Evaluate2Int<G1,int_t>::Value << endl; 
// cout<< Evaluate2Int<D,int_t>::value << endl; 
// cout<< Evaluate2Float<M,double>::value() << endl; 
// cout<< Evaluate2Float<D,double>::value() << endl; 
// cout << (int)NL::Compare<BI1, BI2>::value << endl;
// cout << (int)NL::Compare<BI2, S>::value << endl;
// cout << (int)NL::Compare<BI2, SS>::value << endl;
//cout << (int)GCD<BI2, BI2>::Result::value << endl;
//typedef Pi<13>::Result TPi;
// typedef GCD<TPi::Numer, TPi::Denom>::Result G1;
// cout << ">>>>>>> " << Evaluate2Int<G1,RetType>::Value << endl;
// typedef Div<TPi::Numer, TPi::Denom>::ModResult M1;
// cout << "1) mod("<< Evaluate2Int<TPi::Numer,RetType>::Value << "," << Evaluate2Int<TPi::Denom,RetType>::Value << ") = " << Evaluate2Int<M1,RetType>::Value << endl;
// typedef Div<TPi::Denom, M1>::ModResult M2;
// cout << "2) mod("<< Evaluate2Int<TPi::Denom,RetType>::Value << "," << Evaluate2Int<M1,RetType>::Value << ") = " << Evaluate2Int<M2,RetType>::Value << endl;
//typedef Print<TPi::Numer>::Result TT2;
// //typedef Print<M2>::Result TT1;
// typedef Div<M1, M2>::ModResult M3;
// cout << "3) mod("<< Evaluate2Int<M1,RetType>::Value << "," << Evaluate2Int<M2,RetType>::Value << ") = " << Evaluate2Int<M3,RetType>::Value << endl;
// typedef Div<M2, M3>::ModResult M4;
// cout << "4) mod("<< Evaluate2Int<M2,RetType>::Value << "," << Evaluate2Int<M3,RetType>::Value << ") = " << Evaluate2Int<M4,RetType>::Value << endl;
// 
//const RetType numer = Evaluate2Float<TPi::Numer,RetType>::value();
//const RetType denom = Evaluate2Float<TPi::Denom,RetType>::value();
//cout<< numer << "/" << denom << " = " << (double)numer/denom << endl;

// //typedef MetaPi::Result TPi;
// typedef EX::FractionToDecimal<TPi,20,10>::Result TPiDec;
//typedef Print<TPi::Num>::Result TT2;
//cout << " ";
// Cout<Translate<TPiDec::Num,DefaultDecimalBase>::Result>::apply(cout);
//Cout<TPiDec::Num>::apply(cout);
//cout << endl << M_PI << endl;
//cout << Evaluate2Int<SI::Numer,RetType>::value << "/" << Evaluate2Int<SI::Denom,RetType>::value << endl;
// Cout<TPi1::Numer>::apply(cout);
// cout << " / " << endl;
// Cout<TPi1::Denom>::apply(cout);
// cout << endl;
// Cout<TPi::Numer>::apply(cout);
// cout << " / " << endl;
// Cout<TPi::Denom>::apply(cout);
// cout << endl;
// cout << Loki::TL::Length<TPi::Numer::Num>::value << " " << Loki::TL::Length<TPi::Denom::Num>::value << endl;


//    static const int NStart = 5;
//    typedef EX::SinCosAux<X,X,Loki::NullType>::Result Aux;
//    typedef EX::FuncSeries<X,EX::SinFraction,Add,NStart> Sum;
//    typedef typename Sum::Result StartValue;
// // //   typedef typename Simplify<typename Sum::Result>::Result StartValue;
// //    
//    typedef typename Sum::ResultAux Last;
//    typedef typename EX::SinFraction<NStart,X,Last>::Result Step;
// //   typedef typename Add<StartValue,Step>::Result NextValue;
// //    typedef typename Simplify<typename AddOld<StartValue,Step>::Result>::Result NextValue;
//    typedef typename Add<StartValue,Step>::Result NextValue;
// 
// //typedef EX::FractionToDecimal<Last::second,20,10>::Result LastDec;
// typedef EX::FractionToDecimal<Step,20,10>::Result StepDec;
// typedef EX::FractionToDecimal<StartValue,20,10>::Result StartDec;
// typedef EX::FractionToDecimal<NextValue,20,10>::Result NextDec;

//   typedef SFraction<SInt<-1>,SInt<720> > F3;
//   typedef SFraction<SInt<1>,SInt<40320> > F2;
//   typedef SFraction<SInt<-1>,SInt<3628800> > F1;
//   typedef typename Add<F2,F1>::Result F21;
//   typedef typename Add<F3, F21>::Result F;
//   typedef typename F3::Numer N1;
//   typedef typename F3::Denom D1;
//   typedef typename F21::Numer N2;
//   typedef typename F21::Denom D2;
//    typedef typename Mult<N1,D2>::Result T1;
//    typedef typename Mult<N2,D1>::Result T2;
//    typedef typename Add<T1,T2>::Result Num;
//    typedef typename Mult<D1,D2>::Result Den;
  
// typedef EX::PiAcc<2> MetaPi;
// typedef MetaPi::Result TPi;
// typedef EX::FractionToDecimal<TPi,20,10>::Result TPiDec;
//typedef EX::FractionToDecimal<TPi,2,DefaultDecimalBase>::Result TPiDec2;
// // typedef Translate<TPiDec::Num,DefaultDecimalBase>::Result TPiNum;
// // typedef typename IPowBig<10,17>::Result TPiDen;
// // typedef SFraction<TPiNum,TPiDen> TPiShort;
// 
//   
// typedef EX::SinAcc<X4,2> SinA; 
// typedef SinA::Result SinPi;
// typedef EX::FractionToDecimal<SinPi,20,10>::Result SinPiDec;
//typedef EX::FractionToDecimal<X4,20,10>::Result X4Dec;
//typedef Loki::TL::Print<SinPiDec2>::Result PF;

// typedef SFraction<SInt<1>,SInt<4> > F;
// typedef typename EX::PiAcc<2>::Result TPi;
// typedef typename Mult<TPi,F>::Result X;
// typedef typename EX::SinAcc<X,2>::Result SinPi4;
// typedef EX::SinPiFrac<1,6,2>::Result SinPi4;
//typedef EX::SinPiFrac<1,8,2>::Result SinPi8;
//typedef EX::SinPiFrac<1,4,2>::Result CosPi4;
// // typedef EX::FractionToDecimal<TPi,20,10>::Result TPiDec;
// // typedef EX::FractionToDecimal<X,20,10>::Result XDec;

//typedef EX::FractionToDecimal<SinPi4,20,10>::Result SinPi4Dec;
//typedef EX::FractionToDecimal<CosPi4,20,10>::Result CosPi4Dec;
// typedef EX::FractionToDecimal<SinPi8,20,10>::Result SinPi8Dec;

// typedef typename Sub<SInt<1>,typename Mult<SInt<2>,
//         typename Mult<SinPi4,SinPi4>::Result>::Result>::Result WR;

// cout << sin(M_PI/6.) << endl;
// Cout<SinPi4Dec>::apply(cout);
// cout << endl;

static const int_t N = 6;
typedef GenerateRootList<N,1,2>::Result List;
PrintRootList<List,N>::apply();

// cout << sin(M_PI/8.) << endl;
// Cout<SinPi4Dec3::Result>::apply(cout);
// cout << endl;
// Cout<SinPi4s>::apply(cout);
// cout << endl;
// Cout<SinPi4sDec>::apply(cout);
// cout << endl;
// Cout<SinPiDec2::Num>::apply(cout);
// cout << endl;
// Cout<X4Dec>::apply(cout);
// cout << endl;
// cout << EX::Compute<SinPi,2>::value() << endl;

// cout << endl;
// Cout<StartValue>::apply(cout);
// cout << endl;
// Cout<StartDec>::apply(cout);
// cout << endl;
// Cout<Step>::apply(cout);
// cout << endl;
// Cout<StepDec>::apply(cout);
// cout << endl;
// Cout<NextValue>::apply(cout);
// cout << endl;
// Cout<NextDec>::apply(cout);
// cout << endl;
}


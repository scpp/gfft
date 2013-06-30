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
    \brief Sample program to represent %GFFT usage
*/

#include <iostream>

#include "loki/static_check.h"

#include "gfft.h"

using namespace std;

using namespace GFFT;


void dft1(double* output_data, const double* input_data, const unsigned int size, bool inverse)
{
  double pi2 = (inverse) ? 2.0 * M_PI : -2.0 * M_PI;
  double a, ca, sa;
  double invs = 1.0 / size;
  for(unsigned int y = 0; y < size; y++) {
    output_data[2*y] = 0;
    output_data[2*y+1] = 0;
    for(unsigned int x = 0; x < size; x++) {
      a = pi2 * y * x * invs;
      ca = cos(a);
      sa = sin(a);
      output_data[2*y]   += input_data[2*x] * ca - input_data[2*x+1] * sa;
      output_data[2*y+1] += input_data[2*x] * sa + input_data[2*x+1] * ca;
    }
    if(inverse) {
      output_data[2*y]   *= invs;
      output_data[2*y+1] *= invs;
    }
  }
}


typedef DOUBLE ValueType;
//typedef typename GenNumList<2, 3>::Result NList;
typedef TYPELIST_5(SIntID<2>, SIntID<3>, SIntID<4>, SIntID<5>, SIntID<19>) NList;
typedef GenerateTransform<NList, ValueType, TransformTypeGroup::FullList, SIntID<1>, ParallelizationGroup::Default, DecimationGroup::FullList > TransformSet;

int main(int argc, char *argv[])
{
//     unsigned int p = 2;
//     unsigned long i, n = (TransformType::ID == RDFT::ID) ? (1<<(p-1)) : (1<<p);
    int_t i, n = 11;
//    cin >> n;
/*   
    typedef DFT TransformType;

    TransformSet gfft;
    TransformSet::ObjectType* fftobj  = gfft.CreateTransformObject(n, ValueType::ID, TransformType::ID, 1, 
								   ParallelizationGroup::Default::ID, INFREQ::ID);
    TransformSet::ObjectType* ifftobj = gfft.CreateTransformObject(n, ValueType::ID, TransformType::Inverse::ID, 1, 
								   ParallelizationGroup::Default::ID, INFREQ::ID);
    
// create sample data
    ValueType::ValueType* data = new ValueType::ValueType [2*n];
    ValueType::ValueType* dataout = new ValueType::ValueType [2*n];
    ValueType::ValueType* data1 = new ValueType::ValueType [2*n];
    ValueType::ValueType* dataout1 = new ValueType::ValueType [2*n];
    for (i=0; i < n; ++i) {
       data[2*i] = 2*i;
       data[2*i+1] = 2*i+1; 
       data1[2*i] = 2*i;
       data1[2*i+1] = 2*i+1; 
       dataout[2*i] = 0;
       dataout[2*i+1] = 0; 
       dataout1[2*i] = 0;
       dataout1[2*i+1] = 0; 
    }

// print out sample data
    cout<<"Input data:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;

// apply FFT in-place
//     fftobj->fft(data);
//     for (i=0; i < 2*n; ++i)
//       dataout[i] = data[i];

//      fftobj->fft(data);
    fftobj->fft(data, dataout);
    dft1(dataout1, data1, n, false);

// print out transformed data
    cout<<"Result of transform:"<<endl;
//      for (i=0; i < n; ++i)
//        cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
   for (i=0; i < n; ++i)
     cout<<"("<<dataout[2*i]<<","<<dataout[2*i+1]<<")"<<endl;

//    ifftobj->fft(data);
   ifftobj->fft(dataout, data);
   dft1(data1, dataout1, n, true);

// print out transformed data
    cout<<"Result of backward transform:"<<endl;
    for (i=0; i < n; ++i)
      cout<<"("<<data[2*i]<<","<<data[2*i+1]<<")"<<endl;
    
    cout<<"Check against DFT:"<<endl;
    cout.precision(2);
    double mx1(-1), mx2(-1);
    for (i=0; i < n; ++i) {
      cout<<"("<<fabs(data[2*i]-data1[2*i])<<","<<fabs(data[2*i+1]-data1[2*i+1])<<")   \t("
      <<fabs(dataout[2*i]-dataout1[2*i])<<","<<fabs(dataout[2*i+1]-dataout1[2*i+1])<<")"<<endl;
      mx1 = max(mx1, fabs(data[2*i]-data1[2*i]));
      mx1 = max(mx1, fabs(data[2*i+1]-data1[2*i+1]));
      mx2 = max(mx1, fabs(dataout[2*i]-dataout1[2*i]));
      mx2 = max(mx1, fabs(dataout[2*i+1]-dataout1[2*i+1]));
    }
    cout<<"---------------------------------------------"<<endl;
    cout << mx1 << "                 " << mx2 << endl;
*/

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
cout.precision(16);
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

typedef EX::PiLen<1> MetaPi;
typedef Simplify<MetaPi::Result>::Result TPi;

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

//typedef SFraction<SInt<157>,SInt<100> > X;
//typedef SInt<3> X;
typedef TPi X;

typedef Simplify<EX::CosAcc<X,1>::Result>::Result CosPi;
typedef EX::FractionToDecimal<CosPi,20,10>::Result CosPiDec;
cout << cos(M_PI) << endl;
Cout<CosPiDec>::apply(cout);
cout << endl;

}


/***************************************************************************
 *   Copyright (C) 2008 by Volodymyr Myrnyy                                *
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
    \brief Performance benchmarks for GFFT
*/

#include <iostream>
#include <iomanip>

//#include "gfft.h"
#include "gfftconf.h"

#include "timer.h"

#include "boost/date_time/posix_time/posix_time.hpp"

using namespace std;

typedef double ValueType;


template<class GenListResult>
class GFFTbench;

template<class H, class T>
class GFFTbench<Loki::Typelist<H,T> > {
  typedef typename H::ValueType::ValueType Tp;
  GFFTbench<T> next;
  H gfft;
public:
  void apply(const int hardware_id, const int system_id,
             const int compiler_id, const int release_id)
  {
     next.apply(hardware_id,system_id,compiler_id,release_id);

     size_t i,it;
     double t,mt;
     char space = '\t';
     progClock* cl=new progClock();

     it = size_t(1000000./(double)H::Len)+1;

     Tp* data = new Tp [2*H::Len*it];

      // initial data
     for (i=0; i<2*H::Len*it; ++i)
       data[i] = 0;

/*      for (unsigned int j=0; j<it; ++j) {
        for (i=0; i < n; ++i) {
          data[2*(n*j+i)] = (2*i+j)/(double)n;
          data[2*(n*j+i)+1] = (2*i+j+1)/(double)n;
        }
      }*/
      Tp* d=data;
      mt = 1e+100;
      for (i=0; i<3; ++i) {
        d=data;
        cl->start();
        for (size_t j=0; j<it; ++j) {
          gfft.fft(d);
          d+=2*H::Len;
        }
        t = cl->GetTime();
        if (t<mt) mt=t;
      }
      mt /= (double)it;
      cout<<hardware_id<<space<<system_id<<space
          <<compiler_id<<space<<release_id<<space
          <<H::TransformType::ID<<space
          <<H::ValueType::ID<<space
          <<H::DecimationType::ID<<space
          <<H::DirectionType::ID<<space
          <<0<<space
          <<H::PLen-9<<space<<mt<<endl;
      //cout<<k<<"  "<<mt<<"  "<<norm2(d-2*n,2*n)<<"  "<<norminf(d-2*n,2*n)<<endl;

      delete [] data;
      delete cl;

   }
};

template<>
class GFFTbench<Loki::NullType> {
public:
  void apply(const int, const int, const int, const int) { }
};


int main(int argc, char *argv[])
{

//     omp_set_num_threads(4);
//     omp_set_nested(true);

//     unsigned int i,p=2;
//     unsigned int n= 1<<p;

    if (argc<4) {
      cout<<"usage: bench <hardware_id> <system_id> <compiler_id> <release_id>"<<endl;
      exit(1);
    }

// There are three ways to create object to perform FFT of the length 2^p
// 1) Singleton holds the object factory for GFFT
//     DFT::GFFT_Singleton<Min,Max,ValueType,DFT::COMPLEX,DFT::INTIME,DFT::FORWARD>* gfft;
//     DFT::AbstractFFT<ValueType>* fftobj = gfft->Instance().CreateObject(p);
//
//     DFT::GFFT_Singleton<Min,Max,ValueType,DFT::COMPLEX,DFT::INTIME,DFT::BACKWARD>* igfft;
//     DFT::AbstractFFT<ValueType>* ifftobj = igfft->Instance().CreateObject(p);

// 2) Create the object factory without singleton

    typedef DFT::GenList<17,25,DFT::DOUBLE,DFT::COMPLEX,DFT::FORWARD,DFT::OpenMP<4>,DFT::INFREQ> List;
//    typedef DFT::GenList<10,15> List;

//   typedef DFT::Print<List::Result>::Result deb;
//    Loki::Factory<DFT::AbstractFFT<ValueType>,unsigned int> gfft;
//    FactoryInit<List::Result>::apply(gfft);

//    unsigned int id1[5] = {p-1,0,1,1,0};
//    unsigned int id2[5] = {p-1,0,1,1,1};
//    unsigned int p1 = List::trans_id(id1);
//    unsigned int p2 = List::trans_id(id2);
//    cout<<p1<<" "<<p2<<endl;
//    DFT::AbstractFFT<ValueType>* fftobj = gfft.CreateObject(p1);
//    DFT::AbstractFFT<ValueType>* ifftobj = gfft.CreateObject(p2);

   int hardware_id = atoi(argv[1]);
   int system_id = atoi(argv[2]);
   int compiler_id = atoi(argv[3]);
   int release_id = atoi(argv[4]);

   cout<<setprecision(12);

   GFFTbench<List::Result> bench;

   using namespace boost::posix_time;
   using namespace boost::gregorian;

   time_duration td;
   ptime t1,t2;
   td = seconds(0);
   t1 = microsec_clock::universal_time();

   bench.apply(hardware_id,system_id,compiler_id,release_id);

   t2 = microsec_clock::universal_time();
   td = t2 - t1;
   double rt = (td.total_seconds()*1000000+td.fractional_seconds());
   cout<<rt<<"ms"<<endl;

   return 0;
}


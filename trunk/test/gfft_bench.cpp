/***************************************************************************
 *   Copyright (C) 2011 by Vladimir Mirnyy                                 *
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
    \brief Performance benchmarks for %GFFT
*/

#include <iostream>
#include <iomanip>
#include <ctime>

#include "gfft.h"

#include "boost/date_time/posix_time/posix_time.hpp"

using namespace std;

using namespace GFFT;

using namespace boost::posix_time;
using namespace boost::gregorian;

static const unsigned int NumberOfIterations = 5000000;
static const char TransformType_Name[][17] = {" forward", "backward", "    real forward", "   real backward"};
static const char ValueType_Name[][15] = {"    double    ", "     float    ", "complex double", " complex float"};

template<class GenListResult>
class GFFTbench;

template<class H, class T>
class GFFTbench<Loki::Typelist<H,T> > {
  typedef typename H::ValueType::ValueType Tp;
  GFFTbench<T> next;
  H gfft;
public:
  void cputime()
  {
     next.cputime();

     size_t i,it;
     double t,mt;
     char space = '\t';
     clock_t time1, time2;
     
     it = size_t(NumberOfIterations/(double)H::Len)+1;

     Tp* data = new Tp [2*H::Len*it];

      // initial data is taken zero to avoid overflows
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
        // CPU-time
	time1 = clock();
        for (size_t j=0; j<it; ++j) {
          gfft.fft(d);
          d+=2*H::Len;
        }
        time2=clock();
	t=static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
        if (t<mt) mt=t;
      }
      mt /= (double)it;
      cout<< TransformType_Name[H::TransformType::ID] << space
          << ValueType_Name[H::ValueType::ID] << "     " << space
          //<<H::DecimationType::ID<<space  // is always 'in-frequency', since it is slightly faster
          << H::ParallType::ID+1 << space
          << H::PLen << space
          << mt << endl;

      delete [] data;
  }
  
  void realtime()
  {
     next.realtime();

     size_t i,it;
     char space = '\t';

     it = size_t(NumberOfIterations/(double)H::Len)+1;

     Tp* data = new Tp [2*H::Len*it];

      // initial data
     for (i=0; i<2*H::Len*it; ++i)
       data[i] = 0;

      time_duration td;
      ptime t1,t2;
      // real time
      td = seconds(0);
      t1 = microsec_clock::universal_time();

      Tp* d=data;
      for (i=0; i<3; ++i) {
        d=data;
        for (size_t j=0; j<it; ++j) {
          gfft.fft(d);
          d+=2*H::Len;
        }
      }
      t2 = microsec_clock::universal_time();
      td = t2 - t1;
      double rt = (td.total_seconds()*1000000+td.fractional_seconds())/(3.*it*1e+6);
      cout<< TransformType_Name[H::TransformType::ID] << space
          << ValueType_Name[H::ValueType::ID] << "     " << space
          //<<H::DecimationType::ID<<space  // is always 'in-frequency', since it is slightly faster
          << H::ParallType::ID+1 << space
          << H::PLen << space
          << rt << endl;
      //cout<<k<<"  "<<mt<<"  "<<norm2(d-2*n,2*n)<<"  "<<norminf(d-2*n,2*n)<<endl;

      delete [] data;
   }
};

template<>
class GFFTbench<Loki::NullType> {
public:
  void cputime() { }
  void realtime() { }
};

void print_header() 
{
  cout<<"-------------------------------------------------------------------------"<<endl;
  cout<<" Direction        Value type     N-threads  PowerOfTwo       Time [s]"<<endl;
  cout<<"-------------------------------------------------------------------------"<<endl;
}


static const unsigned int MinP = 15;
static const unsigned int MaxP = 15;


int main(int argc, char *argv[])
{

    typedef GenerateTransform<MinP, MaxP, GFFT::DOUBLE, TransformTypeGroup::Default> List_ds;
    typedef GenerateTransform<MinP, MaxP, GFFT::FLOAT,  TransformTypeGroup::Default> List_fs;
    typedef GenerateTransform<MinP, MaxP, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default> List_cds;
    typedef GenerateTransform<MinP, MaxP, GFFT::COMPLEX_FLOAT,  TransformTypeGroup::Default> List_cfs;

    typedef TYPELIST_2(OpenMP<2>, OpenMP<4>) ParallList;
    typedef GenerateTransform<MinP, MaxP, GFFT::DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallList> List_dp;
    typedef GenerateTransform<MinP, MaxP, GFFT::FLOAT,  TransformTypeGroup::Default, SIntID<1>, ParallList> List_fp;
    typedef GenerateTransform<MinP, MaxP, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallList> List_cdp;
    typedef GenerateTransform<MinP, MaxP, GFFT::COMPLEX_FLOAT,  TransformTypeGroup::Default, SIntID<1>, ParallList> List_cfp;

   cout<<setprecision(12)<<setw(12);

//   GFFTbench<List_ds::Result> bench_ds;
//   GFFTbench<List_fs::Result> bench_fs;
   GFFTbench<List_cds::Result> bench_cds;
//   GFFTbench<List_cfs::Result> bench_cfs;
/*
   GFFTbench<List_dp::Result> bench_dp;
   GFFTbench<List_fp::Result> bench_fp;
   GFFTbench<List_cdp::Result> bench_cdp;
   GFFTbench<List_cfp::Result> bench_cfp;

   cout<<endl<<"================= CPU-time, single-threaded ================"<<endl;
   cout<<"Double precision"<<endl;
   print_header();
   bench_ds.cputime();
   
   cout<<endl<<"Single precision"<<endl;
   print_header();
   bench_fs.cputime();
   
   cout<<endl<<"Complex numbers, double precision"<<endl;
   print_header();
   bench_cds.cputime();
   
   cout<<endl<<"Complex numbers, single precision"<<endl;
   print_header();
   bench_cfs.cputime();
*/
/*
   cout<<endl<<"================= Real time, single-threaded ================"<<endl;
   cout<<"Double precision"<<endl;
   print_header();
   bench_ds.realtime();
   
   cout<<endl<<"Single precision"<<endl;
   print_header();
   bench_fs.realtime();
*/   
   cout<<endl<<"Complex numbers, double precision"<<endl;
   print_header();
   bench_cds.realtime();
/*   
   cout<<endl<<"Complex numbers, single precision"<<endl;
   print_header();
   bench_cfs.realtime();

   cout<<endl<<"================= Real time, multi-threaded ================"<<endl;
   cout<<"Double precision"<<endl;
   print_header();
   bench_dp.realtime();
   
   cout<<endl<<"Single precision"<<endl;
   print_header();
   bench_fp.realtime();
   
   cout<<endl<<"Complex numbers, double precision"<<endl;
   print_header();
   bench_cdp.realtime();

   cout<<endl<<"Complex numbers, single precision"<<endl;
   print_header();
   bench_cfp.realtime();
*/
   return 0;
}


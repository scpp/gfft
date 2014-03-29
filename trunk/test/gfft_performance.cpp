/***************************************************************************
 *   Copyright (C) 2009 by Volodymyr Myrnyy                                *
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

const char space = '\t';

static const char TransformType_Name[][17] = {" forward", "backward", "    real forward", "   real backward"};
static const char ValueType_Name[][15] = {"    double    ", "     float    ", "complex double", " complex float"};
static const char Place_Name[][17] = {"    in-place    ", "  out-of-place  "};


template<class T>
class GFFTbenchBase
{
public:
  void init(T* src, const uint_t len, const int_t it)
  {
      // initial data
      for (uint_t i=0; i<2*len*it; ++i)
	src[i] = 0;

//       for (unsigned int j=0; j<it; ++j) {
// 	for (uint_t i=0; i < n; ++i) {
// 	  src[2*(n*j+i)] = (2*i+j)/(double)n;
// 	  src[2*(n*j+i)+1] = (2*i+j+1)/(double)n;
// 	}
//       }
  }
  void print_line(const uint_t TransformTypeID, const uint_t ValueTypeID, 
		  const uint_t PlaceTypeID, const uint_t ParallTypeID, 
		  const uint_t len, const double t)
  {
     cout<<TransformType_Name[TransformTypeID]<<space
         <<ValueType_Name[ValueTypeID]<<space
         <<Place_Name[PlaceTypeID]<<space
         <<ParallTypeID+1<<space
         <<len<<space
         <<t<<endl;
      //cout<<k<<"  "<<mt<<"  "<<norm2(d-2*n,2*n)<<"  "<<norminf(d-2*n,2*n)<<endl;
  }
  
};

template<class NList, class Place,
int Counter = Loki::TL::Length<NList>::value>
class GFFTbench;

template<class H, class T, int Counter>
class GFFTbench<Loki::Typelist<H,T>,IN_PLACE,Counter> 
: public GFFTbenchBase<typename H::ValueType::ValueType> 
{
  typedef typename H::ValueType::ValueType Tp;
  typedef GFFTbenchBase<Tp> Base;
  GFFTbench<T,IN_PLACE,Counter-1> next;
  H gfft;
public:
  void cputime()
  {
     next.cputime();

     uint_t i,it;
     double t,mt;
     clock_t time1, time2;

     it = size_t(2000000./(double)H::Len)+1;
     Tp* data    = new Tp [2*H::Len*it];
     Base::init(data, H::Len, it);
     
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
	t=static_cast<double>(time2-time1)/static_cast<double>(CLOCKS_PER_SEC);
        if (t<mt) mt=t;
     }
     
     mt /= (double)it;
     Base::print_line(H::TransformType::ID,H::ValueType::ID,H::PlaceType::ID,H::ParallType::ID,H::Len,mt);

     delete [] data;
  }
  void realtime()
  {
     next.realtime();

     uint_t i,it;
     
     it = size_t(2000000./(double)H::Len)+1;
     Tp* data    = new Tp [2*H::Len*it];
     Base::init(data, H::Len, it);
 
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
     Base::print_line(H::TransformType::ID,H::ValueType::ID,H::PlaceType::ID,H::ParallType::ID,H::Len,rt);

     delete [] data;
   }
};

template<class H, class T, int Counter>
class GFFTbench<Loki::Typelist<H,T>,OUT_OF_PLACE,Counter> 
: public GFFTbenchBase<typename H::ValueType::ValueType> 
{
  typedef typename H::ValueType::ValueType Tp;
  typedef GFFTbenchBase<Tp> Base;
  GFFTbench<T,OUT_OF_PLACE,Counter-1> next;
  H gfft;
public:
  void cputime()
  {
     next.cputime();

     uint_t i,it;
     double t,mt;
     clock_t time1, time2;

     it = size_t(2000000./(double)H::Len)+1;
     Tp* data    = new Tp [2*H::Len*it];
     Tp* dataout = new Tp [2*H::Len*it];
     Base::init(data, H::Len, it);
     
     Tp* d=data;
     mt = 1e+100;
     for (i=0; i<3; ++i) {
        d=data;
        // CPU-time
	time1 = clock();
        for (size_t j=0; j<it; ++j) {
          gfft.fft(d, dataout);
//          gfft.fft(d);
          d+=2*H::Len;
        }
        time2=clock();
	t=static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
        if (t<mt) mt=t;
     }
     
     mt /= (double)it;
     Base::print_line(H::TransformType::ID,H::ValueType::ID,H::PlaceType::ID,H::ParallType::ID,H::Len,mt);

     delete [] data;
     delete [] dataout;
  }
  void realtime()
  {
     next.realtime();

     uint_t i,it;
     
     it = size_t(10000000./(double)H::Len)+1;
     Tp* data    = new Tp [2*H::Len*it];
     Tp* dataout = new Tp [2*H::Len*it];
     Base::init(data, H::Len, it);
 
     time_duration td;
     ptime t1,t2;
     // real time
     td = seconds(0);
     t1 = microsec_clock::universal_time();

     Tp* d=data;
     for (i=0; i<3; ++i) {
        d=data;
        for (size_t j=0; j<it; ++j) {
          gfft.fft(d, dataout);
//          gfft.fft(d);
          d+=2*H::Len;
        }
     }
     t2 = microsec_clock::universal_time();
     td = t2 - t1;
     double rt = (td.total_seconds()*1000000+td.fractional_seconds())/(3.*it*1e+6);
     Base::print_line(H::TransformType::ID,H::ValueType::ID,H::PlaceType::ID,H::ParallType::ID,H::Len,rt);

     delete [] data;
     delete [] dataout;
   }
};

template<>
class GFFTbench<Loki::NullType,OUT_OF_PLACE> {
public:
  void cputime() { }
  void realtime() { }
};


void print_header() 
{
  cout<<"---------------------------------------------------------------------------------"<<endl;
  cout<<" Direction        Value type          Place       N-threads    Length   Time [s] "<<endl;
  cout<<"---------------------------------------------------------------------------------"<<endl;
}

static const unsigned int MinP = 1;
static const unsigned int MaxP = 8;

int main(int argc, char *argv[])
{

   //typedef GenNumList<20000, 20010, SIntID>::Result NList;
   typedef GenPowerList<MinP, MaxP, 2>::Result NList;
    
   // Single-threaded out-of-place transforms
   typedef GenerateTransform<NList, GFFT::DOUBLE, TransformTypeGroup::Default> List_ds;
   typedef GenerateTransform<NList, GFFT::FLOAT, TransformTypeGroup::Default>  List_fs;
   typedef GenerateTransform<NList, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default> List_cds;
   typedef GenerateTransform<NList, GFFT::COMPLEX_FLOAT, TransformTypeGroup::Default>  List_cfs;

   // Single- and multi-threaded out-of-place transforms
   typedef GenerateTransform<NList, GFFT::DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::Default, OUT_OF_PLACE> List_dp;
   typedef GenerateTransform<NList, GFFT::FLOAT, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::FullList, OUT_OF_PLACE>  List_fp;
   typedef GenerateTransform<NList, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::FullList, OUT_OF_PLACE> List_cdp;
   typedef GenerateTransform<NList, GFFT::COMPLEX_FLOAT, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::FullList, OUT_OF_PLACE>  List_cfp;

   cout.precision(4);
   
   GFFTbench<List_ds::Result,List_ds::PlaceType> bench_ds;
   GFFTbench<List_fs::Result,List_fs::PlaceType> bench_fs;
   GFFTbench<List_cds::Result,List_cds::PlaceType> bench_cds;
   GFFTbench<List_cfs::Result,List_cfs::PlaceType> bench_cfs;

//    print_header();
//    bench_ds.cputime();
//    bench_fs.cputime();
//    bench_cds.cputime();
//    bench_cfs.cputime();

   GFFTbench<List_dp::Result,List_dp::PlaceType> bench_dp;
   GFFTbench<List_fp::Result,List_fp::PlaceType> bench_fp;
   GFFTbench<List_cdp::Result,List_cdp::PlaceType> bench_cdp;
   GFFTbench<List_cfp::Result,List_cfp::PlaceType> bench_cfp;

   print_header();
   bench_dp.realtime();
//    bench_fp.realtime();
//    bench_cdp.realtime();
//    bench_cfp.realtime();

   return 0;
}


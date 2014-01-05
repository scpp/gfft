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

#include "gfft.h"
#include "timer.h"

#include "boost/date_time/posix_time/posix_time.hpp"

using namespace std;

using namespace GFFT;

using namespace boost::posix_time;
using namespace boost::gregorian;

const char space = '\t';


template<class NList>
class GFFTbench;

template<class H, class T>
class GFFTbench<Loki::Typelist<H,T> > {
  typedef typename H::ValueType::ValueType Tp;
  GFFTbench<T> next;
  H gfft;
public:
  void cputime(const int hardware_id, const int system_id,
               const int compiler_id, const int software_id)
  {
     next.cputime(hardware_id,system_id,compiler_id,software_id);

     size_t i,it;
     double t,mt;
     progClock* cl=new progClock();

     it = size_t(2000000./(double)H::Len)+1;

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
        // CPU-time
        cl->start();
        for (size_t j=0; j<it; ++j) {
          gfft.fft(d);
          d+=2*H::Len;
        }
        t = cl->GetTime();
        if (t<mt) mt=t;
      }
      mt /= (double)it;
      cout<<hardware_id<<space
          <<system_id<<space
          <<compiler_id<<space
          <<software_id<<space
          <<H::TransformType::ID<<space
          <<H::ValueType::ID<<space
          <<H::PlaceType::ID<<space
          <<H::ParallType::ID+1<<space
          <<0<<space    // array type
          <<1<<space    // dim
          <<3<<space    // CPU-time
          <<H::Len<<space
          <<mt<<endl;
      //cout<<k<<"  "<<mt<<"  "<<norm2(d-2*n,2*n)<<"  "<<norminf(d-2*n,2*n)<<endl;

      delete [] data;
      delete cl;

  }
  void realtime(const int hardware_id, const int system_id,
             const int compiler_id, const int software_id)
  {
     next.realtime(hardware_id,system_id,compiler_id,software_id);

     size_t i,it;

     it = size_t(2000000./(double)H::Len)+1;

     Tp* data = new Tp [2*H::Len*it];
     Tp* dataout = new Tp [2*H::Len*it];

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
          gfft.fft(d, dataout);
//          gfft.fft(d);
          d+=2*H::Len;
        }
      }
      t2 = microsec_clock::universal_time();
      td = t2 - t1;
      double rt = (td.total_seconds()*1000000+td.fractional_seconds())/(3.*it*1e+6);
      cout<<hardware_id<<space
          <<system_id<<space
          <<compiler_id<<space
          <<software_id<<space
          <<H::TransformType::ID<<space
          <<H::ValueType::ID<<space
          <<H::PlaceType::ID<<space
          <<H::ParallType::ID+1<<space
          <<0<<space    // array type
          <<1<<space    // dim
          <<2<<space    // real time
          <<H::Len<<space
          <<rt<<endl;
      //cout<<k<<"  "<<mt<<"  "<<norm2(d-2*n,2*n)<<"  "<<norminf(d-2*n,2*n)<<endl;

      delete [] dataout;
      delete [] data;
   }
};

template<>
class GFFTbench<Loki::NullType> {
public:
  void cputime(const int, const int, const int, const int) { }
  void realtime(const int, const int, const int, const int) { }
};


static const unsigned int MinP = 1;
static const unsigned int MaxP = 8;


int main(int argc, char *argv[])
{

    if (argc<4) {
      cout<<"usage: bench <hardware_id> <system_id> <compiler_id> <release_id>"<<endl;
      exit(1);
    }

//     typedef GenNumList<MinP, MaxP, Power2holder>::Result NList;
    typedef GenNumList<20000, 20010, SIntID>::Result NList;
    //typedef GenPowerList<20, 26, 2>::Result NList;
    typedef GenerateTransform<NList, GFFT::DOUBLE, TransformTypeGroup::Default, SIntID<1>, 
       OpenMP<1>, OUT_OF_PLACE> List_ds;
//        ParallelizationGroup::Default, OUT_OF_PLACE> List_ds;
//    typedef GeneratePower2Transform<MinP, MaxP, GFFT::FLOAT,  TransformTypeGroup::Default> List_fs;
//     typedef GeneratePower2Transform<MinP, MaxP, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default> List_cds;
//     typedef GeneratePower2Transform<MinP, MaxP, GFFT::COMPLEX_FLOAT,  TransformTypeGroup::Default> List_cfs;

//     typedef TYPELIST_2(OpenMP<2>, OpenMP<4>) ParallList;
//     typedef GenerateTransform<MinP, MaxP, GFFT::DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallList> List_dp;
//     typedef GenerateTransform<MinP, MaxP, GFFT::FLOAT,  TransformTypeGroup::Default, SIntID<1>, ParallList> List_fp;
//     typedef GenerateTransform<MinP, MaxP, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallList> List_cdp;
//     typedef GenerateTransform<MinP, MaxP, GFFT::COMPLEX_FLOAT,  TransformTypeGroup::Default, SIntID<1>, ParallList> List_cfp;

   int hardware_id = atoi(argv[1]);
   int system_id = atoi(argv[2]);
   int compiler_id = atoi(argv[3]);
   int release_id = atoi(argv[4]);

   cout<<setprecision(12);

   GFFTbench<List_ds::Result> bench_ds;
//   GFFTbench<List_fs::Result> bench_fs;
//    GFFTbench<List_cds::Result> bench_cds;
//    GFFTbench<List_cfs::Result> bench_cfs;

  // Header
      cout<<"hdw"<<space
          <<"sys"<<space
          <<"comp"<<space
          <<"soft"<<space
          <<"fwd/bkw"<<space
          <<"dbl/fl"<<space
          <<"inT/inF"<<space
          <<"parall"<<space
          <<0<<space    // array type
          <<1<<space    // dim
          <<2<<space    // real time
          <<"len"<<space
          <<"time"<<endl;
   cout.precision(4);
//   bench_ds.cputime(hardware_id,system_id,compiler_id,release_id);
//   bench_fs.cputime(hardware_id,system_id,compiler_id,release_id);
//    bench_cds.cputime(hardware_id,system_id,compiler_id,release_id);
//    bench_cfs.cputime(hardware_id,system_id,compiler_id,release_id);

   bench_ds.realtime(hardware_id,system_id,compiler_id,release_id);
//   bench_fs.realtime(hardware_id,system_id,compiler_id,release_id);
//   bench_cds.realtime(hardware_id,system_id,compiler_id,release_id);
//   bench_cfs.realtime(hardware_id,system_id,compiler_id,release_id);

/*   bench_dp.realtime(hardware_id,system_id,compiler_id,release_id);
   bench_fp.realtime(hardware_id,system_id,compiler_id,release_id);*/
   //bench_cdp.realtime(hardware_id,system_id,compiler_id,release_id);
//   bench_cfp.realtime(hardware_id,system_id,compiler_id,release_id);

   return 0;
}


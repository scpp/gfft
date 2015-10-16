/***************************************************************************
 *   Copyright (C) 2009-2015 by Vladimir Mirnyy                            *
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

#include "gfft_performance.h"

using namespace std;

using namespace GFFT;


void print_header() 
{
  cout<<"---------------------------------------------------------------------------------"<<endl;
  cout<<" Direction        Value type          Place       N-threads    Length   Time [s] "<<endl;
  cout<<"---------------------------------------------------------------------------------"<<endl;
}

static const unsigned int MinP = 2;
static const unsigned int MaxP = 10;

int main(int argc, char *argv[])
{

   //typedef GenNumList<20000, 20010, SIntID>::Result NList;
   typedef GenPowerList<MinP, MaxP, 2>::Result NList;
    
   // Single-threaded out-of-place transforms
//    typedef GenerateTransform<NList, GFFT::DOUBLE, TransformTypeGroup::Default> List_ds;
//    typedef GenerateTransform<NList, GFFT::FLOAT, TransformTypeGroup::Default>  List_fs;
//    typedef GenerateTransform<NList, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default> List_cds;
//    typedef GenerateTransform<NList, GFFT::COMPLEX_FLOAT, TransformTypeGroup::Default>  List_cfs;

   // Single- and multi-threaded out-of-place transforms
   typedef GenerateTransform<NList, GFFT::DOUBLE, TransformTypeGroup::Default, SIntID<1>, OpenMP<4>, OUT_OF_PLACE> List_dp;
//   typedef GenerateTransform<NList, GFFT::DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::FullList, OUT_OF_PLACE> List_dp;
//    typedef GenerateTransform<NList, GFFT::FLOAT, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::FullList, OUT_OF_PLACE>  List_fp;
//    typedef GenerateTransform<NList, GFFT::COMPLEX_DOUBLE, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::FullList, OUT_OF_PLACE> List_cdp;
//    typedef GenerateTransform<NList, GFFT::COMPLEX_FLOAT, TransformTypeGroup::Default, SIntID<1>, ParallelizationGroup::FullList, OUT_OF_PLACE>  List_cfp;

   cout.precision(4);
   
//    GFFTbench<List_ds::Result,List_ds::PlaceType> bench_ds;
//    GFFTbench<List_fs::Result,List_fs::PlaceType> bench_fs;
//    GFFTbench<List_cds::Result,List_cds::PlaceType> bench_cds;
//    GFFTbench<List_cfs::Result,List_cfs::PlaceType> bench_cfs;

//    print_header();
//    bench_ds.cputime();
//    bench_fs.cputime();
//    bench_cds.cputime();
//    bench_cfs.cputime();

   GFFTbench<List_dp::Result,List_dp::PlaceType> bench_dp;
//    GFFTbench<List_fp::Result,List_fp::PlaceType> bench_fp;
//    GFFTbench<List_cdp::Result,List_cdp::PlaceType> bench_cdp;
//    GFFTbench<List_cfp::Result,List_cfp::PlaceType> bench_cfp;

   //print_header();
   bench_dp.realtime();
//    bench_fp.realtime();
//    bench_cdp.realtime();
//    bench_cfp.realtime();

   return 0;
}


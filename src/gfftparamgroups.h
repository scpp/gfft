/***************************************************************************
 *   Copyright (C) 2012 by Volodymyr Myrnyy                                *
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


#ifndef __gfftparamgroups_h
#define __gfftparamgroups_h

/** \file
    \brief GFFT parameter group classes
*/

#include "gfftpolicy.h"


namespace GFFT {

/*!
\defgroup gr_groups GFFT parameter groups
\brief Classes that represent groups of parameters to define transform
*/

/// \brief Lists all acceptable value types and the default one
/// \ingroup gr_groups
struct ValueTypeGroup
{
  typedef TYPELIST_4(DOUBLE,FLOAT,COMPLEX_DOUBLE,COMPLEX_FLOAT) FullList;
  static const uint Length = 4;
  typedef DOUBLE Default;
  static const uint default_id = DOUBLE::ID;
};

/// \brief Lists all acceptable types of Fast Fourier transform
/// \ingroup gr_groups
struct TransformTypeGroup
{
  typedef TYPELIST_4(DFT,IDFT,RDFT,IRDFT) FullList;
  static const uint Length = 4;
//  typedef TYPELIST_2(DFT,IDFT) Default;
  typedef DFT Default;
  static const uint default_id = DFT::ID;
};

/// \brief Lists all acceptable parallelization methods
/// \ingroup gr_groups
struct ParallelizationGroup
{
  typedef TYPELIST_2(Serial,OpenMP<2>) FullList;
  static const uint Length = 2;
  typedef Serial Default;
  static const uint default_id = Serial::ID;
};

/// \brief Lists all acceptable decimation versions
/// \ingroup gr_groups
struct DecimationGroup
{
  typedef TYPELIST_2(INTIME,INFREQ) FullList;
  static const uint Length = 2;
  typedef INFREQ Default;
  static const uint default_id = INFREQ::ID;
};

/// \brief Lists in-place and out-of-place FFT algorithms
/// \ingroup gr_groups
struct PlaceGroup
{
  typedef TYPELIST_2(IN_PLACE,OUT_OF_PLACE) FullList;
  static const uint Length = 2;
  typedef OUT_OF_PLACE Default;
  static const uint default_id = OUT_OF_PLACE::ID;
};
  
}  //namespace

#endif

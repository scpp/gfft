/***************************************************************************
 *   Copyright (C) 2006-2009 by Volodymyr Myrnyy                           *
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

#ifndef __gfft_h
#define __gfft_h

/** \file
    \brief Main include file, which just include all other headers
*/

// General Doxygen documentation
#ifdef GFFTDOC
#include "gfftdoc.h"
#endif

#include "loki/Typelist.h"
#include "loki/Factory.h"
#include "loki/Singleton.h"

#include "finit.h"
#include "gfftomp.h"
#include "gfftstdalg.h"
#include "gfftpolicy.h"
#include "gfftcaller.h"
#include "gfftgen.h"


#endif /*__gfft_h*/

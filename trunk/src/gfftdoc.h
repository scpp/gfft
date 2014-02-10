/***************************************************************************
 *   Copyright (C) 2009-2014 by Vladimir Mirnyy                            *
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
    \brief General %GFFT documentation
*/

/**
\mainpage

The main feature of GFFT library is that the length of transform 
is defined as a compile-time constant. You may compile FFT of only single length
or easily define and compile a set of transforms, e.g. the powers of 2 from 1 to 25
or a range 10,...,20, etc.

The major advantage is that the first and expensive planning step usually performed 
by the known FFT libraries is fully done at compile-time. Even the twiddle factors 
are partially computed during compilation. Knowing the length of transform, 
the compiler can orginize and optimize the code better. As result,
at run-time GFFT performs the transform in one step avoiding expensive planning.

%GFFT is a header-only library with a couple of sample programs 
and some files from Loki-lib library by Andrei Alexandrescu.
It doesn't need any installation or additional packages.

\section news What's new in the version 0.3

-# Out-of-place transforms. In the version 0.2, the only in-place transforms were possible.
-# Complex transforms of any N. Although, transforms for prime factors are computed as DFT
   without recursion. Therefore, N with big prime factors will be significantly slower.
-# Better accuracy due to compile-time computation of the roots-of-unity as rationals. Their precision 
   is 18 decimal digits.
-# The syntax changed not significantly: GenerateTransform template class receives the list of transform lengths
   as the first parameter instead of the min. and max. powers of 2. The last new parameter should define
   either IN_PLACE or OUT_OF_PLACE
-# Real transforms remain unchanged, only in-place of the size 2^N. They will be reimplemented in the next release
   according to the new architecture of the complex transforms.

\section compilation Compilation

Starting with version 0.2 %GFFT has 
<a href="http://www.cmake.org/" target="_blank">cmake</a>-based build management. 
That means, you can take a look and make changes in CMakeList.txt file,
where all the compiler options and targets are listed.\n
To compile the project, run cmake in %GFFT directory,
then run make:
\verbatim
cmake .
make
\endverbatim

Compilation of %GFFT is compiler-challenging process, because it intensively applies
template class recursion, which must be completely resolved during compilation.
For instance, it is known that gcc 3.x hangs up, when trying to compile %GFFT in optimized mode.
Newer compilers can handle template classes much better and faster.
GFFT 0.3 is succesfully tested with
- GNU gcc 4.7.x, 4.8.x
- MS Visual Studio 2012

Other compilers will be tested soon.

The compilation should take from a few seconds to a few minutes depending on
amount of transforms you would like to compile within your code.


\section usage Basic usage

To start using %GFFT please take a look into simple example programs (gfft*.cpp and cgfft*.cpp).
The first one defines the data as a C-like array, where each even element represents the real and
odd - the imaginary part of a complex number. The second one uses array of std::complex type.

Since all the parameters of %GFFT are static constants and defined as template parameters,
you have to deside before compilation, which kind of transforms of which length you might need
and declare them in instantiation of template class GFFT::GenerateTransform (see its documentation for details).
Its template parameters may be given as a single type or as a typelist with multiple options.

An object of this instantiated template class contains then object factory of all the needed transforms.
Each of them can be obtained from the object factory on demand. Following example declares 
forward and backward complex transform of double precision in single-threaded mode and creates 
a transform object for forward transform of length 2^10:
\code
using namespace GFFT;
typedef TYPELIST_2(DFT, IDFT) ComplexTransforms;
typedef TYPELIST_3(SIntID<8>,SIntID<16>,SIntID<32>) NList;
typedef GenerateTransform<NList, DOUBLE, ComplexTransforms, SIntID<1>, Serial, OUT_OF_PLACE> TransformSet;
TransformSet gfft;
TransformSet::ObjectType* fftobj = gfft.CreateTransformObject(10, DOUBLE::ID, DFT::ID, 1, Serial::ID, OUT_OF_PLACE::ID);
\endcode

If you need only single transform type of fixed length, then you can use directly template class 
GFFT::Transform without object factory.
*/

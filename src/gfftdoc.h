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
    \brief General %GFFT documentation
*/

/**
\mainpage


%GFFT consists from only a few headers, two sample programs (gfft.cpp and cgfft.cpp)
and some files from Loki-lib library by Andrei Alexandrescu.
It doesn't need any installation or additional packages.

\section compilation Compilation

Starting with version 0.2 %GFFT has 
<a href="http://www.cmake.org/" target="_blank">cmake</a>-based build management. 
That means, you can take a look and made changes in CMakeList.txt file,
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
Succesfully tested compilers are:
- GNU gcc 4.x
- Intel C++ 11.x
- MS Visual Studio 8 and 9

The compilation should take from a few seconds to a few minutes depending on
amount of transforms you would like to compile within you code.


\section usage Basic usage

To start using %GFFT please take a look into very simple example programs (gfft.cpp and cgfft.cpp).
The first one defines the data as a C-like array, where every pair represents real and
imaginary part of a complex number. The second one uses array of std::complex type.

Since all the parameters of %GFFT are static constants and defined as template parameters,
you have to deside before compilation, which kind of transforms of which length you might need
and declare them in instantiation of template class GFFT::GenerateTransform (see its documentation for details).
Its template parameters may be given as a single type or as a typelist of certain options of necessary transforms.

An object of this instantiated template class contains then object factory of all the needed transforms.
Each of them can be obtained from the object factory on demand. Following example declares 
forward and backward complex transform of double precision in single-threaded mode and creates 
a transform object for forward transform of length 2^10:
\code
using namespace GFFT;
typedef TYPELIST_2(DFT, IDFT) ComplexTransforms;
typedef GenerateTransform<5, 15, DOUBLE, ComplexTransforms, SIntID<1>, Serial> TransformSet;
TransformSet gfft;
TransformSet::ObjectType* fftobj = gfft.CreateTransformObject(10, DOUBLE::ID, DFT::ID);
\endcode

If you need only single transform type of fixed length, then you can use directly template class 
GFFT::Transform without object factory.
*/

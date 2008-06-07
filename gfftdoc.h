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

/**
\mainpage
    Generic simple and efficient Fast Fourier Transforms (FFT) implementation
    using policy-based design and template metaprogramming. \n
    Features:
    - Transforms in-place
    - Complex FFT of power-of-two length. \n
      See Numerical recipes in C, Chap. 12.2 for theory \n
      http://www.nrbook.com/a/bookcpdf/c12-2.pdf
    - Single real function FFT of power-of-two length \n
      See Numerical recipes in C, Chap. 12.3 for theory and data format \n
      http://www.nrbook.com/a/bookcpdf/c12-3.pdf
    - High and cache-independent performance
    - No additional data is stored. \n
      You can use all available RAM for your transformed data
    - One-step transform. \n
      Many known FFT implementation perform to steps: initialization and transform.
      Initialization for a given length is usually computationally expensive,
      while transform is fast. GFFT needs only to create an object instance that
      includes FFT-algorithm, but no additional data or computation.

\section start Getting started

GFFT consists now from only 5 headers, one sample programm gfft.cpp.
and some files from Loki-lib library by Andrei Alexandrescu.
It doesn't need any installation or additional packages.

\section license Licensing

  GFFT is open source software distributed under terms of GPL.

\section download Download

  Download the source code from project web site at SourceForge \n
  http://sourceforge.net/projects/gfft/download

\section refs References

- Myrnyy, V. A Simple and Efficient FFT Implementation in C++\n
  http://www.ddj.com/cpp/199500857
- Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T.
  Numerical recipes in C, 2th edition, 1992
- Alexandrescu, A. Modern C++ Design. Addison-Wesley, 2001
*/

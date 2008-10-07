//
// C++ Interface: sint
//
// Description:
//
//
// Author: Volodymyr Myrnyy <vm@scientificcpp.com>, (C) 2008
//

#ifndef __sint_h
#define __sint_h

/// \brief Integer number metacontainer.
///        Integer N is wrapped into container class to handle integers and other
///        compile-time number-classes using the same operation classes specializing
///        them for particular number-container.
/// \param N an integer number
/////////////////////////////////////////////////////////////////////////

template<int N>
struct SInt {
   enum { Value = N };
};

#endif

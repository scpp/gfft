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

#define STATIC_INTEGER_CLASS(Type, Name) \
template<Type N>                         \
struct s_##Name {                        \
   typedef Type value_type;              \
   static const Type value = N;          \
};

STATIC_INTEGER_CLASS(int, int)
STATIC_INTEGER_CLASS(unsigned int, uint)
STATIC_INTEGER_CLASS(long, long)
STATIC_INTEGER_CLASS(unsigned long, ulong)

#undef STATIC_INTEGER_CLASS

#endif

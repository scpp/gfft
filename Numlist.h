/*! \file
    \brief
*/
/***************************************************************************
                          Numlist.h  -  description
                             -------------------
    begin                : Fri Feb 13 2004
    copyright            : (C) 2004 by Volodymyr Myrnyy (Vladimir Mirnyi)
    email                : myrnyy@math.tu-cottbus.de
 ***************************************************************************/

#ifdef GCC296
#include <limits.h>
#else
#include <limits>
#endif

#ifndef __numlist_h
#define __numlist_h


#ifndef INT_MAX
#define INT_MAX __INT_MAX__
#endif

#ifndef IntT
#define IntT int
#endif


/// \brief Metadefinition of integer compile-time arrays
/// \ingroup gr_meta
namespace NL {

/// \brief An empty class defines end of a numlist
class NullType {};

/// \brief Integer number metacontainer to build compile-time arrays
/// \param Number an integer element of an array
/// \param T continuation of an array: either a next numlist or NullType
/////////////////////////////////////////////////////////////////////////
template<IntT Number, class T>
struct Numlist {
   enum { value = Number };
   typedef T Tail;
};

/// \class Print
/// \brief prints a numlist in compiler messages
///
/// This metaalgorithm is developed for the debugging purposes only.
/// If you type Print<SomeNumlist>, you can see the whole hierarchy of
/// SomeNumlist in the compiler messages.
/// Specialization of the end of a Numlist generates the error: no type named 'Result' in struct ...
/// The compilation succeeds only if SomeNumlist is empty.
/////////////////////////////////////////////////////////////////////////
template<class NList> struct Print;

template<> struct Print<NullType> { };

template<IntT Num, class Tail>
struct Print<Numlist<Num,Tail> > {
   typedef typename Print<Tail>::Result Result;
};

/// \class Length
/// \brief Computes the length of a numlist
/// \param NList a numlist
/// \return a compile-time constant containing the length of NList, not counting
///     the end terminator (which is NullType by convention) \n
/// Length<NList>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Length;
        template <> struct Length<NullType>
        {
            enum { value = 0 };
        };

        template <IntT T, class U>
        struct Length< Numlist<T, U> >
        {
            enum { value = 1 + Length<U>::value };
        };

/// \class NumAt
/// \brief Finds an integer at a given index in a numlist
/// \param NList a numlist
/// \param index an integer constant
/// \return an integer in position 'index' in NList:
/// NumAt<NList, index>::value \n
/// If you pass an out-of-bounds index, the result is a compile-time error
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index> struct NumAt;

        template <IntT Num, class Tail>
        struct NumAt<Numlist<Num, Tail>, 0>
        {
            enum { value = Num };
        };

        template <IntT Num, class Tail, unsigned int i>
        struct NumAt<Numlist<Num, Tail>, i>
        {
            enum { value = NumAt<Tail, i - 1>::value };
        };

/// \class NumAtNonStrict
/// \brief Finds an integer at a given index in a numlist
/// \param NList a numlist
/// \param index an integer constant
/// \param DefaultNum returning value, if index is out-of-bounds
/// \return the integer in position 'index' in NList, or DefaultNum if index is out-of-bounds \n
/// NumAtNonStrict<NList, index, D>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, IntT DefaultNum = 0>
        struct NumAtNonStrict
        {
            enum { value = DefaultNum };
        };

        template <IntT Num, class Tail, IntT DefaultNum>
        struct NumAtNonStrict<Numlist<Num, Tail>, 0, DefaultNum>
        {
            enum { value = Num };
        };

        template <IntT Num, class Tail, unsigned int i, IntT DefaultNum>
        struct NumAtNonStrict<Numlist<Num, Tail>, i, DefaultNum>
        {
            enum { value = NumAtNonStrict<Tail, i-1, DefaultNum>::value };
        };

/// \class IndexOf
/// \brief Finds the index of an integer in a numlist
/// \param NList a numlist
/// \param Num an integer
/// \return the position of Num in NList, or -1 if Num is not found in NList \n
/// IndexOf<NList, T>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT Num> struct IndexOf;

        template <IntT Num>
        struct IndexOf<NullType, Num>
        {
            enum { value = -1 };
        };

        template <IntT Num, class Tail>
        struct IndexOf<Numlist<Num, Tail>, Num>
        {
            enum { value = 0 };
        };

        template <IntT Head, class Tail, IntT Num>
        struct IndexOf<Numlist<Head, Tail>, Num>
        {
        private:
            enum { temp = IndexOf<Tail, Num>::value };
        public:
            enum { value = (temp == -1 ? -1 : 1 + temp) };
        };

/// \class Append
/// \brief Appends an integer to a numlist
/// \param NList a numlist
/// \param Num an integer to be added
/// \return a numlist that is NList followed by Num and NullType-terminated \n
/// Append<NList, Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT Num> struct Append;

        template <IntT Num> struct Append<NullType, Num>
        {
            typedef Numlist<Num,NullType> Result;
        };

        template <IntT Head, class Tail, IntT Num>
        struct Append<Numlist<Head, Tail>, Num>
        {
            typedef Numlist<Head,
                    typename Append<Tail, Num>::Result> Result;
        };

/// \class AppendList
/// \brief Appends a numlist to another
/// \param NList a numlist
/// \param List another numlist to be added
/// \return a numlist that is NList followed by List and NullType-terminated \n
/// AppendList<NList, List>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, class List> struct AppendList;

        template <> struct AppendList<NullType, NullType>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct AppendList<NullType, Numlist<Head, Tail> >
        {
            typedef Numlist<Head, Tail> Result;
        };

        template <IntT Head, class Tail, IntT Head1, class Tail1>
        struct AppendList<Numlist<Head, Tail>, Numlist<Head1,Tail1> >
        {
            typedef Numlist<Head,
                    typename AppendList<Tail, Numlist<Head1,Tail1> >::Result> Result;
        };

/// \class Erase
/// \brief Erases the first occurence, if any, of an integer in a numlist
/// \param NList a numlist
/// \param Num an integer
/// \return a numlist that is NList without the first occurence of Num \n
/// Erase<NList, Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT Num> struct Erase;

        template <IntT Num>                         // Specialization 1
        struct Erase<NullType, Num>
        {
            typedef NullType Result;
        };

        template <IntT Num, class Tail>             // Specialization 2
        struct Erase<Numlist<Num, Tail>, Num>
        {
            typedef Tail Result;
        };

        template <IntT Head, class Tail, IntT Num> // Specialization 3
        struct Erase<Numlist<Head, Tail>, Num>
        {
            typedef Numlist<Head,
                    typename Erase<Tail, Num>::Result>
                Result;
        };

/// \class EraseAll
/// \brief Erases all occurences, if any, of an integer in a numlist
/// \param NList a numlist
/// \param Num an integer
/// \return a numlist that is NList without any occurence of Num \n
/// EraseAll<NList, Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT Num> struct EraseAll;
        template <IntT Num>
        struct EraseAll<NullType, Num>
        {
            typedef NullType Result;
        };
        template <IntT Num, class Tail>
        struct EraseAll<Numlist<Num, Tail>, Num>
        {
            // Go all the way down the list removing the type
            typedef typename EraseAll<Tail, Num>::Result Result;
        };
        template <IntT Head, class Tail, IntT Num>
        struct EraseAll<Numlist<Head, Tail>, Num>
        {
            // Go all the way down the list removing the type
            typedef Numlist<Head,
                    typename EraseAll<Tail, Num>::Result>
                Result;
        };

/// \class EraseAt
/// \brief Erases an element of a numlist given by the index
/// \param NList a numlist
/// \param index element index (starting with 0)
/// \return a numlist that is NList without an element in the place index \n
/// EraseAt<NList, Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index> struct EraseAt;

        template <unsigned int index>
        struct EraseAt<NullType, index>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct EraseAt<Numlist<Head, Tail>, 0>
        {
            typedef Tail Result;
        };

        template <IntT Head, class Tail, unsigned int index>
        struct EraseAt<Numlist<Head, Tail>, index>
        {
            typedef Numlist<Head,
                    typename EraseAt<Tail, index-1>::Result> Result;
        };

/// \class NoDuplicates
/// \brief Removes all duplicate integers in a numlist
/// \param NList is a numlist
/// \return NoDuplicates<NList>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct NoDuplicates;

        template <> struct NoDuplicates<NullType>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct NoDuplicates< Numlist<Head, Tail> >
        {
        private:
            typedef typename NoDuplicates<Tail>::Result L1;
            typedef typename Erase<L1, Head>::Result L2;
        public:
            typedef Numlist<Head, L2> Result;
        };

/// \class Replace
/// \brief Replaces the first occurence of an integer in a numlist, with another integer
/// \param NList a numlist
/// \param T an integer
/// \param U an integer
/// \return a numlist in which the first occurence of T is replaced with U \n
/// Replace<NList, T, U>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT T, IntT U> struct Replace;

        template <IntT T, IntT U>
        struct Replace<NullType, T, U>
        {
            typedef NullType Result;
        };

        template <IntT T, class Tail, IntT U>
        struct Replace<Numlist<T, Tail>, T, U>
        {
            typedef Numlist<U, Tail> Result;
        };

        template <IntT Head, class Tail, IntT T, IntT U>
        struct Replace<Numlist<Head, Tail>, T, U>
        {
            typedef Numlist<Head,
                    typename Replace<Tail, T, U>::Result>
                Result;
        };

/// \class ReplaceAt
/// \brief Replaces an integer with the given index in a numlist, with another integer
/// \param NList a numlist
/// \param index non-negative integer
/// \param U an integer
/// \return a numlist in which the i-th element (starting with 0) is replaced with U \n
/// ReplaceAt<NList, index, U>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, IntT U> struct ReplaceAt;

        template <unsigned int i, IntT U>
        struct ReplaceAt<NullType, i, U>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail, IntT U>
        struct ReplaceAt<Numlist<Head, Tail>, 0, U>
        {
            typedef Numlist<U,Tail> Result;
        };

        template <IntT Head, class Tail, unsigned int i, IntT U>
        struct ReplaceAt<Numlist<Head, Tail>, i, U>
        {
            typedef Numlist<Head,
                    typename ReplaceAt<Tail, i-1, U>::Result>
                Result;
        };

/// \class ReplaceList
/// \brief Replaces the first occurence of an integer in a numlist, with another numlist
/// \param NList a numlist
/// \param T an integer
/// \param List another numlist
/// \return a numlist in which the first occurence of T is replaced with List \n
/// ReplaceList<NList, T, List>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT T, class List> struct ReplaceList;

        template <IntT T, class List>
        struct ReplaceList<NullType, T, List>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail, IntT T, class List>
        struct ReplaceList<Numlist<Head, Tail>, T, List>
        {
            typedef Numlist<Head,
                    typename ReplaceList<Tail, T, List>::Result>
                Result;
        };

        template <IntT T, class Tail, IntT Head1, class Tail1>
        struct ReplaceList<Numlist<T, Tail>, T, Numlist<Head1,Tail1> >
        {
            typedef Numlist<Head1,
                    typename AppendList<Tail1,Tail>::Result> Result;
        };

/// \class ReplaceAll
/// \brief Replaces all occurences of an integer in a numlist, with another integer
/// \param NList a numlist
/// \param T an integer
/// \param U an integer
/// \return a numlist in which all occurences of T is replaced with U \n
/// ReplaceAll<NList, T, U>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT T, IntT U> struct ReplaceAll;

        template <IntT T, IntT U>
        struct ReplaceAll<NullType, T, U>
        {
            typedef NullType Result;
        };

        template <IntT T, class Tail, IntT U>
        struct ReplaceAll<Numlist<T, Tail>, T, U>
        {
            typedef Numlist<U, typename ReplaceAll<Tail, T, U>::Result> Result;
        };

        template <IntT Head, class Tail, IntT T, IntT U>
        struct ReplaceAll<Numlist<Head, Tail>, T, U>
        {
            typedef Numlist<Head,
                    typename ReplaceAll<Tail, T, U>::Result>
                Result;
        };

/// \class ReplaceListAll
/// \brief Replaces all occurences of an integer in a numlist, with another numlist
/// \param NList a numlist
/// \param T an integer
/// \param List another numlist
/// \return a numlist in which all occurences of T is replaced with List \n
/// ReplaceListAll<NList, T, List>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT T, class List> struct ReplaceListAll;

        template <IntT T, class List>
        struct ReplaceListAll<NullType, T, List>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail, IntT T, class List>
        struct ReplaceListAll<Numlist<Head, Tail>, T, List>
        {
            typedef Numlist<Head,
                    typename ReplaceListAll<Tail, T, List>::Result>
                Result;
        };

        template <IntT T, class Tail, IntT Head1, class Tail1>
        struct ReplaceListAll<Numlist<T, Tail>, T, Numlist<Head1,Tail1> >
        {
            typedef Numlist<Head1,
                    typename ReplaceListAll<
                    typename AppendList<Tail1,Tail>::Result,T,
                             Numlist<Head1,Tail1> >::Result> Result;
        };

/// \class Reverse
/// \brief Reverses a numlist
/// \param NList a numlist
/// \return a numlist that is NList reversed:
/// Reverse<NList>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Reverse;

        template <IntT Num>
        struct Reverse< Numlist<Num,NullType> >
        {
            typedef Numlist<Num,NullType> Result;
        };

        template <IntT Head, class Tail>
        struct Reverse< Numlist<Head, Tail> >
        {
            typedef typename Append<
                typename Reverse<Tail>::Result, Head>::Result Result;
        };

/// \class Range
/// \brief Returns a range of a numlist given by start and end index
/// \param NList a numlist
/// \param Start index of the first element to be included in the range
/// \param End index of the last element to be included in the range
/// \return a numlist that is the (start-end) range of NList:
/// Range<NList,Start,End>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int Start, unsigned int End,
           unsigned int I=0, bool C=((I>=Start) && (I<=End))>
        struct Range;

        template <IntT H, class T, unsigned int Start,
           unsigned int End, unsigned int I>
        struct Range<Numlist<H,T>,Start,End,I,true>
        {
            typedef Numlist<H,typename Range<T,Start,End,I+1>::Result> Result;
        };

        template <IntT H, class T, unsigned int Start,
           unsigned int End, unsigned int I>
        struct Range<Numlist<H,T>,Start,End,I,false>
        {
            typedef typename Range<T,Start,End,I+1>::Result Result;
        };

        template <unsigned int Start,
           unsigned int End, unsigned int I, bool C>
        struct Range<NullType,Start,End,I,C>
        {
            typedef NullType Result;
        };

/// \class Max
/// \brief Calculates the maximum value of a numlist
/// \param NList a numlist
/// \return an integer that is the maximum value of NList \n
/// Max<NList>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Max;
        template <> struct Max<NullType>
        {
            enum { value = -INT_MAX-1 };
        };

        template <IntT Num, class Tail>
        struct Max< Numlist<Num, Tail> >
        {
        private:
            enum { temp = Max<Tail>::value };
        public:
            enum { value = temp > Num ? temp : Num };
        };

/// \class Min
/// \brief Calculates the minimum value of a numlist
/// \param NList a numlist
/// \return an integer that is the minimum value of NList \n
/// Min<NList>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Min;
        template <> struct Min<NullType>
        {
            enum { value = INT_MAX };
        };

        template <IntT Num, class Tail>
        struct Min< Numlist<Num, Tail> >
        {
        private:
            enum { temp = Min<Tail>::value };
        public:
            enum { value = temp < Num ? temp : Num };
        };

/// \class AddConst
/// \brief Adds a constant integer to all elements of a numlist
/// \param NList a numlist
/// \param Num a constant integer
/// \return a numlist which elements increased on Num comparing to NList \n
/// AddConst<NList,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT Num> struct AddConst;
        template <IntT Num>
        struct AddConst<NullType,Num>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail, IntT Num>
        struct AddConst<Numlist<Head,Tail>,Num>
        {
            typedef Numlist<Head+Num,
                    typename AddConst<Tail,Num>::Result> Result;
        };

/// \class AddAt
/// \brief Adds an integer to the element with given index of a numlist
/// \param NList a numlist
/// \param index a non-negative integer
/// \param Num an integer
/// \return a numlist in which the i-th element increased on Num comparing to NList
/// An out-of-bounds index results in a compile-time error \n
/// AddAt<NList,index,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, IntT Num> struct AddAt;

        template <IntT Head, class Tail, IntT Num>
        struct AddAt<Numlist<Head, Tail>, 0, Num>
        {
            typedef Numlist<Head + Num, Tail> Result;
        };

        template <IntT Head, class Tail, unsigned int i, IntT Num>
        struct AddAt<Numlist<Head, Tail>, i, Num>
        {
            typedef Numlist<Head,
                    typename AddAt<Tail, i-1, Num>::Result>
                Result;
        };

/// \class Add
/// \brief Adds two numlists even then if their lengths are different
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return a numlist which elements are sums of correspondent in NList1 and NList2
/// the tail of the longest numlist is appended \n
/// Add<NList1,NList2>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Add;

        template <>
        struct Add<NullType,NullType>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct Add<NullType,Numlist<Head,Tail> >
        {
            typedef Numlist<Head,Tail> Result;
        };

        template <IntT Head, class Tail>
        struct Add<Numlist<Head,Tail>,NullType>
        {
            typedef Numlist<Head,Tail> Result;
        };

        template <IntT Head1, class Tail1, IntT Head2, class Tail2>
        struct Add<Numlist<Head1,Tail1>,Numlist<Head2,Tail2> >
        {
            typedef Numlist<Head1+Head2,
                    typename Add<Tail1,Tail2>::Result> Result;
        };


/// \class SubConst
/// \brief Subtracts a constant integer to all elements of a numlist
/// \param NList a numlist
/// \param Num a constant integer
/// \return a numlist which elements decreased on Num comparing to NList \n
/// SubConst<NList,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT Num> struct SubConst;

        template <IntT Num>
        struct SubConst<NullType,Num>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail, IntT Num>
        struct SubConst<Numlist<Head,Tail>,Num>
        {
            typedef Numlist<Head-Num,
                    typename SubConst<Tail,Num>::Result> Result;
        };

/// \class SubAt
/// \brief Subtracts an integer to the element with given index of a numlist
/// \param NList a numlist
/// \param index a non-negative integer
/// \param Num an integer
/// \return a numlist in which the i-th element is decreased on Num comparing to NList
/// An out-of-bounds index results in a compile-time error \n
/// SubAt<NList,index,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, IntT Num> struct SubAt;

        template <IntT Head, class Tail, IntT Num>
        struct SubAt<Numlist<Head, Tail>, 0, Num>
        {
            typedef Numlist<Head - Num, Tail> Result;
        };

        template <IntT Head, class Tail, unsigned int i, IntT Num>
        struct SubAt<Numlist<Head, Tail>, i, Num>
        {
            typedef Numlist<Head,
                    typename SubAt<Tail, i-1, Num>::Result> Result;
        };

/// \class Sub
/// \brief Subtracts two numlists even then if their lengths are different
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return a numlist which elements are differences of correspondent in NList1 and NList2
/// the tail of the longest numlist is appended \n
/// Sub<NList1,NList2>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Sub;

        template <>
        struct Sub<NullType,NullType>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct Sub<NullType,Numlist<Head,Tail> >
        {
            typedef Numlist<Head,Tail> Result;
        };

        template <IntT Head, class Tail>
        struct Sub<Numlist<Head,Tail>,NullType>
        {
            typedef Numlist<Head,Tail> Result;
        };

        template <IntT Head1, class Tail1, IntT Head2, class Tail2>
        struct Sub<Numlist<Head1,Tail1>,Numlist<Head2,Tail2> >
        {
            typedef Numlist<Head1-Head2,
                    typename Sub<Tail1,Tail2>::Result> Result;
        };

/// \class MultConst
/// \brief Multiplies all elements of a numlist by a constant integer
/// \param NList a numlist
/// \param Num a constant integer
/// \return a numlist which elements multiplied by Num comparing to NList \n
/// MultConst<NList,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, IntT Num> struct MultConst;

        template <IntT Num>
        struct MultConst<NullType,Num>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail, IntT Num>
        struct MultConst<Numlist<Head,Tail>,Num>
        {
            typedef Numlist<Head*Num,
                    typename MultConst<Tail,Num>::Result> Result;
        };

/// \class MultAt
/// \brief Multiplies an element with given index of a numlist by an integer
/// \param NList a numlist
/// \param index a non-negative integer
/// \param Num an integer
/// \return a numlist in which the i-th element multiplied by Num comparing to NList
/// An out-of-bounds index results in a compile-time error \n
/// MultAt<NList,index,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, IntT Num> struct MultAt;

        template <IntT Head, class Tail, IntT Num>
        struct MultAt<Numlist<Head, Tail>, 0, Num>
        {
            typedef Numlist<Head * Num, Tail> Result;
        };

        template <IntT Head, class Tail, unsigned int i, IntT Num>
        struct MultAt<Numlist<Head, Tail>, i, Num>
        {
            typedef Numlist<Head,
                    typename MultAt<Tail, i-1, Num>::Result>
                Result;
        };

/// \class Mult
/// \brief Multiplies element by element of two numlists even with different lengths
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return a numlist which elements are products of correspondent in NList1 and NList2
/// The tail of the longest numlist is truncated \n
/// Mult<NList1,NList2>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Mult;
        template <>
        struct Mult<NullType,NullType>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct Mult<NullType,Numlist<Head,Tail> >
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct Mult<Numlist<Head,Tail>,NullType>
        {
            typedef Numlist<Head,Tail> Result;
        };

        template <IntT Head1, class Tail1, IntT Head2, class Tail2>
        struct Mult<Numlist<Head1,Tail1>,Numlist<Head2,Tail2> >
        {
            typedef Numlist<Head1*Head2,
                    typename Add<Tail1,Tail2>::Result> Result;
        };

/// \class Sum
/// \brief Computes the sum of elements in a numlist
/// \param NList a numlist
/// \return an integer that is the sum of NList \n
/// Sum<NList>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Sum;
        template <> struct Sum<NullType>
        {
            enum { value = 0 };
        };

        template <IntT Num, class Tail>
        struct Sum< Numlist<Num, Tail> >
        {
            enum { value = Num + Sum<Tail>::value };
        };

/// \class Compare
/// \brief Compares elements in the numlist NList1 with NList2
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return Positive value, if the numlist NList1 greater than NList2,
///         negative value otherwise. Returns zero, if NList1 and NList2 are equal.
///         (last element is the most significant) \n
/// Compare<NList1,NList2>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Compare;

        template <IntT H1, class T1, IntT H2, class T2>
        struct Compare<Numlist<H1,T1>,Numlist<H2,T2> >
        {
            enum { v = Compare<T1,T2>::value };
            enum { value = (v==0) ? (H1-H2) : v };
        };

        template <IntT H, class T>
        struct Compare<Numlist<H,T>,NullType>
        {
            enum { value = 1 };
        };

        template <IntT H, class T>
        struct Compare<NullType,Numlist<H,T> >
        {
            enum { value = -1 };
        };

        template <IntT H1, IntT H2>
        struct Compare<Numlist<H1,NullType>,Numlist<H2,NullType> >
        {
            enum { value = (H1-H2) };
        };

/// \class Sort
/// \brief Sorts a numlist applying consequent the operation Min
/// \param NList a numlist
/// \return a numlist which is sorted NList \n
/// Sort<NList>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Sort;
        template <> struct Sort<NullType>
        {
            typedef NullType Result;
        };

        template <IntT Head, class Tail>
        struct Sort<Numlist<Head,Tail> >
        {
        private:
            enum { min = Min<Numlist<Head,Tail> >::value };
            typedef typename Replace<Tail,min,Head>::Result temp;

        public:
            typedef Numlist<min,typename Sort<temp>::Result> Result;
        };

/// \class InitNumlist
/// \brief Generates a numlist of a given length with default elements
/// \param NList a numlist
/// \param Len a positive length of the resulted numlist
/// \param Value a default numlist element value
/// \return a numlist which contains Len numbers Value \n
/// InitNumlist<Len,Value>::Result
////////////////////////////////////////////////////////////////////////////////

        template<unsigned int Len, IntT Value=0>
        struct InitNumlist
        {
           typedef Numlist<Value, typename InitNumlist<Len-1,Value>::Result> Result;
        };

        template<IntT Value>
        struct InitNumlist<0,Value>
        {
           typedef NullType Result;
        };

/// \class FillNumlist
/// \brief Fills a numlist up to a given length with default elements
/// \param NList a numlist
/// \param Len a positive length of the resulted numlist
/// \param Value a default numlist element value
/// \return a numlist which contains Len numbers Value \n
/// InitNumlist<Len,Value>::Result
////////////////////////////////////////////////////////////////////////////////

        template<class NList, unsigned int Len, IntT Value=0,
                 bool C=(Len>Length<NList>::value)>
        struct FillNumlist;

        template<class NList, unsigned int Len, IntT Value>
        struct FillNumlist<NList,Len,Value,true>
        {
        private:
           typedef typename InitNumlist<Len-Length<NList>::value,Value>::Result VList;
        public:
           typedef typename AppendList<NList,VList>::Result Result;
        };

        template<class NList, unsigned int Len, IntT Value>
        struct FillNumlist<NList,Len,Value,false>
        {
           typedef NList Result;
        };

/// \class ReturnNumList
/// \brief Returns a numlist to an array of integers
/// \param NList a numlist to be assigned to an array
/// \param array the resulted array of integers
/// \return an integer array which contains elements of NList \n
/// ReturnNumList<NList>::apply(array)
////////////////////////////////////////////////////////////////////////////////
template<class NList> struct ReturnNumList;

template<IntT N, class Tail>
struct ReturnNumList<Numlist<N,Tail> >
{
   inline static void apply(IntT* array)
   {
      *array = N;
      ReturnNumList<Tail>::apply(array+1);
   }
};

template<>
struct ReturnNumList<NullType>
{
   inline static void apply(IntT*) { }
};


}  // namespace


// macros NUMLIST_1, NUMLIST_2, ... NUMLIST_50
// Each takes a number of arguments equal to its numeric suffix
// The arguments are type names. NUMLIST_NN generates a numlist containing
//     all types passed as arguments, in that order.
// Example: NUMLIST_2(-5, 2) generates an array containing -5 and 2.

#define NUMLIST_1(T1) ::NL::Numlist<T1, ::NL::NullType>

#define NUMLIST_2(T1, T2) ::NL::Numlist<T1, NUMLIST_1(T2) >

#define NUMLIST_3(T1, T2, T3) ::NL::Numlist<T1, NUMLIST_2(T2, T3) >

#define NUMLIST_4(T1, T2, T3, T4) \
    ::NL::Numlist<T1, NUMLIST_3(T2, T3, T4) >

#define NUMLIST_5(T1, T2, T3, T4, T5) \
    ::NL::Numlist<T1, NUMLIST_4(T2, T3, T4, T5) >

#define NUMLIST_6(T1, T2, T3, T4, T5, T6) \
    ::NL::Numlist<T1, NUMLIST_5(T2, T3, T4, T5, T6) >

#define NUMLIST_7(T1, T2, T3, T4, T5, T6, T7) \
    ::NL::Numlist<T1, NUMLIST_6(T2, T3, T4, T5, T6, T7) >

#define NUMLIST_8(T1, T2, T3, T4, T5, T6, T7, T8) \
    ::NL::Numlist<T1, NUMLIST_7(T2, T3, T4, T5, T6, T7, T8) >

#define NUMLIST_9(T1, T2, T3, T4, T5, T6, T7, T8, T9) \
    ::NL::Numlist<T1, NUMLIST_8(T2, T3, T4, T5, T6, T7, T8, T9) >

#define NUMLIST_10(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10) \
    ::NL::Numlist<T1, NUMLIST_9(T2, T3, T4, T5, T6, T7, T8, T9, T10) >

#define NUMLIST_11(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11) \
    ::NL::Numlist<T1, NUMLIST_10(T2, T3, T4, T5, T6, T7, T8, T9, T10, T11) >

#define NUMLIST_12(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12) \
    ::NL::Numlist<T1, NUMLIST_11(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12) >

#define NUMLIST_13(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13) \
    ::NL::Numlist<T1, NUMLIST_12(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13) >

#define NUMLIST_14(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14) \
    ::NL::Numlist<T1, NUMLIST_13(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14) >

#define NUMLIST_15(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15) \
    ::NL::Numlist<T1, NUMLIST_14(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15) >

#define NUMLIST_16(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16) \
    ::NL::Numlist<T1, NUMLIST_15(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16) >

#define NUMLIST_17(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17) \
    ::NL::Numlist<T1, NUMLIST_16(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17) >

#define NUMLIST_18(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18) \
    ::NL::Numlist<T1, NUMLIST_17(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18) >

#define NUMLIST_19(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19) \
    ::NL::Numlist<T1, NUMLIST_18(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19) >

#define NUMLIST_20(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20) \
    ::NL::Numlist<T1, NUMLIST_19(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20) >

#define NUMLIST_21(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21) \
    ::NL::Numlist<T1, NUMLIST_20(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21) >

#define NUMLIST_22(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22) \
    ::NL::Numlist<T1, NUMLIST_21(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22) >

#define NUMLIST_23(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23) \
    ::NL::Numlist<T1, NUMLIST_22(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23) >

#define NUMLIST_24(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24) \
    ::NL::Numlist<T1, NUMLIST_23(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24) >

#define NUMLIST_25(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25) \
    ::NL::Numlist<T1, NUMLIST_24(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25) >

#define NUMLIST_26(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26) \
    ::NL::Numlist<T1, NUMLIST_25(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26) >

#define NUMLIST_27(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27) \
    ::NL::Numlist<T1, NUMLIST_26(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27) >

#define NUMLIST_28(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28) \
    ::NL::Numlist<T1, NUMLIST_27(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28) >

#define NUMLIST_29(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29) \
    ::NL::Numlist<T1, NUMLIST_28(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29) >

#define NUMLIST_30(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30) \
    ::NL::Numlist<T1, NUMLIST_29(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30) >

#define NUMLIST_31(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31) \
    ::NL::Numlist<T1, NUMLIST_30(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31) >

#define NUMLIST_32(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32) \
    ::NL::Numlist<T1, NUMLIST_31(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32) >

#define NUMLIST_33(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33) \
    ::NL::Numlist<T1, NUMLIST_32(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33) >

#define NUMLIST_34(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34) \
    ::NL::Numlist<T1, NUMLIST_33(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34) >

#define NUMLIST_35(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35) \
    ::NL::Numlist<T1, NUMLIST_34(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35) >

#define NUMLIST_36(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36) \
    ::NL::Numlist<T1, NUMLIST_35(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36) >

#define NUMLIST_37(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37) \
    ::NL::Numlist<T1, NUMLIST_36(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37) >

#define NUMLIST_38(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38) \
    ::NL::Numlist<T1, NUMLIST_37(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38) >

#define NUMLIST_39(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39) \
    ::NL::Numlist<T1, NUMLIST_38(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39) >

#define NUMLIST_40(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40) \
    ::NL::Numlist<T1, NUMLIST_39(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40) >

#define NUMLIST_41(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41) \
    ::NL::Numlist<T1, NUMLIST_40(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41) >

#define NUMLIST_42(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41, T42) \
    ::NL::Numlist<T1, NUMLIST_41(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41, T42) >

#define NUMLIST_43(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41, T42, T43) \
    ::NL::Numlist<T1, NUMLIST_42(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41, T42, T43) >

#define NUMLIST_44(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41, T42, T43, T44) \
    ::NL::Numlist<T1, NUMLIST_43(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, T41, T42, T43, T44) >

#define NUMLIST_45(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45) \
    ::NL::Numlist<T1, NUMLIST_44(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45) >

#define NUMLIST_46(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46) \
    ::NL::Numlist<T1, NUMLIST_45(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45) >

#define NUMLIST_47(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47) \
    ::NL::Numlist<T1, NUMLIST_46(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47) >

#define NUMLIST_48(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47, T48) \
    ::NL::Numlist<T1, NUMLIST_47(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47, T48) >

#define NUMLIST_49(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47, T48, T49) \
    ::NL::Numlist<T1, NUMLIST_48(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47, T48, T49) >

#define NUMLIST_50(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47, T48, T49, T50) \
    ::NL::Numlist<T1, NUMLIST_49(T2, T3, T4, T5, T6, T7, T8, T9, T10, \
        T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, \
        T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, \
        T31, T32, T33, T34, T35, T36, T37, T38, T39, T40, \
        T41, T42, T43, T44, T45, T46, T47, T48, T49, T50) >


#endif /*__numlist_h*/

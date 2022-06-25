#define assert(x) { if (!(x)) dbgPoint(); }
#include "MandelMath.hpp"
#include <math.h>

//#include <cassert>
#include <cmath>
//#define assert(x) { }
#define assert(x) { if (!(x)) dbgPoint(); }

void doNothing(int &x)
{
  x++;
}

void dbgPoint()
{
  int x=3;
  doNothing(x);
}

template <typename Tout, typename Tin> //typechecked hard cast
Tout specific_cast(Tin in) { return (Tout)in; }


namespace MandelMath {

int gcd(int m, int n)
{
  if (m==n)
    return m;
  if (m==0)
    return n;
  if (n==0)
    return m;
  int c=0;
  while (((m|n)&1)==0)
  {
    m>>=1; n>>=1; c++;
  }

  /* Dividing n by 2 until a becomes odd */
  while ((n & 1) == 0)
    n >>= 1;

  /* From here on, 'n' is always odd. */
  do
  {
    /* If m is even, remove all factor of 2 in m */
    while ((m & 1) == 0)
      m >>= 1;

    /* Now m and n are both odd.
       Swap if necessary so n <= m,
       then set m = m - n (which is even).*/
    if (n > m)
    { int x=n; n=m; m=x; }

    m = (m - n)>>1;
  } while (m != 0);

  /* restore common factors of 2 */
  return n << c;
}

int ctz16(int x)
{
  int ctzidx=(((0xF65*(x&-x))&0x7800)>>10);
  int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);
  return ctz1;
  /* shift, discard top bit and append i, rotate to be max; is i followed by 0 only? -> negate else not negate the discarded bit -> append
      i is considered 1-epsilon (1 except when a tie with another 1)
  ctz 8 bit
  00011101
  000-00i-i00-y   1
   001-01i-1i0-y   1
    011-11i-11i-y   1
     111-11i-11i-y   0
      110-10i-i10-n   1
       101-01i-1i0-y   0
        010-10i-i10-n   0
         100-00i-i00-y   0
          000
   0x1D * 0=0x0000 >>4&7=0
   0x1D * 1=0x001D >>4&7=1
   0x1D * 2=0x003A >>4&7=3   MAGIC[3]:=ctz(2)
   0x1D * 4=0x0074 >>4&7=7   MAGIC[7]:=ctz(4)
   0x1D * 8=0x00E8 >>4&7=6
   0x1D *16=0x01D0 >>4&7=5
   0x1D *32=0x03A0 >>4&7=2
   0x1D *64=0x0740 >>4&7=4
   0x1D*128=0x0E80 >>4&7=0   MAGIC[0]=ctz(128)
   0x1D*256=0x1D00 >>4&7=0
   (0x23461507 >> (((0x1D*(i&-i))&0x70)>>2))&0x07 ?= ctz(0..255)

   ctz 16 bit
   0000-000i-i000-y 1
    0001-001i-1i00-y 1
     0011-011i-11i0-y 1
      0111-111i-111i-y 1
       1111-111i-111i-y 0
        1110-110i-i110-n 1
         1101-101i-1i10-n 1
          1011-011i-11i0-y 0
           0110-110i-i110-n 0
            1100-100i-i100-n 1
             1001-001i-1i00-y 0
              0010-010i-(10i0 or i010)-(y or n)-(1 or 0)  (10i0 correct)
               0101-101i-1i10-n 0
                1010-010i-10i0 or i010-y or n-0 or 1 (10i0 correct)
                 0100-100i-i100-n 0
                  1000-000i-i000-y 0
                   0000
   0000111101100101 = 0xF65
   (0x... >> (((0xF65*(i&-i))&0x7800)>>10))&0x0F ?= ctz(0..255)
         FEDCBA9876543210
   MAGIC=34586C9E27BD1A0F
   MAGIC23=0x59EC6C8C; //00 0101 1001 1110 1100 0110 1100 1000 1100
   MAGIC01=0xC486BD63; // 1100 0100 1000 0110 1011 1101 0110 0011
    int ctzidx=(((0xF65*(i&-i))&0x7800)>>10);
    int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);

   ctz 32 bit
   00000-0000i-i0000-y 1
    00001-0001i-1i000-y 1                                        0
     00011-0011i-11i00-y 1                                       1
      00111-0111i-111i0-y 1                                      2
       01111-1111i-1111i-y 1                                     3
        11111-1111i-1111i-y 0                                    4
         11110-1110i-i1110-n 1                                   5
          11101-1101i-1i110-n 1                                  6
           11011-1011i-11i10-n 1                                 7
            10111-0111i-111i0-y 0                                8
             01110-1110i-i1110-n 0                               9
              11100-1100i-i1100-n 1                              a
               11001-1001i-1i100-n 1                             b
                10011-0011i-11i00-y 0                            c
                 00110-0110i-110i0-y 1
                  01101-1101i-1i110-n 0
                   11010-1010i-i1010-n 1
                    10101-0101i-1i010-n 1
                     01011-1011i-11i10-n 0
                      10110-0110i-110i0-y 0
                       01100-1100i-i1100-n 0
                        11000-1000i-i1000-n 1
                         10001-0001i-1i000-y 0
                          00010-0010i-10i00-y 1
                           00101-0101i-1i010-n 0
                            01010-1010i-i1010-n 0
                             10100-0100i-i0100-n 1
                              01001-1001i-1i100-n 0
                               10010-0010i-10i00-y 0
                                00100-0100i-i0100-n 0
                                 01000-1000i-i1000-n 0
                                  10000-0000i-i0000-y 0
                                   00000
   00000111110111001101011000101001 = 0000 0111 1101 1100 1101 0110 0010 1001 = 0x7DCD629
         1f 1e 1d 1c 1b 1a 19 18 17 16 15 14 13 12 11 10 0f 0e 0d 0c 0b 0a 09 08 07 06 05 04 03 02 01 00
   MAGIC= 4  5  6  a  7  f  b 14  8 12 10 19  c 1b 15 1e  3  9  e 13 11 18 1a 1d  2  d 17 1c  1 16  0 1f
   MAGIC4=00000001011101110001111100110101 = 0000 0001 0111 0111 0001 1111 0011 0101 0000 = 0x1771f350
   MAGIC3=00010110100111010110011101010001 = 000 1011 0100 1110 1011 0011 1010 1000 1000 =  0xb4eb3a88
   MAGIC2=11101101000010110010000101110101 = .. 1110 1101 0000 1011 0010 0001 0111 0101 = 0xed0b2175; 0x3b42c85d4 won't fit
   MAGIC1=00111110010001011011001010100101 = 0 0111 1100 1000 1011 0110 0101 0100 1010 = 0x7c8b654a
   MAGIC0=01001110000101101101100101101001 = 0100 1010 0001 0110 1101 1001 0110 1001 = 0x4e16d969

    int ctzidx=(((0x7DCD629*(i&-i))&0x7C000000)>>26);
    int ctz1=((0x1771f350 >> ctzidx)&0x10) |
             ((0xb4eb3a88 >> ctzidx)&0x08) |
             ((0xed0b2175 >> ctzidx)&1)<<2 |
             ((0x7c8b654a >> ctzidx)&0x02) |
             ((0x4e16d969 >> ctzidx)&0x01);
   */

  /*for (int i=0; i<256; i++) //0 returns 7
  {
    int ctz1=(0x23461507 >> (((0x1D*(i&-i))&0x70)>>2))&0x07;
    int ctz2=0;
    if (i&1) ctz2=0;
    else if (i&0x03) ctz2=1;
    else if (i&0x07) ctz2=2;
    else if (i&0x0F) ctz2=3;
    else if (i&0x1F) ctz2=4;
    else if (i&0x3F) ctz2=5;
    else if (i&0x7F) ctz2=6;
    else if (i&0xFF) ctz2=7;
    else ctz2=7;
    if (ctz2!=ctz1)
      dbgPoint();
  }

  for (int i=0; i<65536; i++) //0 returns 15
  {
    int ctzidx=(((0xF65*(i&-i))&0x7800)>>10);
    int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);
    int ctz2=0;
    if (i&1) ctz2=0;
    else if (i&0x03) ctz2=1;
    else if (i&0x07) ctz2=2;
    else if (i&0x0F) ctz2=3;
    else if (i&0x1F) ctz2=4;
    else if (i&0x3F) ctz2=5;
    else if (i&0x7F) ctz2=6;
    else if (i&0xFF) ctz2=7;
    else if (i&0x1FF) ctz2=8;
    else if (i&0x3FF) ctz2=9;
    else if (i&0x7FF) ctz2=10;
    else if (i&0xFFF) ctz2=11;
    else if (i&0x1FFF) ctz2=12;
    else if (i&0x3FFF) ctz2=13;
    else if (i&0x7FFF) ctz2=14;
    else if (i&0xFFFF) ctz2=15;
    else ctz2=15;
    if (ctz2!=ctz1)
      dbgPoint();
  }

  for (int ii=0; ii<32; ii++)
  {
    unsigned int i=(1u<<ii);
    int ctzidx=(((0x7DCD629*(i&-i))&0x7C000000)>>26);
    int ctz1=((0x1771f350 >> ctzidx)&0x10) |
             ((0xb4eb3a88 >> ctzidx)&0x08) |
             ((0xed0b2175 >> ctzidx)&1)<<2 |
             ((0x7c8b654a >> ctzidx)&0x02) |
             ((0x4e16d969 >> ctzidx)&0x01);
    if (ii==32)
    {
      if (ctz1!=31)
        dbgPoint();
    }
    else if (ctz1!=ii)
      dbgPoint();
  }*/
}


/*double two_pow_n(unsigned int n)
{
  assert(n<4096);
  double result=1;
  while (n>=31)
  {
    result*=(1u<<31);
    n-=31;
  }
  if (n>0)
    result*=(1u<<n);
  return result;
}*/

worker_multi::Allocator::Allocator(int capacity):
  lowest(0), allocUpTo(0), capacity(capacity), above(nullptr), worker(nullptr)
{
}

worker_multi::Allocator::Allocator(worker_multi *worker, int capacity):
  lowest(0), allocUpTo(0), capacity(capacity), above(nullptr), worker(worker)
{
}

worker_multi::Allocator::Allocator(Allocator *allo, int capacity):
  lowest(allo->allocUpTo), allocUpTo(allo->allocUpTo), capacity(capacity), above(allo), worker(allo->worker)
{
  allo->allocUpTo+=capacity;
  assert(allo->allocUpTo<=allo->lowest+allo->capacity);
}

worker_multi::Allocator::Allocator(Allocator *allo, int lowest, int count):
  lowest(allo->lowest+lowest), allocUpTo(allo->lowest+lowest), capacity(count), above(nullptr), worker(allo->worker)
{
  assert(lowest>=allo->lowest);
  assert(lowest+capacity<=allo->allocUpTo);
}

worker_multi::Allocator::~Allocator()
{
   assert(allocUpTo==lowest);
   if (above!=nullptr)
   {
     assert(above->allocUpTo==lowest+capacity);
     assert(above->allocUpTo-above->lowest>=capacity);
     above->allocUpTo-=capacity;
   };
}

number_pointer worker_multi::Allocator::alloc()
{
  if (worker==nullptr) //support for Qt that requires parameter-less constructors for ShareableViewInfo
    return number_pointer();
  assert(allocUpTo<lowest+capacity);
  allocUpTo++;
  return worker->getNumber(allocUpTo-1);
}

void worker_multi::Allocator::dealloc(number_pointer ptr)
{
  if (worker==nullptr) //support for Qt that requires parameter-less constructors for ShareableViewInfo
    return;
  assert(allocUpTo>lowest);
  number_pointer expected=worker->getNumber(allocUpTo-1);
  allocUpTo--;
  assert(expected.asf64==ptr.asf64);
}



#if 0
number_store::number_store(): dbgType(number_worker::Type::typeEmpty)
{
  static_assert(sizeof(dd_real)==2*sizeof(double), "dd_real suspiciously big");
  as.doubl=0.0;
}

number_store::~number_store()
{
  assert(dbgType==number_worker::Type::typeEmpty);
}

void number_store::cleanup(number_worker::Type ntype)
{
  switch (ntype)
  {
#if NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDouble:
      assert(dbgType==number_worker::Type::typeDouble);
      dbgType=number_worker::Type::typeEmpty;
      as.doubl=0.0;
      break;
#endif //NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDDouble:
      assert(dbgType==number_worker::Type::typeDDouble);
      dbgType=number_worker::Type::typeEmpty;
      as.ddouble_.deinit_();
      break;
    case number_worker::Type::typeMulti:
      assert(dbgType==number_worker::Type::typeMulti);
      dbgType=number_worker::Type::typeEmpty;
      as.multi_.deinit_();
      break;
    case number_worker::Type::typeEmpty:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeEmpty;
      break;
  }
}

void number_store::init_(number_worker::Type ntype, void *placement, double val)
{
  switch (ntype)
  {
#if NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDouble:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeDouble;
      as.doubl=val;
      break;
#endif //NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDDouble:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeDDouble;
      as.ddouble_.init_((dd_real *)placement);
      as.ddouble_.dd->hi=val;
      as.ddouble_.dd->lo_=0.0;
      break;
    case number_worker::Type::typeMulti:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeMulti;
      as.multi_.init_((multiprec *)placement);
      as.multi_.bytes->set(val);
      break;
    case number_worker::Type::typeEmpty:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeEmpty;
      break;
  }
}

void number_store::zero(number_worker::Type ntype, double val)
{
  switch (ntype)
  {
#if NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDouble:
      assert(dbgType==number_worker::Type::typeDouble);
      as.doubl=val;
      break;
#endif //NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDDouble:
      assert(dbgType==number_worker::Type::typeDDouble);
      as.ddouble_.dd->hi=val;
      as.ddouble_.dd->lo_=0.0;
      break;
    case number_worker::Type::typeMulti:
      assert(dbgType==number_worker::Type::typeMulti);
      as.multi_.bytes->set(val);
      break;
    case number_worker::Type::typeEmpty: ;
  }
}

void number_store::promote_(number_worker::Type oldType, number_worker::Type newType, void *placement, const number_store *src)
{
  if (src==nullptr)
    src=this;
  assert(src->dbgType==oldType);
  switch (oldType)
  {
    case number_worker::Type::typeEmpty:
    {
      //done by init() assert(dbgType==number_worker::Type::typeEmpty);
      switch (newType)
      {
        case number_worker::Type::typeEmpty:
        {
          init_(newType, placement);
        } break;
    #if NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDouble:
          init_(newType, placement);
          break;
    #endif //NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDDouble:
          init_(newType, placement);
          break;
        case number_worker::Type::typeMulti:
          init_(newType, placement);
          break;
      }
    } break;
#if NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDouble:
      assert(dbgType==number_worker::Type::typeDouble);
      switch (newType)
      {
        case number_worker::Type::typeEmpty:
        {
          cleanup(oldType);
        } break;
    #if NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDouble:
        {
          double val=src->as.doubl;
          as.doubl=val;
        } break;
    #endif //NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDDouble:
        {
          double val=src->as.doubl;
          cleanup(oldType);
          init_(newType, placement, val);
        } break;
        case number_worker::Type::typeMulti:
        {
          double val=src->as.doubl;
          cleanup(oldType);
          init_(newType, placement, val);
        } break;
      }
      break;
#endif //NUMBER_DOUBLE_EXISTS
    case number_worker::Type::typeDDouble:
      assert(dbgType==number_worker::Type::typeDDouble);
      switch (newType)
      {
        case number_worker::Type::typeEmpty:
        {
          cleanup(oldType);
        } break;
    #if NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDouble:
        {
          double val=src->as.ddouble_.dd->hi;
          cleanup(oldType);
          init_(newType, placement, val);
        } break;
    #endif //NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDDouble:
          as.ddouble_.dd->assign(*src->as.ddouble_.dd);
          break;
        case number_worker::Type::typeMulti:
        {
          double val=src->as.ddouble_.dd->hi; //``` dd->multi not implemented
          cleanup(oldType);
          init_(newType, placement, val);
        } break;
      }
      break;
    case number_worker::Type::typeMulti:
      assert(dbgType==number_worker::Type::typeMulti);
      switch (newType)
      {
        case number_worker::Type::typeEmpty:
        {
          cleanup(oldType);
        } break;
    #if NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDouble:
        {
          double val=src->as.multi_.bytes->toDouble();
          cleanup(oldType);
          init_(newType, placement, val);
        } break;
    #endif //NUMBER_DOUBLE_EXISTS
        case number_worker::Type::typeDDouble:
        {
          double val=src->as.multi_.bytes->toDouble(); //``` multi->dd not implemented
          cleanup(oldType);
          init_(newType, placement, val);
        } break;
        case number_worker::Type::typeMulti:
          as.multi_.bytes->assign(*src->as.multi_.bytes);
          break;
      }
      break;
  }
}

/*
template <class T>
void number_store::assign(const number_store &other)
{
  assert(other.dbgType==number_to_type<T>::ntype);
  assert(dbgType==number_to_type<T>::ntype);
  //dbgType=number::Type::typeDouble;
  as.doubl=other.as.doubl;
}

template <>
void number_store::assign<number_double>(const number_store &other)
{
  assert(other.dbgType==number::Type::typeDouble);
  assert(dbgType==number::Type::typeDouble);
  //dbgType=number::Type::typeDouble;
  as.doubl=other.as.doubl;
}

template <>
void number_store::assign<number_ddouble>(const number_store &other)
{
  assert(other.dbgType==number::Type::typeDDouble);
  assert(dbgType==number::Type::typeDDouble);
  //dbgType=number::Type::typeDDouble;
  //as.ddouble_.dd->assign(*other.as.ddouble_.dd);
  *as.ddouble_.dd=*other.as.ddouble_.dd;
}

template <>
void number_store::assign<number_multi>(const number_store &other)
{
  assert(other.dbgType==number::Type::typeMulti);
  assert(dbgType==number::Type::typeMulti);
  //dbgType=number::Type::typeMulti;
  *as.multi_.bytes=*other.as.multi_.bytes;
}

void number_store::assignTo_double(number_store &other)
{
  assert(other.dbgType==number::Type::typeDouble);
  assert(dbgType==number::Type::typeDouble);
  other.as.doubl=as.doubl;
}

void number_store::assignTo_ddouble(number_store &other)
{
  assert(other.dbgType==number::Type::typeDDouble);
  assert(dbgType==number::Type::typeDDouble);
  //dbgType=number::Type::typeDDouble;
  //as.ddouble_.dd->assign(*other.as.ddouble_.dd);
  *other.as.ddouble_.dd=*as.ddouble_.dd;
}

void number_store::assignTo_multi(number_store &other)
{
  assert(other.dbgType==number::Type::typeMulti);
  assert(dbgType==number::Type::typeMulti);
  //dbgType=number::Type::typeMulti;
  *other.as.multi_.bytes=*as.multi_.bytes;
}
*/
#endif

worker_multi *worker_multi::create(Type ntype, int capacity)
{
  switch (ntype)
  {
    default://case Type::typeEmpty:
      dbgPoint();
      goto lolwut;
#if NUMBER_DOUBLE_EXISTS
    case Type::typeDouble: lolwut:
      return new worker_multi_double(capacity);
#endif
#if !ONLY_DOUBLE_WORKER
    case Type::typeFloat128:
      return new worker_multi_float128(capacity);
    case Type::typeDDouble:
      return new worker_multi_ddouble(capacity);
    case Type::typeQDouble:
      return new worker_multi_qdouble(capacity);
    case Type::typeReal642:
      return new worker_multi_real642(capacity);
#endif
  }
}

worker_multi *worker_multi::create(Type ntype, Allocator *allocator)
{
  switch (ntype)
  {
    default://case Type::typeEmpty:
      dbgPoint();
      goto lolwut;
#if NUMBER_DOUBLE_EXISTS
    case Type::typeDouble: lolwut:
      return new worker_multi_double(allocator);
#endif
#if !ONLY_DOUBLE_WORKER
    case Type::typeFloat128:
      return new worker_multi_float128(allocator);
    case Type::typeDDouble:
      return new worker_multi_ddouble(allocator);
    case Type::typeQDouble:
      return new worker_multi_qdouble(allocator);
    case Type::typeReal642:
      return new worker_multi_real642(allocator);
#endif
  }
}


#if NUMBER_DOUBLE_EXISTS
//worker_multi_double::worker_multi_double(worker_multi *source):
worker_multi_double::worker_multi_double(Allocator *source):
  worker_multi(Type::typeDouble, specific_cast<worker_multi_double *, worker_multi *>(source->worker)->capacity)
{
  int oldfirst, oldlast;
  source->_getRange(oldfirst, oldlast);
  storage=new double[oldlast-oldfirst];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage[i]=0;
    } break;
    case Type::typeDouble:
    {
      double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->storage;
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
        storage[i]=src_storage[oldfirst+i];
    } break;
#if !ONLY_DOUBLE_WORKER
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i]=src_storage[oldfirst+i];
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++) //segfault
      {
        storage[i]=src_storage[oldfirst+i].hi;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i]=src_storage[oldfirst+i].hi;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i]=src_storage[oldfirst+i].val1;
      }
    } break;
#endif
  }
}

worker_multi_double::~worker_multi_double()
{
  delete[] storage;
  storage=nullptr;
}

number_pointer worker_multi_double::getNumber(int index)
{
  assert(allocator.checkIndex(index));
  assert(index>=0 && index<capacity);
  return number_pointer(&storage[index]);
}

void worker_multi_double::getTmp12(number_pointer &t1, number_pointer &t2)
{
  t1.asf64=&tmp1;
  t2.asf64=&tmp2;
}

void worker_multi_double::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  assert(src_worker->ntype()==Type::typeDouble);
  worker_multi_double *access=(worker_multi_double *)src_worker;
  assert(allocator.checkIndex(dst));
  assert(allocator.checkIndex(dst+len-1));
  assert(access->allocator.checkIndex(src));
  assert(access->allocator.checkIndex(src+len-1));
  assert(dst>=0 && dst+len<=capacity);
  assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage[dst+i]=access->storage[src+i];
}

void worker_multi_double::zero(const number_pointer store, double val)
{
  *store.asf64=val;
}

void worker_multi_double::assign(const number_pointer store, const number_pointer_c src)
{
  *store.asf64=*src.asf64;
}

void worker_multi_double::chs(const number_pointer store)
{
  *store.asf64 = -*store.asf64;
}

void worker_multi_double::lshift(const number_pointer store, int shoft)
{
  *store.asf64=ldexp(*store.asf64, shoft);
}

void worker_multi_double::round(const number_pointer store)
{
  *store.asf64=std::round(*store.asf64);
}

void worker_multi_double::frac(const number_pointer store)
{
  if (*store.asf64<0)
    *store.asf64 -= ceil(*store.asf64);
  else
    *store.asf64 -= floor(*store.asf64);
}

void worker_multi_double::mod1(const number_pointer store)
{
  *store.asf64 -= floor(*store.asf64);
}

void worker_multi_double::add_double(const number_pointer store, double x)
{
  *store.asf64 += x;
}

void worker_multi_double::mul_double(const number_pointer store, double x)
{
  *store.asf64 *= x;
}

void worker_multi_double::add(const number_pointer store, const number_pointer_c other)
{
  *store.asf64 += *other.asf64;
}

void worker_multi_double::sub(const number_pointer store, const number_pointer_c other)
{
  *store.asf64 -= *other.asf64;
}

void worker_multi_double::rsub(const number_pointer store, const number_pointer_c other)
{
  *store.asf64 = *other.asf64-*store.asf64;
}

void worker_multi_double::mul(const number_pointer store, const number_pointer_c other)
{
  *store.asf64 *= *other.asf64;
}

void worker_multi_double::sqr(const number_pointer store)
{
  *store.asf64 *= *store.asf64;
}

double worker_multi_double::radixfloor(const number_pointer_c store1, const number_pointer_c store2)
{
  int ilog1=std::ilogb(*store1.asf64);
  int ilog2=std::ilogb(*store2.asf64);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

void worker_multi_double::recip(const number_pointer store)
{
  *store.asf64 = 1/(*store.asf64);
}

void worker_multi_double::sqrt(const number_pointer store)
{
  *store.asf64 = std::sqrt(*store.asf64);
}

int worker_multi_double::compare(const number_pointer_c store, const number_pointer_c other)
{
  if (*store.asf64 == *other.asf64)
    return 0;
  else if (*store.asf64 < *other.asf64)
    return -1;
  else
    return +1;
}

bool worker_multi_double::isequal(const number_pointer_c store, const number_pointer_c other)
{
  return *store.asf64 == *other.asf64;
}

bool worker_multi_double::is0(const number_pointer_c store)
{
  return *store.asf64 == 0;
}

bool worker_multi_double::isle(const number_pointer_c store, const number_pointer_c other)
{
  return *store.asf64 <= *other.asf64;
}

bool worker_multi_double::isle0(const number_pointer_c store)
{
  return *store.asf64 <= 0;
}

bool worker_multi_double::isl0(const number_pointer_c store)
{
  return *store.asf64 < 0;
}

bool worker_multi_double::isl1(const number_pointer_c store)
{
  return *store.asf64 < 1;
}

QString worker_multi_double::toString(const number_pointer_c store)
{
  return QString::number(*store.asf64, 'f', 16);
}

int worker_multi_double::toRound(const number_pointer_c store)
{
  return qRound(*store.asf64);
}

double worker_multi_double::toDouble(const number_pointer_c store)
{
  return *store.asf64;
}
#endif //NUMBER_DOUBLE_EXISTS


#if !ONLY_DOUBLE_WORKER

worker_multi_float128::worker_multi_float128(Allocator *source):
  worker_multi(Type::typeFloat128, specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->capacity)
{
  int oldfirst, oldlast;
  source->_getRange(oldfirst, oldlast);
  storage=new __float128[oldlast-oldfirst];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage[i]=0;
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage();
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
        storage[i]=src_storage[oldfirst+i];
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i]=src_storage[oldfirst+i];
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i]=src_storage[oldfirst+i].hi+(__float128)src_storage[oldfirst+i].lo_;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i]=src_storage[oldfirst+i].hi+(__float128)src_storage[oldfirst+i].lo_;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i]=src_storage[oldfirst+i].val1;
      }
    } break;
  }
}

worker_multi_float128::~worker_multi_float128()
{
  delete[] storage;
  storage=nullptr;
}

number_pointer worker_multi_float128::getNumber(int index)
{
  assert(allocator.checkIndex(index));
  assert(index>=0 && index<capacity);
  return number_pointer(&storage[index]);
}

void worker_multi_float128::getTmp12(number_pointer &t1, number_pointer &t2)
{
  t1.asf128=&tmp1;
  t2.asf128=&tmp2;
}

void worker_multi_float128::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  assert(src_worker->ntype()==Type::typeFloat128);
  worker_multi_float128 *access=(worker_multi_float128 *)src_worker;
  assert(allocator.checkIndex(dst));
  assert(allocator.checkIndex(dst+len-1));
  assert(access->allocator.checkIndex(src));
  assert(access->allocator.checkIndex(src+len-1));
  assert(dst>=0 && dst+len<=capacity);
  assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage[dst+i]=access->storage[src+i];
}

void worker_multi_float128::zero(const number_pointer store, double val)
{
  *store.asf128=val;
}

void worker_multi_float128::assign(const number_pointer store, const number_pointer_c src)
{
  *store.asf128=*src.asf128;
}

void worker_multi_float128::chs(const number_pointer store)
{
  *store.asf128 = -*store.asf128;
}

void worker_multi_float128::lshift(const number_pointer store, int shoft)
{
  //*store.asf128=ldexp(*store.asf128, shoft);
  //could try ldexp only on upper half but...
  if (*store.asf128 != 0) //don't play with 0's exponent or you get NaN
    ((uint16_t *)store.asf128)[7]+=shoft;
}

void worker_multi_float128::round(const number_pointer store)
{
  __float128 remainder=*store.asf128;
  *store.asf128=0;
  for (;;)
  {
    double some=std::round((double)remainder);
    if (some==0)
      break;
    *store.asf128+=some;
    remainder-=some;
  }
}

void worker_multi_float128::frac(const number_pointer store)
{
  if (*store.asf128<0)
  {
    double some=ceil((double)*store.asf128);
    while (some)
    {
      *store.asf128 -= some;
      some=ceil((double)*store.asf128);
    }
  }
  else
  {
    double some=floor((double)*store.asf128);
    while (some)
    {
      *store.asf128 -= some;
      some=floor((double)*store.asf128);
    }
  }
}

void worker_multi_float128::mod1(const number_pointer store)
{
  double some=floor((double)*store.asf128);
  while (some)
  {
    *store.asf128 -= some;
    some=floor((double)*store.asf128);
  }
}

void worker_multi_float128::add_double(const number_pointer store, double x)
{
  *store.asf128 += x;
}

void worker_multi_float128::mul_double(const number_pointer store, double x)
{
  *store.asf128 *= x;
}

void worker_multi_float128::add(const number_pointer store, const number_pointer_c other)
{
  *store.asf128 += *other.asf128;
}

void worker_multi_float128::sub(const number_pointer store, const number_pointer_c other)
{
  *store.asf128 -= *other.asf128;
}

void worker_multi_float128::rsub(const number_pointer store, const number_pointer_c other)
{
  *store.asf128 = *other.asf128-*store.asf128;
}

void worker_multi_float128::mul(const number_pointer store, const number_pointer_c other)
{
  *store.asf128 *= *other.asf128;
}

void worker_multi_float128::sqr(const number_pointer store)
{
  *store.asf128 *= *store.asf128;
}

double worker_multi_float128::radixfloor(const number_pointer_c store1, const number_pointer_c store2)
{
  int ilog1=std::ilogb((double)*store1.asf128);
  int ilog2=std::ilogb((double)*store2.asf128);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

void worker_multi_float128::recip(const number_pointer store)
{
  *store.asf128 = 1/(*store.asf128);
}

void worker_multi_float128::sqrt(const number_pointer store)
{
  double x=1/std::sqrt((double)*store.asf128);
  double ax=x*(double)*store.asf128;
  *store.asf128 = ax + double(*store.asf128 - __float128(ax)*__float128(ax))*x/2;
  //x=1/(sqrt + e)
  //ax=a/(sqrt + e)+f
  //result=a/(sqrt + e)+f + (a - (a/(sqrt + e)+f)*(a/(sqrt + e)+f) + g)*x/2
  //result=a/(sqrt + e)+f + (a - (a*a/(sqrt + e)/(sqrt + e)+2*f*a/(sqrt + e)+f*f) + g)*x/2
  //E=sqrt/(sqrt+e), a=sqrt*sqrt
  //result=sqrt*E+f + (a - a*E*E-2*f*sqrt*E-f*f + g)*1/(sqrt + e)/2
  //result ~ sqrt*E+f + (a*(1-E*E)-2*f*sqrt + g)*1/(sqrt + e)/2
  //result ~ sqrt*E+f + sqrt*(1-E*E)*E/2-f*E + g/sqrt*E/2
  //1-E*E(1-E)*(1+E)
  //result ~ sqrt*E+f + sqrt*(1+E)*(1-E)*E/2-f*E + g/sqrt*E/2
  //E=1+h
  //result ~ sqrt*(1+h)+f + sqrt*(2+h)*(-h)*(1+h)/2-f*(1+h) + g/sqrt*(1+h)/2
  //result ~ sqrt*(1+h) + sqrt*(-h)*(1+h)+sqrt*(-h)*(1+h)*h/2 + g/sqrt*(1+h)/2 -f*h
  //result ~ sqrt + g/sqrt/2 + second order terms
}

int worker_multi_float128::compare(const number_pointer_c store, const number_pointer_c other)
{
  if (*store.asf128 == *other.asf128)
    return 0;
  else if (*store.asf128 < *other.asf128)
    return -1;
  else
    return +1;
}

bool worker_multi_float128::isequal(const number_pointer_c store, const number_pointer_c other)
{
  return *store.asf128 == *other.asf128;
}

bool worker_multi_float128::is0(const number_pointer_c store)
{
  return *store.asf128 == 0;
}

bool worker_multi_float128::isle(const number_pointer_c store, const number_pointer_c other)
{
  return *store.asf128 <= *other.asf128;
}

bool worker_multi_float128::isle0(const number_pointer_c store)
{
  return *store.asf128 <= 0;
}

bool worker_multi_float128::isl0(const number_pointer_c store)
{
  return *store.asf128 < 0;
}

bool worker_multi_float128::isl1(const number_pointer_c store)
{
  return *store.asf128 < 1; //enough to check top 16-17 bits?
}

QString worker_multi_float128::toString(const number_pointer_c store)
{
  return QString::number(*store.asf128, 'f', 16);
}

int worker_multi_float128::toRound(const number_pointer_c store)
{
  return qRound((double)*store.asf128);
}

double worker_multi_float128::toDouble(const number_pointer_c store)
{
  return (double)*store.asf128;
}







number_pointer worker_multi_ddouble::getNumber(int index)
{
  assert(allocator.checkIndex(index));
  assert(index>=0 && index<capacity);
  return number_pointer(&storage[index]);
}

void worker_multi_ddouble::getTmp12(number_pointer &t1, number_pointer &t2)
{
  t1.asdd=&tmp1;
  t2.asdd=&tmp2;
}

worker_multi_ddouble::worker_multi_ddouble(Allocator *source):
  worker_multi(Type::typeDDouble, specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->capacity)
{
  int oldfirst, oldlast;
  source->_getRange(oldfirst, oldlast);
  storage=new dd_real[oldlast-oldfirst];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage[i].zero(0);
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage();
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
        storage[i].zero(src_storage[oldfirst+i]);
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].hi=(double)src_storage[oldfirst+i];
        storage[i].lo_=src_storage[oldfirst+i]-storage[i].hi;
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].assign(src_storage[oldfirst+i]);
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].assign(src_storage[oldfirst+i]);
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].hi=src_storage[oldfirst+i].val1;
        storage[i].lo_=0;
      }
    } break;
  }
}

worker_multi_ddouble::~worker_multi_ddouble()
{
  delete[] storage;
  storage=nullptr;
}

void worker_multi_ddouble::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  assert(src_worker->ntype()==Type::typeDDouble);
  worker_multi_ddouble *access=(worker_multi_ddouble *)src_worker;
  assert(allocator.checkIndex(dst));
  assert(allocator.checkIndex(dst+len-1));
  assert(access->allocator.checkIndex(src));
  assert(access->allocator.checkIndex(src+len-1));
  assert(dst>=0 && dst+len<=capacity);
  assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage[dst+i].assign(access->storage[src+i]);
}


void worker_multi_ddouble::zero(const number_pointer store, double val)
{
  store.asdd->zero(val);
}

/*void number_worker_ddouble::swap(number_store *store, number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeEmpty);
  assert(src->dbgType==Type::typeDDouble);
  assert(store->as.ddouble_.dd==nullptr);
  store->as.ddouble_.dd=src->as.ddouble_.dd;
  src->as.ddouble_.dd=nullptr;
  store->dbgType=Type::typeDDouble;
  src->dbgType=Type::typeEmpty;
}*/

void worker_multi_ddouble::assign(const number_pointer store, const number_pointer_c src)
{
  *store.asdd = *src.asdd;
}

void worker_multi_ddouble::chs(const number_pointer store)
{
  store.asdd->chs();
}

void worker_multi_ddouble::lshift(const number_pointer store, int shoft)
{
  store.asdd->lshift(shoft);
}

void worker_multi_ddouble::round(const number_pointer store)
{
  store.asdd->round();
}

void worker_multi_ddouble::frac(const number_pointer store)
{
  store.asdd->frac();
}

void worker_multi_ddouble::mod1(const number_pointer store)
{
  store.asdd->mod1();
}

void worker_multi_ddouble::add_double(const number_pointer store, double x)
{
  store.asdd->add_double(x);
}

void worker_multi_ddouble::mul_double(const number_pointer store, double x)
{
  store.asdd->mul_double(x);
}

void worker_multi_ddouble::add(const number_pointer store, const number_pointer_c other)
{
  store.asdd->add(other.asdd->hi, other.asdd->lo_);
}

void worker_multi_ddouble::sub(const number_pointer store, const number_pointer_c other)
{
  store.asdd->add(-other.asdd->hi, -other.asdd->lo_);
}
void worker_multi_ddouble::rsub(const number_pointer store, const number_pointer_c other)
{
  store.asdd->chs();
  store.asdd->add(other.asdd->hi, other.asdd->lo_);
}

void worker_multi_ddouble::mul(const number_pointer store, const number_pointer_c other)
{
  store.asdd->mul(other.asdd->hi, other.asdd->lo_);
}

void worker_multi_ddouble::sqr(const number_pointer store)
{
  store.asdd->sqr();
}

double worker_multi_ddouble::radixfloor(const number_pointer_c store1, const number_pointer_c store2)
{
  double rf1=store1.asdd->radixfloor();
  double rf2=store2.asdd->radixfloor();
  if (rf1<rf2)
    return rf2;
  return rf1;
}

void worker_multi_ddouble::recip(const number_pointer store)
{
  store.asdd->recip();
}

void worker_multi_ddouble::sqrt(const number_pointer store)
{
  store.asdd->sqrt();
}

int worker_multi_ddouble::compare(const number_pointer_c store, const number_pointer_c other)
{
  return store.asdd->compare(other.asdd);
}

bool worker_multi_ddouble::isequal(const number_pointer_c store, const number_pointer_c other)
{
  return store.asdd->isequal(other.asdd);
}

bool worker_multi_ddouble::is0(const number_pointer_c store)
{
  return store.asdd->is0();
}

bool worker_multi_ddouble::isle(const number_pointer_c store, const number_pointer_c other)
{
  return store.asdd->isle(other.asdd);
}

bool worker_multi_ddouble::isle0(const number_pointer_c store)
{
  return store.asdd->isle0();
}

bool worker_multi_ddouble::isl0(const number_pointer_c store)
{
  return store.asdd->isl0();
}

bool worker_multi_ddouble::isl1(const number_pointer_c store)
{
  return store.asdd->isl1();
}

QString worker_multi_ddouble::toString(const number_pointer_c store)
{
  return QString("dd(%1,%2)").arg(store.asdd->hi, 0, 'f', 16).arg(store.asdd->lo_, 0, 'g', 16);
}

int worker_multi_ddouble::toRound(const number_pointer_c store)
{
  return floor(store.asdd->hi+0.5)+floor(store.asdd->lo_+0.5);
}

double worker_multi_ddouble::toDouble(const number_pointer_c store)
{
  return store.asdd->hi;
}


number_pointer worker_multi_qdouble::getNumber(int index)
{
  assert(allocator.checkIndex(index));
  assert(index>=0 && index<capacity);
  return number_pointer(&storage[index]);
}

void worker_multi_qdouble::getTmp12(number_pointer &t1, number_pointer &t2)
{
  t1.asqd=&tmp1;
  t2.asqd=&tmp2;
}

worker_multi_qdouble::worker_multi_qdouble(Allocator *source):
  worker_multi(Type::typeDDouble, specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->capacity)
{
  int oldfirst, oldlast;
  source->_getRange(oldfirst, oldlast);
  storage=new dq_real[oldlast-oldfirst];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage[i].zero(0);
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage();
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
        storage[i].zero(src_storage[oldfirst+i]);
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].hi=(double)src_storage[oldfirst+i];
        storage[i].lo_=src_storage[oldfirst+i]-storage[i].hi;
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].assign(src_storage[oldfirst+i]);
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].assign(src_storage[oldfirst+i]);
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].hi=src_storage[oldfirst+i].val1;
        storage[i].lo_=0;
      }
    } break;
  }
}

worker_multi_qdouble::~worker_multi_qdouble()
{
  delete[] storage;
  storage=nullptr;
}

void worker_multi_qdouble::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  assert(src_worker->ntype()==Type::typeQDouble);
  worker_multi_qdouble *access=(worker_multi_qdouble *)src_worker;
  assert(allocator.checkIndex(dst));
  assert(allocator.checkIndex(dst+len-1));
  assert(access->allocator.checkIndex(src));
  assert(access->allocator.checkIndex(src+len-1));
  assert(dst>=0 && dst+len<=capacity);
  assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage[dst+i].assign(access->storage[src+i]);
}

void worker_multi_qdouble::zero(const number_pointer store, double val)
{
  store.asqd->zero(val);
}

/*void number_worker_ddouble::swap(number_store *store, number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeEmpty);
  assert(src->dbgType==Type::typeDDouble);
  assert(store->as.ddouble_.dd==nullptr);
  store->as.ddouble_.dd=src->as.ddouble_.dd;
  src->as.ddouble_.dd=nullptr;
  store->dbgType=Type::typeDDouble;
  src->dbgType=Type::typeEmpty;
}*/

void worker_multi_qdouble::assign(const number_pointer store, const number_pointer_c src)
{
  *store.asqd = *src.asqd;
}

void worker_multi_qdouble::chs(const number_pointer store)
{
  store.asqd->chs();
}

void worker_multi_qdouble::lshift(const number_pointer store, int shoft)
{
  store.asqd->lshift(shoft);
}

void worker_multi_qdouble::round(const number_pointer store)
{
  store.asqd->round();
}

void worker_multi_qdouble::frac(const number_pointer store)
{
  store.asqd->frac();
}

void worker_multi_qdouble::mod1(const number_pointer store)
{
  store.asqd->mod1();
}

void worker_multi_qdouble::add_double(const number_pointer store, double x)
{
  store.asqd->add_double(x);
}

void worker_multi_qdouble::mul_double(const number_pointer store, double x)
{
  store.asqd->mul_double(x);
}

void worker_multi_qdouble::add(const number_pointer store, const number_pointer_c other)
{
  store.asqd->add(other.asqd->hi, other.asqd->lo_);
}

void worker_multi_qdouble::sub(const number_pointer store, const number_pointer_c other)
{
  store.asqd->add(-other.asqd->hi, -other.asqd->lo_);
}

void worker_multi_qdouble::rsub(const number_pointer store, const number_pointer_c other)
{
  store.asqd->chs();
  store.asqd->add(other.asqd->hi, other.asqd->lo_);
}

void worker_multi_qdouble::mul(const number_pointer store, const number_pointer_c other)
{
  store.asqd->mul(other.asqd->hi, other.asqd->lo_);
}

void worker_multi_qdouble::sqr(const number_pointer store)
{
  store.asqd->sqr();
}

double worker_multi_qdouble::radixfloor(const number_pointer_c store1, const number_pointer_c store2)
{
  double rf1=store1.asqd->radixfloor();
  double rf2=store2.asqd->radixfloor();
  if (rf1<rf2)
    return rf2;
  return rf1;
}

void worker_multi_qdouble::recip(const number_pointer store)
{
  store.asqd->recip();
}

void worker_multi_qdouble::sqrt(const number_pointer store)
{
  store.asqd->sqrt();
}

int worker_multi_qdouble::compare(const number_pointer_c store, const number_pointer_c other)
{
  return store.asqd->compare(other.asqd);
}

bool worker_multi_qdouble::isequal(const number_pointer_c store, const number_pointer_c other)
{
  return store.asqd->isequal(other.asqd);
}

bool worker_multi_qdouble::is0(const number_pointer_c store)
{
  return store.asqd->is0();
}

bool worker_multi_qdouble::isle(const number_pointer_c store, const number_pointer_c other)
{
  return store.asqd->isle(other.asqd);
}

bool worker_multi_qdouble::isle0(const number_pointer_c store)
{
  return store.asqd->isle0();
}

bool worker_multi_qdouble::isl0(const number_pointer_c store)
{
  return store.asqd->isl0();
}

bool worker_multi_qdouble::isl1(const number_pointer_c store)
{
  return store.asqd->isl1();
}

QString worker_multi_qdouble::toString(const number_pointer_c store)
{
  return QString("dd(%1,%2)").arg(store.asqd->hi, 0, 'f', 16).arg(store.asqd->lo_, 0, 'g', 16);
}

int worker_multi_qdouble::toRound(const number_pointer_c store)
{
  return floor(store.asqd->hi+0.5)+floor(store.asqd->lo_+0.5);
}

double worker_multi_qdouble::toDouble(const number_pointer_c store)
{
  return store.asqd->hi;
}





worker_multi_real642::worker_multi_real642(Allocator *source):
  worker_multi(Type::typeReal642, specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->capacity)
{
  int oldfirst, oldlast;
  source->_getRange(oldfirst, oldlast);
  storage=new real642[oldlast-oldfirst];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
      {
        storage[i].val1=0;
        storage[i].val2=0;
      }
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage();
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
      {
        storage[i].val1=src_storage[oldfirst+i];
        storage[i].val2=src_storage[oldfirst+i];
      }
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].val1=src_storage[oldfirst+i];
        storage[i].val2=src_storage[oldfirst+i];
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].val1=src_storage[oldfirst+i].hi;
        storage[i].val2=src_storage[oldfirst+i].hi+(__float128)src_storage[oldfirst+i].lo_;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].val1=src_storage[oldfirst+i].hi;
        storage[i].val2=src_storage[oldfirst+i].hi+(__float128)src_storage[oldfirst+i].lo_;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage();
      for (int i=0; i<capacity; i++)
      {
        storage[i].val1=src_storage[oldfirst+i].val1;
        storage[i].val2=src_storage[oldfirst+i].val2;
      }
    } break;
  }
}

worker_multi_real642::~worker_multi_real642()
{
  delete[] storage;
  storage=nullptr;
}

number_pointer worker_multi_real642::getNumber(int index)
{
  assert(allocator.checkIndex(index));
  assert(index>=0 && index<capacity);
  return number_pointer(&storage[index]);
}

void worker_multi_real642::getTmp12(number_pointer &t1, number_pointer &t2)
{
  t1.as642=&tmp1;
  t2.as642=&tmp2;
}

void worker_multi_real642::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  assert(src_worker->ntype()==Type::typeReal642);
  worker_multi_real642 *access=(worker_multi_real642 *)src_worker;
  assert(allocator.checkIndex(dst));
  assert(allocator.checkIndex(dst+len-1));
  assert(access->allocator.checkIndex(src));
  assert(access->allocator.checkIndex(src+len-1));
  assert(dst>=0 && dst+len<=capacity);
  assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage[dst+i]=access->storage[src+i];
}

void worker_multi_real642::zero(const number_pointer store, double val)
{
  store.as642->val1=val;
  store.as642->val2=val;
}

void worker_multi_real642::assign(const number_pointer store, const number_pointer_c src)
{
  *store.as642=*src.as642;
}

void worker_multi_real642::chs(const number_pointer store)
{
  store.as642->val1 = -store.as642->val1;
  store.as642->val2 = -store.as642->val2;
}

void worker_multi_real642::lshift(const number_pointer store, int shoft)
{
  //*store.asf128=ldexp(*store.asf128, shoft);
  //could try ldexp only on upper half but...
  store.as642->val1=ldexp(store.as642->val1, shoft);
  if (store.as642->val2 != 0) //don't play with 0's exponent or you get NaN
    ((uint16_t *)&store.as642->val2)[7]+=shoft;
}

void worker_multi_real642::round(const number_pointer store)
{
  store.as642->val1=std::round(store.as642->val1);
  __float128 remainder=store.as642->val2;
  store.as642->val2=0;
  for (;;)
  {
    double some=std::round((double)remainder);
    if (some==0)
      break;
    store.as642->val2+=some;
    remainder-=some;
  }
}

void worker_multi_real642::frac(const number_pointer store)
{
  if (store.as642->val1<0)
    store.as642->val1 -= ceil(store.as642->val1);
  else
    store.as642->val1 -= floor(store.as642->val1);
  if (store.as642->val2<0)
  {
    double some=ceil((double)store.as642->val2);
    while (some)
    {
      store.as642->val2 -= some;
      some=ceil((double)store.as642->val2);
    }
  }
  else
  {
    double some=floor((double)store.as642->val2);
    while (some)
    {
      store.as642->val2 -= some;
      some=floor((double)store.as642->val2);
    }
  }
}

void worker_multi_real642::mod1(const number_pointer store)
{
  store.as642->val1 -= floor(store.as642->val1);
  double some=floor((double)store.as642->val2);
  while (some)
  {
    store.as642->val2 -= some;
    some=floor((double)store.as642->val2);
  }
}

void worker_multi_real642::add_double(const number_pointer store, double x)
{
  store.as642->val1 += x;
  store.as642->val2 += x;
}

void worker_multi_real642::mul_double(const number_pointer store, double x)
{
  store.as642->val1 *= x;
  store.as642->val2 *= x;
}

void worker_multi_real642::add(const number_pointer store, const number_pointer_c other)
{
  store.as642->val1 += other.as642->val1;
  store.as642->val2 += other.as642->val2;
}

void worker_multi_real642::sub(const number_pointer store, const number_pointer_c other)
{
  store.as642->val1 -= other.as642->val1;
  store.as642->val2 -= other.as642->val2;
}

void worker_multi_real642::rsub(const number_pointer store, const number_pointer_c other)
{
  store.as642->val1 = other.as642->val1-store.as642->val1;
  store.as642->val2 = other.as642->val2-store.as642->val2;
}

void worker_multi_real642::mul(const number_pointer store, const number_pointer_c other)
{
  store.as642->val1 *= other.as642->val1;
  store.as642->val2 *= other.as642->val2;
}

void worker_multi_real642::sqr(const number_pointer store)
{
  store.as642->val1 *= store.as642->val1;
  store.as642->val2 *= store.as642->val2;
}

double worker_multi_real642::radixfloor(const number_pointer_c store1, const number_pointer_c store2)
{
  int ilog1=std::ilogb(store1.as642->val1);
  int ilog2=std::ilogb(store2.as642->val1);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

void worker_multi_real642::recip(const number_pointer store)
{
  store.as642->val1 = 1/store.as642->val1;
  store.as642->val2 = 1/store.as642->val2;
}

void worker_multi_real642::sqrt(const number_pointer store)
{
  store.as642->val1 = std::sqrt(store.as642->val1);
  double x=1/std::sqrt((double)store.as642->val2);
  double ax=x*(double)store.as642->val2;
  store.as642->val2 = ax + double(store.as642->val2 - __float128(ax)*__float128(ax))*x/2;
}

int worker_multi_real642::compare(const number_pointer_c store, const number_pointer_c other)
{
  if (store.as642->val1 == other.as642->val1)
    return 0;
  else if (store.as642->val1 < other.as642->val1)
    return -1;
  else
    return +1;
}

bool worker_multi_real642::isequal(const number_pointer_c store, const number_pointer_c other)
{
  return store.as642->val1 == other.as642->val1;
}

bool worker_multi_real642::is0(const number_pointer_c store)
{
  return store.as642->val1 == 0;
}

bool worker_multi_real642::isle(const number_pointer_c store, const number_pointer_c other)
{
  return store.as642->val1 <= other.as642->val1;
}

bool worker_multi_real642::isle0(const number_pointer_c store)
{
  return store.as642->val1 <= 0;
}

bool worker_multi_real642::isl0(const number_pointer_c store)
{
  return store.as642->val1 < 0;
}

bool worker_multi_real642::isl1(const number_pointer_c store)
{
  return store.as642->val1 < 1;
}

QString worker_multi_real642::toString(const number_pointer_c store)
{
  return QString::number(store.as642->val1, 'f', 16);
}

int worker_multi_real642::toRound(const number_pointer_c store)
{
  return qRound(store.as642->val1);
}

double worker_multi_real642::toDouble(const number_pointer_c store)
{
  return store.as642->val1;
}

#endif


/*
number::number(number_pointer *src):
  worker(src->worker), allocator(src->allocator),
  free_storage_on_destroy(src->free_storage_on_destroy),
  ptr(src->ptr)
{
  src->ptr.asf64=nullptr;
}*/

void number::zero(double v)
{
  worker->zero(ptr, v);
}

void number::assign(const number_pointer_c src)
{
  worker->assign(ptr, src);
}

void number::lshift(int shoft)
{
  worker->lshift(ptr, shoft);
}

void number::add(const number_pointer_c other)
{
  worker->add(ptr, other);
}

void number::sub(const number_pointer_c other)
{
  worker->sub(ptr, other);
}

void number::sqr()
{
  worker->sqr(ptr);
}

void number::mul(const number_pointer_c other)
{
  worker->mul(ptr, other);
}

void number::recip()
{
  worker->recip(ptr);
}

void number::sqrt()
{
  worker->sqrt(ptr);
}

void number::add_double(double x)
{
  worker->add_double(ptr, x);
}

double number::toDouble() const
{
  return worker->toDouble(ptr);
}


void complex::zero(double r, double i)
{
  worker->zero(re, r);
  worker->zero(im, i);
}

void complex::assign(const complex *other)
{
  worker->assign(re, other->re);
  worker->assign(im, other->im);
}

void complex::lshift(int shoft)
{
  worker->lshift(re, shoft);
  worker->lshift(im, shoft);
}

void complex::add(const complex *other)
{
  worker->add(re, other->re);
  worker->add(im, other->im);
}

void complex::sub(const complex *other)
{
  worker->sub(re, other->re);
  worker->sub(im, other->im);
}

void complex::rsub(const complex *other)
{
  worker->rsub(re, other->re);
  worker->rsub(im, other->im);
}

void complex::mul(const complex *other)
{
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r1*r2-i1*i2
  //i:=r1*i2+i1*r2
  worker->assign(tmp1, re);
  worker->mul(tmp1, other->re);
  worker->assign(tmp2, im);
  worker->mul(tmp2, other->im);
  worker->sub(tmp1, tmp2); //r1*r2-i1*i2
  worker->mul(re, other->im);
  worker->mul(im, other->re);
  worker->add(im, re); //i1*r2+r1*i2
  worker->assign(re, tmp1);
}

void complex::sqr()
{
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r*r-i*i
  //i=2*r*i
  worker->assign(tmp1, im);
  worker->sqr(tmp1);
  worker->mul(im, re);
  worker->lshift(im, 1);
  worker->sqr(re);
  worker->sub(re, tmp1);
}

void complex::recip()
{
  getMag_tmp_();
  recip_prepared();
}

void complex::recip_prepared()
{ // 1/(re+i*im) = (re-i*im)/((re+i*im)*(re-i*im)) = (re-i*im)/(re*re+im*im)
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->recip(tmp1);
  worker->mul(re, tmp1);
  worker->chs(im);
  worker->mul(im, tmp1);
}

void complex::sqrt()
{
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->assign(tmp2, im);
  worker->sqr(tmp1);
  worker->sqr(tmp2);
  worker->add(tmp1, tmp2); //re*re+im*im
  worker->sqrt(tmp1);
  if (!worker->isle0(re))
  {
    worker->add(tmp1, re);
    worker->lshift(tmp1, -1);
    worker->sqrt(tmp1); //t1=sqrt((sqrt(re*re+im*im)+re)/2);
    worker->assign(re, tmp1); //re=t1
    /*if (t1==0)
      *res_im=0;
    else*/
    worker->lshift(tmp1, 1);
    worker->recip(tmp1);
    worker->mul(im, tmp1); //im=im/(2*t1);
  }
  else
  {
    worker->sub(tmp1, re);
    worker->lshift(tmp1, -1);
    worker->sqrt(tmp1); //t1=sqrt((sqrt(re*re+im*im)-re)/2);
    if (worker->isle0(tmp1)) //t1==0
    {
      worker->zero(re);
      worker->zero(im);
    }
    else
    {
      worker->assign(re, im);
      worker->assign(im, tmp1); //new im=t1
      worker->lshift(tmp1, 1);
      worker->recip(tmp1);
      worker->mul(re, tmp1); //re=old im/(2*t1);
      if (worker->isl0(re))
      {
        worker->chs(re);
        worker->chs(im);
      };
    }
  };
}

double complex::getMag_double() const
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->sqr(tmp1);
  worker->assign(tmp2, im);
  worker->sqr(tmp2);
  worker->add(tmp1, tmp2);
  return worker->toDouble(tmp1);
}

const number_pointer_c complex::getMag_tmp_() const
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->sqr(tmp1);
  worker->assign(tmp2, im);
  worker->sqr(tmp2);
  worker->add(tmp1, tmp2);
  return tmp1;
}

const number_pointer_c complex::getMag1_tmp() const
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->sqr(tmp1);
  worker->assign(tmp2, im);
  worker->sqr(tmp2);
  worker->add(tmp1, tmp2);
  worker->add_double(tmp1, -1);
  return tmp1;
}

const number_pointer_c complex::getDist1_tmp() const
{
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->add_double(tmp1, -1);
  worker->sqr(tmp1);
  worker->assign(tmp2, im);
  worker->sqr(tmp2);
  worker->add(tmp1, tmp2);
  return tmp1;
}

const number_pointer_c complex::mulreT_tmp(const complex *other) const //Re(this*conjugate(other))
{ //re*ore+im*oim
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->assign(tmp2, im);
  worker->mul(tmp1, other->re);
  worker->mul(tmp2, other->im);
  worker->add(tmp1, tmp2);
  return tmp1;
}

double complex::dist2_double(const complex *other) const
{
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->assign(tmp2, im);
  worker->sub(tmp1, other->re);
  worker->sub(tmp2, other->im);
  worker->sqr(tmp1);
  worker->sqr(tmp2);
  worker->add(tmp1, tmp2);
  return worker->toDouble(tmp1);
}

const number_pointer_c complex::dist2_tmp_(const complex *other) const
{
  number_pointer tmp1, tmp2;
  worker->getTmp12(tmp1, tmp2);
  worker->assign(tmp1, re);
  worker->assign(tmp2, im);
  worker->sub(tmp1, other->re);
  worker->sub(tmp2, other->im);
  worker->sqr(tmp1);
  worker->sqr(tmp2);
  worker->add(tmp1, tmp2);
  return tmp1;
}

bool complex::isequal(const complex *other) const
{
  return worker->isequal(re, other->re) &&
         worker->isequal(im, other->im);
}

bool complex::is0() const
{
  return worker->is0(re) && worker->is0(im);
}

QString complex::toString()
{
  return worker->toString(re)+" +i* "+worker->toString(im);
}

double sqr_double(double x)
{
  return x*x;
}

void complex_double_sqrt(double *res_re, double *res_im, double in_re, double in_im)
{
  double t1;
  if (in_re>=0)
  {
    t1=sqrt((sqrt(in_re*in_re+in_im*in_im)+in_re)/2);
    *res_re=t1;
    if (t1==0)
      *res_im=0;
    else
      *res_im=in_im/(2*t1);
  }
  else
  {
    t1=sqrt((sqrt(in_re*in_re+in_im*in_im)-in_re)/2);
    if (in_im>=0)
    {
      *res_re=in_im/(2*t1);
      *res_im=t1;
    }
    else
    {
      *res_re=-in_im/(2*t1);
      *res_im=-t1;
    };
  };
}

/*
 finds smaller (in abs) root of ax^2+2bx+c=0
*/
void complex_double_quadratic(double *res_re, double *res_im,
                              double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im)
{
  /*
  a x^2 + 2 b x + c = 0
  x1,x2= -(b +- sqrt(b^2-a*c))/a
  for small a
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/a
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/(a*c/b^2)*c/b^2
  HELP(x)=(1+-sqrt(1-x))/x =1/(1-+sqrt(1-x))
  x1,x2= -c/b*HELP(a*c/b^2)   or rather, for a<b, c<b
  for b<a
  (F1)  x1,x2= -(b/a +- sqrt(b^2/a^2-c/a))     good for a>b, c>b until ~ c>b^2/a
  (F1') x1,x2= -(1/a)*(b - sqrt(b^2-a*c))      good for b^2<<a*c
  (F2)  x1,x2= -(b/a)*(1 +- sqrt(1-a*c/b^2))   good for a>b, c<b^2/a
  for b>a  b^2/a>b
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/a/c*b^2*c/b^2
        x1,x2= -(c/b)*HELP(a*c/b^2)            good for c<b, until c<b^2/a
  (F3)  x1=    -(c/b)/(1+sqrt(1-a*c/b^2))      good for b>0 (cancellation for x2)
  (F3')        -(c)/(b+-sqrt(b^2-a*c))         good except both b,c small e.g. 0

  Muller's method: B=2b
  x1,x2= c/(-b+-sqrt(b^2-ac))     -c/(b+sqrt(b^2-ac))
  x1,x2= c/b/(-1+-sqrt(1-ac/b^2))
  */
  double bb_re=b2_re*b2_re-b2_im*b2_im;
  double bb_im=2*b2_re*b2_im;
  double bbac_re=bb_re-a_re*c_re+a_im*c_im;
  double bbac_im=bb_im-a_im*c_re-a_re*c_im; //b^2-ac
  double d_re, d_im;
  complex_double_sqrt(&d_re, &d_im, bbac_re, bbac_im);
  double t1_re, t1_im;
  if (b2_re*d_re+b2_im*d_im<0)//if (b2_re<0)
  {
    t1_re=b2_re-d_re;
    t1_im=b2_im-d_im;
  }
  else
  {
    t1_re=b2_re+d_re;
    t1_im=b2_im+d_im;
  }
  double am=a_re*a_re+a_im*a_im;
  double t1m=t1_re*t1_re+t1_im*t1_im;
  if (t1m>am)
  { //F3'
    *res_re=-(c_re*t1_re+c_im*t1_im)/t1m;
    *res_im=-(c_im*t1_re-c_re*t1_im)/t1m; //-c/(b+sqrt(b^2-a*c))
  }
  else
  {
    t1_re=2*b2_re-t1_re;//b2_re-d_re;
    t1_im=2*b2_im-t1_im;//b2_im-d_im;
    *res_re=-(t1_re*a_re+t1_im*a_im)/am;
    *res_im=-(t1_im*a_re-t1_re*a_im)/am; //-(b - sqrt(b^2-a*c))/a
  }
}

/*
 finds both roots of ax^2+2bx+c=0
*/
void complex_double_quadratic2(double *res1_re, double *res1_im,
                               double *res2_re, double *res2_im,
                               double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im)
{
  //see complex_double_quadratic()
  double bb_re=b2_re*b2_re-b2_im*b2_im;
  double bb_im=2*b2_re*b2_im;
  double bbac_re=bb_re-a_re*c_re+a_im*c_im;
  double bbac_im=bb_im-a_im*c_re-a_re*c_im; //b^2-ac
  double d_re, d_im;
  complex_double_sqrt(&d_re, &d_im, bbac_re, bbac_im);
  double t1_re, t1_im;
  double t2_re, t2_im;
  if (b2_re*d_re+b2_im*d_im<0)//if (b2_re<0)
  {
    t1_re=b2_re-d_re;
    t1_im=b2_im-d_im;
    t2_re=b2_re+d_re;
    t2_im=b2_im+d_im;
  }
  else
  {
    t1_re=b2_re+d_re;
    t1_im=b2_im+d_im;
    t2_re=b2_re-d_re;
    t2_im=b2_im-d_im;
  }
  double am=a_re*a_re+a_im*a_im;
  double t1m=t1_re*t1_re+t1_im*t1_im;
  if (t1m>am)
  { //F3'
    *res1_re=-(c_re*t1_re+c_im*t1_im)/t1m;
    *res1_im=-(c_im*t1_re-c_re*t1_im)/t1m; //-c/(b+sqrt(b^2-a*c))
  }
  else
  {
    *res1_re=-(t2_re*a_re+t2_im*a_im)/am;
    *res1_im=-(t2_im*a_re-t2_re*a_im)/am; //-(b - sqrt(b^2-a*c))/a
  }
  double t2m=t2_re*t2_re+t2_im*t2_im;
  if (t2m>am)
  { //F3'
    *res2_re=-(c_re*t2_re+c_im*t2_im)/t2m;
    *res2_im=-(c_im*t2_re-c_re*t2_im)/t2m; //-c/(b+sqrt(b^2-a*c))
  }
  else
  {
    *res2_re=-(t1_re*a_re+t1_im*a_im)/am;
    *res2_im=-(t1_im*a_re-t1_re*a_im)/am; //-(b - sqrt(b^2-a*c))/a
  }
}


} // namespace MandelMath

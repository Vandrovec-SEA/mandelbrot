#define assert(x) { if (!(x)) dbgPoint(); }
#include "MandelMath.hpp"

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
  /* shift, append i, rotate to be max; is i followed by 0 only? -> negate else not negate
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
              0010-010i-10i0 or i010-y or n-1 or 0  (10i0 correct)
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
    case number_worker::Type::typeDouble:
      assert(dbgType==number_worker::Type::typeDouble);
      dbgType=number_worker::Type::typeEmpty;
      as.doubl=0.0;
      break;
    case number_worker::Type::typeDDouble:
      assert(dbgType==number_worker::Type::typeDDouble);
      dbgType=number_worker::Type::typeEmpty;
      as.ddouble_.deinit();
      break;
    case number_worker::Type::typeMulti:
      assert(dbgType==number_worker::Type::typeMulti);
      dbgType=number_worker::Type::typeEmpty;
      as.multi_.deinit();
      break;
    case number_worker::Type::typeEmpty: ;
  }
}

void number_store::init(number_worker::Type ntype, double val)
{
  switch (ntype)
  {
    case number_worker::Type::typeDouble:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeDouble;
      as.doubl=val;
      break;
    case number_worker::Type::typeDDouble:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeDDouble;
      as.ddouble_.init();
      as.ddouble_.dd->hi=val;
      as.ddouble_.dd->lo=0.0;
      break;
    case number_worker::Type::typeMulti:
      assert(dbgType==number_worker::Type::typeEmpty);
      dbgType=number_worker::Type::typeMulti;
      as.multi_.init();
      as.multi_.bytes->set(val);
      break;
    case number_worker::Type::typeEmpty: ;
  }
}

void number_store::zero(number_worker::Type ntype, double val)
{
  switch (ntype)
  {
    case number_worker::Type::typeDouble:
      assert(dbgType==number_worker::Type::typeDouble);
      as.doubl=val;
      break;
    case number_worker::Type::typeDDouble:
      assert(dbgType==number_worker::Type::typeDDouble);
      as.ddouble_.dd->hi=val;
      as.ddouble_.dd->lo=0.0;
      break;
    case number_worker::Type::typeMulti:
      assert(dbgType==number_worker::Type::typeMulti);
      as.multi_.bytes->set(val);
      break;
    case number_worker::Type::typeEmpty: ;
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




void number_worker_double::init(number_store *store, double val)
{
  store->init(number_worker::Type::typeDouble, val);
}

void number_worker_double::zero(number_store *store, double val)
{
  store->zero(number_worker::Type::typeDouble, val);
}

void number_worker_double::swap(number_store *store, number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeEmpty);
  assert(src->dbgType==Type::typeDouble);
  store->as.doubl=src->as.doubl;
  store->dbgType=Type::typeDouble;
  src->dbgType=Type::typeEmpty;
}

void number_worker_double::assign(number_store *store, const number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeDouble);
  assert(src->dbgType==Type::typeDouble);
  store->as.doubl=src->as.doubl;
}

void number_worker_double::chs(number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  store->as.doubl=-store->as.doubl;
}

void number_worker_double::lshift(number_store *store, int shoft)
{
  assert(store->dbgType==number_worker::Type::typeDouble);
  store->as.doubl=ldexp(store->as.doubl, shoft);
  /*if (shoft>0)
    store->as.doubl*=two_pow_n(shoft);
  else if (shoft<0)
    store->as.doubl/=two_pow_n(-shoft);*/
}

void number_worker_double::frac_pos(number_store *store)
{
  assert(store->dbgType==number_worker::Type::typeDouble);
  store->as.doubl-=floor(store->as.doubl);
}

void number_worker_double::add_double(number_store *store, double x)
{
  assert(store->dbgType==number_worker::Type::typeDouble);
  store->as.doubl+=x;
}

void number_worker_double::add(number_store *store, const number_store *other)
{
  assert(store->dbgType==number_worker::Type::typeDouble);
  assert(other->dbgType==number_worker::Type::typeDouble);
  store->as.doubl+=other->as.doubl;
}

void number_worker_double::sub(number_store *store, const number_store *other)
{
  assert(store->dbgType==number_worker::Type::typeDouble);
  assert(other->dbgType==number_worker::Type::typeDouble);
  store->as.doubl-=other->as.doubl;
}

void number_worker_double::rsub(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDouble);
  assert(other->dbgType==Type::typeDouble);
  store->as.doubl=other->as.doubl-store->as.doubl;
}

void number_worker_double::mul(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDouble);
  assert(other->dbgType==Type::typeDouble);
  store->as.doubl*=other->as.doubl;
}

void number_worker_double::sqr(number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  store->as.doubl*=store->as.doubl;
}

double number_worker_double::radixfloor(number_store *store1, number_store *store2)
{
  assert(store1->dbgType==Type::typeDouble);
  assert(store2->dbgType==Type::typeDouble);
  int ilog1=std::ilogb(store1->as.doubl);
  int ilog2=std::ilogb(store2->as.doubl);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

void number_worker_double::recip(number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  store->as.doubl=1/store->as.doubl;
}

void number_worker_double::sqrt(number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  store->as.doubl=std::sqrt(store->as.doubl);
}

int number_worker_double::compare(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDouble);
  assert(other->dbgType==Type::typeDouble);
  if (store->as.doubl==other->as.doubl)
    return 0;
  else if (store->as.doubl<other->as.doubl)
    return -1;
  else
    return +1;
}

bool number_worker_double::isequal(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDouble);
  assert(other->dbgType==Type::typeDouble);
  return store->as.doubl==other->as.doubl;
}

bool number_worker_double::is0(const number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  return store->as.doubl==0;
}

bool number_worker_double::isle(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDouble);
  assert(other->dbgType==Type::typeDouble);
  return store->as.doubl<=other->as.doubl;
}

bool number_worker_double::isle0(const number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  return store->as.doubl<=0;
}

bool number_worker_double::isl0(const number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  return store->as.doubl<0;
}

bool number_worker_double::isl1(const number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  return store->as.doubl<1;
}

QString number_worker_double::toString(const number_store *store)
{
  assert(store->dbgType==number_worker::Type::typeDouble);
  return QString::number(store->as.doubl, 'f', 16);
}

int number_worker_double::toRound(const number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  return qRound(store->as.doubl);
}

double number_worker_double::toDouble(const number_store *store)
{
  assert(store->dbgType==Type::typeDouble);
  return store->as.doubl;
}



void number_worker_ddouble::init(number_store *store, double val)
{
  store->init(number_worker::Type::typeDDouble, val);
}

void number_worker_ddouble::zero(number_store *store, double val)
{
  store->zero(number_worker::Type::typeDDouble, val);
}

void number_worker_ddouble::swap(number_store *store, number_store *src)
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
}

void number_worker_ddouble::assign(number_store *store, const number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeDouble);
  assert(src->dbgType==Type::typeDouble);
  *store->as.ddouble_.dd=*src->as.ddouble_.dd;
}

void number_worker_ddouble::chs(number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->chs();
}

void number_worker_ddouble::lshift(number_store *store, int shoft)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->lshift(shoft);
}

void number_worker_ddouble::frac_pos(number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->hi-=floor(store->as.ddouble_.dd->hi);
  store->as.ddouble_.dd->lo-=floor(store->as.ddouble_.dd->lo);
}

void number_worker_ddouble::add_double(number_store *store, double x)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->add_double(x);
}

void number_worker_ddouble::add(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->add(other->as.ddouble_.dd->hi, other->as.ddouble_.dd->lo);
}

void number_worker_ddouble::sub(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->add(-other->as.ddouble_.dd->hi, -other->as.ddouble_.dd->lo);
}
void number_worker_ddouble::rsub(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->chs();
  store->as.ddouble_.dd->add(other->as.ddouble_.dd->hi, other->as.ddouble_.dd->lo);
}

void number_worker_ddouble::mul(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->mul(other->as.ddouble_.dd->hi, other->as.ddouble_.dd->lo);
}

void number_worker_ddouble::sqr(number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->sqr();
}

double number_worker_ddouble::radixfloor(number_store *store1, number_store *store2)
{
  assert(store1->dbgType==Type::typeDDouble);
  assert(store2->dbgType==Type::typeDDouble);
  double rf1=store1->as.ddouble_.dd->radixfloor();
  double rf2=store2->as.ddouble_.dd->radixfloor();
  if (rf1<rf2)
    return rf2;
  return rf1;
}

void number_worker_ddouble::recip(number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->recip();
}

void number_worker_ddouble::sqrt(number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->sqrt();
}

int number_worker_ddouble::compare(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->compare(other->as.ddouble_.dd);
}

bool number_worker_ddouble::isequal(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->isequal(other->as.ddouble_.dd);
}

bool number_worker_ddouble::is0(const number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->is0();
}

bool number_worker_ddouble::isle(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->isle(other->as.ddouble_.dd);
}

bool number_worker_ddouble::isle0(const number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->isle0();
}

bool number_worker_ddouble::isl0(const number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->isl0();
}

bool number_worker_ddouble::isl1(const number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->isl1();
}

QString number_worker_ddouble::toString(const number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  return QString("dd(%1,%2)").arg(store->as.ddouble_.dd->hi).arg(store->as.ddouble_.dd->lo);
}

int number_worker_ddouble::toRound(const number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  return floor(store->as.ddouble_.dd->hi+0.5)+floor(store->as.ddouble_.dd->lo+0.5);
}

double number_worker_ddouble::toDouble(const number_store *store)
{
  assert(store->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->hi;
}




void number_worker_multi::init(number_store *store, double val)
{
  store->init(Type::typeMulti, val);
}

void number_worker_multi::zero(number_store *store, double val)
{
  store->zero(Type::typeMulti, val);
}

void number_worker_multi::swap(number_store *store, number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeEmpty);
  assert(src->dbgType==Type::typeMulti);
  assert(store->as.multi_.bytes==nullptr);
  store->as.multi_.bytes=src->as.multi_.bytes;
  src->as.multi_.bytes=nullptr;
  store->dbgType=Type::typeMulti;
  src->dbgType=Type::typeEmpty;
}

void number_worker_multi::assign(number_store *store, const number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeDouble);
  assert(src->dbgType==Type::typeDouble);
  *store->as.multi_.bytes=*src->as.multi_.bytes;
}

void number_worker_multi::chs(number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->chs();
}

void number_worker_multi::lshift(number_store *store, int shoft)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->lshift(shoft);
}

void number_worker_multi::frac_pos(number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->frac_pos();
}

void number_worker_multi::add_double(number_store *store, double x)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->add_double(x);
}

void number_worker_multi::add(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->add(other->as.multi_.bytes);
}

void number_worker_multi::sub(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->chs();
  store->as.multi_.bytes->add(other->as.multi_.bytes);
  store->as.multi_.bytes->chs();
}
void number_worker_multi::rsub(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->chs();
  store->as.multi_.bytes->add(other->as.multi_.bytes);
}

void number_worker_multi::mul(number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->mul(other->as.multi_.bytes);
}

void number_worker_multi::sqr(number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->sqr();
}

double number_worker_multi::radixfloor(number_store *store1, number_store *store2)
{
  assert(store1->dbgType==Type::typeMulti);
  assert(store2->dbgType==Type::typeMulti);
  double rf1=store1->as.multi_.bytes->radixfloor();
  double rf2=store2->as.multi_.bytes->radixfloor();
  if (rf1<rf2)
    return rf2;
  return rf1;
}

void number_worker_multi::recip(number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->recip();
}

void number_worker_multi::sqrt(number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->sqrt();
}

int number_worker_multi::compare(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->compare(other->as.multi_.bytes);
}

bool number_worker_multi::isequal(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->isequal(other->as.multi_.bytes);
}

bool number_worker_multi::is0(const number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->is0();
}

bool number_worker_multi::isle(const number_store *store, const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->isle(other->as.multi_.bytes);
}

bool number_worker_multi::isle0(const number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->isle0();
}

bool number_worker_multi::isl0(const number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->isl0();
}

bool number_worker_multi::isl1(const number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->isl1();
}

QString number_worker_multi::toString(const number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  return QString("multi");
}

int number_worker_multi::toRound(const number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->round();
}

double number_worker_multi::toDouble(const number_store *store)
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->toDouble();
}


#if !COMPLEX_IS_TEMPLATE

const number_store *complex::getMagTmp()
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  worker->assign(&tmp1_s, re_s);
  worker->sqr(&tmp1_s);
  worker->assign(&tmp2_s, im_s);
  worker->sqr(&tmp2_s);
  worker->add(&tmp1_s, &tmp2_s);
  return &tmp1_s;
}

const number_store *complex::getMag1Tmp()
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  worker->assign(&tmp1_s, re_s);
  worker->sqr(&tmp1_s);
  worker->assign(&tmp2_s, im_s);
  worker->sqr(&tmp2_s);
  worker->add(&tmp1_s, &tmp2_s);
  worker->add_double(&tmp1_s, -1);
  return &tmp1_s;
}

void complex::add(const complex *other)
{
  worker->add(re_s, other->re_s);
  worker->add(im_s, other->im_s);
}

void complex::mul(const complex *other)
{
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r1*r2-i1*i2
  //i:=r1*i2+i1*r2
  worker->assign(&tmp1_s, re_s);
  worker->mul(&tmp1_s, other->re_s);
  worker->assign(&tmp2_s, im_s);
  worker->mul(&tmp2_s, other->im_s);
  worker->sub(&tmp1_s, &tmp2_s);
  worker->assign(&tmp2_s, other->re_s);
  worker->mul(re_s, other->im_s);
  worker->mul(im_s, &tmp2_s);
  worker->add(im_s, re_s);
  worker->assign(re_s, &tmp1_s);
}

void complex::sqr()
{
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r*r-i*i
  //i=2*r*i
  worker->assign(&tmp1_s, im_s);
  worker->sqr(&tmp1_s);
  worker->mul(im_s, re_s);
  worker->lshift(im_s, 1);
  worker->sqr(re_s);
  worker->sub(re_s, &tmp1_s);
}

void complex::recip()
{
  getMagTmp();
  recip_prepared();
}

void complex::recip_prepared()
{ // 1/(re+i*im) = (re-i*im)/((re+i*im)*(re-i*im)) = (re-i*im)/(re*re+im*im)
  worker->recip(&tmp1_s);
  worker->mul(re_s, &tmp1_s);
  worker->chs(im_s);
  worker->mul(im_s, &tmp1_s);
}

void complex::sqrt()
{
  worker->assign(&tmp1_s, re_s);
  worker->assign(&tmp2_s, im_s);
  worker->sqr(&tmp1_s);
  worker->sqr(&tmp2_s);
  worker->add(&tmp1_s, &tmp2_s); //re*re+im*im
  worker->sqrt(&tmp1_s);
  if (!worker->isle0(re_s))
  {
    worker->add(&tmp1_s, re_s);
    worker->lshift(&tmp1_s, -1);
    worker->sqrt(&tmp1_s); //t1=sqrt((sqrt(re*re+im*im)+re)/2);
    worker->assign(re_s, &tmp1_s);
    /*if (t1==0)
      *res_im=0;
    else*/
    worker->lshift(&tmp1_s, 1);
    worker->recip(&tmp1_s);
    worker->mul(im_s, &tmp1_s); //im=im/(2*t1);
  }
  else
  {
    worker->sub(&tmp1_s, re_s);
    worker->lshift(&tmp1_s, -1);
    worker->sqrt(&tmp1_s); //t1=sqrt((sqrt(re*re+im*im)-re)/2);
    if (worker->isle0(&tmp1_s)) //t1==0
    {
      worker->zero(re_s);
      worker->zero(im_s);
    }
    else
    {
      worker->assign(re_s, im_s);
      worker->assign(im_s, &tmp1_s); //new im=t1
      worker->lshift(&tmp1_s, 1);
      worker->recip(&tmp1_s);
      worker->mul(re_s, &tmp1_s); //re=old im/(2*t1);
      if (worker->isle0(re_s)) //isl0 would be nicer here
      {
        worker->chs(re_s);
        worker->chs(im_s);
      };
    }
  };
}

const number_store *complex::mulreT(const complex *other) //Re(this*conjugate(other))
{ //re*ore+im*oim
  worker->assign(&tmp1_s, re_s);
  worker->assign(&tmp2_s, im_s);
  worker->mul(&tmp1_s, other->re_s);
  worker->mul(&tmp2_s, other->im_s);
  worker->add(&tmp1_s, &tmp2_s);
  return &tmp1_s;
}

const number_store *complex::dist2_tmp(const complex *other)
{
  worker->assign(&tmp1_s, re_s);
  worker->assign(&tmp2_s, im_s);
  worker->sub(&tmp1_s, other->re_s);
  worker->sub(&tmp2_s, other->im_s);
  worker->sqr(&tmp1_s);
  worker->sqr(&tmp2_s);
  worker->add(&tmp1_s, &tmp2_s);
  return &tmp1_s;
}

bool complex::isequal(const complex *other)
{
  return worker->isequal(re_s, other->re_s) &&
      worker->isequal(im_s, other->im_s);
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

#else //COMPLEX_IS_TEMPLATE

template <class NW>
const number_store *complex<NW>::getMagTmp()
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  worker.assign(&tmp1_s, re_s);
  worker.sqr(&tmp1_s);
  worker.assign(&tmp2_s, im_s);
  worker.sqr(&tmp2_s);
  worker.add(&tmp1_s, &tmp2_s);
  return &tmp1_s;
}

template <class NW>
void complex<NW>::add(const complex<NW> *other)
{
  assert(other->re_s);
  assert(other->im_s);
  worker.add(re_s, other->re_s);
  worker.add(im_s, other->im_s);
}

template <class NW>
void complex<NW>::mul(const complex<NW> *other)
{
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r1*r2-i1*i2
  //i:=r1*i2+i1*r2
  worker.assign(&tmp1_s, re_s);
  worker.mul(&tmp1_s, other->re_s);
  worker.assign(&tmp2_s, im_s);
  worker.mul(&tmp2_s, other->im_s);
  worker.sub(&tmp1_s, &tmp2_s);
  worker.assign(&tmp2_s, other->re_s);
  worker.mul(re_s, other->im_s);
  worker.mul(im_s, &tmp2_s);
  worker.add(im_s, re_s);
  worker.assign(re_s, &tmp1_s);
}

template <class NW>
void complex<NW>::sqr()
{
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r*r-i*i
  //i=2*r*i
  worker.assign(&tmp1_s, im_s);
  worker.sqr(&tmp1_s);
  worker.mul(im_s, re_s);
  worker.lshift_(im_s, 1);
  worker.sqr(re_s);
  worker.sub(re_s, &tmp1_s);
}
#endif //COMPLEX_IS_TEMPLATE

} // namespace MandelMath

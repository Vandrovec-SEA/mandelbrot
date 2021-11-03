#include "MandelMath.hpp"

#include <cassert>
#include <cmath>

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

number_store::number_store(): dbgType(number::Type::typeEmpty)
{
  static_assert(sizeof(dd_real)==2*sizeof(double), "dd_real suspiciously big");
  as.doubl=0.0;
}

number_store::~number_store()
{
  assert(dbgType==number::Type::typeEmpty);
}

void number_store::cleanup(number::Type ntype)
{
  switch (ntype)
  {
    case number::Type::typeDouble: cleanup_double_(); break;
    case number::Type::typeDDouble: cleanup_ddouble_(); break;
    case number::Type::typeMulti: cleanup_multi_(); break;
    case number::Type::typeEmpty: ;
  }
}

void number_store::cleanup_double_()
{
  assert(dbgType==number::Type::typeDouble);
  dbgType=number::Type::typeEmpty;
  as.doubl=0.0;
}

void number_store::cleanup_ddouble_()
{
  assert(dbgType==number::Type::typeDDouble);
  dbgType=number::Type::typeEmpty;
  as.ddouble_.deinit();
}

void number_store::cleanup_multi_()
{
  assert(dbgType==number::Type::typeMulti);
  dbgType=number::Type::typeEmpty;
  as.multi_.deinit();
}

void number_store::init(number::Type ntype, double val)
{
  switch (ntype)
  {
    case number::Type::typeDouble: init_double_(val); break;
    case number::Type::typeDDouble: init_ddouble_(val); break;
    case number::Type::typeMulti: init_multi_(val); break;
    case number::Type::typeEmpty: ;
  }
}

void number_store::init_double_(double val)
{
  assert(dbgType==number::Type::typeEmpty);
  dbgType=number::Type::typeDouble;
  as.doubl=val;
}

void number_store::init_ddouble_(double val)
{
  assert(dbgType==number::Type::typeEmpty);
  dbgType=number::Type::typeDDouble;
  as.ddouble_.init();
  as.ddouble_.dd->hi=val;
  as.ddouble_.dd->lo=0.0;
}

void number_store::init_multi_(double val)
{
  assert(dbgType==number::Type::typeEmpty);
  dbgType=number::Type::typeMulti;
  as.multi_.init();
  as.multi_.bytes->set(val);
}

void number_store::zero(number::Type ntype, double val)
{
  switch (ntype)
  {
    case number::Type::typeDouble: zero_double_(val); break;
    case number::Type::typeDDouble: zero_ddouble_(val); break;
    case number::Type::typeMulti: zero_multi_(val); break;
    case number::Type::typeEmpty: ;
  }
}

void number_store::zero_double_(double val)
{
  assert(dbgType==number::Type::typeDouble);
  //dbgType=number::Type::typeDouble;
  as.doubl=val;
}

void number_store::zero_ddouble_(double val)
{
  assert(dbgType==number::Type::typeDDouble);
  //dbgType=number::Type::typeDDouble;
  as.ddouble_.dd->hi=val;
  as.ddouble_.dd->lo=0.0;
}

void number_store::zero_multi_(double val)
{
  assert(dbgType==number::Type::typeMulti);
  //dbgType=number::Type::typeMulti;
  as.multi_.bytes->set(val);
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

void number_store::assign_double(const number_store &other)
{
  assert(other.dbgType==number::Type::typeDouble);
  assert(dbgType==number::Type::typeDouble);
  //dbgType=number::Type::typeDouble;
  as.doubl=other.as.doubl;
}

void number_store::assign_ddouble(const number_store &other)
{
  assert(other.dbgType==number::Type::typeDDouble);
  assert(dbgType==number::Type::typeDDouble);
  //dbgType=number::Type::typeDDouble;
  //as.ddouble_.dd->assign(*other.as.ddouble_.dd);
  *as.ddouble_.dd=*other.as.ddouble_.dd;
}

void number_store::assign_multi(const number_store &other)
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




QString number_double::toString()
{
  assert(store->dbgType==number::Type::typeDouble);
  return QString::number(store->as.doubl, 'f', 16);
}

void number_double::init(double val)
{
  store->init_double_(val);
}

void number_double::zero(double val)
{
  store->zero_double_(val);
}

void number_double::lshift_(int shoft)
{
  assert(store->dbgType==number::Type::typeDouble);
  store->as.doubl=ldexp(store->as.doubl, shoft);
  /*if (shoft>0)
    store->as.doubl*=two_pow_n(shoft);
  else if (shoft<0)
    store->as.doubl/=two_pow_n(-shoft);*/
}

void number_double::frac_pos()
{
  assert(store->dbgType==number::Type::typeDouble);
  store->as.doubl-=floor(store->as.doubl);
}

void number_double::add_double(double x)
{
  assert(store->dbgType==number::Type::typeDouble);
  store->as.doubl+=x;
}

void number_double::add(const number_store *other)
{
  assert(store->dbgType==number::Type::typeDouble);
  assert(other->dbgType==number::Type::typeDouble);
  store->as.doubl+=other->as.doubl;
}

void number_double::sub(const number_store *other)
{
  assert(store->dbgType==number::Type::typeDouble);
  assert(other->dbgType==number::Type::typeDouble);
  store->as.doubl-=other->as.doubl;
}

void number_double::rsub(const number_store *other)
{
  assert(store->dbgType==Type::typeDouble);
  assert(other->dbgType==Type::typeDouble);
  store->as.doubl=other->as.doubl-store->as.doubl;
}

void number_double::mul(const number_store *other)
{
  assert(store->dbgType==Type::typeDouble);
  assert(other->dbgType==Type::typeDouble);
  store->as.doubl*=other->as.doubl;
}

void number_double::sqr()
{
  assert(store->dbgType==Type::typeDouble);
  store->as.doubl*=store->as.doubl;
}

int number_double::toRound()
{
  assert(store->dbgType==Type::typeDouble);
  return qRound(store->as.doubl);
}

double number_double::toDouble()
{
  assert(store->dbgType==Type::typeDouble);
  return store->as.doubl;
}



QString number_ddouble::toString()
{
  assert(store->dbgType==Type::typeDDouble);
  return QString("dd(%1,%2)").arg(store->as.ddouble_.dd->hi).arg(store->as.ddouble_.dd->lo);
}

void number_ddouble::init(double val)
{
  store->init_ddouble_(val);
}

void number_ddouble::zero(double val)
{
  store->zero_ddouble_(val);
}

void number_ddouble::lshift_(int shoft)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->lshift(shoft);
}

void number_ddouble::frac_pos()
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->hi-=floor(store->as.ddouble_.dd->hi);
  store->as.ddouble_.dd->lo-=floor(store->as.ddouble_.dd->lo);
}

void number_ddouble::add_double(double x)
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->add_double(x);
}

void number_ddouble::add(const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->add(other->as.ddouble_.dd->hi, other->as.ddouble_.dd->lo);
}

void number_ddouble::sub(const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->add(-other->as.ddouble_.dd->hi, -other->as.ddouble_.dd->lo);
}
void number_ddouble::rsub(const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->chs();
  store->as.ddouble_.dd->add(other->as.ddouble_.dd->hi, other->as.ddouble_.dd->lo);
}

void number_ddouble::mul(const number_store *other)
{
  assert(store->dbgType==Type::typeDDouble);
  assert(other->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->mul(other->as.ddouble_.dd->hi, other->as.ddouble_.dd->lo);
}

void number_ddouble::sqr()
{
  assert(store->dbgType==Type::typeDDouble);
  store->as.ddouble_.dd->sqr();
}

int number_ddouble::toRound()
{
  assert(store->dbgType==Type::typeDDouble);
  return floor(store->as.ddouble_.dd->hi+0.5)+floor(store->as.ddouble_.dd->lo+0.5);
}

double number_ddouble::toDouble()
{
  assert(store->dbgType==Type::typeDDouble);
  return store->as.ddouble_.dd->hi;
}




QString number_multi::toString()
{
  assert(store->dbgType==Type::typeMulti);
  return QString("multi");
}

void number_multi::init(double val)
{
  store->init_multi_(val);
}

void number_multi::zero(double val)
{
  store->zero_multi_(val);
}

void number_multi::lshift_(int shoft)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->lshift(shoft);
}

void number_multi::frac_pos()
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->frac_pos();
}

void number_multi::add_double(double x)
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->add_double(x);
}

void number_multi::add(const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->add(other->as.multi_.bytes);
}

void number_multi::sub(const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->chs();
  store->as.multi_.bytes->add(other->as.multi_.bytes);
  store->as.multi_.bytes->chs();
}
void number_multi::rsub(const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->chs();
  store->as.multi_.bytes->add(other->as.multi_.bytes);
}

void number_multi::mul(const number_store *other)
{
  assert(store->dbgType==Type::typeMulti);
  assert(other->dbgType==Type::typeMulti);
  store->as.multi_.bytes->mul(other->as.multi_.bytes);
}

void number_multi::sqr()
{
  assert(store->dbgType==Type::typeMulti);
  store->as.multi_.bytes->sqr();
}

int number_multi::toRound()
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->round();
}

double number_multi::toDouble()
{
  assert(store->dbgType==Type::typeMulti);
  return store->as.multi_.bytes->toDouble();
}




template <class T>
T *complex<T>::getMagTmp()
{
  tmp1.assign(re.store);
  tmp1.sqr();
  tmp2.assign(im.store);
  tmp2.sqr();
  tmp1.add(tmp2.store);
  return &tmp1;
}

template <class T>
void complex<T>::add(const complex<T> *other)
{
  re.add(other->re.store);
  im.add(other->im.store);
}

template <class T>
void complex<T>::mul(const complex<T> *other)
{
  if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
    dbgPoint();
  assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r1*r2-i1*i2
  //i:=r1*i2+i1*r2
  tmp1.assign(re.store);
  tmp1.mul(other->re.store);
  tmp2.assign(im.store);
  tmp2.mul(other->im.store);
  tmp1.sub(tmp2.store);
  tmp2.assign(other->re.store);
  re.mul(other->im.store);
  im.mul(tmp2.store);
  im.add(re.store);
  re.assign(tmp1.store);
}

template <class T>
void complex<T>::sqr()
{
  if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
    dbgPoint();
  assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r*r-i*i
  //i=2*r*i
  tmp1.assign(im.store);
  tmp1.sqr();
  im.mul(re.store);
  im.lshift_(1);
  re.sqr();
  re.sub(tmp1.store);
}

} // namespace MandelMath

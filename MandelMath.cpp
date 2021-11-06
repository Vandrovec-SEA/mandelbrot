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

void number_worker_double::assign(number_store *store, const number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeDouble);
  assert(src->dbgType==Type::typeDouble);
  store->as.doubl=src->as.doubl;
}

void number_worker_double::lshift_(number_store *store, int shoft)
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

void number_worker_ddouble::assign(number_store *store, const number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeDouble);
  assert(src->dbgType==Type::typeDouble);
  *store->as.ddouble_.dd=*src->as.ddouble_.dd;
}

void number_worker_ddouble::lshift_(number_store *store, int shoft)
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

void number_worker_multi::assign(number_store *store, const number_store *src)
{
  assert(store);
  assert(src);
  assert(store->dbgType==Type::typeDouble);
  assert(src->dbgType==Type::typeDouble);
  *store->as.multi_.bytes=*src->as.multi_.bytes;
}

void number_worker_multi::lshift_(number_store *store, int shoft)
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
  worker->lshift_(im_s, 1);
  worker->sqr(re_s);
  worker->sub(re_s, &tmp1_s);
}

} // namespace MandelMath

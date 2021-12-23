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

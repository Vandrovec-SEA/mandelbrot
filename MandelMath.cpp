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

double two_pow_n(unsigned int n)
{
  assert(n<4096);
  double result=1;
  while (n>=31)
  {
    result*=(1<<31);
    n-=31;
  }
  if (n>0)
    result*=(1<<n);
  return result;
}

number_store::number_store(): dbgType(DbgType::typeEmpty)
{
  static_assert(sizeof(double_double)==2*sizeof(double), "double_double suspiciously big");
  as.doubl=0.0;
}

number_store::~number_store()
{
  assert(dbgType==number_store::DbgType::typeEmpty);
}

void number_store::cleanup_double()
{
  assert(dbgType==number_store::DbgType::typeDouble);
  dbgType=DbgType::typeEmpty;
  as.doubl=0.0;
}

void number_store::cleanup_ddouble()
{
  assert(dbgType==number_store::DbgType::typeDDouble);
  dbgType=DbgType::typeEmpty;
  as.ddouble.hi=0.0;
  as.ddouble.lo=0.0;
}

void number_store::cleanup_multi()
{
  assert(dbgType==number_store::DbgType::typeMulti);
  assert(as.multi._filler==nullptr);
  dbgType=DbgType::typeEmpty;
  delete as.multi.bytes;
  as.multi.bytes=nullptr;
  as.multi._filler=nullptr;
}

void number_store::init_double(double val)
{
  assert(dbgType==number_store::DbgType::typeEmpty);
  dbgType=DbgType::typeDouble;
  as.doubl=val;
}

void number_store::init_ddouble(double val)
{
  assert(dbgType==number_store::DbgType::typeEmpty);
  dbgType=DbgType::typeDDouble;
  as.ddouble.hi=val;
  as.ddouble.lo=0.0;
}

void number_store::init_multi(double val)
{
  assert(dbgType==number_store::DbgType::typeEmpty);
  dbgType=DbgType::typeMulti;
  as.multi.bytes->set(val);
}

void number_store::zero_double(double val)
{
  assert(dbgType==number_store::DbgType::typeDouble);
  dbgType=DbgType::typeDouble;
  as.doubl=val;
}

void number_store::zero_ddouble(double val)
{
  assert(dbgType==number_store::DbgType::typeDDouble);
  dbgType=DbgType::typeDDouble;
  as.ddouble.hi=val;
  as.ddouble.lo=0.0;
}

void number_store::zero_multi(double val)
{
  assert(dbgType==number_store::DbgType::typeMulti);
  dbgType=DbgType::typeMulti;
  as.multi.bytes->set(val);
}

void number_store::assign_double(const number_store &other)
{
  assert(other.dbgType==DbgType::typeDouble);
  assert(dbgType==DbgType::typeDouble);
  dbgType=DbgType::typeDouble;
  as.doubl=other.as.doubl;
}

void number_store::assign_ddouble(const number_store &other)
{
  assert(other.dbgType==DbgType::typeDDouble);
  assert(dbgType==DbgType::typeDDouble);
  dbgType=DbgType::typeDDouble;
  as.ddouble=other.as.ddouble;
}

void number_store::assign_multi(const number_store &other)
{
  assert(other.dbgType==DbgType::typeMulti);
  assert(dbgType==DbgType::typeMulti);
  dbgType=DbgType::typeMulti;
  *as.multi.bytes=*other.as.multi.bytes;
}




QString number_double::toString()
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  return QString::number(store->as.doubl, 'f', 16);
}

void number_double::init(double val)
{
  store->init_double(val);
}

void number_double::zero(double val)
{
  store->zero_double(val);
}

void number_double::lshift_(int shoft)
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  if (shoft>0)
    store->as.doubl*=two_pow_n(shoft);
  else if (shoft<0)
    store->as.doubl/=two_pow_n(-shoft);
}

void number_double::frac_pos()
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  store->as.doubl-=floor(store->as.doubl);
}

void number_double::add_double(double x)
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  store->as.doubl+=x;
}

void number_double::add(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  assert(other->dbgType==number_store::DbgType::typeDouble);
  store->as.doubl+=other->as.doubl;
}

void number_double::sub(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  assert(other->dbgType==number_store::DbgType::typeDouble);
  store->as.doubl-=other->as.doubl;
}

void number_double::rsub(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  assert(other->dbgType==number_store::DbgType::typeDouble);
  store->as.doubl=other->as.doubl-store->as.doubl;
}

void number_double::mul(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  assert(other->dbgType==number_store::DbgType::typeDouble);
  store->as.doubl*=other->as.doubl;
}

int number_double::toRound()
{
  assert(store->dbgType==number_store::DbgType::typeDouble);
  return qRound(store->as.doubl);
}



QString number_ddouble::toString()
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  return QString("dd(%1,%2)").arg(store->as.ddouble.hi).arg(store->as.ddouble.lo);
}

void number_ddouble::init(double val)
{
  store->init_ddouble(val);
}

void number_ddouble::zero(double val)
{
  store->zero_ddouble(val);
}

void number_ddouble::lshift_(int shoft)
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  if (shoft>0)
  {
    double coeff=two_pow_n(shoft);
    store->as.ddouble.hi*=coeff;
    store->as.ddouble.lo*=coeff;
  }
  else if (shoft<0)
  {
    double coeff=two_pow_n(shoft);
    store->as.ddouble.hi/=coeff;
    store->as.ddouble.lo/=coeff;
  }
}

void number_ddouble::frac_pos()
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  store->as.ddouble.hi-=floor(store->as.ddouble.hi);
  store->as.ddouble.lo-=floor(store->as.ddouble.lo);
}

void number_ddouble::add_double(double x)
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  store->as.ddouble.add_double(x);
}

void number_ddouble::add(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  assert(other->dbgType==number_store::DbgType::typeDDouble);
  store->as.ddouble.add(other->as.ddouble.hi, other->as.ddouble.lo);
}

void number_ddouble::sub(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  assert(other->dbgType==number_store::DbgType::typeDDouble);
  store->as.ddouble.add(-other->as.ddouble.hi, -other->as.ddouble.lo);
}
void number_ddouble::rsub(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  assert(other->dbgType==number_store::DbgType::typeDDouble);
  store->as.ddouble.chs();
  store->as.ddouble.add(other->as.ddouble.hi, other->as.ddouble.lo);
}

void number_ddouble::mul(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  assert(other->dbgType==number_store::DbgType::typeDDouble);
  store->as.ddouble.mul(other->as.ddouble.hi, other->as.ddouble.lo);
}

int number_ddouble::toRound()
{
  assert(store->dbgType==number_store::DbgType::typeDDouble);
  return floor(store->as.ddouble.hi+0.5)+floor(store->as.ddouble.lo+0.5);
}




QString number_multi::toString()
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  return QString("multi");
}

void number_multi::init(double val)
{
  store->init_multi(val);
}

void number_multi::zero(double val)
{
  store->zero_multi(val);
}

void number_multi::lshift_(int shoft)
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  if (shoft>0)
  {
    double coeff=two_pow_n(shoft);
    store->as.multi.bytes->mul_double(coeff);
  }
  else if (shoft<0)
  {
    double coeff=two_pow_n(shoft);
    store->as.multi.bytes->mul_double(1/coeff);
  }
}

void number_multi::frac_pos()
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  store->as.multi.bytes->frac_pos();
}

void number_multi::add_double(double x)
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  store->as.multi.bytes->add_double(x);
}

void number_multi::add(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  assert(other->dbgType==number_store::DbgType::typeMulti);
  store->as.multi.bytes->add(other->as.multi.bytes);
}

void number_multi::sub(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  assert(other->dbgType==number_store::DbgType::typeMulti);
  store->as.multi.bytes->chs();
  store->as.multi.bytes->add(other->as.multi.bytes);
  store->as.multi.bytes->chs();
}
void number_multi::rsub(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  assert(other->dbgType==number_store::DbgType::typeMulti);
  store->as.multi.bytes->chs();
  store->as.multi.bytes->add(other->as.multi.bytes);
}

void number_multi::mul(const number_store *other)
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  assert(other->dbgType==number_store::DbgType::typeMulti);
  store->as.multi.bytes->mul(other->as.multi.bytes);
}

int number_multi::toRound()
{
  assert(store->dbgType==number_store::DbgType::typeMulti);
  return store->as.multi.bytes->round();
}





void complex_double::add(const complex *other)
{
  assert(re->dbgType==number_store::DbgType::typeDouble);
  assert(im->dbgType==number_store::DbgType::typeDouble);
  assert(other->re->dbgType==number_store::DbgType::typeDouble);
  assert(other->im->dbgType==number_store::DbgType::typeDouble);
  re->as.doubl+=other->re->as.doubl;
  im->as.doubl+=other->im->as.doubl;
}

void complex_double::mul(const complex *other)
{
  assert(re->dbgType==number_store::DbgType::typeDouble);
  assert(im->dbgType==number_store::DbgType::typeDouble);
  assert(other->re->dbgType==number_store::DbgType::typeDouble);
  assert(other->im->dbgType==number_store::DbgType::typeDouble);
  double tmp=re->as.doubl*other->im->as.doubl + im->as.doubl*other->re->as.doubl;
  im->as.doubl=re->as.doubl*other->re->as.doubl - im->as.doubl*other->im->as.doubl;
  re->as.doubl=tmp;
}

} // namespace MandelMath

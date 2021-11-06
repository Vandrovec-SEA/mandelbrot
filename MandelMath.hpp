#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <QString>

#include "double_double.hpp"
#include "multiprec.hpp"

void dbgPoint();

namespace MandelMath {

//double two_pow_n(unsigned int n);

struct number_store;

class number
{
public:
  enum Type { typeEmpty, typeDouble, typeDDouble, typeMulti };

  number(number_store *store): store(store) { }
  number_store *store;
  virtual QString toString()=0;
  virtual void init(double val=0)=0;
  virtual void zero(double val=0)=0;
  virtual void assign(const number_store *src)=0;
  virtual void assignTo(number_store *src)=0;
  virtual void cleanup()=0;
  virtual void lshift_(int shoft)=0; // self <<= shoft
  virtual void frac_pos()=0; //0<=result<1
  virtual void add_double(double x)=0;
  virtual void add(const number_store *other)=0;
  virtual void sub(const number_store *other)=0;
  virtual void rsub(const number_store *other)=0;
  virtual void mul(const number_store *other)=0;
  virtual void sqr()=0;

  virtual int toRound()=0;
  virtual double toDouble()=0;
};

struct number_store
{
protected:
  class dd_real_managed
  {
  public:
    dd_real_managed(): dd(nullptr)
    {
      static_assert (sizeof(dd_real_managed)<=sizeof(double), "Try to keep multiprec itself small");
      dbgPoint(); //should not actually be constructed or destructed since it lives within the union number_store::As
    }
    void init() { dd=new dd_real(); }
    void deinit() { delete dd; dd=nullptr; }
    ~dd_real_managed() { dbgPoint(); delete dd; dd=nullptr; }
    dd_real *dd;
  };

  class multiprec_managed
  {
  public:
    multiprec_managed(): bytes(nullptr)
    {
      static_assert (sizeof(multiprec_managed)<=sizeof(double), "Try to keep multiprec itself small");
      dbgPoint(); //should not actually be constructed or destructed since it lives within the union number_store::As
    }
    void init() { bytes=new multiprec(); }
    void deinit() { delete bytes; bytes=nullptr; }
    ~multiprec_managed() { dbgPoint(); delete bytes; bytes=nullptr; }
    multiprec *bytes;
  };

public:
  number::Type dbgType;
  number_store();
  ~number_store();
  void cleanup(number::Type ntype);
  void init(number::Type ntype, double val=0);
  void zero(number::Type ntype, double val=0);
  template <class T> void assign(const number_store &other);
  void assignTo_double(number_store &other);
  void assignTo_ddouble(number_store &other);
  void assignTo_multi(number_store &other);

  union As
  {
    As() {}
    ~As() {}
    double doubl;
    dd_real_managed ddouble_;
    multiprec_managed multi_;
  } as;
};

class number_double: public number
{
public:
  number_double(number_store *store): number(store)
  { assert((store->dbgType==Type::typeDouble) ||
           (store->dbgType==Type::typeEmpty)); }
  QString toString() override;
  void init(double val=0) override;
  void zero(double val=0) override;
  void assign(const number_store *src) override { store->assign<number_double>(*src); };
  void assignTo(number_store *src) override { store->assignTo_double(*src); };
  void cleanup() override { store->cleanup(Type::typeDouble); }
  void lshift_(int shoft) override;
  void frac_pos() override;
  void add_double(double x) override;
  void add(const number_store *other) override;
  void sub(const number_store *other) override;
  void rsub(const number_store *other) override;
  void mul(const number_store *other) override;
  void sqr() override;

  int toRound() override;
  double toDouble() override;
};

class number_ddouble: public number
{
public:
  number_ddouble(number_store *store): number(store)
  { assert((store->dbgType==Type::typeDDouble) ||
           (store->dbgType==Type::typeEmpty)); }
  QString toString() override;
  void init(double val=0) override;
  void zero(double val=0) override;
  void assign(const number_store *src) override { store->assign<number_ddouble>(*src); };
  void assignTo(number_store *src) override { store->assignTo_ddouble(*src); };
  void cleanup() override { store->cleanup(Type::typeDDouble); }
  void lshift_(int shoft) override;
  void frac_pos() override;
  void add_double(double x) override;
  void add(const number_store *other) override;
  void sub(const number_store *other) override;
  void rsub(const number_store *other) override;
  void mul(const number_store *other) override;
  void sqr() override;

  int toRound() override;
  double toDouble() override;
};

class number_multi: public number
{
public:
  number_multi(number_store *store): number(store)
  { assert((store->dbgType==Type::typeMulti) ||
           (store->dbgType==Type::typeEmpty)); }
  QString toString() override;
  void init(double val=0) override;
  void zero(double val=0) override;
  void assign(const number_store *src) override { store->assign<number_multi>(*src); };;
  void assignTo(number_store *src) override { store->assignTo_multi(*src); };;
  void cleanup() override { store->cleanup(Type::typeMulti); }
  void lshift_(int shoft) override;
  void frac_pos() override;
  void add_double(double x) override;
  void add(const number_store *other) override;
  void sub(const number_store *other) override;
  void rsub(const number_store *other) override;
  void mul(const number_store *other) override;
  void sqr() override;

  int toRound() override;
  double toDouble() override;
};

template <class T> struct number_to_type
{
public:
  static const number::Type ntype;
};

template <> struct number_to_type<number_double>
{ static const number::Type ntype=number::Type::typeDouble; };
template <> struct number_to_type<number_ddouble>
{ static const number::Type ntype=number::Type::typeDDouble; };
template <> struct number_to_type<number_multi>
{ static const number::Type ntype=number::Type::typeMulti; };

template <class T>
class complex
{
  bool external_stores;
  number_store tmp1_s;
  number_store tmp2_s;
  T tmp1;
  T tmp2;
public:
  //complex(number *re, number *im, number *tmp1, number *tmp2): tmp1(tmp1), tmp2(tmp2), re(re), im(im) { }
  complex(number_store *re, number_store *im, bool external_stores=false):
    external_stores(external_stores), tmp1(&tmp1_s), tmp2(&tmp2_s), re(re), im(im)
  {
    tmp1.init(number_to_type<T>::ntype);
    tmp2.init(number_to_type<T>::ntype);
  }
  ~complex()
  {
    tmp2.cleanup();
    tmp1.cleanup();
    if (!external_stores)
    {
      im.cleanup();
      re.cleanup();
    };
  }
  T re;
  T im;
  T *getMagTmp();
  void add(const complex *other);
  void mul(const complex *other);
  void sqr();
};

template class complex<number_double>;
template class complex<number_ddouble>;
template class complex<number_multi>;

} // namespace MandelMath

#endif // MANDELMATH_NUMBER_HPP

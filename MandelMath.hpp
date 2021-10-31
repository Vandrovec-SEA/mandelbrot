#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <QString>

#include "double_double.hpp"
#include "multiprec.hpp"

void dbgPoint();

namespace MandelMath {

double two_pow_n(unsigned int n);

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
  virtual void cleanup()=0;
  virtual void lshift_(int shoft)=0; // self <<= shoft
  virtual void frac_pos()=0; //0<=result<1
  virtual void add_double(double x)=0;
  virtual void add(const number_store *other)=0;
  virtual void sub(const number_store *other)=0;
  virtual void rsub(const number_store *other)=0;
  virtual void mul(const number_store *other)=0;

  virtual int toRound()=0;
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
  void cleanup_double_();
  void cleanup_ddouble_();
  void cleanup_multi_();
  void init(number::Type ntype, double val=0);
  void init_double_(double val=0); //should already be switched to empty
  void init_ddouble_(double val=0); //should already be switched to empty
  void init_multi_(double val=0); //should already be switched to empty
  void zero(number::Type ntype, double val=0);
  void zero_double_(double val=0); //should already be switched to double
  void zero_ddouble_(double val=0); //should already be switched to ddouble
  void zero_multi_(double val=0); //should already be switched to multi
  void assign_double(const number_store &other);
  void assign_ddouble(const number_store &other);
  void assign_multi(const number_store &other);

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
  void assign(const number_store *src) override { store->assign_double(*src); };
  void cleanup() override { store->cleanup_double_(); }
  void lshift_(int shoft) override;
  void frac_pos() override;
  void add_double(double x) override;
  void add(const number_store *other) override;
  void sub(const number_store *other) override;
  void rsub(const number_store *other) override;
  void mul(const number_store *other) override;

  int toRound() override;
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
  void assign(const number_store *src) override { store->assign_ddouble(*src); };
  void cleanup() override { store->cleanup_ddouble_(); }
  void lshift_(int shoft) override;
  void frac_pos() override;
  void add_double(double x) override;
  void add(const number_store *other) override;
  void sub(const number_store *other) override;
  void rsub(const number_store *other) override;
  void mul(const number_store *other) override;

  int toRound() override;
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
  void assign(const number_store *src) override { store->assign_multi(*src); };;
  void cleanup() override { store->cleanup_multi_(); }
  void lshift_(int shoft) override;
  void frac_pos() override;
  void add_double(double x) override;
  void add(const number_store *other) override;
  void sub(const number_store *other) override;
  void rsub(const number_store *other) override;
  void mul(const number_store *other) override;

  int toRound() override;
};

class complex
{
public:
  complex(number_store *re, number_store *im): re(re), im(im) { }
  number_store *re;
  number_store *im;
  virtual void add(const complex *other)=0;
  virtual void mul(const complex *other)=0;
};

class complex_double: public complex
{
public:
  complex_double(number_store *re, number_store *im): complex(re, im) { }
  void add(const complex *other) override;
  void mul(const complex *other) override;
};

//complex_multi might have some temporaries to prevent malloc all the time

} // namespace MandelMath

#endif // MANDELMATH_NUMBER_HPP

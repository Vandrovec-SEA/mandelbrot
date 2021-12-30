#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <QString>

#include "double_double.hpp"
#include "multiprec.hpp"

void dbgPoint();

namespace MandelMath {

//double two_pow_n(unsigned int n);
int gcd(int m, int n);
int ctz16(int x);

struct number_store;

class number_worker
{
public:
  enum Type { typeEmpty, typeDouble, typeDDouble, typeMulti };
  virtual double eps2() { return 1.23e-32; }

  number_worker() { }
  virtual Type ntype() { return typeEmpty; }
  virtual void init(number_store *store, double val=0)=0;
  virtual void zero(number_store *store, double val=0)=0;
  virtual void swap(number_store *store, number_store *src)=0; //dst uninit + src init -> dst init + src uninit
  virtual void assign(number_store *store, const number_store *src)=0;
  //virtual void assignTo(number_store *src)=0;
  virtual void cleanup(number_store *store)=0;
  virtual void chs(number_store *store)=0;
  virtual void lshift(number_store *store, int shoft)=0; // self <<= shoft; 1 lshift -10000 = 0 not error
  virtual void frac_pos(number_store *store)=0; //0<=result<1
  virtual void add_double(number_store *store, double x)=0;
  virtual void add(number_store *store, const number_store *other)=0;
  virtual void sub(number_store *store, const number_store *other)=0;
  virtual void rsub(number_store *store, const number_store *other)=0;
  virtual void mul(number_store *store, const number_store *other)=0;
  virtual void sqr(number_store *store)=0;
  virtual void recip(number_store *store)=0;
  virtual void sqrt(number_store *store)=0;
  virtual bool isequal(const number_store *store, const number_store *other)=0; //return store==other
  virtual bool is0(const number_store *store)=0;
  virtual bool isle(const number_store *store, const number_store *other)=0; //return store<=other
  virtual bool isle0(const number_store *store)=0; //return store<=0
  virtual bool isl1(const number_store *store)=0; //return store<1

  virtual QString toString(const number_store *store)=0;
  virtual int toRound(const number_store *store)=0;
  virtual double toDouble(const number_store *store)=0;
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
  number_worker::Type dbgType;
  number_store();
  ~number_store();
  void cleanup(number_worker::Type ntype);
  void init(number_worker::Type ntype, double val=0);
  void zero(number_worker::Type ntype, double val=0);
  //template <class T> void assign(const number_store &other);
  //void assignTo_double(number_store &other);
  //void assignTo_ddouble(number_store &other);
  //void assignTo_multi(number_store &other);

  union As
  {
    As() {}
    ~As() {}
    double doubl;
    dd_real_managed ddouble_;
    multiprec_managed multi_;
  } as;
};

class number_worker_double: public number_worker
{
public:
  double eps2() override { return 1.23e-32;  /* 2^-(2*53) */ }

  number_worker_double() {}
  //{ assert((store->dbgType==Type::typeDouble) ||
  //         (store->dbgType==Type::typeEmpty)); }
  virtual Type ntype() override { return typeDouble; }
  void init(number_store *store, double val=0) override;
  void zero(number_store *store, double val=0) override;
  void swap(number_store *store, number_store *src) override; //dst uninit + src init -> dst init + src uninit
  void assign(number_store *store, const number_store *src) override;// { store->assign<number_double>(*src); };
  //void assignTo(number_store *src) override { store->assignTo_double(*src); };
  void cleanup(number_store *store) override { store->cleanup(Type::typeDouble); }
  void chs(number_store *store) override;
  void lshift(number_store *store, int shoft) override;
  void frac_pos(number_store *store) override;
  void add_double(number_store *store, double x) override;
  void add(number_store *store, const number_store *other) override;
  void sub(number_store *store, const number_store *other) override;
  void rsub(number_store *store, const number_store *other) override;
  void mul(number_store *store, const number_store *other) override;
  void sqr(number_store *store) override;
  void recip(number_store *store) override;
  void sqrt(number_store *store) override;
  bool isequal(const number_store *store, const number_store *other) override;
  bool is0(const number_store *store) override;
  bool isle(const number_store *store, const number_store *other) override;
  bool isle0(const number_store *store) override;
  bool isl1(const number_store *store) override;

  QString toString(const number_store *store) override;
  int toRound(const number_store *store) override;
  double toDouble(const number_store *store) override;
};

class number_worker_ddouble: public number_worker
{
public:
  double eps2() override { return 6.1e-64; /* 2^-(2*(53+52)) */ }
  number_worker_ddouble() {}
  //{ assert((store->dbgType==Type::typeDDouble) ||
  //         (store->dbgType==Type::typeEmpty)); }
  virtual Type ntype() override { return typeDDouble; }
  void init(number_store *store, double val=0) override;
  void zero(number_store *store, double val=0) override;
  void swap(number_store *store, number_store *src) override; //dst uninit + src init -> dst init + src uninit
  void assign(number_store *store, const number_store *src) override;// { store->assign<number_ddouble>(*src); };
  //void assignTo(number_store *src) override { store->assignTo_ddouble(*src); };
  void cleanup(number_store *store) override { store->cleanup(Type::typeDDouble); }
  void chs(number_store *store) override;
  void lshift(number_store *store, int shoft) override;
  void frac_pos(number_store *store) override;
  void add_double(number_store *store, double x) override;
  void add(number_store *store, const number_store *other) override;
  void sub(number_store *store, const number_store *other) override;
  void rsub(number_store *store, const number_store *other) override;
  void mul(number_store *store, const number_store *other) override;
  void sqr(number_store *store) override;
  void recip(number_store *store) override;
  void sqrt(number_store *store) override;
  bool isequal(const number_store *store, const number_store *other) override;
  bool is0(const number_store *store) override;
  bool isle(const number_store *store, const number_store *other) override;
  bool isle0(const number_store *store) override;
  bool isl1(const number_store *store) override;

  QString toString(const number_store *store) override;
  int toRound(const number_store *store) override;
  double toDouble(const number_store *store) override;
};

class number_worker_multi: public number_worker
{
public:
  double eps2() override { return 6.1e-64; /* 2^-(2*(53+52)) */ }
  number_worker_multi() {}
  //{ assert((store->dbgType==Type::typeMulti) ||
  //         (store->dbgType==Type::typeEmpty)); }
  virtual Type ntype() override { return typeMulti; }
  void init(number_store *store, double val=0) override;
  void zero(number_store *store, double val=0) override;
  void swap(number_store *store, number_store *src) override; //dst uninit + src init -> dst init + src uninit
  void assign(number_store *store, const number_store *src) override;// { store->assign<number_multi>(*src); };;
  //void assignTo(number_store *src) override { store->assignTo_multi(*src); };;
  void cleanup(number_store *store) override { store->cleanup(Type::typeMulti); }
  void chs(number_store *store) override;
  void lshift(number_store *store, int shoft) override;
  void frac_pos(number_store *store) override;
  void add_double(number_store *store, double x) override;
  void add(number_store *store, const number_store *other) override;
  void sub(number_store *store, const number_store *other) override;
  void rsub(number_store *store, const number_store *other) override;
  void mul(number_store *store, const number_store *other) override;
  void sqr(number_store *store) override;
  void recip(number_store *store) override;
  void sqrt(number_store *store) override;
  bool isequal(const number_store *store, const number_store *other) override;
  bool is0(const number_store *store) override;
  bool isle(const number_store *store, const number_store *other) override;
  bool isle0(const number_store *store) override;
  bool isl1(const number_store *store) override;

  QString toString(const number_store *store) override;
  int toRound(const number_store *store) override;
  double toDouble(const number_store *store) override;
};
/*
template <class T> struct number_to_type
{
public:
  static const number_worker::Type ntype;
};

template <> struct number_to_type<number_worker_double>
{ static const number_worker::Type ntype=number_worker::Type::typeDouble; };
template <> struct number_to_type<number_worker_ddouble>
{ static const number_worker::Type ntype=number_worker::Type::typeDDouble; };
template <> struct number_to_type<number_worker_multi>
{ static const number_worker::Type ntype=number_worker::Type::typeMulti; };
*/

void complex_double_sqrt(double *res_re, double *res_im, double in_re, double in_im); //res_re>=0
void complex_double_quadratic(double *res_re, double *res_im,
                              double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);
void complex_double_quadratic2(double *res1_re, double *res1_im,
                               double *res2_re, double *res2_im,
                               double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);

#define COMPLEX_IS_TEMPLATE 0 //virtual or templated class complex
#if !COMPLEX_IS_TEMPLATE
class complex
{
  number_worker *worker;
  bool external_stores;
  number_store tmp1_s;
  number_store tmp2_s;
public:
  //complex(number *re, number *im, number *tmp1, number *tmp2): tmp1(tmp1), tmp2(tmp2), re(re), im(im) { }
  complex(number_worker *worker, number_store *re_s, number_store *im_s, bool external_stores=false):
    worker(worker), external_stores(external_stores), re_s(re_s), im_s(im_s)
  {
    worker->init(&tmp1_s);
    worker->init(&tmp2_s);
  }
  ~complex()
  {
    worker->cleanup(&tmp2_s);
    worker->cleanup(&tmp1_s);
    if (!external_stores)
    {
      worker->cleanup(im_s);
      worker->cleanup(re_s);
    };
  }
  number_store *re_s;
  number_store *im_s;
  //could use assign(const complex *)
  const number_store *getMagTmp();
  //could use double toMagDouble();, assuming toDouble is cheaper than mul
  //could use lshift
  void add(const complex *other);
  //could use sub(const complex *)
  void mul(const complex *other);
  void sqr();
  void recip();
  void recip_prepared();
  void sqrt();
  const number_store *mulreT(const complex *other); //Re(this*conjugate(other)) = re*o->re+im*o->im
  const number_store *dist2_tmp(const complex *other);
  bool isequal(const complex *other);
};

#else //COMPLEX_IS_TEMPLATE
template <class NW>
class complex
{
  NW worker;
  bool external_stores;
  number_store tmp1_s;
  number_store tmp2_s;
public:
  //complex(number *re, number *im, number *tmp1, number *tmp2): tmp1(tmp1), tmp2(tmp2), re(re), im(im) { }
  complex(number_store *re_s, number_store *im_s, bool external_stores=false):
    external_stores(external_stores), re_s(re_s), im_s(im_s)
  {
    worker.init(&tmp1_s);
    worker.init(&tmp2_s);
  }
  ~complex()
  {
    worker.cleanup(&tmp2_s);
    worker.cleanup(&tmp1_s);
    if (!external_stores)
    {
      worker.cleanup(im_s);
      worker.cleanup(re_s);
    };
  }
  number_store *re_s;
  number_store *im_s;
  const number_store *getMagTmp();
  void add(const complex *other);
  void mul(const complex *other);
  void sqr();
};

template class complex<number_worker_double>;
template class complex<number_worker_ddouble>;
template class complex<number_worker_multi>;
#endif //COMPLEX_IS_TEMPLATE

} // namespace MandelMath
#endif // MANDELMATH_NUMBER_HPP

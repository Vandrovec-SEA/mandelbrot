#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <QString>

#include "double_double.hpp"
#include "multiprec.hpp"

#define NUMBER_DOUBLE_EXISTS 1
#define ONLY_DOUBLE_WORKER 1

#ifndef M_PI
#define M_PI 3.14159265358979323846
//constexpr double M_PI=4*atan2(1., 0);//acos(-1);
#endif

void dbgPoint();

namespace MandelMath {

int gcd(int m, int n);
int ctz16(int x);

/*struct number_store;

union number_place //might have to switch to struct to allow inplace promote()
{
  dd_real dd;
  multiprec multi;
  number_place() { memset(this, 0x00, sizeof(number_place)); }
};*/

typedef dd_real dq_real;

union number_pointer
{
  double *asf64;
  __float128 *asf128;
  dd_real *asdd;
  dq_real *asdq;
  number_pointer(): asf64(nullptr) {}
  number_pointer(double *asf64): asf64(asf64) {}
};

union number_pointer_c
{
  const double *asf64;
  const __float128 *asf128;
  const dd_real *asdd;
  const dq_real *asdq;
  number_pointer_c(): asf64(nullptr) {}
  number_pointer_c(double *asf64): asf64(asf64) {}
  number_pointer_c(const number_pointer &src): asf64(src.asf64) {}
};

class complex; //forward declaration for friend

class upgrademe
{

};

class worker_multi
{
public:
  enum Type { typeEmpty
#if NUMBER_DOUBLE_EXISTS
              , typeDouble
#endif
#if !ONLY_DOUBLE_WORKER
              , typeFloat128
              , typeDDouble
              , typeQDouble
#endif
            };
  virtual double eps2() { return 1.23e-32; }

  class Allocator
  {
  protected:
    int lowest;
    int allocUpTo;
    int capacity;
    Allocator *above;
  public:
    worker_multi *worker;
    Allocator(int capacity); //fake constructor
    Allocator(worker_multi *worker, int capacity, MandelMath::upgrademe *promote);
    Allocator(Allocator *allo, int capacity); //alloc a block
    Allocator(Allocator *allo, int lowest, int count, upgrademe *relative); //reuse allocated block
    ~Allocator();
    bool checkFill() { return (allocUpTo==lowest+capacity); }
    bool checkIndex(int index) { return (index>=lowest && index<lowest+allocUpTo); }
    number_pointer alloc();
    void dealloc(number_pointer ptr);
    void _getRange(int &first, int &last) { first=lowest; last=lowest+capacity; }
  };

protected:
  Type _ntype;
  virtual number_pointer getNumber(int index)=0; //support for Allocator
  Allocator allocator;
  int capacity;

  virtual void getTmp12(number_pointer &t1, number_pointer &t2)=0;
  friend class MandelMath::complex;

public:
  worker_multi(Type ntype, int capacity): _ntype(ntype), allocator(this, capacity, nullptr), capacity(capacity) { }
  virtual ~worker_multi() { capacity=0; }
  virtual Type ntype() { return this->_ntype; }
  Allocator *getAllocator() { return &allocator; }
  /*virtual number_pointer alloc()=0;
  virtual void alloc_array(int count)=0;
  virtual void dealloc(number_pointer store)=0;
  virtual void dealloc_array(int count)=0;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len)=0; //for now, assert(this.ntype==src.ntype)

  //virtual void init(number_store *store, void *placement, double val=0)=0;
  virtual void zero(const number_pointer store, double val=0)=0;
  //virtual void swap(number_store *store, number_store *src)=0; //dst uninit + src init -> dst init + src uninit
  virtual void assign(number_pointer store, const number_pointer_c src)=0;
  //virtual void cleanup(number_store *store)=0;
  //void promote_(number_worker::Type oldType, number_worker::Type newType, void *placement, const number_store *src=nullptr/*this*/);
  virtual void chs(const number_pointer store)=0;
  virtual void lshift(const number_pointer store, int shoft)=0; // self <<= shoft; 1 lshift -10000 = 0 not error
  virtual void round(const number_pointer store)=0;
  virtual void frac(const number_pointer store)=0; //-1<result<1
  virtual void mod1(const number_pointer store)=0; //0<=result<1
  virtual void add_double(const number_pointer store, double x)=0;
  virtual void mul_double(const number_pointer store, double x)=0;
  virtual void add(const number_pointer store, const number_pointer_c other)=0;
  virtual void sub(const number_pointer store, const number_pointer_c other)=0;
  virtual void rsub(const number_pointer store, const number_pointer_c other)=0;
  virtual void mul(const number_pointer store, const number_pointer_c other)=0;
  virtual void sqr(const number_pointer store)=0;
  virtual double radixfloor(const number_pointer_c store1, number_pointer_c store2)=0; //nearest smaller power of 2 (1.5->1->1)
  virtual void recip(const number_pointer store)=0;
  virtual void sqrt(const number_pointer store)=0;
  virtual int compare(const number_pointer_c store, const number_pointer_c other)=0; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  virtual bool isequal(const number_pointer_c store, const number_pointer_c other)=0; //return store==other
  virtual bool is0(const number_pointer_c store)=0;
  virtual bool isle(const number_pointer_c store, const number_pointer_c other)=0; //return store<=other
  virtual bool isle0(const number_pointer_c store)=0; //return store<=0
  virtual bool isl0(const number_pointer_c store)=0; //return store<0
  virtual bool isl1(const number_pointer_c store)=0; //return store<1

  virtual QString toString(const number_pointer_c store)=0;
  virtual int toRound(const number_pointer_c store)=0;
  virtual double toDouble(const number_pointer_c store)=0;
};

//TODO: combine worker_multi_xxx into template<dd_real, Type::typeDDouble, 1e-64>
//will need to write a wrapper around double and float128

#if NUMBER_DOUBLE_EXISTS
class worker_multi_double: public worker_multi
{
  double *storage;
  double tmp1;
  double tmp2;
protected:
  virtual number_pointer getNumber(int index) override; //support for Allocator
  virtual void getTmp12(number_pointer &t1, number_pointer &t2) override;
public:
  double eps2() override { return 1.23e-32;  /* 2^-(2*53) */ }

  worker_multi_double(int capacity): worker_multi(Type::typeDouble, capacity),
      storage(new double[capacity]) {}
  //worker_multi_double(worker_multi *source);
  worker_multi_double(Allocator *source);
  virtual ~worker_multi_double() override;
  virtual Type ntype() override { return typeDouble; }
  /*virtual number_pointer alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer store) override;
  virtual void dealloc_array(int count) override;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len) override; //for now, assert(this.ntype==src.ntype)

  //void init_(number_store *store, void *placement, double val=0) override;
  void zero(const number_pointer store, double val=0) override;
  //void swap(const number_pointer store, const number_store src) override; //dst uninit + src init -> dst init + src uninit
  void assign(const number_pointer store, const number_pointer_c src) override;// { store->assign<number_double>(*src); };
  //void assignTo(number_store *src) override { store->assignTo_double(*src); };
  //void cleanup(number_store *store) override { store->cleanup(Type::typeDouble); }
  void chs(const number_pointer store) override;
  void lshift(const number_pointer store, int shoft) override;
  void round(const number_pointer store) override;
  void frac(const number_pointer store) override; //-1<result<1
  void mod1(const number_pointer store) override; //0<=result<1
  void add_double(const number_pointer store, double x) override;
  void mul_double(const number_pointer store, double x) override;
  void add(const number_pointer store, const number_pointer_c other) override;
  void sub(const number_pointer store, const number_pointer_c other) override;
  void rsub(const number_pointer store, const number_pointer_c other) override;
  void mul(const number_pointer store, const number_pointer_c other) override;
  void sqr(const number_pointer store) override;
  double radixfloor(const number_pointer_c store1, number_pointer_c store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_pointer store) override;
  void sqrt(const number_pointer store) override;
  int compare(const number_pointer_c store, const number_pointer_c other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_pointer_c store, const number_pointer_c other) override;
  bool is0(const number_pointer_c store) override;
  bool isle(const number_pointer_c store, const number_pointer_c other) override;
  bool isle0(const number_pointer_c store) override;
  bool isl0(const number_pointer_c store) override;
  bool isl1(const number_pointer_c store) override;

  QString toString(const number_pointer_c store) override;
  int toRound(const number_pointer_c store) override;
  double toDouble(const number_pointer_c store) override;
};
#endif //NUMBER_DOUBLE_EXISTS

#if !ONLY_DOUBLE_WORKER
class worker_multi_float128: public worker_multi
{
  __float128 *storage;
  __float128 tmp1;
  __float128 tmp2;
protected:
  virtual number_pointer getNumber(int index) override; //support for Allocator
  virtual void getTmp12(number_pointer &t1, number_pointer &t2) override;
public:
  double eps2() override { return 1.23e-32;  /* 2^-(2*53) */ }

  worker_multi_float128(int capacity): worker_multi(Type::typeFloat128, capacity),
      storage(new __float128[capacity+4]) {}
  worker_multi_float128(worker_multi *source);
  ~worker_multi_float128();
  virtual Type ntype() override { return typeDouble; }
  /*virtual number_pointer alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer store) override;
  virtual void dealloc_array(int count) override;*/

  //void init_(number_store *store, void *placement, double val=0) override;
  void zero(const number_pointer store, double val=0) override;
  //void swap(const number_pointer store, const number_store src) override; //dst uninit + src init -> dst init + src uninit
  void assign(const number_pointer store, const number_pointer_c src) override;// { store->assign<number_double>(*src); };
  //void assignTo(number_store *src) override { store->assignTo_double(*src); };
  //void cleanup(number_store *store) override { store->cleanup(Type::typeDouble); }
  void chs(const number_pointer store) override;
  void lshift(const number_pointer store, int shoft) override;
  void round(const number_pointer store) override;
  void frac(const number_pointer store) override; //-1<result<1
  void mod1(const number_pointer store) override; //0<=result<1
  void add_double(const number_pointer store, double x) override;
  void mul_double(const number_pointer store, double x) override;
  void add(const number_pointer store, const number_pointer_c other) override;
  void sub(const number_pointer store, const number_pointer_c other) override;
  void rsub(const number_pointer store, const number_pointer_c other) override;
  void mul(const number_pointer store, const number_pointer_c other) override;
  void sqr(const number_pointer store) override;
  double radixfloor(const number_pointer_c store1, number_pointer_c store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_pointer store) override;
  void sqrt(const number_pointer store) override;
  int compare(const number_pointer_c store, const number_pointer_c other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_pointer_c store, const number_pointer_c other) override;
  bool is0(const number_pointer_c store) override;
  bool isle(const number_pointer_c store, const number_pointer_c other) override;
  bool isle0(const number_pointer_c store) override;
  bool isl0(const number_pointer_c store) override;
  bool isl1(const number_pointer_c store) override;

  QString toString(const number_pointer_c store) override;
  int toRound(const number_pointer_c store) override;
  double toDouble(const number_pointer_c store) override;
};

class worker_multi_ddouble: public worker_multi
{
protected:
  dd_real *storage;
  dd_real tmp1;
  dd_real tmp2;
protected:
  virtual number_pointer getNumber(int index) override; //support for Allocator
  virtual void getTmp12(number_pointer &t1, number_pointer &t2) override;
public:
  double eps2() override { return 6.1e-64; /* 2^-(2*(53+52)) */ }
  worker_multi_ddouble(int capacity): worker_multi(Type::typeDDouble, capacity),
    storage(new dd_real[capacity+4]) {}
  worker_multi_ddouble(worker_multi *source);
  ~worker_multi_ddouble();
  //{ assert((store->dbgType==Type::typeDDouble) ||
  //         (store->dbgType==Type::typeEmpty)); }
  virtual Type ntype() override { return typeDDouble; }
  /*virtual number_pointer alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer store) override;
  virtual void dealloc_array(int count) override;*/

  //void init_(number_store *store, void *placement, double val=0) override;
  void zero(const number_pointer store, double val=0) override;
  //void swap(const number_pointer store, const number_store src) override; //dst uninit + src init -> dst init + src uninit
  void assign(const number_pointer store, const number_pointer_c src) override;// { store->assign<number_double>(*src); };
  //void assignTo(number_store *src) override { store->assignTo_double(*src); };
  //void cleanup(number_store *store) override { store->cleanup(Type::typeDDouble); }
  void chs(const number_pointer store) override;
  void lshift(const number_pointer store, int shoft) override;
  void round(const number_pointer store) override;
  void frac(const number_pointer store) override; //-1<result<1
  void mod1(const number_pointer store) override; //0<=result<1
  void add_double(const number_pointer store, double x) override;
  void mul_double(const number_pointer store, double x) override;
  void add(const number_pointer store, const number_pointer_c other) override;
  void sub(const number_pointer store, const number_pointer_c other) override;
  void rsub(const number_pointer store, const number_pointer_c other) override;
  void mul(const number_pointer store, const number_pointer_c other) override;
  void sqr(const number_pointer store) override;
  double radixfloor(const number_pointer_c store1, number_pointer_c store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_pointer store) override;
  void sqrt(const number_pointer store) override;
  int compare(const number_pointer_c store, const number_pointer_c other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_pointer_c store, const number_pointer_c other) override;
  bool is0(const number_pointer_c store) override;
  bool isle(const number_pointer_c store, const number_pointer_c other) override;
  bool isle0(const number_pointer_c store) override;
  bool isl0(const number_pointer_c store) override;
  bool isl1(const number_pointer_c store) override;

  QString toString(const number_pointer_c store) override;
  int toRound(const number_pointer_c store) override;
  double toDouble(const number_pointer_c store) override;
};

class worker_multi_qdouble: public worker_multi
{
protected:
  dq_real *storage;
  dq_real tmp1;
  dq_real tmp2;
protected:
  virtual number_pointer getNumber(int index) override; //support for Allocator
  virtual void getTmp12(number_pointer &t1, number_pointer &t2) override;
public:
  double eps2() override { return 6.1e-64; /* 2^-(2*(53+52)) */ }
  worker_multi_qdouble(int capacity): worker_multi(Type::typeQDouble, capacity),
    storage(new dq_real[capacity+4]) {}
  worker_multi_qdouble(worker_multi *source);
  ~worker_multi_qdouble();
  //{ assert((store->dbgType==Type::typeDDouble) ||
  //         (store->dbgType==Type::typeEmpty)); }
  virtual Type ntype() override { return typeQDouble; }
  /*virtual number_pointer alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer store) override;
  virtual void dealloc_array(int count) override;*/

  //void init_(number_store *store, void *placement, double val=0) override;
  void zero(const number_pointer store, double val=0) override;
  //void swap(const number_pointer store, const number_store src) override; //dst uninit + src init -> dst init + src uninit
  void assign(const number_pointer store, const number_pointer_c src) override;// { store->assign<number_double>(*src); };
  //void assignTo(number_store *src) override { store->assignTo_double(*src); };
  //void cleanup(number_store *store) override { store->cleanup(Type::typeDDouble); }
  void chs(const number_pointer store) override;
  void lshift(const number_pointer store, int shoft) override;
  void round(const number_pointer store) override;
  void frac(const number_pointer store) override; //-1<result<1
  void mod1(const number_pointer store) override; //0<=result<1
  void add_double(const number_pointer store, double x) override;
  void mul_double(const number_pointer store, double x) override;
  void add(const number_pointer store, const number_pointer_c other) override;
  void sub(const number_pointer store, const number_pointer_c other) override;
  void rsub(const number_pointer store, const number_pointer_c other) override;
  void mul(const number_pointer store, const number_pointer_c other) override;
  void sqr(const number_pointer store) override;
  double radixfloor(const number_pointer_c store1, number_pointer_c store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_pointer store) override;
  void sqrt(const number_pointer store) override;
  int compare(const number_pointer_c store, const number_pointer_c other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_pointer_c store, const number_pointer_c other) override;
  bool is0(const number_pointer_c store) override;
  bool isle(const number_pointer_c store, const number_pointer_c other) override;
  bool isle0(const number_pointer_c store) override;
  bool isl0(const number_pointer_c store) override;
  bool isl1(const number_pointer_c store) override;

  QString toString(const number_pointer_c store) override;
  int toRound(const number_pointer_c store) override;
  double toDouble(const number_pointer_c store) override;
};
#endif //ONLY_DOUBLE_WORKER

double sqr_double(double x); //no one ever needed this function before year 2022, right?
void complex_double_sqrt(double *res_re, double *res_im, double in_re, double in_im); //res_re>=0
void complex_double_quadratic(double *res_re, double *res_im,
                              double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);
void complex_double_quadratic2(double *res1_re, double *res1_im,
                               double *res2_re, double *res2_im,
                               double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);

class number
{
protected:
  worker_multi *worker;
  worker_multi::Allocator *allocator;
  //bool free_storage_on_destroy;
public:
  number_pointer ptr;
  //or maybe implement cast to number_pointer_c but I don't like that
  number(worker_multi *worker): worker(worker), allocator(worker->getAllocator()), //free_storage_on_destroy(true),
    ptr(allocator->alloc()) {}
  number(worker_multi::Allocator *allocator): worker(allocator->worker), allocator(allocator), //free_storage_on_destroy(true),
    ptr(allocator->alloc()) {}
  //number(number_pointer *src); //swap into me, and clear src, src->freeond=false because already freed
  ~number()
  {
    //if (free_storage_on_destroy)
    {
      allocator->dealloc(ptr);
    };
  }
  void zero(double v=0);
  void assign(const number_pointer_c src);
  //assign(number_pointer) to read complex.re
  void lshift(int shoft);
  void add(const number_pointer_c other);
  //void chs
  void sub(const number_pointer_c other);
  void sqr();
  void mul(const number_pointer_c other);
  void recip();
  //void recip_prepared();
  void sqrt();
  void add_double(double x);
  //toDouble()
};

class complex
{
protected:
  worker_multi *worker;
  worker_multi::Allocator *allocator;
  //number_pointer tmp1;
  //number_pointer tmp2;

public:
  //bool free_storage_on_destroy; //TODO: private
  number_pointer re;
  number_pointer im;
  //complex(number *re, number *im, number *tmp1, number *tmp2): tmp1(tmp1), tmp2(tmp2), re(re), im(im) { }
  complex(worker_multi *worker): worker(worker), allocator(worker->getAllocator()), //free_storage_on_destroy(true),
    re(allocator->alloc()), im(allocator->alloc()) {}
  complex(worker_multi::Allocator *allocator): worker(allocator->worker), allocator(allocator), //free_storage_on_destroy(true),
    re(allocator->alloc()), im(allocator->alloc()) {}
  //complex(number_pointer *src_re, number_pointer *src_im); //swap into me, and clear src, src->freeond=false because already freed
  ~complex()
  {
    //if (free_storage_on_destroy)
    {
      allocator->dealloc(im);
      allocator->dealloc(re);
    };
  }
  void zero(double r=0, double i=0);
  void assign(const complex *other);
  //assign(number_pointer_c re, im:=0)
  void lshift(int shoft);
  void add(const complex *other);
  //chs()
  void sub(const complex *other); //this=this-other
  void rsub(const complex *other); //this=other-this
  //add_double(r), add_double(r, i)
  void sqr();
  void mul(const complex *other);
  //mul(number or number_pointer_c)
  //or mul_double, div_double using tmp (usually mul_int, div_int but also *(1/m-1/n))
  void recip();
  void recip_prepared();
  void sqrt();
  double getMag_double() const;
  const number_pointer_c  getMag_tmp_() const;
  const number_pointer_c getMag1_tmp() const;
  const number_pointer_c getDist1_tmp() const;
  const number_pointer_c mulreT_tmp(const complex *other) const; //Re(this*conjugate(other)) = re*o->re+im*o->im
  double dist2_double(const complex *other) const;
  const number_pointer_c dist2_tmp_(const complex *other) const;
  bool isequal(const complex *other) const;
  //assign_re(number_pointer) or assign_re(number) & assign_re_re(complex) & assign_re_im(complex)
  //is0 used a few times (re==0 && im==0)
};


} // namespace MandelMath
#endif // MANDELMATH_NUMBER_HPP

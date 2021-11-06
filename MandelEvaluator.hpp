#ifndef MANDELEVALUATOR_HPP
#define MANDELEVALUATOR_HPP

#include <QThread>

#include "MandelMath.hpp"
#include "double_double.hpp"

namespace MandelMath
{

struct number_any
{
protected:
  number_store my_store;
  /*number_worker_double d;
  number_worker_ddouble dd;
  number_worker_multi m;*/
public:
  //number_any(MandelMath::number_store::DbgType ntype, MandelMath::number_store *src);
  number_any();
  number_any(number_store *store);
  number_any(number_any *src);
  ~number_any();
  number_worker *impl;
  number_store *store;
  void reinit(number_worker *worker);
  number_worker *ntype_();
};

} //namespace MandelMath

struct MandelPoint
{
  MandelPoint();
  enum State { stUnknown, stOutside, stMaxIter } state;
  MandelMath::number_store zr_, zi_;
  int iter;
  MandelPoint &operator =(MandelPoint &src) = delete;
  void assign(MandelMath::number_worker *worker, const MandelPoint &src);
  void init(MandelMath::number_worker *worker);
  void zero(MandelMath::number_worker *worker);
  void cleanup(MandelMath::number_worker *worker);
};

class MandelEvaluator: public QThread
{
  Q_OBJECT
public:
  MandelEvaluator();
  ~MandelEvaluator();
  static void simple_double(double cr, double ci, MandelPoint &data, int maxiter);
  static void simple_ddouble(MandelMath::dd_real *cr, MandelMath::dd_real *ci, MandelPoint &data, int maxiter);
  static void simple_multi(MandelMath::multiprec *cr, MandelMath::multiprec *ci, MandelPoint &data, int maxiter);

  //MandelMath::number_worker::Type currentType;
  MandelMath::number_worker *currentWorker;
  void switchType(MandelMath::number_worker *worker);
  bool wantStop;
  int pointsComputed;
  QElapsedTimer timeOuter;
  qint64 timeOuterTotal;
  QElapsedTimer timeInner;
  qint64 timeInnerTotal;
  QElapsedTimer timeInvoke;
  qint64 timeInvokePostTotal;
  qint64 timeInvokeSwitchTotal;

  struct ComputeParams
  {
    MandelMath::number_store cr_s;
    MandelMath::number_store ci_s;
    //MandelMath::number_any cr_n;
    //MandelMath::number_any ci_n;
    int epoch;
    int pixelIndex;
    int maxiter;
    ComputeParams();
  } currentParams;

  MandelPoint currentData;
  MandelMath::number_store data_zr_s;
  MandelMath::number_store data_zi_s;
  //MandelMath::number_any data_zr_n;
  //MandelMath::number_any data_zi_n;
  //MandelMath::number_any data_z_tmp1;
  //MandelMath::number_any data_z_tmp2;

  bool startCompute(const MandelPoint *data, bool no_quick_route);
protected:
  void evaluate();
protected slots:
  void doCompute();
signals:
  void doneCompute(MandelEvaluator *me);
};

#endif // MANDELEVALUATOR_HPP

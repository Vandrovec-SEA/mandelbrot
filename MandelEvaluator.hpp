#ifndef MANDELEVALUATOR_HPP
#define MANDELEVALUATOR_HPP

#include <QThread>

#include "MandelMath.hpp"
#include "double_double.hpp"

namespace MandelMath
{
/*
struct number_any
{
protected:
  number_store my_store;
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
*/
} //namespace MandelMath

struct MandelPoint
{
  MandelPoint();
  enum State { stUnknown, stOutside, stOutAngle, stBoundary, stMisiur, stDiverge, stPeriod2, stPeriod3, stMaxIter } state;
  MandelMath::number_store f_re, f_im;
  MandelMath::number_store fc_c_re, fc_c_im; //fc_c, or fz_r if stPeriod2 or stPeriod3
  MandelMath::number_store fz_c_mag;
  int near0iter;
  MandelMath::number_store near0f_re, near0f_im;
  int period;
  MandelMath::number_store root_re, root_im;
  int iter;
  double exterior_hits, exterior_avoids; //upper and lower bound
  double interior;
  MandelPoint &operator =(MandelPoint &src) = delete;
  void assign(MandelMath::number_worker *worker, const MandelPoint &src);
  void init(MandelMath::number_worker *worker);
  void zero(MandelMath::number_worker *worker, const MandelMath::number_store *c_re, const MandelMath::number_store *c_im);
  void cleanup(MandelMath::number_worker *worker);
};

class MandelEvaluator: public QThread
{
  Q_OBJECT
public:
  constexpr static double LARGE_FLOAT2=1e60;
  constexpr static double MAGIC_MIN_SHRINK=1.5;
  constexpr static int MAX_PERIOD=8000;
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
    MandelMath::number_store c_re, c_im;
    int epoch;
    int pixelIndex;
    int maxiter;
    ComputeParams();
  } currentParams;
  MandelPoint currentData;

  bool startCompute_(const MandelPoint *data, int quick_route); //qr: -1..never 0..auto 1..always
protected:
  typedef MandelMath::complex complex;
  struct
  {
    MandelMath::number_store fz_r_re, fz_r_im;
    MandelMath::number_store near0fmag;
    MandelMath::number_store lookper_startf_re, lookper_startf_im;
    MandelMath::number_store lookper_bestf_re, lookper_bestf_im;
    MandelMath::number_store lookper_near0f_re, lookper_near0f_im; //TODO: near0f->nearc
    MandelMath::number_store lookper_near0f_dist;
    MandelMath::number_store lookper_dist2;
    int lookper_startiter, lookper_prevGuess, lookper_lastGuess;
    MandelMath::number_store lookper_totalFzmag;
  } eval;
  struct
  {
    MandelMath::number_store bestr_re, bestr_im;
    MandelMath::number_store f_r_re, f_r_im;
    MandelMath::number_store fz_r_re, fz_r_im;
    MandelMath::number_store fzz_r_re, fzz_r_im;
    MandelMath::number_store tmp1_re, tmp1_im;
    MandelMath::number_store fzfix_re, fzfix_im;
    MandelMath::number_store laguH_re, laguH_im;
    MandelMath::number_store laguG_re, laguG_im;
    MandelMath::number_store laguG2_re, laguG2_im;
    MandelMath::number_store laguX_re, laguX_im;
    MandelMath::number_store newtX_re, newtX_im;
    MandelMath::number_store fzzf_re, fzzf_im;
    MandelMath::number_store tmp2;
  } newt;
  struct InteriorInfo
  {
    MandelMath::number_store inte_re, inte_im, inte_abs;
    MandelMath::number_store fz_re, fz_im, fz_mag;
  } interior;

  int newton(int period, const complex *c, complex *r, const bool fastHoming, const int suggestedMultiplicity);
  int periodCheck(int period, const complex *c);
  int estimateInterior(int period, const complex *c, const complex *root);//, InteriorInfo *interior);
  void eval_until_bailout(complex *c, complex *f, complex *fc_c);
#if COMPLEX_IS_TEMPLATE
  template <class NW>
#endif
  void evaluate();
protected slots:
  void doCompute();
signals:
  void doneCompute(MandelEvaluator *me);
};

#endif // MANDELEVALUATOR_HPP

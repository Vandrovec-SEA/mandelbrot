#ifndef MANDELEVALUATOR_HPP
#define MANDELEVALUATOR_HPP

#include <QThread>

#include "MandelMath.hpp"
#include "double_double.hpp"

enum NewtonNaiveChoice { nc03, nc05, nc08, ncWide, nc100, nc90_, nc80, nc60, ncClose };

struct LaguerrePoint
{
  LaguerrePoint();
  enum State { stUnknown, stResolved, stFail } state;
  MandelMath::number_store f_re, f_im;
  MandelMath::number_store fz_r_re, fz_r_im;
  double firstM;
  NewtonNaiveChoice naiveChoice;
  //MandelMath::number_store fz_c_mag;
  //int period;
  //MandelMath::number_store root_re, root_im;
  int iter;
  LaguerrePoint &operator =(LaguerrePoint &src) = delete;
  void assign(MandelMath::number_worker *worker, const LaguerrePoint &src);
  void init(MandelMath::number_worker *worker);
  void zero(MandelMath::number_worker *worker, const MandelMath::number_store *c_re, const MandelMath::number_store *c_im);
  void cleanup(MandelMath::number_worker *worker);
  void promote(MandelMath::number_worker::Type oldType, MandelMath::number_worker::Type newType);
  union Place
  {
    static constexpr size_t LEN=4;
    MandelMath::dd_real (*dd)[LEN];
    MandelMath::multiprec (*multi)[LEN];
    Place(): dd(nullptr) { }
  } place;
};

struct MandelPoint
{
  MandelPoint();
  enum State { stUnknown, stOutside, stOutAngle, stBoundary, stMisiur, stDiverge, stPeriod2, stPeriod3, stMaxIter } state;
  MandelMath::number_store f_re, f_im;
  MandelMath::number_store fc_c_re_, fc_c_im_; //fc_c, or fz_r if stPeriod2 or stPeriod3
  bool has_fc_r;
  MandelMath::number_store fz_r_re, fz_r_im;
  MandelMath::number_store fz_c_mag_;
  int lookper_startiter, lookper_prevGuess_, lookper_lastGuess;
  MandelMath::number_store lookper_startf_re, lookper_startf_im;
  MandelMath::number_store lookper_nearr_dist_;
  bool lookper_nearr_dist_touched; //check if equal but only once; in theory should only happen at dist=0
  MandelMath::number_store lookper_totalFzmag;
  int near0iter;
  MandelMath::number_store near0f_re, near0f_im; //TODO: delete? useless
  int period;
  MandelMath::number_store root_re, root_im;
  int iter, newton_iter;
  double exterior_hits, exterior_avoids; //upper and lower bound
  double interior;
  MandelPoint &operator =(MandelPoint &src) = delete;
  void assign(MandelMath::number_worker *worker, const MandelPoint &src);
  void init(MandelMath::number_worker *worker);
  void zero(MandelMath::number_worker *worker, const MandelMath::number_store *c_re, const MandelMath::number_store *c_im);
  void cleanup(MandelMath::number_worker *worker);
  void promote(MandelMath::number_worker::Type oldType, MandelMath::number_worker::Type newType);
  union Place
  {
    static constexpr size_t LEN=13;
    MandelMath::dd_real (*dd)[LEN];
    MandelMath::multiprec (*multi)[LEN];
    Place(): dd(nullptr) { }
  } place;
};

class ShareableViewInfo: public QObject
{
  Q_OBJECT
public:
  ShareableViewInfo();
  ShareableViewInfo(ShareableViewInfo &src);
  ShareableViewInfo(const ShareableViewInfo &src);
  ShareableViewInfo(ShareableViewInfo &&src); //important
  ShareableViewInfo &operator=(ShareableViewInfo &src);
  ShareableViewInfo &operator=(ShareableViewInfo &&src);
  MandelMath::number_worker *worker;
  MandelMath::number_store re_, im;
  MandelMath::number_store root_re, root_im;
  MandelMath::number_place re_p, im_p, rre_p, rim_p;
  double scale;
  int period;
};
Q_DECLARE_METATYPE(ShareableViewInfo);

class LaguerreStep
{
public:
  LaguerreStep();
  union Place
  {
    static constexpr size_t LEN=20;
    MandelMath::dd_real (*dd)[LEN];
    MandelMath::multiprec (*multi)[LEN];
    Place(): dd(nullptr) { }
  } place;
  MandelMath::number_worker *currentWorker;
  void switchType(MandelMath::number_worker *worker);
  //do one Laguerre step
  bool eval(int lg2_degree, const MandelMath::complex *f, const MandelMath::complex *f_z, const MandelMath::complex *f_zz); //->step_re, step_im

  //results
  MandelMath::number_store step_re, step_im;
  struct
  {
    double mum_re, mum_im;
    double mu_re, mu_im;
    int lastm;
  } dbg;

  //temporary
protected:
  MandelMath::number_store s1_re, s1_im;
  MandelMath::number_store s2_re, s2_im;
  MandelMath::number_store tmp1_re, tmp1_im, tmp2;
  MandelMath::number_store laguG_re, laguG_im;
  MandelMath::number_store laguG2_re, laguG2_im;
  MandelMath::number_store laguH_re, laguH_im;
  MandelMath::number_store laguX_re, laguX_im;
  //->step_re, step_im MandelMath::number_store newtX_re, newtX_im;
  MandelMath::number_store fzzf_re, fzzf_im;
};

class MandelLoopEvaluator
{
public:
  MandelLoopEvaluator();
  union Place
  {
    static constexpr size_t LEN=16;
    MandelMath::dd_real (*dd)[LEN];
    MandelMath::multiprec (*multi)[LEN];
    Place(): dd(nullptr) { }
  } place;
  MandelMath::number_worker *currentWorker;
  void switchType(MandelMath::number_worker *worker);
  //eval F^o(period)_c(c)
  bool evalg(int period, const MandelMath::complex *c); //->f, f_c, f_cc
  //eval F^o(period)_c(z)
  bool eval2(int period, const MandelMath::complex *c, const MandelMath::complex *z); //->f, f_z, f_c, f_zz, f_zc, f_cc
  bool eval_zz(int period, const MandelMath::complex *c, const MandelMath::complex *z); //->f, f_z, f_zz
  bool eval_multi(int period, const MandelMath::complex *c, const MandelMath::complex *z, const MandelMath::complex *f_z_target); //->f, f_z, f_zz, multi, first_multi
  bool eval2zzc(int period, const MandelMath::complex *c, const MandelMath::complex *z); //->f, f_z, f_c, f_zz, f_zc, f_cc, f_zzc

  //inputs
  //MandelMath::number_store z_re, z_im;
  //MandelMath::number_store c_re, c_im;

  //results
  MandelMath::number_store f_re, f_im;
  MandelMath::number_store f_z_re, f_z_im;
  MandelMath::number_store f_c_re, f_c_im;
  MandelMath::number_store f_zz_re, f_zz_im;
  MandelMath::number_store f_zc_re, f_zc_im;
  MandelMath::number_store f_cc_re, f_cc_im;
  MandelMath::number_store f_zzc_re, f_zzc_im;
  int multi;
  MandelMath::number_store first_multi_re, first_multi_im;
  MandelMath::number_store sumA_re, sumA_im;

  /*MandelMath::number_store g_re, g_im;
  MandelMath::number_store g_c_re, g_c_im;//, g_c_mag;
  MandelMath::number_store g_cc_re, g_cc_im;
  MandelMath::number_store g_c2_re, g_c2_im;*/

  //temporary
protected:
  MandelMath::number_store s1_re, s1_im;
  MandelMath::number_store s2_re, s2_im;
  //MandelMath::number_store rbx_re, rbx_im;
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
#if NUMBER_DOUBLE_EXISTS
  static void simple_double(double cr, double ci, MandelPoint &data, int maxiter);
#endif //NUMBER_DOUBLE_EXISTS
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
    int maxiter_;
    bool breakOnNewNearest;
    bool want_fc_r;
    ComputeParams();
  } currentParams;
  MandelPoint currentData;

  bool startCompute(const MandelPoint *data, int quick_route); //qr: -1..never 0..auto 1..always
  void startNewton(int period, const MandelMath::complex *c /*, currentData.f const *root, */);
  int newton(int period, const MandelMath::complex *c, MandelMath::complex *r, const bool fastHoming, const int suggestedMultiplicity);
  struct NewtRes
  {
    int cyclesNeeded;
    //now firstMum_re double firstM;
    MandelMath::number_store fz_r_re_, fz_r_im_;
    MandelMath::number_store first_guess_lagu_re, first_guess_lagu_im;
    MandelMath::number_store first_guess_newt_re, first_guess_newt_im;
    double first_fejer_re, first_fejer_im;
    double first_naive1_re_, first_naive1_im, first_naive2_re, first_naive2_im, first_naive_re, first_naive_im;
    NewtonNaiveChoice naiveChoice;
    double first_neumaier1_re_, first_neumaier1_im_; //circle enclosing 1 but not 2 roots (approx)
    double first_neumaier2_re, first_neumaier2_im;   //circle enclosing 2 but not 3 roots (approx) (never valid without  f''')
    //double first_lagum_re, first_lagum_im;
    double first_lagu1_re, first_lagu1_im, first_lagu1o_re, first_lagu1o_im;
    double firstMu_re_, firstMu_im, firstMum_re_, firstMum_im_;
    double accy_tostop, accy_multiplier; //in units of eps2()
  } newtres_;
protected:
  typedef MandelMath::complex complex;
  struct
  {
    MandelMath::number_store fz_r_re, fz_r_im;
    MandelMath::number_store near0fmag;
    //MandelMath::number_store lookper_bestf_re, lookper_bestf_im;
    MandelMath::number_store lookper_nearr_re, lookper_nearr_im;//, lookper_nearr_dist;
    //MandelMath::number_store lookper_dist2;
    //int lookper_startiter, lookper_prevGuess;
    //int lookper_lastGuess;
    //MandelMath::number_store lookper_totalFzmag;
  } eval;
  struct
  {
    MandelMath::number_store bestr_re, bestr_im;
    MandelMath::number_store f_r_re, f_r_im;
    //MandelMath::number_store fz_r_re, fz_r_im;
    MandelMath::number_store fzz_r_re, fzz_r_im;
    MandelMath::number_store tmp1_re, tmp1_im;
    //MandelMath::number_store fzfix_re, fzfix_im;
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
  /*struct //Bulb
  {
    MandelMath::number_store rb_re, rb_im;
    MandelMath::number_store cb_re, cb_im;
    MandelMath::number_store xc_re, xc_im;
    MandelMath::number_store baseZC_re, baseZC_im;
    MandelMath::number_store baseCC_re, baseCC_im;
    MandelMath::number_store s1_re, s1_im;
    MandelMath::number_store s2_re, s2_im;
    MandelMath::number_store cbx_re, cbx_im;
    MandelMath::number_store rbx_re, rbx_im;

    MandelMath::number_store g_re, g_im;
    MandelMath::number_store g_c_re, g_c_im;//, g_c_mag;
    MandelMath::number_store g_cc_re, g_cc_im;
    MandelMath::number_store g_c2_re, g_c2_im;

    MandelMath::number_store f_re, f_im;
    MandelMath::number_store f_z_re, f_z_im;
    MandelMath::number_store f_c_re, f_c_im;
    MandelMath::number_store f_zz_re, f_zz_im;
    MandelMath::number_store f_zc_re, f_zc_im;
    MandelMath::number_store f_cc_re, f_cc_im;
  } bulb;*/
public:
  struct
  {
    MandelLoopEvaluator bulbe;
    MandelMath::number_store rb_re_, rb_im;
    MandelMath::number_store cb_re, cb_im;
    MandelMath::number_store xc_re, xc_im;
    MandelMath::number_store cbx_re, cbx_im;
    MandelMath::number_store rbx_re, rbx_im;
    MandelMath::number_store baseZC_re, baseZC_im;
    MandelMath::number_store baseCC_re, baseCC_im;
    MandelMath::number_store s1_re, s1_im;
    MandelMath::number_store s2_re_, s2_im_;
    MandelMath::number_store s3_re, s3_im_;
    MandelMath::number_store deltac_re, deltac_im;
    MandelMath::number_store deltar_re, deltar_im;
    MandelMath::number_store B_re, B_im;
    MandelMath::number_store C_re, C_im;
    MandelMath::number_store dbg_first_cb_re, dbg_first_cb_im;
    MandelMath::number_store dbg_first_rb_re, dbg_first_rb_im;
    MandelMath::number_store target_f_z_re, target_f_z_im;
    //MandelMath::number_store test_x0_re, test_x0_im;
    //MandelMath::number_store test_xn_re, test_xn_im;
    double dbg_guessmult;
    LaguerreStep lagu;
  } bulb;
protected:
  union Place
  {
    static constexpr size_t LEN=78;
    MandelMath::dd_real (*dd)[LEN];
    MandelMath::multiprec (*multi)[LEN];
    Place(): dd(nullptr) { }
  } place;

  public:
  bool findBulbBase(int period2, const MandelMath::complex *c, MandelMath::complex *cb, MandelMath::complex *rb, MandelMath::complex *xc, MandelMath::complex *baseZC, MandelMath::complex *baseCC, bool *is_card, int *foundMult);
  protected:
  void fixRnearBase(MandelMath::complex *r, const MandelMath::complex *c, int period, int *mult);
  int periodCheck(int period, const complex *c);
  int estimateInterior(int period, const complex *c, const complex *root);//, InteriorInfo *interior);
  void eval_until_bailout(complex *c, complex *f, complex *fc_c);
#if COMPLEX_IS_TEMPLATE
  template <class NW>
#endif
  void evaluate();
protected slots:
  void doCompute();
  void doNewton();
signals:
  void doneCompute(MandelEvaluator *me);
  void doneNewton(MandelEvaluator *me, int result);
};

#endif // MANDELEVALUATOR_HPP

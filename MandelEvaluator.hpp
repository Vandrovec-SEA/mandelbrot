#ifndef MANDELEVALUATOR_HPP
#define MANDELEVALUATOR_HPP

#include <QThread>

#include "MandelMath.hpp"
#include "double_double.hpp"

enum NewtonNaiveChoice { nc03, nc05, nc08, ncWide, nc100, nc90_, nc80, nc60, ncClose };

struct LaguerrePointStore
{
  enum State { stUnknown, stResolved, stFail } state;
  double firstM;
  NewtonNaiveChoice naiveChoice;
  int iter;
  LaguerrePointStore();
  void assign(const LaguerrePointStore *src);
};

struct LaguerrePoint
{
  static constexpr size_t LEN=4;
  MandelMath::worker_multi::Allocator self_allocator;
  LaguerrePointStore *store;
  LaguerrePoint(LaguerrePointStore *store, MandelMath::worker_multi::Allocator *allocator);
  MandelMath::complex f;
  MandelMath::complex fz_r;
  LaguerrePoint &operator =(LaguerrePoint &src) = delete;
  void assign(const LaguerrePoint &src);
  //void init(MandelMath::worker_multi *worker, int index_into_worker);
  void zero(const MandelMath::complex *c);
  //void cleanup(MandelMath::worker_multi *worker, int index_into_worker);
  //void promote(MandelMath::number_worker::Type oldType, MandelMath::number_worker::Type newType);
};

struct MandelPointStore
{
  //TODO: split into state_working { idle, working, resolved }, result { unknown, period2, ... }
  //because we still work on known result for fc_r
  enum WorkState { stIdle, stWorking, stDone };
  std::atomic<WorkState> wstate;
  enum ResultState { stUnknown_, stOutside, stOutAngle, stBoundary, stMisiur, stDiverge, stPeriod2, stPeriod3, stMaxIter };
  ResultState rstate;
  bool has_fc_r;
  int lookper_startiter, lookper_prevGuess_, lookper_lastGuess;
  bool lookper_nearr_dist_touched; //check if equal but only once; in theory should only happen at dist=0
  int near0iter;
  int period;
  int iter, newton_iter;
  double exterior_hits, exterior_avoids; //upper and lower bound
  double interior;

  MandelPointStore();
  void assign(const MandelPointStore *src);
};

struct MandelPoint
{
  static constexpr size_t LEN=15;
  MandelMath::worker_multi::Allocator self_allocator;
  MandelPointStore *store;
  MandelPoint(MandelPointStore *store, MandelMath::worker_multi::Allocator *allocator);
  MandelMath::complex f;
  MandelMath::complex fc_c; //fc_c, or fz_r if stPeriod2 or stPeriod3
  MandelMath::complex fz_r;
  MandelMath::number fz_c_mag;
  MandelMath::complex lookper_startf;
  MandelMath::number lookper_nearr_dist;
  MandelMath::number lookper_totalFzmag;
  MandelMath::complex near0f; //TODO: delete? useless
  MandelMath::complex root;
  MandelPoint &operator =(MandelPoint &src) = delete;
  void assign(const MandelPoint &src);
  //void init(MandelMath::worker_multi *worker, int index_into_worker);
  void zero(/*MandelMath::worker_multi *worker, int index_into_worker, */const MandelMath::complex *c);
  //void cleanup(MandelMath::worker_multi *worker, int index_into_worker);
  //void promote(MandelMath::number_worker::Type oldType, MandelMath::number_worker::Type newType);
};

class ShareableViewInfo: public QObject
{
  Q_OBJECT
protected:
  MandelMath::worker_multi::Allocator selfAllocator;
  int refcount;
public:
  MandelMath::worker_multi::Allocator *originalAllocator; //need it when reading out
  static constexpr int LEN=4;
  ShareableViewInfo(): selfAllocator(LEN), c(&selfAllocator), root(&selfAllocator) {}
  ShareableViewInfo(MandelMath::worker_multi::Allocator *allocator);
  ShareableViewInfo(ShareableViewInfo &src);
  ShareableViewInfo(const ShareableViewInfo &src);
  ShareableViewInfo(ShareableViewInfo &&src); //important
  //~ShareableViewInfo();
  ShareableViewInfo &operator=(ShareableViewInfo &src);
  ShareableViewInfo &operator=(ShareableViewInfo &&src);
  MandelMath::complex c;
  MandelMath::complex root;
  double scale;
  int period;
};
Q_DECLARE_METATYPE(ShareableViewInfo);

class LaguerreStep
{
public:
  MandelMath::worker_multi::Allocator self_allocator;
  static constexpr size_t LEN=19;
  MandelMath::worker_multi *currentWorker;
  LaguerreStep(MandelMath::worker_multi::Allocator *allocator);
  //~LaguerreStep();
  //void switchType(MandelMath::number_worker *worker);
  //do one Laguerre step
  bool eval(int lg2_degree, const MandelMath::complex *const f, const MandelMath::complex *const f_z, const MandelMath::complex *const f_zz); //->step_re, step_im

  //results
  MandelMath::complex step;
  struct
  {
    double mum_re, mum_im;
    double mu_re, mu_im;
    int lastm;
  } dbg;

  //temporary
protected:
  MandelMath::complex s1;
  MandelMath::complex s2;
  MandelMath::complex tmp1;
  MandelMath::number tmp2;
  MandelMath::complex laguG;
  MandelMath::complex laguG2;
  MandelMath::complex laguH;
  MandelMath::complex laguX;
  //->step_re, step_im MandelMath::number_store newtX_re, newtX_im;
  MandelMath::complex fzzf;
};

class MandelLoopEvaluator
{
public:
  MandelMath::worker_multi::Allocator self_allocator;
  static constexpr size_t LEN=23;
  MandelMath::worker_multi *currentWorker;
  MandelLoopEvaluator(MandelMath::worker_multi::Allocator *allocator);
  //~MandelLoopEvaluator();
  /*union Place
  {
    MandelMath::dd_real (*dd)[LEN];
    MandelMath::multiprec (*multi)[LEN];
    Place(): dd(nullptr) { }
  } place;*/
  //void switchType(MandelMath::number_worker *worker);
  //eval F_c(c)
  bool evalg(int period, const MandelMath::complex *const c); //->f, f_c, f_cc
  //eval F_c(z)
  bool eval2(int period, const MandelMath::complex *const c, const MandelMath::complex *const z);     //->f, f_z, f_zz, f_c, f_zc, f_cc
  bool eval2_mag(int period, const MandelMath::complex *const c, const MandelMath::complex *const z); //->f, f_z, f_zz, f_c, f_zc, f_z_mag
  bool eval_zz(int period, const MandelMath::complex *const c, const MandelMath::complex *const z, bool minusR);   //->f, f_z, f_zz
  bool eval2zzc(int period, const MandelMath::complex *const c, const MandelMath::complex *const z);  //->f, f_z, f_zz, f_c, f_zc, f_cc, f_zzc
  bool eval_multi(int period, const MandelMath::complex *const c, const MandelMath::complex *const z, const MandelMath::complex *const f_z_target); //->f, f_z, f_zz, multi, first_multi

  //inputs
  //MandelMath::number_store z_re, z_im;
  //MandelMath::number_store c_re, c_im;

  //results
  MandelMath::complex f;
  MandelMath::complex f_z;
  MandelMath::complex f_c;
  MandelMath::complex f_zz;
  MandelMath::complex f_zc;
  MandelMath::complex f_cc;
  MandelMath::complex f_zzc;
  int multi;
  MandelMath::complex first_multi;
  MandelMath::complex sumA;
  MandelMath::number f_z_mag;

  //temporary
protected:
  MandelMath::complex s1;
  MandelMath::complex s2;
};


class MandelEvaluator: public QThread
{
  Q_OBJECT
protected:
  //MandelMath::worker_multi::Allocator self_allocator;
public:
  constexpr static double LARGE_FLOAT2=1e60;
  constexpr static double MAGIC_MIN_SHRINK=1.5;
  constexpr static int MAX_PERIOD=8000;
  int busyEpoch;
  MandelMath::worker_multi *currentWorker;
  MandelEvaluator(MandelMath::worker_multi::Type ntype, bool dontRun);
  ~MandelEvaluator();
#if NUMBER_DOUBLE_EXISTS
  static void simple_double(double cr, double ci, MandelPoint &data, int maxiter);
#endif //NUMBER_DOUBLE_EXISTS
  static void simple_ddouble(MandelMath::dd_real *cr, MandelMath::dd_real *ci, MandelPoint &data, int maxiter);
  static void simple_multi(MandelMath::multiprec *cr, MandelMath::multiprec *ci, MandelPoint &data, int maxiter);

  //MandelMath::number_worker::Type currentType;
  //MandelMath::worker_multi::Allocator *currentAllocator;
  //void switchType(MandelMath::number_worker *worker);
  int workIfEpoch;
  int pointsComputed;
  qint64 totalNewtonIterations;
  QElapsedTimer timeThreaded;
  qint64 timeThreadedTotal;
  QElapsedTimer timeOuter_;
  qint64 timeOuterTotal_;
  QElapsedTimer timeInner_;
  qint64 timeInnerTotal_;
  QElapsedTimer timeInvoke_;
  qint64 timeInvokePostTotal_;
  qint64 timeInvokeSwitchTotal_;

  MandelMath::worker_multi::Allocator params_allocator;
  struct ComputeParams
  {
    static constexpr int LEN=2;
    MandelMath::complex c;
    int epoch;
    int pixelIndex;
    int maxiter_;
    bool breakOnNewNearest;
    bool want_fc_r;
    ComputeParams(MandelMath::worker_multi::Allocator *allocator);
  } currentParams;
  //MandelMath::worker_multi::Allocator currentDataAllocator;
  MandelPointStore currentDataStore;
  MandelPoint currentData;
  LaguerrePoint tmpLaguerrePoint;

  bool startCompute(/*const MandelPoint *data,*/ int quick_route); //works directly on currentData //qr: -1..never 0..auto 1..always
  void startNewton(int period, const MandelMath::complex *c /*, currentData.f const *root, */);
  int newton(int period, const MandelMath::complex *c, MandelMath::complex *r, const bool fastHoming, const int suggestedMultiplicity);
  struct NewtRes
  {
    MandelMath::worker_multi::Allocator self_allocator;
    static constexpr int LEN=6;
    int cyclesNeeded;
    //now firstMum_re double firstM;
    MandelMath::complex fz_r; //TODO: would be nice to move into newt and move newton() to newt
    MandelMath::complex first_guess_lagu;
    MandelMath::complex first_guess_newt;
    double first_fejer_re, first_fejer_im;
    double first_naive1_re_, first_naive1_im, first_naive2_re, first_naive2_im, first_naive_re, first_naive_im;
    NewtonNaiveChoice naiveChoice;
    double first_neumaier1_re_, first_neumaier1_im_; //circle enclosing 1 but not 2 roots (approx)
    double first_neumaier2_re, first_neumaier2_im;   //circle enclosing 2 but not 3 roots (approx) (never valid without  f''')
    //double first_lagum_re, first_lagum_im;
    double first_lagu1_re, first_lagu1_im, first_lagu1o_re, first_lagu1o_im;
    double firstMu_re_, firstMu_im, firstMum_re_, firstMum_im_;
    double accy_tostop, accy_multiplier; //in units of eps2()
    NewtRes(MandelMath::worker_multi::Allocator *allocator);
  } newtres_;
protected:
  struct Eval
  {
    MandelMath::worker_multi::Allocator self_allocator;
    static constexpr int LEN=6;
    MandelMath::complex fz_r;
    MandelMath::number fz_mag1;
    MandelMath::number near0fmag;
    //MandelMath::number_store lookper_bestf_re, lookper_bestf_im;
    MandelMath::complex lookper_nearr;//, lookper_nearr_dist;
    //MandelMath::number_store lookper_dist2;
    //int lookper_startiter, lookper_prevGuess;
    //int lookper_lastGuess;
    //MandelMath::number_store lookper_totalFzmag;
    Eval(MandelMath::worker_multi::Allocator *allocator);
  } eval;
  struct Newt
  {
    MandelMath::worker_multi::Allocator self_allocator;
    static constexpr int LEN=25;
    MandelMath::complex bestr;
    MandelMath::complex f_r;
    //MandelMath::number_store fz_r_re, fz_r_im;
    MandelMath::complex fzz_r;
    MandelMath::complex tmp1;
    //MandelMath::number_store fzfix_re, fzfix_im;
    MandelMath::complex laguH;
    MandelMath::complex laguG;
    MandelMath::complex laguG2;
    MandelMath::complex laguX;
    MandelMath::complex newtX;
    MandelMath::complex prevR;
    MandelMath::complex prevGz;
    MandelMath::complex fzzf;
    MandelMath::number tmp2;
    Newt(MandelMath::worker_multi::Allocator *allocator);
  } newt;
  struct InteriorInfo
  {
    MandelMath::worker_multi::Allocator self_allocator;
    static constexpr int LEN=6;
    MandelMath::complex inte;
    MandelMath::number inte_abs;
    MandelMath::complex fz;
    MandelMath::number fz_mag;
    InteriorInfo(MandelMath::worker_multi::Allocator *allocator);
  } interior;
public:
  struct Bulb
  {
    MandelMath::worker_multi::Allocator self_allocator;
    MandelMath::worker_multi *currentWorker;
    MandelLoopEvaluator bulbe;
    MandelMath::complex rb;
    MandelMath::complex cb;
    MandelMath::complex xc;
    MandelMath::complex cbx; //TODO: unused?
    MandelMath::complex rbx; //TODO: unused?
    MandelMath::complex baseZC;
    MandelMath::complex baseCC;
    MandelMath::complex s1;
    MandelMath::complex s2;
    MandelMath::complex s3;
    MandelMath::complex deltac;
    MandelMath::complex deltar;
    MandelMath::complex B; //TODO: unused?
    MandelMath::complex C; //TODO: unused?
    MandelMath::complex dbg_first_cb;
    MandelMath::complex dbg_first_rb;
    MandelMath::complex target_f_z;
    //MandelMath::number_store test_x0_re, test_x0_im;
    //MandelMath::number_store test_xn_re, test_xn_im;
    double dbg_guessmult;
    LaguerreStep lagu;
    Bulb(MandelMath::worker_multi::Allocator *allocator);
  protected:
    void fixRnearBase(MandelMath::complex *r, const MandelMath::complex *c, int period, int *mult);
  public:
    bool findBulbBase(int period2, const MandelMath::complex *c, MandelMath::complex *cb, MandelMath::complex *rb, MandelMath::complex *xc, MandelMath::complex *baseZC, MandelMath::complex *baseCC, bool *is_card, int *foundMult);
    static constexpr int LEN=MandelLoopEvaluator::LEN+34+LaguerreStep::LEN;
  } bulb;
  static constexpr size_t LEN=ComputeParams::LEN+MandelPoint::LEN+LaguerrePoint::LEN+NewtRes::LEN+Eval::LEN+Newt::LEN+InteriorInfo::LEN+Bulb::LEN;
public:
  struct
  {
    std::function<bool (MandelEvaluator *me)> give;
    std::function<bool (MandelEvaluator *me, bool giveWork)> done;
  } threaded;
protected:

  protected:
  int periodCheck(int period, const MandelMath::complex *c, const MandelMath::complex *root_seed);
  int estimateInterior(int period, const MandelMath::complex *c, const MandelMath::complex *root);//, InteriorInfo *interior);
  void eval_until_bailout(const MandelMath::complex *c, MandelMath::complex *f, MandelMath::complex *fc_c);
  void evaluate();
protected slots:
  void doCompute();
  void doNewton();
public slots:
  void doComputeThreaded(int epoch);
signals:
  void doneCompute(MandelEvaluator *me);
  void doneComputeThreaded(MandelEvaluator *me);
  void doneNewton(MandelEvaluator *me, int result);
};

#endif // MANDELEVALUATOR_HPP

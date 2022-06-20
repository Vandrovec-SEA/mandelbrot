#ifndef LAGUERREMODEL_H
#define LAGUERREMODEL_H

#include <QObject>
#include <QImage>

#include "ShareableImageWrapper.hpp"
#include "MandelEvaluator.hpp"

class LaguerreModel: public QObject
{
  Q_OBJECT
public:
  LaguerreModel();
  ~LaguerreModel();
  void transformStore(MandelMath::worker_multi *old_worker, MandelMath::worker_multi *old_sworker, LaguerrePointStore *old_store, int old_width, int old_height, const MandelMath::complex *old_c,
                      MandelMath::worker_multi *new_worker, MandelMath::worker_multi *new_sworker, LaguerrePointStore *new_store, int new_width, int new_height, const MandelMath::complex *new_c,
                      int inlog, int new_step_log);
  Q_INVOKABLE void setParams(ShareableViewInfo viewInfo);
  void setView_(const MandelMath::complex *c, double scale);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  //void setWorker(MandelMath::worker_multi *newWorker);
  void startNewEpoch();
  void giveWorkAll();
  Q_INVOKABLE int writeToImage(ShareableImageWrapper img);
  void reimToPixel(int *circ_x, int *circ_y, const MandelMath::complex *c, MandelMath::number *tmp);
  Q_INVOKABLE void paintOrbit(ShareableImageWrapper image, int x, int y);
  Q_INVOKABLE QString pixelXtoRE_str(int x);
  Q_INVOKABLE QString pixelYtoIM_str(int y);
  Q_INVOKABLE QString getTimes();
  Q_INVOKABLE QString getTextXY();
  Q_INVOKABLE QString getTextInfoGen();
  Q_INVOKABLE QString getTextInfoSpec();
  MandelMath::worker_multi::Allocator *paramsAllocator;
  struct Params
  {
    static constexpr int LEN=4;
    int period;
    MandelMath::complex base;
    MandelMath::complex root;
    Params(MandelMath::worker_multi::Allocator *allocator);
    void assign(Params *src);
  } *params_;


  enum paintStyle
  {
    paintStyleCls=0
  };
  Q_ENUM(paintStyle);
  paintStyle _selectedPaintStyle;
  Q_PROPERTY(paintStyle selectedPaintStyle READ getselectedPaintStyle WRITE setselectedPaintStyle NOTIFY selectedPaintStyleChanged)
  paintStyle getselectedPaintStyle() { return _selectedPaintStyle; }
  void setselectedPaintStyle(paintStyle ps) { _selectedPaintStyle=ps; }
  int _threadsWorking;
  Q_PROPERTY(int threadsWorking READ getThreadsWorking CONSTANT)
  int getThreadsWorking() { return _threadsWorking; }
  Q_PROPERTY(int threadsMax READ getThreadCount CONSTANT)
  int getThreadCount() { return threadCount; }

  enum precision
  {
    precisionDouble=0,
#if !ONLY_DOUBLE_WORKER
    precisionFloat128=1,
    precisionDDouble=2,
    precisionQDouble=3
#endif
  };
  Q_ENUM(precision);
  precision _selectedPrecision;
  Q_PROPERTY(precision selectedPrecision READ getselectedPrecision WRITE setselectedPrecision NOTIFY selectedPrecisionChange)
  precision getselectedPrecision() { return _selectedPrecision; }
  void setselectedPrecision(precision ps) { _selectedPrecision=ps; emit selectedPrecisionChange(); }

  void giveWork(MandelEvaluator *worker);
  void donePixel1(MandelEvaluator *me, int result);
public slots:
  void donePixel(MandelEvaluator *me, int result);
  void selectedPrecisionChanged();
signals:
  void selectedPaintStyleChanged();
  void selectedPrecisionChange();
protected:
  MandelMath::worker_multi *currentWorker; //for Position and Orbit
  MandelMath::worker_multi::Allocator *storeAllocator;
  MandelMath::worker_multi *storeWorker; //pointStore
  LaguerrePointStore *pointStore_;
  //constexpr static int MAX_ZOOM_IN_DOUBLE=55;//53;
  //MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  //LaguerrePoint *pointStore;
  int lastGivenPointIndex_;
  int effortBonus;
  //constexpr static int MAX_EFFORT=17;//18;
  int threadCount;
  MandelEvaluator **threads;
  QElapsedTimer timerWriteToImage;

  LaguerrePoint *wtiPoint;
  struct Position
  {
    static constexpr int LEN=2;
    MandelMath::worker_multi *worker;
    MandelMath::complex center;
    int step_log;
    double step_size__; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position(MandelMath::worker_multi::Allocator *allocator);
    ~Position();
    void assign(Position *src);
    //void setNumberType(MandelMath::worker_multi::Type ntype);
    void setView(const MandelMath::complex *c, double scale);
    void move(int delta_x, int delta_y);
    void scale(int inlog, int center_x, int center_y);
    void updateCachedDepth();
    void pixelXtoRE(int x, MandelMath::number_pointer result);
    void pixelYtoIM(int y, MandelMath::number_pointer result);
  } *position_;
  struct Orbit
  {
    MandelMath::worker_multi *currentWorker;
    //MandelMath::worker_multi::Allocator evaluatorAllocator;
    MandelEvaluator evaluator;
    MandelMath::worker_multi::Allocator pointAllocator;
    LaguerrePointStore pointDataStore;
    LaguerrePoint pointData;
    double first_mu_re_, first_mu_im, first_mum_re_, first_mum_im_;
    Orbit(MandelMath::worker_multi::Allocator *allocator);
    ~Orbit();
    constexpr static int LEN=MandelEvaluator::LEN+LaguerrePoint::LEN;
  } *orbit_;
  constexpr static int LEN=Params::LEN+LaguerrePoint::LEN+Position::LEN+Orbit::LEN  +4;
signals:
  void doneWork(MandelEvaluator *evaluator);
protected slots:
  void giveWork1(MandelEvaluator *me) { giveWork(me); }
};

#endif // LAGUERREMODEL_H

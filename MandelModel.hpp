#ifndef MANDELMODEL_H
#define MANDELMODEL_H

#include <QObject>
#include <QImage>

#include "ShareableImageWrapper.hpp"
#include "MandelEvaluator.hpp"

class MandelModel: public QObject
{
  Q_OBJECT
public:
  MandelModel();
  ~MandelModel();
  void transformStore(MandelMath::worker_multi *old_worker, MandelMath::worker_multi *old_sworker, MandelPointStore *old_store, int old_width, int old_height, const MandelMath::complex *old_c,
                      MandelMath::worker_multi *new_worker, MandelMath::worker_multi *new_sworker, MandelPointStore *new_store, int new_width, int new_height, const MandelMath::complex *new_c,
                      int inlog, int new_step_log);
  Q_INVOKABLE void setView_double(double c_re, double c_im, double scale);
  void setView_(const MandelMath::complex *c, double scale);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  //volano pouze ze setPrecision void setWorker(MandelMath::worker_multi *newWorker);
  void startNewEpoch();
  Q_INVOKABLE int writeToImage(ShareableImageWrapper img);
  void reimToPixel(int *circ_x, int *circ_y, const MandelMath::complex *c, MandelMath::number *tmp);
  Q_INVOKABLE void paintOrbit(ShareableImageWrapper image, int x, int y);
  Q_INVOKABLE QString pixelXtoRE_str(int x);
  Q_INVOKABLE QString pixelYtoIM_str(int y);
  Q_INVOKABLE QString getTimes();
  Q_INVOKABLE QString getTextXY();
  Q_INVOKABLE QString getTextInfoGen();
  Q_INVOKABLE QString getTextInfoSpec();
  MandelMath::worker_multi::Allocator *shareableViewInfoAllocator;
  ShareableViewInfo getViewInfo();
  Q_PROPERTY(ShareableViewInfo viewInfo READ getViewInfo CONSTANT)// WRITE setViewInfo NOTIFY viewInfoChanged)

  enum paintStyle
  {
    paintStyleKind=0,
    paintStyleCls=1,
    paintStyleExter=2,
    paintStyleInter=3,
    paintStyleNear=4,
    paintStyleFZ=5,
    paintStyleFC=6
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

  QVector<int> periodToIndexCache;
  int periodToIndex(int period);
  void giveWork(MandelEvaluator *worker);
  void donePixel1(MandelEvaluator *me);
public slots:
  void donePixel(MandelEvaluator *me);
  void selectedPrecisionChanged();
signals:
  void selectedPaintStyleChanged();
  void selectedPrecisionChange();
protected:
  MandelMath::worker_multi *currentWorker_; //for Position and Orbit
  MandelMath::worker_multi::Allocator *storeAllocator;
  MandelMath::worker_multi *storeWorker; //pointStore
  MandelPointStore *pointStore_;
  //constexpr static int MAX_ZOOM_IN_DOUBLE=55;//53;
  //MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  //MandelPoint *pointStore;
  int lastGivenPointIndex_;
  int effortBonus_;
  //constexpr static int MAX_EFFORT=17;//131072 iters;
  static constexpr int MAX_EFFORT=22;//
  int threadCount;
  MandelEvaluator **threads;
  QElapsedTimer timerWriteToImage;

  struct Position
  {
    static constexpr int LEN=2;
    MandelMath::worker_multi *worker;
#if NUMBER_DOUBLE_EXISTS
    //MandelMath::worker_multi_double number_worker_double_template;
#endif //NUMBER_DOUBLE_EXISTS
    //MandelMath::number_worker_ddouble number_worker_ddouble_template;
    //MandelMath::number_worker_multi number_worker_multi_template;
    MandelMath::complex center;
    int step_log;
    double step_size__; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position(MandelMath::worker_multi::Allocator *allocator);//: worker(worker), center(worker), ;
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
    MandelPointStore pointDataStore;
    MandelPoint pointData;
    MandelMath::complex lagu_c;
    MandelMath::complex lagu_r;
    MandelMath::number tmp;
    struct Bulb
    {
      bool valid;
      MandelMath::complex cb;
      MandelMath::complex rb;
      MandelMath::complex xc;
      MandelMath::complex baseZC;
      MandelMath::complex baseCC;
      int foundMult;
      bool is_card;
      Bulb(MandelMath::worker_multi::Allocator *allocator);
      ~Bulb();
      constexpr static int LEN=10;
    } bulb;
    Orbit(MandelMath::worker_multi::Allocator *allocator);
    ~Orbit();
    constexpr static int LEN=MandelEvaluator::LEN+MandelPoint::LEN+5+Bulb::LEN;
  } *orbit_;
  //constexpr static int INDEX_OF_POINTDATA=Position::LEN+MandelEvaluator::LEN;
  constexpr static int LEN=ShareableViewInfo::LEN+Position::LEN+Orbit::LEN  +4;
};



#endif // MANDELMODEL_H

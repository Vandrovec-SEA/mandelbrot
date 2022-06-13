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
  void transformStore(MandelPoint *old_store, int old_width, int old_height, MandelMath::number_store *old_cre, MandelMath::number_store *old_cim,
                      MandelPoint *new_store, int new_width, int new_height, const MandelMath::number_store *new_cre, const MandelMath::number_store *new_cim,
                      int inlog, int new_step_log);
  Q_INVOKABLE void setView_double(double c_re, double c_im, double scale);
  void setView_(const MandelMath::number_store *c_re, const MandelMath::number_store *c_im, double scale);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  void setWorker(MandelMath::number_worker *newWorker);
  void startNewEpoch();
  Q_INVOKABLE void writeToImage(ShareableImageWrapper img);
  void reimToPixel(int *circ_x, int *circ_y, const MandelMath::number_store *re, const MandelMath::number_store *im, MandelMath::number_store *tmp);
  Q_INVOKABLE void paintOrbit(ShareableImageWrapper image, int x, int y);
  Q_INVOKABLE QString pixelXtoRE_str(int x);
  Q_INVOKABLE QString pixelYtoIM_str(int y);
  Q_INVOKABLE QString getTimes();
  Q_INVOKABLE QString getTextXY();
  Q_INVOKABLE QString getTextInfoGen();
  Q_INVOKABLE QString getTextInfoSpec();
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
    precisionDDouble=1,
    precisionQDouble=2,
    precisionMulti=3
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
  //constexpr static int MAX_ZOOM_IN_DOUBLE=55;//53;
  //MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  MandelPoint *pointStore;
  int lastGivenPointIndex_;
  int effortBonus_;
  //constexpr static int MAX_EFFORT=17;//131072 iters;
  constexpr static int MAX_EFFORT=22;//
  int threadCount;
  MandelEvaluator *threads;

  struct Position
  {
    MandelMath::number_worker *worker;
#if NUMBER_DOUBLE_EXISTS
    MandelMath::number_worker_double number_worker_double_template;
#endif //NUMBER_DOUBLE_EXISTS
    MandelMath::number_worker_ddouble number_worker_ddouble_template;
    MandelMath::number_worker_multi number_worker_multi_template;
    MandelMath::number_store center_re_s;
    MandelMath::number_store center_im_s;
    MandelMath::number_place center_re_p;
    MandelMath::number_place center_im_p;
    int step_log;
    double step_size__; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position();
    ~Position();
    void setNumberType(MandelMath::number_worker::Type ntype);
    void setView(const MandelMath::number_store *c_re, const MandelMath::number_store *c_im, double scale);
    void move(int delta_x, int delta_y);
    void scale(int inlog, int center_x, int center_y);
    void updateCachedDepth();
    void pixelXtoRE(int x, MandelMath::number_store *result);
    void pixelYtoIM(int y, MandelMath::number_store *result);
  } position;
  struct
  {
    MandelMath::number_worker *worker;
    MandelEvaluator evaluator;
    MandelPoint pointData;
    MandelMath::number_store lagu_c_re_, lagu_c_im;
    MandelMath::number_store lagu_r_re, lagu_r_im;
    MandelMath::number_store tmp;
    MandelMath::number_place lagu_c_re_p, lagu_c_im_p, lagu_r_re_p, lagu_r_im_p, tmp_p;
    struct
    {
      bool valid;
      MandelMath::number_store cb_re, cb_im;
      MandelMath::number_store rb_re_, rb_im;
      MandelMath::number_store xc_re, xc_im;
      MandelMath::number_store baseZC_re, baseZC_im;
      MandelMath::number_store baseCC_re, baseCC_im;
      int foundMult;
      bool is_card;
      MandelMath::number_place cb_re_p, cb_im_p, rb_re_p, rb_im_p, xc_re_p, xc_im_p;
      MandelMath::number_place baseZC_re_p, baseZC_im_p, baseCC_re_p, baseCC_im_p;
    } bulb;
  } orbit;
};



#endif // MANDELMODEL_H

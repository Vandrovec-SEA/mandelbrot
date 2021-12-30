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
  Q_INVOKABLE void setView(double c_re, double c_im, double scale);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  void startNewEpoch();
  Q_INVOKABLE void writeToImage(ShareableImageWrapper img);
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
    paintStyleNear=4
  };
  Q_ENUM(paintStyle);
  paintStyle _selectedPaintStyle;
  Q_PROPERTY(paintStyle selectedPaintStyle READ getselectedPaintStyle WRITE setselectedPaintStyle NOTIFY selectedPaintStyleChanged)
  paintStyle getselectedPaintStyle() { return _selectedPaintStyle; }
  void setselectedPaintStyle(paintStyle ps) { _selectedPaintStyle=ps; }

  QVector<int> periodToIndexCache;
  int periodToIndex(int period);
  void giveWork(MandelEvaluator *worker);
  void donePixel1(MandelEvaluator *me);
public slots:
  void donePixel(MandelEvaluator *me);
signals:
  void selectedPaintStyleChanged();
protected:
  constexpr static int MAX_ZOOM_IN_DOUBLE=55;//53;
  //MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  MandelPoint *pointStore;
  int lastGivenPointIndex_;
  int effortBonus_;
  constexpr static int MAX_EFFORT=17;//18;
  int threadCount;
  MandelEvaluator *threads;

  struct Position
  {
    //MandelMath::number_store center_re_s;
    //MandelMath::number_store center_im_s;
    MandelMath::number_worker *worker;
    MandelMath::number_worker_double number_worker_double_template;
    MandelMath::number_worker_ddouble number_worker_ddouble_template;
    MandelMath::number_worker_multi number_worker_multi_template;
    MandelMath::number_store center_re_s;
    MandelMath::number_store center_im_s;
    int step_log;
    double step_size__; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position();
    ~Position();
    void setNumberType(MandelMath::number_worker::Type ntype);
    void setView(double c_re, double c_im, double scale);
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
  } orbit;
};



#endif // MANDELMODEL_H

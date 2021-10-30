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
  void transformStore(MandelPoint *old_store, int old_width, int old_height, MandelMath::number *old_cre, MandelMath::number *old_cim,
                      MandelPoint *new_store, int new_width, int new_height, const MandelMath::number_store *new_cre, const MandelMath::number_store *new_cim,
                      int inlog, int new_step_log);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  void startNewEpoch();
  Q_INVOKABLE void writeToImage(ShareableImageWrapper img);
  Q_INVOKABLE QString pixelXtoRE_str(int x);
  Q_INVOKABLE QString pixelYtoIM_str(int y);
  Q_INVOKABLE QString getTimes();

  void giveWork(MandelEvaluator *worker);
  void donePixel1(MandelEvaluator *me);
public slots:
  void donePixel(MandelEvaluator *me);
protected:
  MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  MandelPoint *pointStore;
  int lastGivenPointIndex_;
  int effortBonus;
  constexpr static int MAX_EFFORT=18;
  int threadCount;
  MandelEvaluator *threads;

  struct Position
  {
    //MandelMath::number_store center_re_s;
    //MandelMath::number_store center_im_s;
    MandelMath::number_any center_re_n;
    MandelMath::number_any center_im_n;
    int step_log;
    double step_size__; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position();
    void setNumberType(MandelMath::number_store::DbgType ntype);
    void move(int delta_x, int delta_y);
    void scale(int inlog, int center_x, int center_y);
    void updateCachedDepth();
    void pixelXtoRE(int x, MandelMath::number_any *result);
    void pixelYtoIM(int y, MandelMath::number_any *result);
  } position;
};

#endif // MANDELMODEL_H

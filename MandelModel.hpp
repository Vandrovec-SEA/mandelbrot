#ifndef MANDELMODEL_H
#define MANDELMODEL_H

#include <QObject>
#include <QImage>

#include "MandelEvaluator.hpp"
#include "ShareableImageWrapper.hpp"

class MandelModel: public QObject
{
  Q_OBJECT
public:
  MandelModel();
  ~MandelModel();
  void transformStore(MandelPoint *old_store, int old_width, int old_height, double old_cre, double old_cim,
                      MandelPoint *new_store, int new_width, int new_height, double new_cre, double new_cim,
                      int inlog, double new_step);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  void startNewEpoch();
  Q_INVOKABLE void writeToImage(ShareableImageWrapper img);
  Q_INVOKABLE double pixelXtoRE(int x);
  Q_INVOKABLE double pixelYtoIM(int y);
  Q_INVOKABLE QString getTimes();

  void giveWork(MandelEvaluator *worker);
  void donePixel1(MandelEvaluator *me);
public slots:
  void donePixel(MandelEvaluator *me);
protected:
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
    double center_re;
    double center_im;
    int step_log;
    double step_size;
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position();
    void move(int delta_x, int delta_y);
    void scale(int inlog, int center_x, int center_y);
    void updateCachedDepth();
  } position;
};

#endif // MANDELMODEL_H

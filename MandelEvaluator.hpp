#ifndef MANDELEVALUATOR_HPP
#define MANDELEVALUATOR_HPP

#include <QThread>

void dbgPoint();

struct MandelPoint
{
  MandelPoint();
  enum State { stUnknown, stOutside, stMaxIter } state;
  double zr, zi;
  int iter;
  //MandelPoint &operator =(MandelPoint &src) = default;
  void reset();
};

class MandelEvaluator: public QThread
{
  Q_OBJECT
public:
  MandelEvaluator();
  ~MandelEvaluator();
  void setHint(int hint);
  static void simple(double cr, double ci, MandelPoint &data, int maxiter);

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
    double cr;
    double ci;
    int epoch;
    int pixelIndex;
    int maxiter;
    ComputeParams();
  } currentParams;
  MandelPoint currentData;

  bool startCompute(const ComputeParams &params, const MandelPoint *data, bool no_quick_route);
protected slots:
  void doCompute();
signals:
  void doneCompute(MandelEvaluator *me);
};

#endif // MANDELEVALUATOR_HPP

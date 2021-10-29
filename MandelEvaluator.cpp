#include "MandelEvaluator.hpp"

void doNothing(int &x)
{
  x++;
}

void dbgPoint()
{
  int x=3;
  doNothing(x);
}



MandelPoint::MandelPoint()
{
  reset();
}

void MandelPoint::reset()
{
  zr=0;
  zi=0;
  state=State::stUnknown;
  iter=0;
}



MandelEvaluator::MandelEvaluator(): QThread(nullptr), currentParams()
{
  QThread::start(QThread::Priority::LowestPriority);
  wantStop=false;
  pointsComputed=0;
  timeOuterTotal=0;
  timeInnerTotal=0;
  timeInvokePostTotal=0;
  timeInvokeSwitchTotal=0;
  QObject::moveToThread(this);
}

MandelEvaluator::~MandelEvaluator()
{
}

void MandelEvaluator::setHint(int hint)
{
  currentParams.cr=hint;
  currentParams.ci=0;
  currentData.iter=-2;
}

void MandelEvaluator::simple(double cr, double ci, MandelPoint &data, int maxiter)
{
  double zr=0;
  double zi=0;
  for (int iter=0; iter<maxiter; iter++)
  {
    if (zr*zr+zi*zi>4)
    {
      data.state=MandelPoint::State::stOutside;
      data.iter=iter;
      data.zr=zr;
      data.zi=zi;
      return;
    };
    double tmp=zr*zr-zi*zi+cr;
    zi=2*zr*zi+ci;
    zr=tmp;
  }
  //data.state=MandelPoint::State::stMaxIter;
  data.iter=maxiter;
  data.zr=zr;
  data.zi=zi;
}

bool MandelEvaluator::startCompute(const ComputeParams &params, const MandelPoint *data, bool no_quick_route)
{
  /*currentParams.cr=cr;
  currentParams.ci=ci;
  currentParams.epoch=epoch;
  currentParams.pixelIndex=pixelIndex;*/
  currentParams=params;
  currentData=*data;
  if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
  {
    simple(currentParams.cr, currentParams.ci, currentData, currentParams.maxiter);
    pointsComputed++;
    return false;
  };
  timeInvoke.start();
  QMetaObject::invokeMethod(this,
                            //[cr, ci, pointData, worker]() { worker->startCompute(cr, ci, *pointData); },
                            &MandelEvaluator::doCompute,
                            Qt::ConnectionType::QueuedConnection);
  timeInvokePostTotal+=timeInvoke.nsecsElapsed();
  return true;
}

void MandelEvaluator::doCompute()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  simple(currentParams.cr, currentParams.ci, currentData, currentParams.maxiter);
  pointsComputed++;
  //currentData.state=MandelPoint::State::stOutside;
  //currentData.iter=3;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}


MandelEvaluator::ComputeParams::ComputeParams()
{
  cr=0; ci=0;
  epoch=-1;
  pixelIndex=-1;
  maxiter=1;
}

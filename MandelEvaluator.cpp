#include "MandelEvaluator.hpp"


MandelMath::number_any::number_any():
  my_store(), impl(nullptr), store(&my_store)
{

}

MandelMath::number_any::number_any(number_store *store):
  my_store(), impl(nullptr), store(store)
{

}

MandelMath::number_any::number_any(number_any *src):
  my_store(), impl(nullptr), store(&my_store)
{
  reinit(src->ntype_());
  if (impl==nullptr)
    dbgPoint();
  else
    impl->assign(store, src->store);
}

MandelMath::number_any::~number_any()
{
  if (store==&my_store)
  {
    if (impl==nullptr)
      dbgPoint();
    else
      impl->cleanup(store);
  }
  impl=nullptr;
}

/*
MandelMath::number_any::number_any(MandelMath::number_store::DbgType ntype, MandelMath::number_store *src):
  d(src), dd(src), m(src)
{
  switch (ntype)
  {
    case MandelMath::number_store::DbgType::typeDouble:
    {
      impl=&d;
    } break;
    case MandelMath::number_store::DbgType::typeDDouble:
    {
      impl=&dd;
    } break;
    case MandelMath::number_store::DbgType::typeMulti:
    {
      impl=&m;
    } break;
    case MandelMath::number_store::DbgType::typeEmpty:
    {
      impl=nullptr;
    } break;
  }
}
*/

void MandelMath::number_any::reinit(MandelMath::number_worker *worker)
{ //TODO: should try to convert old value to new type
  if (worker==impl)
    return;
  if (impl)
    impl->cleanup(store);
  impl=worker;
  if (impl)
    impl->init(store);
}

MandelMath::number_worker *MandelMath::number_any::ntype_()
{
  return impl;
}



MandelPoint::MandelPoint(): zr_(), zi_()
{
  //reset();
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::assign(MandelMath::number_worker *worker, const MandelPoint &src)
{
  worker->assign(&zr_, &src.zr_);
  worker->assign(&zi_, &src.zi_);
  state=src.state;
  iter=src.iter;
}

void MandelPoint::init(MandelMath::number_worker *worker)
{
  worker->init(&zr_, 0.0);
  worker->init(&zi_, 0.0);
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::zero(MandelMath::number_worker *worker)
{
  worker->zero(&zr_, 0.0);
  worker->zero(&zi_, 0.0);
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::cleanup(MandelMath::number_worker *worker)
{
  if (worker==nullptr)
    dbgPoint();
  worker->cleanup(&zr_);
  worker->cleanup(&zi_);
}





MandelEvaluator::MandelEvaluator(): QThread(nullptr),
  currentWorker(nullptr),
  currentParams()//,
  //data_zr_n(&currentData.zr_),
  //data_zi_n(&currentData.zi_)
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
  switchType(nullptr);
}

void MandelEvaluator::simple_double(double cr, double ci, MandelPoint &data, int maxiter)
{
  double zr=0;
  double zi=0;
  for (int iter=0; iter<maxiter; iter++)
  {
    if (zr*zr+zi*zi>4)
    {
      data.state=MandelPoint::State::stOutside;
      data.iter=iter;
      data.zr_.as.doubl=zr;
      data.zi_.as.doubl=zi;
      return;
    };
    double tmp=zr*zr-zi*zi+cr;
    zi=2*zr*zi+ci;
    zr=tmp;
  }
  //data.state=MandelPoint::State::stMaxIter;
  data.iter=maxiter;
  data.zr_.as.doubl=zr;
  data.zi_.as.doubl=zi;
}

void MandelEvaluator::simple_ddouble(MandelMath::dd_real *cr, MandelMath::dd_real *ci, MandelPoint &data, int maxiter)
{
  MandelMath::dd_real zr;
  MandelMath::dd_real zi;
  MandelMath::dd_real r2;
  MandelMath::dd_real i2;
  MandelMath::dd_real t;
  for (int iter=0; iter<maxiter; iter++)
  {
    r2.assign(zr); r2.sqr();
    i2.assign(zi); i2.sqr();
    t.assign(r2); t.add(i2.hi, i2.lo);
    if (t.hi>4)
    {
      data.state=MandelPoint::State::stOutside;
      data.iter=iter;
      data.zr_.as.ddouble_.dd->assign(zr);
      data.zi_.as.ddouble_.dd->assign(zi);
      return;
    };
    t.assign(r2); t.add(-i2.hi, -i2.lo); t.add(cr->hi, cr->lo); //double tmp=zr*zr-zi*zi+cr;
    zi.mul(2*zr.hi, 2*zr.lo); zi.add(ci->hi, ci->lo);
    zr.assign(t);
  }
  data.iter=maxiter;
  data.zr_.as.ddouble_.dd->assign(zr);
  data.zi_.as.ddouble_.dd->assign(zi);
}

void MandelEvaluator::simple_multi(MandelMath::multiprec *cr, MandelMath::multiprec *ci, MandelPoint &data, int maxiter)
{
  (void)cr;
  (void)ci;
  data.state=MandelPoint::State::stMaxIter;
  data.iter=maxiter;
}

void MandelEvaluator::switchType(MandelMath::number_worker *worker)
{
  if (worker==currentWorker)
    return;
  if (currentWorker!=nullptr)
  {
    currentWorker->cleanup(&currentParams.ci_s);
    currentWorker->cleanup(&currentParams.cr_s);
    currentData.cleanup(currentWorker);
    currentWorker->cleanup(&data_zr_s);
    currentWorker->cleanup(&data_zi_s);
  }
  if (worker)
  if (worker)
  {
    worker->init(&data_zr_s);
    worker->init(&data_zi_s);
    currentData.init(worker);
    worker->init(&currentParams.ci_s);
    worker->init(&currentParams.cr_s);
  }
  currentWorker=worker;
}

bool MandelEvaluator::startCompute(const MandelPoint *data, bool no_quick_route)
{
  //currentParams=params;
  /*data_zr_n.reinit(currentParams.cr_n.ntype());
  data_zi_n.reinit(currentParams.ci_n.ntype());
  data_z_tmp1.reinit(currentParams.cr_n.ntype());
  data_z_tmp2.reinit(currentParams.cr_n.ntype());*/
  if (currentWorker==nullptr)
  {
    dbgPoint();
    currentData.state=MandelPoint::State::stMaxIter;
    return false;
  }
  currentData.assign(currentWorker, *data);
  if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
  {
    //simple_double(currentParams.cr_n.impl->store->as.doubl, currentParams.ci_n.impl->store->as.doubl, currentData, currentParams.maxiter);
    evaluate();
    pointsComputed++;
    return false;
  };
  timeInvoke.start();
  QMetaObject::invokeMethod(this,
                            &MandelEvaluator::doCompute,
                            Qt::ConnectionType::QueuedConnection);
  timeInvokePostTotal+=timeInvoke.nsecsElapsed();
  return true;
}

void MandelEvaluator::evaluate()
{
  MandelMath::complex c(currentWorker, &currentParams.cr_s, &currentParams.ci_s, true);
  MandelMath::complex z(currentWorker, &this->data_zr_s, &this->data_zi_s, true);
  currentWorker->assign(z.re_s, &currentData.zr_);
  currentWorker->assign(z.im_s, &currentData.zi_);
  for (int iter=currentData.iter; iter<currentParams.maxiter; iter++)
  {
    const MandelMath::number_store *magtmp=z.getMagTmp();
    if (currentWorker->toDouble(magtmp)>4)
    {
      currentData.state=MandelPoint::State::stOutside;
      currentData.iter=iter;
      currentWorker->assign(&currentData.zr_, z.re_s);
      currentWorker->assign(&currentData.zi_, z.im_s);
      return;
    };
    z.sqr();
    z.add(&c);
  }
  //data.state=MandelPoint::State::stMaxIter;
  currentData.iter=currentParams.maxiter;
  currentWorker->assign(&currentData.zr_, z.re_s);
  currentWorker->assign(&currentData.zi_, z.im_s);
}

void MandelEvaluator::doCompute()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  //simple_double(currentParams.cr_n.impl->store->as.doubl, currentParams.ci_n.impl->store->as.doubl, currentData, currentParams.maxiter);
  evaluate();
  pointsComputed++;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}


MandelEvaluator::ComputeParams::ComputeParams():
  cr_s(),
  ci_s()
{
  epoch=-1;
  pixelIndex=-1;
  maxiter=1;
}

#include "MandelEvaluator.hpp"


MandelMath::number_any::number_any():
  my_store(), d(&my_store), dd(&my_store), m(&my_store), impl(nullptr)
{

}

MandelMath::number_any::number_any(number_store *store):
  my_store(), d(store), dd(store), m(store), impl(nullptr)
{

}

MandelMath::number_any::number_any(number_any *src):
  d(&my_store), dd(&my_store), m(&my_store)
{
  switch (src->ntype())
  {
    case MandelMath::number_store::DbgType::typeDouble:
    {
      my_store.init_double();
      my_store.assign_double(*src->impl->store);
      impl=&d;
    } break;
    case MandelMath::number_store::DbgType::typeDDouble:
    {
      my_store.init_ddouble();
      my_store.assign_ddouble(*src->impl->store);
      impl=&dd;
    } break;
    case MandelMath::number_store::DbgType::typeMulti:
    {
      my_store.init_multi();
      my_store.assign_multi(*src->impl->store);
      impl=&m;
    } break;
    case MandelMath::number_store::DbgType::typeEmpty:
    {
      dbgPoint();
      impl=nullptr;
    } break;
  }
}

MandelMath::number_any::~number_any()
{
  if (impl && (impl->store==&my_store))
  {
    impl->cleanup();
    impl=nullptr;
  }
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

void MandelMath::number_any::reinit(MandelMath::number_store::DbgType ntype)
{ //TODO: should try to convert old value to new type
  switch (ntype)
  {
    case MandelMath::number_store::DbgType::typeDouble:
    {
      if (impl==&d)
        return;
      if (impl==&dd)
        dd.store->cleanup_ddouble();
      else if (impl==&m)
        m.store->cleanup_multi();
      else if (impl!=nullptr)
        dbgPoint();
      impl=&d;
      d.store->init_double();
    } break;
    case MandelMath::number_store::DbgType::typeDDouble:
    {
      if (impl==&dd)
        return;
      if (impl==&d)
        dd.store->cleanup_double();
      else if (impl==&m)
        m.store->cleanup_multi();
      else if (impl!=nullptr)
        dbgPoint();
      impl=&dd;
      dd.store->init_ddouble();
    } break;
    case MandelMath::number_store::DbgType::typeMulti:
    {
      if (impl==&m)
        return;
      if (impl==&dd)
        dd.store->cleanup_ddouble();
      else if (impl==&d)
        d.store->cleanup_double();
      else if (impl!=nullptr)
        dbgPoint();
      impl=&m;
      d.store->init_multi();
    } break;
    case MandelMath::number_store::DbgType::typeEmpty:
    {
      if (impl==nullptr)
        return;
      if (impl==&d)
        d.store->cleanup_double();
      else if (impl==&dd)
        dd.store->cleanup_ddouble();
      else if (impl==&m)
        m.store->cleanup_multi();
      else
        dbgPoint();
      impl=nullptr;
    } break;
  }
}


MandelPoint::MandelPoint(): zr_(), zi_()
{
  //reset();
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::assign_double(const MandelPoint &src)
{
  zr_.assign_double(src.zr_);
  zi_.assign_double(src.zi_);
  state=src.state;
  iter=src.iter;
}

void MandelPoint::assign_ddouble(const MandelPoint &src)
{
  zr_.assign_ddouble(src.zr_);
  zi_.assign_ddouble(src.zi_);
  state=src.state;
  iter=src.iter;
}

void MandelPoint::assign_multi(const MandelPoint &src)
{
  zr_.assign_multi(src.zr_);
  zi_.assign_multi(src.zi_);
  state=src.state;
  iter=src.iter;
}

void MandelPoint::init_double()
{
  //zr_.cleanup_double();
  zr_.init_double(0.0);
  //zi_.cleanup_double();
  zi_.init_double(0.0);
  state=State::stUnknown;
  iter=0;
}



MandelEvaluator::MandelEvaluator(): QThread(nullptr), currentParams(),
  data_zr_n(&currentData.zr_),
  data_zi_n(&currentData.zi_)
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
  //data_zi_n won't autoclean because the store is external
  data_zi_n.reinit(MandelMath::number_store::DbgType::typeEmpty);
  data_zr_n.reinit(MandelMath::number_store::DbgType::typeEmpty);
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

void MandelEvaluator::simple_ddouble(MandelMath::double_double cr, MandelMath::double_double ci, MandelPoint &data, int maxiter)
{
  (void)cr;
  (void)ci;
  data.state=MandelPoint::State::stMaxIter;
  data.iter=maxiter;
}

void MandelEvaluator::simple_multi(MandelMath::multiprec *cr, MandelMath::multiprec *ci, MandelPoint &data, int maxiter)
{
  (void)cr;
  (void)ci;
  data.state=MandelPoint::State::stMaxIter;
  data.iter=maxiter;
}

bool MandelEvaluator::startCompute(const MandelPoint *data, bool no_quick_route)
{
  //currentParams=params;
  data_zr_n.reinit(currentParams.cr_n.ntype());
  data_zi_n.reinit(currentParams.ci_n.ntype());
  switch (currentParams.cr_n.ntype())
  {
    case MandelMath::number_store::DbgType::typeDouble:
    {
      currentData.assign_double(*data);
      if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
      {
        simple_double(currentParams.cr_n.impl->store->as.doubl, currentParams.ci_n.impl->store->as.doubl, currentData, currentParams.maxiter);
        pointsComputed++;
        return false;
      };
      timeInvoke.start();
      QMetaObject::invokeMethod(this,
                                &MandelEvaluator::doCompute_double,
                                Qt::ConnectionType::QueuedConnection);
      timeInvokePostTotal+=timeInvoke.nsecsElapsed();
      return true;
    } break;
    case MandelMath::number_store::DbgType::typeDDouble:
    {
      currentData.assign_ddouble(*data);
      if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
      {
        simple_ddouble(currentParams.cr_n.impl->store->as.ddouble, currentParams.ci_n.impl->store->as.ddouble, currentData, currentParams.maxiter);
        pointsComputed++;
        return false;
      };
      timeInvoke.start();
      QMetaObject::invokeMethod(this,
                                &MandelEvaluator::doCompute_ddouble,
                                Qt::ConnectionType::QueuedConnection);
      timeInvokePostTotal+=timeInvoke.nsecsElapsed();
      return true;
    } break;
    case MandelMath::number_store::DbgType::typeMulti:
    {
      currentData.assign_multi(*data);
      if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
      {
        simple_multi(currentParams.cr_n.impl->store->as.multi.bytes, currentParams.ci_n.impl->store->as.multi.bytes, currentData, currentParams.maxiter);
        pointsComputed++;
        return false;
      };
      timeInvoke.start();
      QMetaObject::invokeMethod(this,
                                &MandelEvaluator::doCompute_multi,
                                Qt::ConnectionType::QueuedConnection);
      timeInvokePostTotal+=timeInvoke.nsecsElapsed();
      return true;
    } break;
   case MandelMath::number_store::DbgType::typeEmpty:
      dbgPoint();
  }

  dbgPoint();
  currentData.state=MandelPoint::State::stMaxIter;
  return false;
}

void MandelEvaluator::doCompute_double()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  simple_double(currentParams.cr_n.impl->store->as.doubl, currentParams.ci_n.impl->store->as.doubl, currentData, currentParams.maxiter);
  pointsComputed++;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}

void MandelEvaluator::doCompute_ddouble()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  simple_ddouble(currentParams.cr_n.impl->store->as.ddouble, currentParams.ci_n.impl->store->as.ddouble, currentData, currentParams.maxiter);
  pointsComputed++;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}

void MandelEvaluator::doCompute_multi()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  simple_multi(currentParams.cr_n.impl->store->as.multi.bytes, currentParams.ci_n.impl->store->as.multi.bytes, currentData, currentParams.maxiter);
  pointsComputed++;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}


MandelEvaluator::ComputeParams::ComputeParams():
  cr_n(),
  ci_n()
{
  epoch=-1;
  pixelIndex=-1;
  maxiter=1;
}

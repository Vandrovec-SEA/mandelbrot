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
  d(&my_store), dd(&my_store), m(&my_store), impl(nullptr)
{
  reinit(src->ntype());
  switch (src->ntype())
  {
    case MandelMath::number::Type::typeDouble:
    {
      my_store.assign_double(*src->impl->store);
    } break;
    case MandelMath::number::Type::typeDDouble:
    {
      my_store.assign_ddouble(*src->impl->store);
    } break;
    case MandelMath::number::Type::typeMulti:
    {
      my_store.assign_multi(*src->impl->store);
    } break;
    case MandelMath::number::Type::typeEmpty:
    {
      dbgPoint();
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

void MandelMath::number_any::reinit(MandelMath::number::Type ntype)
{ //TODO: should try to convert old value to new type
  MandelMath::number::Type old_ntype=this->ntype();
  if (ntype==old_ntype)
    return;
  if (impl)
    impl->store->cleanup(old_ntype);
  switch (ntype)
  {
    case MandelMath::number::Type::typeDouble:
    {
      impl=&d;
    } break;
    case MandelMath::number::Type::typeDDouble:
    {
      impl=&dd;
    } break;
    case MandelMath::number::Type::typeMulti:
    {
      impl=&m;
    } break;
    case MandelMath::number::Type::typeEmpty:
    {
      impl=nullptr;
    } break;
  }
  if (impl)
    impl->store->init(ntype);

}

MandelMath::number::Type MandelMath::number_any::ntype()
{
  if (impl==nullptr)
    return MandelMath::number::Type::typeEmpty;
  else if (impl==&d)
    return MandelMath::number::Type::typeDouble;
  else if (impl==&dd)
    return MandelMath::number::Type::typeDDouble;
  else if (impl==&m)
    return MandelMath::number::Type::typeMulti;
  else
    dbgPoint();
  return MandelMath::number::Type::typeEmpty;
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

void MandelPoint::init(MandelMath::number::Type ntype)
{
  zr_.init(ntype, 0.0);
  zi_.init(ntype, 0.0);
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::zero(MandelMath::number::Type ntype)
{
  zr_.zero(ntype, 0.0);
  zi_.zero(ntype, 0.0);
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::cleanup(MandelMath::number::Type ntype)
{
  zr_.cleanup(ntype);
  zi_.cleanup(ntype);
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
  data_zi_n.reinit(MandelMath::number::Type::typeEmpty);
  data_zr_n.reinit(MandelMath::number::Type::typeEmpty);
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

bool MandelEvaluator::startCompute(const MandelPoint *data, bool no_quick_route)
{
  //currentParams=params;
  data_zr_n.reinit(currentParams.cr_n.ntype());
  data_zi_n.reinit(currentParams.ci_n.ntype());
  data_z_tmp1.reinit(currentParams.cr_n.ntype());
  data_z_tmp2.reinit(currentParams.cr_n.ntype());
  switch (currentParams.cr_n.ntype())
  {
    case MandelMath::number::Type::typeDouble:
    {
      currentData.assign_double(*data);
      if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
      {
        //simple_double(currentParams.cr_n.impl->store->as.doubl, currentParams.ci_n.impl->store->as.doubl, currentData, currentParams.maxiter);
        evaluate();
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
    case MandelMath::number::Type::typeDDouble:
    {
      currentData.assign_ddouble(*data);
      if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
      {
        //simple_ddouble(currentParams.cr_n.impl->store->as.ddouble_.dd, currentParams.ci_n.impl->store->as.ddouble_.dd, currentData, currentParams.maxiter);
        evaluate();
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
    case MandelMath::number::Type::typeMulti:
    {
      currentData.assign_multi(*data);
      if (!no_quick_route && (currentParams.maxiter-currentData.iter<=1000))
      {
        //simple_multi(currentParams.cr_n.impl->store->as.multi_.bytes, currentParams.ci_n.impl->store->as.multi_.bytes, currentData, currentParams.maxiter);
        evaluate();
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
   case MandelMath::number::Type::typeEmpty:
      dbgPoint();
  }

  dbgPoint();
  currentData.state=MandelPoint::State::stMaxIter;
  return false;
}

void MandelEvaluator::evaluate()
{
  MandelMath::complex c(currentParams.cr_n.impl, currentParams.ci_n.impl, nullptr, nullptr);
  MandelMath::complex z(this->data_zr_n.impl, this->data_zi_n.impl, data_z_tmp1.impl, data_z_tmp2.impl);
  for (int iter=currentData.iter; iter<currentParams.maxiter; iter++)
  {
    MandelMath::number *magtmp=z.getMagTmp();
    if (magtmp->toDouble()>4)
    {
      currentData.state=MandelPoint::State::stOutside;
      currentData.iter=iter;
      z.re->assignTo(&currentData.zr_);
      z.im->assignTo(&currentData.zi_);
      return;
    };
    z.sqr();
    z.add(&c);
  }
  //data.state=MandelPoint::State::stMaxIter;
  currentData.iter=currentParams.maxiter;
  z.re->assignTo(&currentData.zr_);
  z.im->assignTo(&currentData.zi_);
}

void MandelEvaluator::doCompute_double()
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

void MandelEvaluator::doCompute_ddouble()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  //simple_ddouble(currentParams.cr_n.impl->store->as.ddouble_.dd, currentParams.ci_n.impl->store->as.ddouble_.dd, currentData, currentParams.maxiter);
  evaluate();
  pointsComputed++;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}

void MandelEvaluator::doCompute_multi()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  //simple_multi(currentParams.cr_n.impl->store->as.multi_.bytes, currentParams.ci_n.impl->store->as.multi_.bytes, currentData, currentParams.maxiter);
  evaluate();
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

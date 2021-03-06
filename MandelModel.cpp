#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "MandelModel.hpp"
#include "MandelEvaluator.hpp"

#define CURRENT_STORE_DIRECT 0 //either both or none of store and numbers must be direct, can't mix
#define UPDATE_CACHED_MOD 0 //works until precision doesn't allow

MandelModel::MandelModel(): QObject(),
  storeAllocator(nullptr), storeWorker(nullptr), pointStore(nullptr),
  precisionRecord(nullptr)
{
  unsigned int oldcw; //524319 = 0x8001F = mask all interrupts, 80bit precision
  MandelMath::fpu_fix_start(&oldcw);

   /*__float128 p=1.4;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;*/

  //selftest of dd_real.split / x87 rounding
  MandelMath::dd_real dd1, dd2;
  /*dd2.hi=0x1fffffff;
  ((uint64_t &)dd2.hi)&=0xfffffffffc000000ull;
  dd2.lo=0x1fffffff-dd2.hi;
  double dd2hi=0x1fffffff;
  ((uint64_t &)dd2hi)&=0xfffffffffc000000ull;
  dd2hi+=3;
  dd1.split(dd2hi);*/
  dd1.split(0x1fffffff); //-> 0x20000000 - 0x1
  dd1.split(0x20000000); //-> 0x20000000 + 0x0
  dd1.split(0x20000001); //-> 0x20000000 + 0x1
  dd1.split(0x20000002); //-> 0x20000000 + 0x2
  dd1.split(0x20000003); //-> 0x20000000 + 0x3
  dd1.split(0x20000004); //-> 0x20000000 + 0x4
  dd1.split(0x20000005); //-> 0x20000000 + 0x5
  dd1.split(0x20000006); //-> 0x20000000 + 0x6
  dd1.split(0x20000007); //-> 0x20000000 + 0x7
  dd1.split(0x20000008); //-> 0x20000000 + 0x8                       100 0000 0000 002                                     80 0000 0000 0010
  dd1.split(0x20000009); //-> 0x20000010 - 0x7  0x20000009*0x8000001=100 0000 6800 0009 @52b= 100 0000 0000 0000 -20000009=FF FFFF DFFF FFF7 rounds to FF FFFF DFFF FFF0
  dd1.split(0x2000000a); //-> 0x20000010 - 0x6
  dd1.split(0x2000000b); //-> 0x20000010 - 0x5
  dd1.split(0x2000000c); //-> 0x20000010 - 0x4
  dd1.split(0x2000000d); //-> 0x20000010 - 0x3
  dd1.split(0x2000000e); //-> 0x20000010 - 0x2
  dd1.split(0x2000000f); //-> 0x20000010 - 0x1
  dd1.split(0x20000010); //-> 0x20000010 + 0x0
  dd1.split(0x20000011); //-> 0x20000010 + 0x1
  dd1.split(-0.00043541188577694845);

  _selectedPaintStyle=paintStyleCls;//Kind;
  _selectedPrecision=precisionDouble;
  epoch=1;
  imageWidth=0;
  imageHeight=0;
  //pointStore=nullptr;
  nextGivenPointIndex=0;
  effortBonus=0;
  //threadCount=4;
  _threadsWorking=0;
  QObject::connect(this, &MandelModel::selectedPrecisionChange,
                   this, &MandelModel::selectedPrecisionChanged);
  selectedPrecisionChanged();
}

MandelModel::~MandelModel()
{
  epoch=(epoch%2000000000)+1;
  delete precisionRecord;
  delete storeAllocator;
  delete storeWorker;
  storeWorker=nullptr;

  delete[] pointStore;
  pointStore=nullptr;
  imageWidth=0;
  imageHeight=0;
}

QString MandelModel::pixelXtoRE_str(int x)
{
  MandelMath::number num(precisionRecord->currentWorker.get());
  num.assign(precisionRecord->position.center.re);
  num.add_double((x - imageWidth/2)*precisionRecord->position.step_size__);
  QString result=precisionRecord->currentWorker->toString(num.ptr);
  return result;
}

QString MandelModel::pixelYtoIM_str(int y)
{
//  return (y - imageHeight/2)*position.step_size+position.center_im;
  MandelMath::number num(precisionRecord->currentWorker.get());
  num.assign(precisionRecord->position.center.im);
  num.add_double((y - imageHeight/2)*precisionRecord->position.step_size__);
  QString result=precisionRecord->currentWorker->toString(num.ptr);
  return result;
}

QString MandelModel::getTimes()
{
  QString result;
  /*
  for (int t=0; t<precisionRecord->threadCount; t++)
    result+=QString("%1-%2[%3,%4],").
        arg((precisionRecord->threads[t]->timeOuterTotal)/1000000000.0, 0, 'f', 3).
        arg((precisionRecord->threads[t]->timeInnerTotal)/1000000000.0, 0, 'f', 3).
        arg((precisionRecord->threads[t]->timeInvokePostTotal)/1000000000.0, 0, 'f', 3).
        arg((precisionRecord->threads[t]->timeInvokeSwitchTotal)/1000000000.0, 0, 'f', 3);
  */
  qint64 outer=0, inner=0, invokepost=0, invokeswitch=0, threaded=0;
  for (int t=0; t<precisionRecord->threadCount; t++)
  {
    outer+=precisionRecord->threads[t]->timeOuterTotal_;
    inner+=precisionRecord->threads[t]->timeInnerTotal_;
    invokepost+=precisionRecord->threads[t]->timeInvokePostTotal_;
    invokeswitch+=precisionRecord->threads[t]->timeInvokeSwitchTotal_;
    threaded+=precisionRecord->threads[t]->timeThreadedTotal;
  }
  result=QString("%1-%2,%3-%4 =%5,").
      arg((inner)/1000000000.0, 0, 'f', 3).
      arg((invokeswitch)/1000000000.0, 0, 'f', 3).
      arg((invokepost)/1000000000.0, 0, 'f', 3).
      arg((outer)/1000000000.0, 0, 'f', 3).
      arg((threaded)/1000000000.0, 0, 'f', 3);
  return result;
}

QString MandelModel::getTextXY()
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
    return "-";
  return precisionRecord->currentWorker->toString(precisionRecord->orbit.evaluator.currentParams.c.re)+" +i* "+
         precisionRecord->currentWorker->toString(precisionRecord->orbit.evaluator.currentParams.c.im);
}

static QString mandDoubleToString(double x)
{
  if (x>0)
  {
    if (x>10000)
      return QString("+%1").arg(x, 10, 'g');
    else
      return QString("+%1").arg(x, 10, 'f');
  }
  else
  {
    if (x<-10000)
      return QString("%1").arg(x, 10, 'g');
    else
      return QString("%1").arg(x, 10, 'f');//returns like 50 digits QString::number(x, 10, 'f');
  }
}

QString MandelModel::getTextInfoGen()
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number tmp(precisionRecord->currentWorker.get());
    reimToPixel(&orbit_x, &orbit_y, &precisionRecord->orbit.evaluator.currentParams.c, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  //MandelPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  precisionRecord->wtiPoint.store=&pointStore[orbit_x+imageWidth*orbit_y];
  int wtiIndexFirst, wtiIndexLast;
  precisionRecord->wtiPoint.self_allocator._getRange(wtiIndexFirst, wtiIndexLast);
  precisionRecord->currentWorker->assign_block(wtiIndexFirst, storeWorker, (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN);
  MandelPoint *data=&precisionRecord->wtiPoint;

  QString state;
  switch (data->store->rstate)
  {
    case MandelPointStore::ResultState::stUnknown_:
      if (data->store->wstate.load()==MandelPointStore::WorkState::stIdle)
        state="Unk";
      else
        state="Working...";
      break;
    case MandelPointStore::ResultState::stOutside:
      state="Out"; break;
    case MandelPointStore::ResultState::stOutAngle:
      state="OutA"; break;
    case MandelPointStore::ResultState::stBoundary:
      state="Bound"; break;
    case MandelPointStore::ResultState::stDiverge:
      state="Diver"; break;
    case MandelPointStore::ResultState::stMisiur:
      state="Misiur"; break;
    case MandelPointStore::ResultState::stPeriod2:
      state="Per2"; break;
    case MandelPointStore::ResultState::stPeriod3:
      state="Per3"; break;
    case MandelPointStore::ResultState::stMaxIter:
      state="Max"; break;
  }

  return state+" iter="+QString::number(data->store->iter)+" near="+QString::number(data->store->near0iter)+
      " fc="+mandDoubleToString(storeWorker->toDouble(data->fc_c.re))+mandDoubleToString(storeWorker->toDouble(data->fc_c.im))+"i";
}

QString MandelModel::getTextInfoSpec()
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number tmp(precisionRecord->currentWorker.get());
    reimToPixel(&orbit_x, &orbit_y, &precisionRecord->orbit.evaluator.currentParams.c, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  //MandelPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  precisionRecord->wtiPoint.store=&pointStore[orbit_x+imageWidth*orbit_y];
  int wtiIndexFirst, wtiIndexLast;
  precisionRecord->wtiPoint.self_allocator._getRange(wtiIndexFirst, wtiIndexLast);
  precisionRecord->currentWorker->assign_block(wtiIndexFirst, storeWorker, (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN);
  MandelPoint *data=&precisionRecord->wtiPoint;

  switch (data->store->rstate)
  {
    case MandelPointStore::ResultState::stUnknown_:
      if (data->store->wstate.load()==MandelPointStore::WorkState::stIdle)
        return " ";
      else
        return "Working...";
      break;
    case MandelPointStore::ResultState::stOutside:
      return QString("ext=")+QString::number(data->store->exterior_hits); break;
    case MandelPointStore::ResultState::stOutAngle:
      return QString("ext=")+QString::number(data->store->exterior_hits); break;
    case MandelPointStore::ResultState::stBoundary:
      return " ";
    case MandelPointStore::ResultState::stDiverge:
      return " ";
    case MandelPointStore::ResultState::stMisiur:
      return " ";
    case MandelPointStore::ResultState::stPeriod2:
      return QString("per=")+QString::number(data->store->period)+" int="+QString::number(data->store->interior)   +
          " mult="+QString::number(precisionRecord->orbit.evaluator.bulb.dbg_guessmult);
    case MandelPointStore::ResultState::stPeriod3:
      return QString("per=")+QString::number(data->store->period)+" int="+QString::number(data->store->interior)   +
          " mult="+QString::number(precisionRecord->orbit.evaluator.bulb.dbg_guessmult);
    case MandelPointStore::ResultState::stMaxIter:
      return " ";
  }
  return "-?-?-";
}

ShareableViewInfo MandelModel::getViewInfo()
{
  ShareableViewInfo result(&precisionRecord->shareableViewInfoAllocator);
  //result.worker=orbit.worker;
  result.period=precisionRecord->orbit.evaluator.currentData.store->near0iter;//evaluator.currentData.lookper_lastGuess;//orbit.pointData.period;
  if (result.period<1)
    result.period=1;
  result.scale=precisionRecord->position.step_size__;
  /*orbit.worker->init_(&result.re_, &result.re_p);
  orbit.worker->init_(&result.im, &result.im_p);
  orbit.worker->init_(&result.root_re, &result.rre_p);
  orbit.worker->init_(&result.root_im, &result.rim_p);*/
  result.c.assign(&precisionRecord->orbit.evaluator.currentParams.c);
  result.root.assign(&precisionRecord->orbit.evaluator.currentData.root);

  //TODO: why here? should be somewhere else
  precisionRecord->lagu_c.assign(&precisionRecord->orbit.evaluator.currentParams.c);
  precisionRecord->lagu_r.assign(&precisionRecord->orbit.evaluator.currentData.root);

  return result;
}

void MandelModel::transformStore(MandelMath::worker_multi *old_worker, MandelMath::worker_multi *old_sworker, MandelPointStore *old_store, int old_width, int old_height, const MandelMath::complex *old_c,
                                 MandelMath::worker_multi *new_worker, MandelMath::worker_multi *new_sworker, MandelPointStore *new_store, int new_width, int new_height, const MandelMath::complex *new_c,
                                 int inlog, int new_step_log)
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
  {
    dbgPoint();
    return;
  };
  int indexOfWtiPoint, wtiIndexLast;
  precisionRecord->wtiPoint.self_allocator._getRange(indexOfWtiPoint, wtiIndexLast);
  assert(wtiIndexLast-indexOfWtiPoint==MandelPoint::LEN);
  //(oldx-width/2)*old_step+old_cre = (newx-width/2)*new_step+new_cre
  //oldx = (width/2) + (newx-width/2+(new_cre-old_cre)/new_step)*new_step/old_step
  int step_scale_n_shift, step_scale_d_shift, step_scale_d_mask;
  if (inlog>=0)
  {
    step_scale_n_shift=0;
    step_scale_d_shift=inlog;
    step_scale_d_mask=(1<<inlog)-1;
  }
  else
  {
    step_scale_n_shift=-inlog;
    step_scale_d_shift=0;
    step_scale_d_mask=0;
  }
  int delta_y_int, delta_x_int;
  {
    //double delta_y=-new_height/2+(new_cim-old_cim)/new_step;
    //delta_y*=(1<<step_scale_n_shift);
    //if (fabs(delta_y-qRound(delta_y))>0.0001)
    //  dbgPoint();
    //int delta_y_int=qRound(delta_y);
    MandelMath::number tmp(old_worker);
    old_worker->assign(tmp.ptr, old_c->im);
    old_worker->sub(tmp.ptr, new_c->im); //and reversing y at the last minute
    old_worker->lshift(tmp.ptr, new_step_log+step_scale_n_shift);
    delta_y_int=old_worker->toRound(tmp.ptr);
    delta_y_int-=(new_height/2)<<step_scale_n_shift;

    //double delta_x=-new_width/2+(new_cre-old_cre)/new_step;
    //delta_x*=(1<<step_scale_n_shift);
    //if (fabs(delta_x-qRound(delta_x))>0.0001)
    //  dbgPoint();
    //int delta_x_int=qRound(delta_x);
    old_worker->assign(tmp.ptr, new_c->re);
    old_worker->sub(tmp.ptr, old_c->re);
    old_worker->lshift(tmp.ptr, new_step_log+step_scale_n_shift);
    delta_x_int=old_worker->toRound(tmp.ptr);
    delta_x_int-=(new_width/2)<<step_scale_n_shift;
  }

  MandelMath::complex c(new_worker);
  for (int newy=0; newy<new_height; newy++)
  {
    int oldy;
    if ((step_scale_d_mask==0) || ((((newy<<step_scale_n_shift)+delta_y_int)&step_scale_d_mask)==0))
      oldy=(old_height/2) + (((newy<<step_scale_n_shift)+delta_y_int)>>step_scale_d_shift);
    else
      oldy=-1;//call reset() in the second loop, we may still need the points  =imageHeight;
    new_worker->assign(c.im, precisionRecord->position.center.im);
    new_worker->add_double(c.im, (new_height/2-newy)*precisionRecord->position.step_size__);
    for (int newx=0; newx<new_width; newx++)
    {
      int oldx;
      if ((step_scale_d_mask==0) || ((((newx<<step_scale_n_shift)+delta_x_int)&step_scale_d_mask)==0))
        oldx=(old_width/2) + (((newx<<step_scale_n_shift)+delta_x_int)>>step_scale_d_shift);
      else
        oldx=-1;//call reset() in the second loop, we may still need the points  =imageWidth;
      if ((oldy>newy) || ((oldy==newy) && (oldx>newx)))
      {
        if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
        { //copy old point to new place
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          if (new_store[newy*new_width+newx].wstate==MandelPointStore::WorkState::stWorking)
            new_store[newy*new_width+newx].wstate=MandelPointStore::WorkState::stIdle; //work will be cancelled because of new epoch
          new_sworker->assign_block((newy*new_width+newx)*MandelPoint::LEN, old_sworker, (oldy*old_width+oldx)*MandelPoint::LEN, MandelPoint::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          new_worker->assign(c.re, precisionRecord->position.center.re);
          new_worker->add_double(c.re, (newx - new_width/2)*precisionRecord->position.step_size__);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(&c);
          new_sworker->assign_block((newy*new_width+newx)*MandelPoint::LEN, precisionRecord->currentWorker.get(), indexOfWtiPoint, MandelPoint::LEN);
        }
      };
    }
  }
  for (int newy=(new_height-1)&0xfffffff; newy>=0; newy--) //avoid Clang warning about newy possibly ~ 2^31
  {
    int oldy;
    if ((step_scale_d_mask==0) || ((((newy<<step_scale_n_shift)+delta_y_int)&step_scale_d_mask)==0))
      oldy=(old_height/2) + (((newy<<step_scale_n_shift)+delta_y_int)>>step_scale_d_shift);
    else
      oldy=-1;
    new_worker->assign(c.im, precisionRecord->position.center.im);
    new_worker->add_double(c.im, (new_height/2-newy)*precisionRecord->position.step_size__);
    for (int newx=new_width-1; newx>=0; newx--)
    {
      int oldx;
      if ((step_scale_d_mask==0) || ((((newx<<step_scale_n_shift)+delta_x_int)&step_scale_d_mask)==0))
        oldx=(old_width/2) + (((newx<<step_scale_n_shift)+delta_x_int)>>step_scale_d_shift);
      else
        oldx=-1;
      if ((oldy<newy) || ((oldy==newy) && (oldx<=newx)))
      {
        if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
        { //copy old to new
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          if (new_store[newy*new_width+newx].wstate==MandelPointStore::WorkState::stWorking)
            new_store[newy*new_width+newx].wstate=MandelPointStore::WorkState::stIdle; //work will be cancelled because of new epoch
          new_sworker->assign_block((newy*new_width+newx)*MandelPoint::LEN, old_sworker, (oldy*old_width+oldx)*MandelPoint::LEN, MandelPoint::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          new_worker->assign(c.re, precisionRecord->position.center.re);
          new_worker->add_double(c.re, (newx - new_width/2)*precisionRecord->position.step_size__);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(&c);
          new_sworker->assign_block((newy*new_width+newx)*MandelPoint::LEN, precisionRecord->currentWorker.get(), indexOfWtiPoint, MandelPoint::LEN);
        }
      }
    }
  }
}

void MandelModel::setView_double(double c_re, double c_im, double scale)
{
  MandelMath::complex c(precisionRecord->currentWorker.get());
  c.zero(c_re, c_im);
  setView(&c, scale);
}

void MandelModel::setView(const MandelMath::complex *c, double scale)
{
  MandelMath::complex old_c(precisionRecord->currentWorker.get());
  old_c.assign(&precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    precisionRecord->position.setView(c, scale);

    transformStore(precisionRecord->currentWorker.get(), storeWorker, pointStore, imageWidth, imageHeight, &old_c,
                   precisionRecord->currentWorker.get(), storeWorker, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

void MandelModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::complex old_c(precisionRecord->currentWorker.get());
  old_c.assign(&precisionRecord->position.center);

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    int dx=qRound(delta_x);
    int dy=qRound(delta_y);
    precisionRecord->position.move(dx, dy);
    //qDebug()<<"new c: re="<<position.worker->toString(&position.center_re_s)<<",im="<<position.worker->toString(&position.center_im_s);

    transformStore(precisionRecord->currentWorker.get(), storeWorker, pointStore, imageWidth, imageHeight, &old_c,
                   precisionRecord->currentWorker.get(), storeWorker, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   0, precisionRecord->position.step_log);

    startNewEpoch();
#if 0
  for (int newy=0; newy<imageHeight; newy++)
  {
    for (int newx=0; newx<imageWidth; newx++)
    {
      if (pointStore_[newy*imageWidth+newx].state==MandelPointStore::State::stWorking)
        dbgPoint();
    }
  }
#endif
  }
}

void MandelModel::zoom(double x, double y, int inlog)
{
  MandelMath::complex old_c(precisionRecord->currentWorker.get());
  old_c.assign(&precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    precisionRecord->position.scale(inlog, qRound(x)-imageWidth/2, imageHeight/2-qRound(y));

    transformStore(precisionRecord->currentWorker.get(), storeWorker, pointStore, imageWidth, imageHeight, &old_c,
                   precisionRecord->currentWorker.get(), storeWorker, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

void MandelModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    int newLength=width*height;
    MandelPointStore *newStore=new MandelPointStore[newLength];
    {
      QString size_as_text=QString::number(newLength*sizeof(MandelPointStore));
      for (int pos=size_as_text.length()-3; pos>0; pos-=3)
        size_as_text.insert(pos, '\'');
      qDebug()<<"pointStore uses"<<size_as_text.toLocal8Bit().constData()<<"B"; //lots of work to skip those quotes... can't skip spaces at all
    }
    MandelMath::worker_multi *newStoreWorker;
    newStoreWorker=MandelMath::worker_multi::create(storeWorker->ntype(), newLength*MandelPoint::LEN);
    MandelMath::worker_multi::Allocator *newStoreAllocator;
    newStoreAllocator=new MandelMath::worker_multi::Allocator(newStoreWorker->getAllocator(), newLength*MandelPoint::LEN);
    MandelMath::complex old_c(precisionRecord->currentWorker.get());
    old_c.assign(&precisionRecord->position.center);

    transformStore(precisionRecord->currentWorker.get(), storeWorker, pointStore, imageWidth, imageHeight, &old_c,
                   precisionRecord->currentWorker.get(), newStoreWorker, newStore, width, height, &precisionRecord->position.center,
                   0, precisionRecord->position.step_log);

    delete storeAllocator;
    delete storeWorker;
    storeWorker=newStoreWorker;
    storeAllocator=newStoreAllocator;
    delete[] pointStore;
    pointStore=newStore;
    imageWidth=width;
    imageHeight=height;

    startNewEpoch();
  }
}

void MandelModel::startNewEpoch()
{
  epoch=(epoch%2000000000)+1;
  nextGivenPointIndex=0;
  effortBonus=0;
#if 0
  for (int t=0; t<precisionRecord->threadCount; t++)
    if (precisionRecord->threads[t]->currentParams.pixelIndex<0)
      giveWork(precisionRecord->threads[t]);
#else
  for (int t=0; t<precisionRecord->threadCount; t++)
  {
    //giveWorkToThread(precisionRecord->threads[t]);
    precisionRecord->threads[t]->workIfEpoch=epoch;
  }
  _threadsWorking+=precisionRecord->threadCount;
  emit triggerComputeThreaded(epoch); //::invokeMethod cannot pass parameters, but ::connect can
#endif
}

void MandelModel::reimToPixel(int *circ_x, int *circ_y, const MandelMath::complex *c, MandelMath::number *tmp)
{
  double scaled;
  tmp->assign(c->re);
  tmp->sub(precisionRecord->position.center.re);
  scaled=tmp->toDouble()/precisionRecord->position.step_size__;
  if (scaled<-10003 || scaled>10003)
    *circ_x=-100;
  else
    *circ_x=qRound(scaled)+imageWidth/2;

  tmp->assign(c->im);
  tmp->sub(precisionRecord->position.center.im);
  scaled=tmp->toDouble()/precisionRecord->position.step_size__;
  if (scaled<-10003 || scaled>10003)
    *circ_y=-100;
  else
    *circ_y=imageHeight/2-qRound(scaled);
}

void MandelModel::paintOrbit(ShareableImageWrapper image, int x, int y)
{
  if ((x<0) || (x>=imageWidth) || (y<0) || (y>=imageHeight))
    return;
  QImage newOverlay(imageWidth, imageHeight, QImage::Format::Format_ARGB32);
  QPainter painter(&newOverlay);
  //QRgb what=newOverlay.pixel(0, 0);
  //if (what!=0) //0xcdcdcdcd in MSVC compiler
  {
    painter.setCompositionMode(QPainter::CompositionMode::CompositionMode_Source);
    painter.fillRect(0, 0, imageWidth, imageHeight, Qt::GlobalColor::transparent);
  };

/* testing how drawLine rounds in Qt
   result: extremely ugly and no way to fix it
  {
    painter.setPen(QColor(0xff, 0xff, 0xff));
    painter.fillRect(0, 0, 20+6*10, 20, Qt::GlobalColor::black);
    for (int i=0; i<6; i++)
    {
      newOverlay.setPixel(10+10*i, 10-5, 0x00ffffff);
      newOverlay.setPixel(10+10*i-1, 10+5, 0x00ffffff);
      newOverlay.setPixel(10+10*i+1, 10+5, 0x00ffffff);

      switch (i)
      {
        case 0: painter.drawLine(10+10*i-2, 10+3, 10+10*i, 10-3); break;
        case 1: painter.drawLine(10+10*i, 10-3, 10+10*i-2, 10+3); break;
        case 2: painter.drawLine(10+10*i+2, 10+3, 10+10*i, 10-3); break;
        case 3: painter.drawLine(10+10*i, 10-3, 10+10*i+2, 10+3); break;
        case 4:
        {
          QLine l3[3]={{10+10*i+1, 10-1, 10+10*i-1, 10-1},
                       {10+10*i-1, 10-1, 10+10*i-1, 10+1},
                       {10+10*i-1, 10+1, 10+10*i+1, 10+1}};
          painter.drawLines(l3, 3);
        } break;
        case 5:
        {
          QLine l2[2]={{10+10*i-2, 10, 10+10*i+2, 10},
                       {10+10*i, 10-2, 10+10*i, 10+2}};
          painter.drawLines(l2, 2);
        } break;
      }
    }
  }*/

  if (precisionRecord==nullptr)
    dbgPoint();
  else if (precisionRecord->currentWorker==nullptr)
    dbgPoint();
  else if (precisionRecord->orbit.currentWorker->ntype()!=precisionRecord->currentWorker->ntype())
    dbgPoint();

  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  //MandelPoint data_(&pointStore_[y*imageWidth+x], &allo);
  //MandelPoint *data=&orbit_->evaluator.currentData;
  MandelPointStore *resultStore=&pointStore[y*imageWidth+x];
  switch (resultStore->rstate)
  {
    case MandelPointStore::ResultState::stOutside:
    case MandelPointStore::ResultState::stOutAngle:
    {
      int exterior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0xff, 0xff, 0));
      exterior=qRound(resultStore->exterior_hits/precisionRecord->position.step_size__);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
      painter.setPen(QColor(0xc0, 0xc0, 0));
      exterior=qRound(resultStore->exterior_avoids/precisionRecord->position.step_size__);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
    } break;
    case MandelPointStore::ResultState::stPeriod2:
    case MandelPointStore::ResultState::stPeriod3:
    {
      int interior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0, 0xff, 0xff));
      interior=qRound(resultStore->interior/precisionRecord->position.step_size__);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
      painter.setPen(QColor(0, 0xc0, 0xc0));
      interior=qRound(resultStore->interior/4/precisionRecord->position.step_size__);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
    } break;
    default: ;
  }
  precisionRecord->position.pixelXtoRE(x-imageWidth/2, precisionRecord->orbit.evaluator.currentParams.c.re);
  precisionRecord->position.pixelYtoIM(imageHeight/2-y, precisionRecord->orbit.evaluator.currentParams.c.im);
  precisionRecord->orbit.evaluator.currentParams.epoch=epoch;
  precisionRecord->orbit.evaluator.workIfEpoch=precisionRecord->orbit.evaluator.busyEpoch;//epoch;
  precisionRecord->orbit.evaluator.currentParams.pixelIndex=0;
  //precisionRecord->orbit.evaluator.currentData.store->rstate=MandelPointStore::ResultState::stUnknown_;
  //precisionRecord->orbit.evaluator.currentData.store->wstate=MandelPointStore::WorkState::stIdle;
  precisionRecord->orbit.evaluator.currentData.zero(&precisionRecord->orbit.evaluator.currentParams.c);
  precisionRecord->orbit.evaluator.currentData.store->wstate=MandelPointStore::WorkState::stWorking;
  precisionRecord->orbit.evaluator.currentParams.breakOnNewNearest=true;
  precisionRecord->orbit.evaluator.currentParams.maxiter_=1<<MAX_EFFORT;
  {
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0xff, 0xff)); //paint path
    int line_sx, line_sy;
    reimToPixel(&line_sx, &line_sy, &precisionRecord->orbit.evaluator.currentData.f, &precisionRecord->tmp);
    while ((precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stUnknown_) &&
           (precisionRecord->orbit.evaluator.currentData.store->iter<(1<<MAX_EFFORT)))
    {
      int line_ex, line_ey;

      if ((resultStore->rstate==MandelPointStore::ResultState::stPeriod2 || resultStore->rstate==MandelPointStore::ResultState::stPeriod3) &&
          precisionRecord->orbit.evaluator.currentData.store->iter<resultStore->period)
      { //paint first period fully
        precisionRecord->orbit.evaluator.currentParams.maxiter_=precisionRecord->orbit.evaluator.currentData.store->iter+1;
      }
      else if (precisionRecord->orbit.evaluator.currentData.store->lookper_lastGuess==0)
        precisionRecord->orbit.evaluator.currentParams.maxiter_=1<<MAX_EFFORT; //dont't know->run fully
      else //stop at multiples of lookper, +1
        precisionRecord->orbit.evaluator.currentParams.maxiter_=1+(precisionRecord->orbit.evaluator.currentData.store->iter/precisionRecord->orbit.evaluator.currentData.store->lookper_lastGuess+1)*precisionRecord->orbit.evaluator.currentData.store->lookper_lastGuess;
      precisionRecord->orbit.evaluator.startCompute(+1);

      reimToPixel(&line_ex, &line_ey, &precisionRecord->orbit.evaluator.currentData.f, &precisionRecord->tmp);
      if (line_ex>=-3 && line_ex<=10003 && line_ey>=-3 && line_ey<=10003)
      {
        if (line_sx>=-3 && line_sx<=10003 && line_sy>=-3 && line_sy<=10003)
          painter.drawLine(line_sx, line_sy, line_ex, line_ey);
        line_sx=line_ex;
        line_sy=line_ey;
      };
    }
  }
  if ((precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod2) ||
      (precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod3))
  {
    int circ_x, circ_y;
    painter.setPen(QColor(0, 0xff, 0xff)); //paint root
    reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.currentData.root, &precisionRecord->tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };

    //orbit.bulbinfo
    /*MandelMath::complex c(orbit.worker, &orbit.evaluator.currentParams.c_re, &orbit.evaluator.currentParams.c_im, true);
    MandelMath::complex cb(orbit.worker, &orbit.bulb.cb_re, &orbit.bulb.cb_im, true);
    MandelMath::complex rb(orbit.worker, &orbit.bulb.rb_re_, &orbit.bulb.rb_im, true);
    MandelMath::complex xc(orbit.worker, &orbit.bulb.xc_re, &orbit.bulb.xc_im, true);
    MandelMath::complex baseZC(orbit.worker, &orbit.bulb.baseZC_re, &orbit.bulb.baseZC_im, true);
    MandelMath::complex baseCC(orbit.worker, &orbit.bulb.baseCC_re, &orbit.bulb.baseCC_im, true);*/
    precisionRecord->orbit.bulb.foundMult=0;
    precisionRecord->orbit.bulb.is_card=false;

    precisionRecord->orbit.bulb.valid=precisionRecord->orbit.evaluator.bulb.findBulbBase(precisionRecord->orbit.evaluator.currentData.store->period,
        &precisionRecord->orbit.evaluator.currentParams.c, &precisionRecord->orbit.bulb.cb, &precisionRecord->orbit.bulb.rb,
        &precisionRecord->orbit.bulb.xc, &precisionRecord->orbit.bulb.baseZC, &precisionRecord->orbit.bulb.baseCC,
        &precisionRecord->orbit.bulb.is_card, &precisionRecord->orbit.bulb.foundMult);
    if (precisionRecord->orbit.bulb.valid)
    {
      painter.setBrush(QBrush(QColor(0, 0xff, 0xff)));
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.bulb.dbg_first_cb, &precisionRecord->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0x80, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.bulb.dbg_first_rb, &precisionRecord->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0x80, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.bulb.xc, &precisionRecord->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        if (precisionRecord->orbit.bulb.is_card)
          painter.setPen(QColor(0, 0xff, 0xff)); //card center
        else
          painter.setPen(QColor(0x80, 0xc0, 0xc0)); //bulb center
        //painter.setBrush(Qt::BrushStyle::SolidPattern);
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.bulb.cb, &precisionRecord->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0xff, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.bulb.rb, &precisionRecord->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0xff, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
    };
  }

  if (!precisionRecord->lagu_c.is0())
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0, 0xff, 0)); //paint c
    reimToPixel(&circ_x, &circ_y, &precisionRecord->lagu_c, &precisionRecord->tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l3[3]={{circ_x+1, circ_y-1, circ_x-1, circ_y-1},
                   {circ_x-1, circ_y-1, circ_x-1, circ_y+1},
                   {circ_x-1, circ_y+1, circ_x+1, circ_y+1}};
      painter.drawLines(l3, 3);
    };

    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0, 0xff, 0)); //paint root as /\    .
    reimToPixel(&circ_x, &circ_y, &precisionRecord->lagu_r, &precisionRecord->tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };

  }
  //donePixel1(&orbit.evaluator);
  //orbit.pointData.cleanup(position.worker);
  image.image->swap(newOverlay);
}

int MandelModel::periodToIndex(int period)
{
  //powers of factors -> color index
  //[1]..1
  //[1,1] .. 2
  //[2] .. 2
  //[1,1,1] .. 2
  //[2,1] .. 3
  //[3] .. 3
  //[1,1,1,1] .. 2
  //[2,1,1] .. 4
  //[2,2] .. 4
  //[3,1] .. 4
  //[4] .. 4
  //[1,1,1,1,1] .. 2
  //[2,1,1,1] .. 5
  //[2,2,1] .. 5
  //[3,1,1] .. 5
  //[3,2] .. 5
  //[4,1] .. 5
  //[5] .. 5
  if (period<=1)
    return 0;
  if (periodToIndexCache.length()<=period)
    periodToIndexCache.resize(period+1);
  if (periodToIndexCache[period]!=0)
    return periodToIndexCache[period];
  int c=1;
  for (int i=1; i<=period/2; i++)
  {
    if (period % i ==0)
    {
      int c2=periodToIndex(i);
      if (c2>=c)
        c=c2+1;
    };
  }
  periodToIndexCache[period]=c;
  return c;
}

int MandelModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return -1;
  timerWriteToImage.start();
  int indexOfWtiPoint, _discard;
  {
    qint64 totalNewtons=0;
    for (int t=0; t<precisionRecord->threadCount; t++)
      totalNewtons+=precisionRecord->threads[t]->totalNewtonIterations;
    indexOfWtiPoint=totalNewtons; //zoom0: 47167, zoom1: 188490, zoom2: 754210  with breaker
    (void)totalNewtons;           //       47238         188929         756264  without breaker
  }
  precisionRecord->wtiPoint.self_allocator._getRange(indexOfWtiPoint, _discard);
  MandelPointStore *wtiStore;
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      bool knownenum=false;
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
      //MandelPoint data_(&pointStore_[y*imageWidth+x], &allo);
      wtiStore=&pointStore[y*imageWidth+x];
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleKind:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown_:
              if (wtiStore->wstate==MandelPointStore::WorkState::stIdle)
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xffff00ff);
                //image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int b=0x9f+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b<<0));
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff8000ff);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffc000);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              /*int r;
              switch (data->period)
              {
                case 1: r=0xc0; break;
                case 2: r=0xff; break;
                case 3: r=0x80; break;
                default: r=0xe0;
              }*/
              if (wtiStore->period>wtiStore->near0iter)
                image.image->setPixel(x, y, 0xffff00ff); //seems to only happen by mistake, not in reality
              else
              {
                int index=periodToIndex(wtiStore->period);
                //reverse bottom 7 bits:
                int rh=0x73516240>>((index&7)<<2); //reverse bits 0..2
                int rl=0x73516240>>((index&0x70)>>2); //reverse bits 4..6
                rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
                image.image->setPixel(x, y, 0xff000000+(rh<<16));
              }
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
              knownenum=true;
            } break;
          }
          if (!knownenum)
            image.image->setPixel(x, y, 0xffffffff);
        } break;
        case paintStyle::paintStyleCls:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown_:
              image.image->setPixel(x, y, 0x00906090);
              knownenum=true;
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              #if 0 //smooth color madness
              int r=128+floor(127*cos(data->iter/10.0*2*3.1415926535));
              int g=128+floor(127*cos((data->iter/10.0+0.333)*2*3.1415926535));
              int b=128+floor(127*cos((data->iter/10.0+0.667)*2*3.1415926535));
              if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                result.setPixel(x, y, 0xffffffff);
              else
                result.setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));*/
              #endif

              #if 0 //by data->iter only
              int iter=data->iter;
              switch (iter % 12)
              {
                case  0: image.image->setPixel(x, y, 0xff0000ff); break;
                case  1: image.image->setPixel(x, y, 0xff0080ff); break;
                case  2: image.image->setPixel(x, y, 0xff00ffff); break;
                case  3: image.image->setPixel(x, y, 0xff00ff80); break;
                case  4: image.image->setPixel(x, y, 0xff00ff00); break;
                case  5: image.image->setPixel(x, y, 0xff80ff00); break;
                case  6: image.image->setPixel(x, y, 0xffffff00); break;
                case  7: image.image->setPixel(x, y, 0xffff8000); break;
                case  8: image.image->setPixel(x, y, 0xffff0000); break;
                case  9: image.image->setPixel(x, y, 0xffff0080); break;
                case 10: image.image->setPixel(x, y, 0xffff00ff); break;
                case 11: image.image->setPixel(x, y, 0xff8000ff); break;
                default: image.image->setPixel(x, y, 0xffffffff);
              }
              #endif

              #if 1 //smooth by iter
              precisionRecord->currentWorker->assign_block(indexOfWtiPoint, storeWorker, (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN);
              double re=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.f.re);
              double im=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.f.im);
              double iter=wtiStore->iter+6-log2(log2(re*re+im*im)); //+6 to match integer coloring
              iter=iter/12;
              iter=(iter-floor(iter))*6;
              int iter_phase=iter;
              int iter_256=(iter-iter_phase)*256;
              int r=0xff, g=0xff, b=0xff;
              switch (iter_phase)
              {
                case  0: r=0x00; g=iter_256; b=0xff; break;
                case  1: r=0x00; g=0xff; b=0xff-iter_256; break;
                case  2: r=iter_256; g=0xff; b=0x00; break;
                case  3: r=0xff; g=0xff-iter_256; b=0x00; break;
                case  4: r=0xff; g=0x00; b=iter_256; break;
                case  5: r=0xff-iter_256; g=0x00; b=0xff; break;
              }
              image.image->setPixel(x, y, 0xff000000|(r<<16)|(g<<8)|b);
              #endif

              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              int index=periodToIndex(wtiStore->period);
              //reverse bottom 7 bits:
              int rh=0x73516240>>((index&7)<<2);
              int rl=0x73516240>>((index&0x70)>>2);
              rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
              /*switch (data->period)
              {
                case 1: r=0xc0; break;
                case 2: r=0xff; break;
                case 3: r=0x80; break;
                default: r=0xe0;
              }*/
              //image.image->setPixel(x, y, 0xffffc0c0);
              image.image->setPixel(x, y, 0xff000000+rh*0x010101);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              //image.image->setPixel(x, y, 0xff808080);
              image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
            } break;
          }
          if (!knownenum)
            image.image->setPixel(x, y, 0xffffffff);
        } break;
        case paintStyle::paintStyleExter:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown_:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              double tf;
              if ((wtiStore->exterior_avoids>10000) || (wtiStore->exterior_avoids<=0))
                tf=0;
              else if (wtiStore->exterior_avoids>=1)
                tf=(1-wtiStore->exterior_avoids)*1;
              else
                tf=sqrt(1-log(wtiStore->exterior_avoids))*2-2;
              int r=0x9f+qRound(0x60*sin(tf*2.828)); //red middle
              int g=0x9f+qRound(0x60*sin(tf*6.928)); //green fastest
              int b=0x9f+qRound(0x60*sin(tf)); //blue slowest
              if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              precisionRecord->currentWorker->assign_block(indexOfWtiPoint, storeWorker, (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN);
              double re=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.fz_r.re);
              double im=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.fz_r.im);
              //double angle=std::atan2(im, re);
              double mag=sqrt(re*re+im*im);
              int r, g, b;
              if (mag<=0)
               { r=0xff; g=0x00; b=0x00; }
              else
              {
                //re/=mag;
                //im/=mag;
                if (mag>=1)
                  mag=0.99;
                //mag=-log(1-mag)/log(2);
                int expo;
                mag=frexp(1-mag, &expo);
                r=0x40+(int)floor(0xc0*mag);
                g=r;//0xc0+qRound(0x3f*re);
                b=r;//0xc0+qRound(0x3f*im);
              }
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+b);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
              knownenum=true;
            } break;
          }
          if (!knownenum)
            image.image->setPixel(x, y, 0xffffffff);
        } break;
        case paintStyle::paintStyleInter:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown_:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              //int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              int b=qRound(log(wtiStore->iter)*100)%256;
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              int ti=30;
              if ((wtiStore->interior>1) || (wtiStore->interior<=0))
                ti=0;
              else
                ti=(qRound(-log(wtiStore->interior/4)*300)+12*0xc0) % (6*0xc0);
              int r, g, b;
              if (ti<0xC0)
              { r=0x3f+ti; g=0xff; b=0x3f; }                           // + H L
              else if (ti<2*0xC0)                                      //
              { r=0xff; g=0xff-(ti-0xc0); b=0x3f; }                    // H - L
              else if (ti<3*0xC0)                                      //
              { r=0xff; g=0x3f; b=0x3f+(ti-2*0xc0); }                  // H L +
              else if (ti<4*0xC0)                                      //
              { r=0xff-(ti-3*0xc0); g=0x3f; b=0xff; }                  // - L H
              else if (ti<5*0xC0)                                      //
              { r=0x3f; g=0x3f+(ti-4*0xc0); b=0xff; }                  // L + H
              else                                                     //
              { r=0x3f; g=0xff; b=0xff-(ti-5*0xc0); }                  // L H -
              if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
              knownenum=true;
            } break;
          }
          if (!knownenum)
            image.image->setPixel(x, y, 0xffffffff);
        } break;
        case paintStyle::paintStyleNear:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown_:
              image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int ti=wtiStore->near0iter;
              if (ti>=0)
              {
                int tj=0;
                while ((ti%2)==0) { tj+=128*2/5; ti/=2; }
                while ((ti%3)==0) { tj+=128*2/3; ti/=3; }
                while ((ti%5)==0) { tj+=30; ti/=5; }
                while ((ti%7)==0) { tj+=40; ti/=7; }
                while ((ti%11)==0) { tj+=50; ti/=11; }
                while ((ti%13)==0) { tj+=60; ti/=13; }
                while ((ti%17)==0) { tj+=70; ti/=17; }
                int b=0x80+(tj%0x80);
                image.image->setPixel(x, y, 0xff000000+(b<<0));
              }
              else
                image.image->setPixel(x, y, 0xff000080);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              /*int ti=data->newton_iter; //very noisy, maybe show <=10, >10, >30, 49
              int r;
              if (ti<=10)
                r=0x60;
              else if (ti<30)
                r=0x90;
              else
                r=0xff;
              image.image->setPixel(x, y, 0xff000000+(r<<16));*/

              int ti=wtiStore->near0iter;
              if (ti>0)
              {
                int tj=0;
                while ((ti%2)==0) { tj+=128*2/5; ti/=2; }
                while ((ti%3)==0) { tj+=128*2/3; ti/=3; }
                while ((ti%5)==0) { tj+=30; ti/=5; }
                while ((ti%7)==0) { tj+=40; ti/=7; }
                while ((ti%11)==0) { tj+=50; ti/=11; }
                while ((ti%13)==0) { tj+=60; ti/=13; }
                while ((ti%17)==0) { tj+=70; ti/=17; }
                int r=0x80+(tj%0x80);
                image.image->setPixel(x, y, 0xff000000+(r<<16));
              }
              else
                image.image->setPixel(x, y, 0xff800000);
              /* we need func(2)!=func(3) here
              int index=periodToIndex(data->near0iter);
              //reverse bottom 7 bits:
              int rh=0x73516240>>((index&7)<<2);
              int rl=0x73516240>>((index&0x70)>>2);
              rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
              image.image->setPixel(x, y, 0xff000000+(rh<<16));*/
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
              knownenum=true;
            } break;
          }
          if (!knownenum)
            image.image->setPixel(x, y, 0xffffffff);
        } break;
        case paintStyle::paintStyleFZ:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown_:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              precisionRecord->currentWorker->assign_block(indexOfWtiPoint, storeWorker, (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN);
              double re=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.fz_r.re);
              double im=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.fz_r.im);
              double mag=sqrt(MandelMath::sqr_double(re)+MandelMath::sqr_double(im));
              if (mag<0) mag=0;
              else if (mag>1) mag=1;
              /*double mag2=mag*127.49+128;
              double phi=std::atan2(position.worker->toDouble(&data->fz_r_im), position.worker->toDouble(&data->fz_r_re))/(2*M_PI);
              if (phi<0) phi+=1;
              int g=qRound(phi*mag2);
              phi+=0.25;
              if (phi>=1) phi-=1;
              int b=qRound(phi*mag2);*/
              int g=mag==0?128:qRound(re/mag*127.49+128);
              //int b=mag==0?128:qRound(im/mag*127.49+128);
              int b=mag==0?128:qRound(im/mag*127.49+256)%256;
              int r=0;
              mag=-log(1-mag)/log(2);
              mag=mag-floor(mag);
              if (mag<0.2)
                r=128;
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
              knownenum=true;
            } break;
          }
          if (!knownenum)
            image.image->setPixel(x, y, 0xffffffff);
        } break;
        case paintStyle::paintStyleFC:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown_:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              if (!wtiStore->has_fc_r)
                image.image->setPixel(x, y, 0xffc0c0c0);
              else
              {
                precisionRecord->currentWorker->assign_block(indexOfWtiPoint, storeWorker, (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN);
                double re=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.fc_c.re);
                double im=precisionRecord->position.worker->toDouble(precisionRecord->wtiPoint.fc_c.im);
                double mag=std::hypot(re, im);//sqrt(MandelMath::sqr_double(re)+MandelMath::sqr_double(im));
                if (mag<0) mag=0;
                else if (mag>1) mag=1;
                //int magi=qRound(mag*127.49)+128;
                double mag2=mag*127.49+128;
                double phi=std::atan2(im, re)/(2*M_PI);
                if (phi<0) phi+=1;
                int g=qRound(phi*mag2);
                phi+=0.25;
                if (phi>=1) phi-=1;
                int b=qRound(phi*mag2);
                int r=0;
                mag=log(mag)/log(2);
                mag=mag-floor(mag);
                //if (mag<0.2)
                  //r=128;
                r=255*mag;
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
              knownenum=true;
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
              knownenum=true;
            } break;
          }
          if (!knownenum)
            image.image->setPixel(x, y, 0xffffffff);
        } break;
      }

    }
  //return result;
  qint64 elapsed=timerWriteToImage.elapsed();
  return elapsed;
}

void MandelModel::doneWorkInThread(MandelEvaluator *)
{
  _threadsWorking--;
}

int MandelModel::giveWorkThreaded(MandelEvaluator *me)
{
  QReadLocker locker(&threading_mutex);
  me->timeInvokeSwitchTotal_+=me->timeInvoke_.nsecsElapsed();
  if (epoch!=me->busyEpoch)
    return 3;
  int retryEffortFrom=0;
  int nextEffortBonus=effortBonus; //don't jump to max instantly once pixels<threads
  int intoEvaluator, _discard;
  me->currentData.self_allocator._getRange(intoEvaluator, _discard);
  while (retryEffortFrom>=0)
  {
    retryEffortFrom=-1;
    for (int pi=0; pi<imageWidth*imageHeight; pi++)
    {
      int pointIndex=(nextGivenPointIndex+pi)%(imageWidth*imageHeight);
      //if ((lastGivenPointIndex_!=0) && (pointIndex==0))
        //dbgPoint();
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), pointIndex*MandelPoint::LEN, MandelPoint::LEN, nullptr);
      //MandelPoint pointData_(&pointStore_[pointIndex], &allo);
      MandelPointStore *storeAtIndex=&pointStore[pointIndex];
      MandelPointStore::WorkState state_expected=MandelPointStore::WorkState::stIdle;
      int extra_effort=0;
      if (!storeAtIndex->wstate.compare_exchange_strong(state_expected, MandelPointStore::WorkState::stWorking))
      {
        if ((_selectedPaintStyle==paintStyleFC) &&
            (!storeAtIndex->has_fc_r) &&
            (storeAtIndex->rstate==MandelPointStore::ResultState::stPeriod2 ||
             storeAtIndex->rstate==MandelPointStore::ResultState::stPeriod3))
        {
          state_expected=MandelPointStore::WorkState::stDone;
          if (!storeAtIndex->wstate.compare_exchange_strong(state_expected, MandelPointStore::WorkState::stWorking))
            continue;
            //should check if it still needs but it should be quite rare, evaluate tests it anyway
          extra_effort=1;
        }
        else
          continue;
      }
      {
        if (me->currentParams.pixelIndex!=-1)
          dbgPoint();
        assert(me->currentParams.pixelIndex==-1);
        {
          int phasex=(pointIndex%imageWidth-imageWidth/2+precisionRecord->position.cached_center_re_mod+32768)%32768;
          int phasey=(pointIndex/imageWidth-imageHeight/2+precisionRecord->position.cached_center_im_mod+32768)%32768;
          //int effort=ctz16(phasex)+ctz16(phasey);
          int effort=MandelMath::ctz16(phasex | phasey);
          if (effort>8)
            effort=8;
          effort+=nextEffortBonus;
          if (effort>=MAX_EFFORT)
            effort=MAX_EFFORT;
          if (effort<2)
            effort=2; //quick run ... at least 4 iterations
          me->currentParams.maxiter_=1<<effort;
          if (storeAtIndex->iter >= me->currentParams.maxiter_+extra_effort)
          {
            if (effort>=MAX_EFFORT)
            {
              storeAtIndex->rstate=MandelPointStore::ResultState::stMaxIter;
              storeAtIndex->wstate=MandelPointStore::WorkState::stDone;
            }
            else
            {
              if (retryEffortFrom<0)
                retryEffortFrom=pointIndex;
              storeAtIndex->wstate=MandelPointStore::WorkState::stIdle;
            }
          }
          else
          {
            if (me->currentWorker->ntype()!=precisionRecord->currentWorker->ntype())
              dbgPoint();
            //evaluator->switchType(position.worker);
            precisionRecord->position.pixelXtoRE(pointIndex%imageWidth - imageWidth/2, me->currentParams.c.re);
            precisionRecord->position.pixelYtoIM(imageHeight/2-pointIndex/imageWidth, me->currentParams.c.im);
            me->currentParams.epoch=me->busyEpoch;
            //storeAtIndex->state=MandelPointStore::State::stWorking;
            me->currentParams.pixelIndex=pointIndex;
            me->currentParams.want_fc_r=(_selectedPaintStyle==paintStyleFC);
#if CURRENT_STORE_DIRECT
            me->currentData.store=storeAtIndex;
#else
            me->currentDataStore.assign(storeAtIndex);
#endif
            me->currentWorker->assign_block(intoEvaluator, storeWorker, pointIndex*MandelPoint::LEN, MandelPoint::LEN);
            nextGivenPointIndex=(pointIndex+1)%(imageWidth*imageHeight);
            effortBonus=nextEffortBonus;
            return 0;
          }
        }
      }
        //MandelEvaluator::simple(cr, ci, pointStore[y*imageWidth+x]);
    }
    if ((retryEffortFrom>=0) && (nextEffortBonus<MAX_EFFORT))
    {
      nextEffortBonus++;
      nextGivenPointIndex=retryEffortFrom;
    }
    else
      retryEffortFrom=-1;
  }
  return 2;
}

int MandelModel::doneWorkThreaded(MandelEvaluator *me, bool giveWork)
{
  QReadLocker locker(&threading_mutex);
  if ((me->currentParams.epoch==epoch) && (me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
  {
    {
      MandelPointStore *dstStore=&pointStore[me->currentParams.pixelIndex];
#if CURRENT_STORE_DIRECT
#else
      if (dstStore->wstate!=MandelPointStore::WorkState::stWorking)
        dbgPoint();
#endif
      if (precisionRecord->position.worker==nullptr)
        dbgPoint();
      else
      {
        int first, last;
        me->currentData.self_allocator._getRange(first, last);
        storeWorker->assign_block(me->currentParams.pixelIndex*MandelPoint::LEN, me->currentWorker, first, last-first);
#if CURRENT_STORE_DIRECT
#else
        dstStore->assign(me->currentData.store);
#endif
      }
      if (dstStore->wstate.load()==MandelPointStore::WorkState::stIdle)
        dbgPoint();
      //else if (dstStore->state==MandelPointStore::State::stWorking)
      //  dstStore->state=MandelPointStore::State::stUnknown;
      if (dstStore->wstate==MandelPointStore::WorkState::stWorking)
      {
         if (dstStore->iter>=(1<<MAX_EFFORT))
         {
           dstStore->rstate=MandelPointStore::ResultState::stMaxIter;
           dstStore->wstate=MandelPointStore::WorkState::stDone;
         }
         else if (dstStore->rstate!=MandelPointStore::ResultState::stUnknown_)
           dstStore->wstate=MandelPointStore::WorkState::stDone;
         else
           dstStore->wstate=MandelPointStore::WorkState::stIdle;
      }
      else
        dbgPoint();
    }
  }
  else if (me->currentParams.epoch!=epoch)
  {
    //qDebug()<<"Old pixel finished";
    /*problem fixed if ((me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
    {
      MandelPointStore *dstStore=&pointStore_[me->currentParams.pixelIndex];
      if (dstStore->state==MandelPointStore::State::stWorking)
        dbgPoint();
    }*/
  }
  else
    qWarning()<<"Invalid pixel finished";
  me->currentParams.pixelIndex=-1;
  if (me->currentParams.epoch==epoch)
  {
    if (giveWork)
      return giveWorkThreaded(me);
    else
      return 0;
  }
  return 1;
}

void MandelModel::selectedPrecisionChanged()
{
  MandelMath::worker_multi *newWorker=nullptr;
  MandelMath::worker_multi *newStoreWorker=nullptr;
  int pointCount;
  if (imageWidth<=0 || imageHeight<=0) //Qt begins with width=0, height=-13
    pointCount=0;
  else
    pointCount=imageWidth*imageHeight;
  MandelMath::worker_multi::Type newPrecision;
  switch (_selectedPrecision)
  {
    using Type=MandelMath::worker_multi::Type;
    case precisionDouble: newPrecision=Type::typeDouble; break;
#if !ONLY_DOUBLE_WORKER
    case precisionFloat128: newPrecision=Type::typeFloat128; break;
    case precisionDDouble: newPrecision=Type::typeDDouble; break;
    case precisionQDouble: newPrecision=Type::typeQDouble; break;
    case precisionReal642: newPrecision=Type::typeReal642; break;
#endif
  }
  if (precisionRecord!=nullptr)
  {
    newWorker=MandelMath::worker_multi::create(newPrecision, precisionRecord->currentWorker->getAllocator());
    newStoreWorker=MandelMath::worker_multi::create(newPrecision, storeWorker->getAllocator());
  }
  else
  {
    newWorker=MandelMath::worker_multi::create(newPrecision, MandelModel::PrecisionRecord::LEN);
    newStoreWorker=MandelMath::worker_multi::create(newPrecision, pointCount*MandelPoint::LEN);
  }
  MandelMath::worker_multi::Allocator *newStoreAllocator=new MandelMath::worker_multi::Allocator(newStoreWorker->getAllocator(), pointCount*MandelPoint::LEN);

  {
    PrecisionRecord *newPrecisionRecord=new PrecisionRecord(newWorker, precisionRecord, this);
    delete precisionRecord;
    precisionRecord=newPrecisionRecord;
  }

  delete storeAllocator;
  delete storeWorker;
  storeWorker=newStoreWorker;
  storeAllocator=newStoreAllocator;
  //pointStore stays

  {
    QWriteLocker locker(&threading_mutex);
    startNewEpoch();
  }
}


MandelModel::Position::Position(MandelMath::worker_multi::Allocator *allocator, const Position *source):
  worker(allocator->worker), center(allocator)
{
  if (source)
  {
    //keep preinitialized values center.zero(-0.5, 0.0);
    step_log=source->step_log;
    step_size__=source->step_size__;
  }
  else
  {
    center.zero(-0.5, 0.0);
    step_log=7;
    step_size__=1.0/128;
  }
  updateCachedDepth();
}

MandelModel::Position::~Position()
{
  /*if (worker!=nullptr)
  {
    worker->cleanup(&center_re_s);
    worker->cleanup(&center_im_s);
  };*/
}

void MandelModel::Position::setView(const MandelMath::complex *c, double scale)
{
  step_log=-ilogb(scale);
  step_size__=ldexp(1.0, -step_log);
  center.assign(c);
  center.lshift(step_log);
  worker->round(center.re);
  worker->round(center.im);
  center.lshift(-step_log);
  updateCachedDepth();
}

void MandelModel::Position::move(int delta_x, int delta_y)
{
  //qDebug()<<"move ("<<delta_x<<","<<delta_y<<")";
  worker->add_double(center.re, -delta_x*step_size__);
  worker->add_double(center.im, +delta_y*step_size__);
  cached_center_re_mod=(cached_center_re_mod-delta_x+32768)%32768;
  cached_center_im_mod=(cached_center_im_mod+delta_y+32768)%32768;
#if UPDATE_CACHED_MOD
  int ccrm=cached_center_re_mod;
  int ccim=cached_center_im_mod;
  updateCachedDepth();
  if ((ccrm!=cached_center_re_mod) || (ccim!=cached_center_im_mod))
    dbgPoint();
#endif
}

void MandelModel::Position::scale(int inlog, int center_x, int center_y)
{
  //qDebug()<<"scale "<<inlog<<" @"<<center_x<<","<<center_y<<"";
  /*
  center_re+center_x*step_size = new_center_re+center_x*new_step_size
  new_step_size=step_size/2^inlog

  center_re+center_x*(step_size-new_step_size) = new_center_re
  */
  if (inlog==0)
    return;
  if (inlog>0)
  {
    int max_zoom=qRound(2-log(worker->eps2())/log(4));
    if (step_log+inlog>max_zoom)//MAX_ZOOM_IN_DOUBLE)
      return;
    double old_step_size=step_size__;
    step_log+=inlog;
    for (int i=0; i<inlog; i++)
    {
      step_size__/=2;
    }
    worker->add_double(center.re, center_x*(old_step_size-step_size__));
    worker->add_double(center.im, center_y*(old_step_size-step_size__));
    //center_re/step_size=center_re/old_step_size*old_step_size/step_size+center_x*(old_step_size-step_size)/step_size;
    cached_center_re_mod=cached_center_re_mod*(1<<inlog)+center_x*((1<<inlog)-1);
    cached_center_re_mod&=0x7fff;//=(cached_center_re_mod+32768)%32768;
    cached_center_im_mod=cached_center_im_mod*(1<<inlog)+center_y*((1<<inlog)-1);
    cached_center_im_mod&=0x7fff;//=(cached_center_im_mod+32768)%32768;
#if UPDATE_CACHED_MOD
    int ccrm=cached_center_re_mod;
    int ccim=cached_center_im_mod;
    updateCachedDepth();
    if ((ccrm!=cached_center_re_mod) || (ccim!=cached_center_im_mod))
      dbgPoint();
#endif
  }
  else
  {
    double old_step_size=step_size__;
    int adjust_x=(cached_center_re_mod+center_x)&((1<<-inlog)-1);
    if (adjust_x&(1<<(1-inlog)))
      adjust_x-=(1<<-inlog);
    int adjust_y=(cached_center_im_mod+center_y)&((1<<-inlog)-1);
    if (adjust_y&(1<<(1-inlog)))
      adjust_y-=(1<<-inlog);
    step_log+=inlog;
    for (int i=0; i<-inlog; i++)
    {
      step_size__*=2;
    }
    worker->add_double(center.re, (center_x-adjust_x)*(old_step_size-step_size__)); //(old_step_size-step_size__)=(1-(1<<-inlog))*old_step_size
    worker->add_double(center.im, (center_y-adjust_y)*(old_step_size-step_size__));
    //need to roll in high bits anyway
    /*cached_center_re_mod=cached_center_re_mod*(old_step_size/step_size)+(center_x-adjust_y)*(old_step_size/step_size-1);
    cached_center_re_mod%=32768;
    cached_center_im_mod=cached_center_im_mod*(old_step_size/step_size)+(center_y-adjust_y)*(old_step_size/step_size-1);
    cached_center_im_mod%=32768;*/
    updateCachedDepth();
  }
}

void MandelModel::Position::updateCachedDepth()
{
  MandelMath::complex d(worker);

  d.assign(&center);
  d.lshift(step_log-15);
  worker->mod1(d.re);
  worker->mod1(d.im);
  d.lshift(15);
  cached_center_re_mod=worker->toRound(d.re);
  cached_center_im_mod=worker->toRound(d.im);
}

void MandelModel::Position::pixelXtoRE(int x, MandelMath::number_pointer result)
{
  worker->assign(result, center.re);
  worker->add_double(result, x*step_size__);
}

void MandelModel::Position::pixelYtoIM(int y, MandelMath::number_pointer result)
{
  worker->assign(result, center.im);
  worker->add_double(result, y*step_size__);
}



MandelModel::Orbit::Orbit(MandelMath::worker_multi::Allocator *allocator):
  currentWorker(allocator->worker),
  evaluator(allocator->worker->ntype(), true),
  bulb(allocator)
{
}

MandelModel::Orbit::~Orbit()
{
  evaluator.workIfEpoch=-1;
  evaluator.quit();
  evaluator.wait(1000);
}

MandelModel::Orbit::Bulb::Bulb(MandelMath::worker_multi::Allocator *allocator):
  cb(allocator), rb(allocator), xc(allocator), baseZC(allocator), baseCC(allocator)
{
}

MandelModel::Orbit::Bulb::~Bulb()
{
}

MandelModel::PrecisionRecord::PrecisionRecord(MandelMath::worker_multi *newWorker, PrecisionRecord *source, MandelModel *doneReceiver):
  currentWorker(newWorker),
  shareableViewInfoAllocator(newWorker->getAllocator(), ShareableViewInfo::LEN),
  shareableVIAuser(&shareableViewInfoAllocator, ShareableViewInfo::LEN),
  wtiPoint(nullptr, newWorker->getAllocator()),
  position(newWorker->getAllocator(), source?&source->position:nullptr),
  orbit(newWorker->getAllocator()),//, source?&source->orbit:nullptr),
  lagu_c(newWorker->getAllocator()), lagu_r(newWorker->getAllocator()), tmp(newWorker->getAllocator()),
  threadCount(source?source->threadCount:0), threads(nullptr)
{
  if (source)
  {
  }
  else
  {
    lagu_c.zero(0, 0);
    lagu_r.zero(0, 0);
  }
  if (threadCount<=0)
  { //first init
    threadCount=QThread::idealThreadCount()-1;
    if (threadCount<1)
      threadCount=1;
  };
  threads=new MandelEvaluator *[threadCount];
  for (int t=0; t<threadCount; t++)
  {
    threads[t]=new MandelEvaluator(newWorker->ntype(), false);
    //Qt not allowing parameters in invokeMethod makes it really easy... pfft
    threads[t]->threaded.give=[doneReceiver](MandelEvaluator *me)
        {
          return doneReceiver->giveWorkThreaded(me);
        };
    threads[t]->threaded.done=[doneReceiver](MandelEvaluator *me, bool giveWork)
        {
          return doneReceiver->doneWorkThreaded(me, giveWork);
        };
    //threads[t].setHint(t);
    //QObject::connect(threads[t], &MandelEvaluator::doneCompute,
    //                 doneReceiver, &MandelModel::donePixel,
    //                 Qt::ConnectionType::QueuedConnection);
    QObject::connect(threads[t], &MandelEvaluator::doneComputeThreaded,
                     doneReceiver, &MandelModel::doneWorkInThread,
                     Qt::ConnectionType::QueuedConnection);
    QObject::connect(doneReceiver, &MandelModel::triggerComputeThreaded,
                     threads[t], &MandelEvaluator::doComputeThreaded,
                     Qt::ConnectionType::QueuedConnection);
  }
}

MandelModel::PrecisionRecord::~PrecisionRecord()
{
  for (int t=threadCount-1; t>=0; t--)
  {
    threads[t]->workIfEpoch=-1;
    threads[t]->quit();
  }
  for (int t=threadCount-1; t>=0; t--)
  {
    threads[t]->wait(1000);
    delete threads[t];
  }
  delete[] threads;
  threadCount=0;
  threads=nullptr;
}

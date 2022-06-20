#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "MandelModel.hpp"
#include "MandelEvaluator.hpp"

MandelModel::MandelModel(): QObject(),
  shareableViewInfoAllocator(nullptr), currentWorker(nullptr),
  storeAllocator(nullptr), storeWorker(nullptr), pointStore_(nullptr),
  wtiPointAllocator(nullptr), wtiPoint(nullptr),
  position_(nullptr), orbit_(nullptr)
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
  epoch=0;
  imageWidth=0;
  imageHeight=0;
  //pointStore=nullptr;
  lastGivenPointIndex_=0;
  effortBonus_=0;
  //threadCount=4;
  _threadsWorking=0;
  threadCount=0;
  threads=nullptr;
  QObject::connect(this, &MandelModel::selectedPrecisionChange,
                   this, &MandelModel::selectedPrecisionChanged);
  selectedPrecisionChanged();
}

MandelModel::~MandelModel()
{
  for (int t=threadCount-1; t>=0; t--)
  {
    threads[t]->wantStop=true;
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

  //orbit.evaluator.switchType(nullptr);
  if (orbit_!=nullptr)
  {
    orbit_->evaluator.wantStop=true;
    orbit_->evaluator.quit();
    orbit_->evaluator.wait(1000);
    delete orbit_;
  }
  delete position_;
  delete wtiPoint;
  delete wtiPointAllocator;
  delete shareableViewInfoAllocator;

  //storeWorker->getAllocator()->dealloc_array(imageWidth*imageHeight);
  delete storeAllocator;
  delete storeWorker;
  delete currentWorker;
  storeWorker=nullptr;
  currentWorker=nullptr;
  /*for (int i=imageWidth*imageHeight-1; i>=0; i--)
  {
    //pointStore[i].cleanup(position.worker);
  }*/

  delete[] pointStore_;
  pointStore_=nullptr;
  imageWidth=0;
  imageHeight=0;
}

QString MandelModel::pixelXtoRE_str(int x)
{
  MandelMath::number num(currentWorker);
  num.assign(position_->center.re);
  num.add_double((x - imageWidth/2)*position_->step_size__);
  QString result=currentWorker->toString(num.ptr);
  return result;
}

QString MandelModel::pixelYtoIM_str(int y)
{
//  return (y - imageHeight/2)*position.step_size+position.center_im;
  MandelMath::number num(currentWorker);
  num.assign(position_->center.im);
  num.add_double((y - imageHeight/2)*position_->step_size__);
  QString result=currentWorker->toString(num.ptr);
  return result;
}

QString MandelModel::getTimes()
{
  QString result;
  for (int t=0; t<threadCount; t++)
    result+=QString("%1-%2[%3,%4],").
        arg((threads[t]->timeOuterTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t]->timeInnerTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t]->timeInvokePostTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t]->timeInvokeSwitchTotal)/1000000000.0, 0, 'f', 3);
  return result;
}

QString MandelModel::getTextXY()
{
  if (currentWorker==nullptr || orbit_==nullptr)
    return "-";
  return currentWorker->toString(orbit_->evaluator.currentParams.c.re)+" +i* "+
         currentWorker->toString(orbit_->evaluator.currentParams.c.im);
}

QString mandDoubleToString(double x)
{
  if (x>0)
    return QString("+%1").arg(x, 10, 'f');
  else
    return QString("%1").arg(x, 10, 'f');//returns like 50 digits QString::number(x, 10, 'f');
}

QString MandelModel::getTextInfoGen()
{
  if (currentWorker==nullptr || orbit_==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number tmp(currentWorker);
    tmp.assign(orbit_->evaluator.currentParams.c.re);
    tmp.sub(position_->center.re);
    orbit_x=qRound(currentWorker->toDouble(tmp.ptr)/position_->step_size__)+imageWidth/2;

    tmp.assign(orbit_->evaluator.currentParams.c.im);
    tmp.sub(position_->center.im);
    orbit_y=imageHeight/2-qRound(currentWorker->toDouble(tmp.ptr)/position_->step_size__);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  MandelPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);

  QString state;
  switch (data_.store->state)
  {
    case MandelPointStore::State::stUnknown:
      state="Unk"; break;
    case MandelPointStore::State::stOutside:
      state="Out"; break;
    case MandelPointStore::State::stOutAngle:
      state="OutA"; break;
    case MandelPointStore::State::stBoundary:
      state="Bound"; break;
    case MandelPointStore::State::stDiverge:
      state="Diver"; break;
    case MandelPointStore::State::stMisiur:
      state="Misiur"; break;
    case MandelPointStore::State::stPeriod2:
      state="Per2"; break;
    case MandelPointStore::State::stPeriod3:
      state="Per3"; break;
    case MandelPointStore::State::stMaxIter:
      state="Max"; break;
  }

  return state+" iter="+QString::number(data_.store->iter)+" near="+QString::number(data_.store->near0iter)+
      " fc="+mandDoubleToString(storeWorker->toDouble(data_.fc_c.re))+mandDoubleToString(storeWorker->toDouble(data_.fc_c.im))+"i";
}

QString MandelModel::getTextInfoSpec()
{
  if (currentWorker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number tmp(currentWorker);
    reimToPixel(&orbit_x, &orbit_y, &orbit_->evaluator.currentParams.c, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  MandelPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);

  switch (data_.store->state)
  {
    case MandelPointStore::State::stUnknown:
      return " ";
      break;
    case MandelPointStore::State::stOutside:
      return QString("ext=")+QString::number(data_.store->exterior_hits); break;
    case MandelPointStore::State::stOutAngle:
      return QString("ext=")+QString::number(data_.store->exterior_hits); break;
    case MandelPointStore::State::stBoundary:
      return " ";
    case MandelPointStore::State::stDiverge:
      return " ";
    case MandelPointStore::State::stMisiur:
      return " ";
    case MandelPointStore::State::stPeriod2:
      return QString("per=")+QString::number(data_.store->period)+" int="+QString::number(data_.store->interior)   +" mult="+QString::number(orbit_->evaluator.bulb.dbg_guessmult); break;
    case MandelPointStore::State::stPeriod3:
      return QString("per=")+QString::number(data_.store->period)+" int="+QString::number(data_.store->interior)   +" mult="+QString::number(orbit_->evaluator.bulb.dbg_guessmult); break;
    case MandelPointStore::State::stMaxIter:
      return " "; break;
  }
  return "-?-?-";
}

ShareableViewInfo MandelModel::getViewInfo()
{
  ShareableViewInfo result(shareableViewInfoAllocator);
  //result.worker=orbit.worker;
  result.period=orbit_->pointData.store->near0iter;//evaluator.currentData.lookper_lastGuess;//orbit.pointData.period;
  if (result.period<1)
    result.period=1;
  result.scale=position_->step_size__;
  /*orbit.worker->init_(&result.re_, &result.re_p);
  orbit.worker->init_(&result.im, &result.im_p);
  orbit.worker->init_(&result.root_re, &result.rre_p);
  orbit.worker->init_(&result.root_im, &result.rim_p);*/
  result.c.assign(&orbit_->evaluator.currentParams.c);
  result.root.assign(&orbit_->evaluator.currentData.root);

  //TODO: why here? should be somewhere else
  orbit_->lagu_c.assign(&orbit_->evaluator.currentParams.c);
  orbit_->lagu_r.assign(&orbit_->evaluator.currentData.root);

  return result;
}

void MandelModel::transformStore(MandelMath::worker_multi *old_worker, MandelMath::worker_multi *old_sworker, MandelPointStore *old_store, int old_width, int old_height, const MandelMath::complex *old_c,
                                 MandelMath::worker_multi *new_worker, MandelMath::worker_multi *new_sworker, MandelPointStore *new_store, int new_width, int new_height, const MandelMath::complex *new_c,
                                 int inlog, int new_step_log)
{
  if (currentWorker==nullptr)
  {
    dbgPoint();
    return;
  };
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
    new_worker->assign(c.im, position_->center.im);
    new_worker->add_double(c.im, (new_height/2-newy)*position_->step_size__);
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
        {
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          new_sworker->assign_block((newy*new_width+newx)*MandelPoint::LEN, old_sworker, (oldy*old_width+oldx)*MandelPoint::LEN, MandelPoint::LEN);
        }
        else
        {
          new_worker->assign(c.re, position_->center.re);
          new_worker->add_double(c.re, (newx - new_width/2)*position_->step_size__);
          MandelMath::worker_multi::Allocator allo(new_sworker->getAllocator(), (newy*new_width+newx)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
          MandelPoint pt(&new_store[newy*new_width+newx], &allo);
          pt.zero(/*new_sworker, (newy*new_width+newx)*MandelPoint::LEN,*/ &c);
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
    new_worker->assign(c.im, position_->center.im);
    new_worker->add_double(c.im, (new_height/2-newy)*position_->step_size__);
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
        {
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          new_sworker->assign_block((newy*new_width+newx)*MandelPoint::LEN, old_sworker, (oldy*old_width+oldx)*MandelPoint::LEN, MandelPoint::LEN);
        }
        else
        {
          new_worker->assign(c.re, position_->center.re);
          new_worker->add_double(c.re, (newx - new_width/2)*position_->step_size__);
          MandelMath::worker_multi::Allocator allo(new_sworker->getAllocator(), (newy*new_width+newx)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
          MandelPoint pt(&new_store[newy*new_width+newx], &allo);
          pt.zero(/*new_worker, (newy*new_width+newx)*MandelPoint::LEN,*/ &c);
        }
      }
    }
  }
}

void MandelModel::setView_double(double c_re, double c_im, double scale)
{
  MandelMath::complex c(currentWorker);
  c.zero(c_re, c_im);
  setView_(&c, scale);
}

void MandelModel::setView_(const MandelMath::complex *c, double scale)
{
  MandelMath::complex old_c(currentWorker);
  old_c.assign(&position_->center);
  int old_step_log=position_->step_log;

  position_->setView(c, scale);

  transformStore(currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &position_->center,
                 position_->step_log-old_step_log, position_->step_log);

  startNewEpoch();
}

void MandelModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::complex old_c(currentWorker);
  old_c.assign(&position_->center);

  int dx=qRound(delta_x);
  int dy=qRound(delta_y);
  position_->move(dx, dy);
  //qDebug()<<"new c: re="<<position.worker->toString(&position.center_re_s)<<",im="<<position.worker->toString(&position.center_im_s);

  transformStore(currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &position_->center,
                 0, position_->step_log);

  startNewEpoch();
}

void MandelModel::zoom(double x, double y, int inlog)
{
  MandelMath::complex old_c(currentWorker);
  old_c.assign(&position_->center);
  int old_step_log=position_->step_log;

  position_->scale(inlog, qRound(x)-imageWidth/2, imageHeight/2-qRound(y));

  transformStore(currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &position_->center,
                 position_->step_log-old_step_log, position_->step_log);

  startNewEpoch();
}

void MandelModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  int newLength=width*height;
  MandelPointStore *newStore_=new MandelPointStore[newLength];
  {
    QString size_as_text=QString::number(newLength*sizeof(MandelPointStore));
    for (int pos=size_as_text.length()-3; pos>0; pos-=3)
      size_as_text.insert(pos, '\'');
    qDebug()<<"pointStore uses"<<size_as_text.toLocal8Bit().constData()<<"B"; //lots of work to skip those quotes... can't skip spaces at all
  }
  MandelMath::worker_multi *newStoreWorker;
  MandelMath::worker_multi::Allocator *newStoreAllocator;
  switch (storeWorker->ntype())
  {
    default://case MandelMath::worker_multi::Type::typeEmpty:
      dbgPoint();
      goto lolwut;
    case MandelMath::worker_multi::Type::typeDouble: lolwut:
      newStoreWorker=new MandelMath::worker_multi_double(newLength*MandelPoint::LEN);
      break;
#if !ONLY_DOUBLE_WORKER
    case MandelMath::worker_multi::Type::typeFloat128:
      newStoreWorker=new MandelMath::worker_multi_float128(newLength*MandelPoint::LEN);
      break;
    case MandelMath::worker_multi::Type::typeDDouble:
      newStoreWorker=new MandelMath::worker_multi_ddouble(newLength*MandelPoint::LEN);
      break;
    case MandelMath::worker_multi::Type::typeQDouble:
      newStoreWorker=new MandelMath::worker_multi_qdouble(newLength*MandelPoint::LEN);
      break;
#endif
  }
  newStoreAllocator=new MandelMath::worker_multi::Allocator(newStoreWorker->getAllocator(), newLength*MandelPoint::LEN);
  MandelMath::complex old_c(currentWorker);
  old_c.assign(&position_->center);

  transformStore(currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 currentWorker, newStoreWorker, newStore_, width, height, &position_->center,
                 0, position_->step_log);

  //for (int i=imageWidth*imageHeight-1; i>=0; i--)
    //pointStore[i].cleanup(position.worker);
  //storeWorker->dealloc_array(imageWidth*imageHeight*MandelPoint::LEN);
  delete storeAllocator;
  delete storeWorker;
  storeWorker=newStoreWorker;
  storeAllocator=newStoreAllocator;
  delete[] pointStore_;
  pointStore_=newStore_;
  imageWidth=width;
  imageHeight=height;

  startNewEpoch();
}

void MandelModel::startNewEpoch()
{
  epoch=(epoch+1)%2000000000;
  lastGivenPointIndex_=0;
  effortBonus_=0;
  for (int t=0; t<threadCount; t++)
    if (threads[t]->currentParams.pixelIndex<0)
      giveWork(threads[t]);
}

void MandelModel::reimToPixel(int *circ_x, int *circ_y, const MandelMath::complex *c, MandelMath::number *tmp)
{
  tmp->assign(c->re);
  tmp->sub(position_->center.re);
  //currentWorker_->lshift(tmp, position_->step_log);
  //*circ_x=currentWorker_->toDouble(tmp)+imageWidth/2;
  *circ_x=qRound(currentWorker->toDouble(tmp->ptr)/position_->step_size__)+imageWidth/2;

  tmp->assign(c->im);
  tmp->sub(position_->center.im);
  //currentWorker_->lshift(tmp, position_->step_log);
  //*circ_y=imageHeight/2-currentWorker_->toDouble(tmp);
  *circ_y=imageHeight/2-qRound(currentWorker->toDouble(tmp->ptr)/position_->step_size__);
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

  if (orbit_==nullptr)
    dbgPoint();
  else if (currentWorker==nullptr)
    dbgPoint();
  else if (orbit_->currentWorker->ntype()!=currentWorker->ntype())
    dbgPoint();

  MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  MandelPoint data_(&pointStore_[y*imageWidth+x], &allo);
  switch (data_.store->state)
  {
    case MandelPointStore::State::stOutside:
    case MandelPointStore::State::stOutAngle:
    {
      int exterior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0xff, 0xff, 0));
      exterior=qRound(data_.store->exterior_hits/position_->step_size__);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
      painter.setPen(QColor(0xc0, 0xc0, 0));
      exterior=qRound(data_.store->exterior_avoids/position_->step_size__);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
    } break;
    case MandelPointStore::State::stPeriod2:
    case MandelPointStore::State::stPeriod3:
    {
      int interior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0, 0xff, 0xff));
      interior=qRound(data_.store->interior/position_->step_size__);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
      painter.setPen(QColor(0, 0xc0, 0xc0));
      interior=qRound(data_.store->interior/4/position_->step_size__);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
    } break;
    default: ;
  }
  position_->pixelXtoRE(x-imageWidth/2, orbit_->evaluator.currentParams.c.re);
  position_->pixelYtoIM(imageHeight/2-y, orbit_->evaluator.currentParams.c.im);
  orbit_->evaluator.currentParams.epoch=epoch;
  orbit_->evaluator.currentParams.pixelIndex=0;
  orbit_->pointData.zero(&orbit_->evaluator.currentParams.c);
  orbit_->evaluator.currentParams.breakOnNewNearest=true;
  orbit_->evaluator.currentParams.maxiter_=1<<MAX_EFFORT;
  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0xff, 0xff)); //paint path
  while ((orbit_->pointData.store->state==MandelPointStore::State::stUnknown) &&
         (orbit_->pointData.store->iter<(1<<MAX_EFFORT)))
  {
    int line_sx, line_sy, line_ex, line_ey;
    reimToPixel(&line_sx, &line_sy, &orbit_->pointData.f, &orbit_->tmp);

    if ((data_.store->state==MandelPointStore::State::stPeriod2 || data_.store->state==MandelPointStore::State::stPeriod3) &&
        orbit_->pointData.store->iter<data_.store->period)
    { //paint first period fully
      orbit_->evaluator.currentParams.maxiter_=orbit_->pointData.store->iter+1;
    }
    else if (orbit_->evaluator.currentData.store->lookper_lastGuess==0)
      orbit_->evaluator.currentParams.maxiter_=1<<MAX_EFFORT; //dont't know->run fully
    else //stop at multiples of lookper
      orbit_->evaluator.currentParams.maxiter_=(orbit_->pointData.store->iter/orbit_->evaluator.currentData.store->lookper_lastGuess+1)*orbit_->evaluator.currentData.store->lookper_lastGuess;
    orbit_->evaluator.startCompute(&orbit_->pointData, +1);
    orbit_->pointData.assign(orbit_->evaluator.currentData);

    reimToPixel(&line_ex, &line_ey, &orbit_->pointData.f, &orbit_->tmp);
    painter.drawLine(line_sx, line_sy, line_ex, line_ey);
  }
  if ((orbit_->pointData.store->state==MandelPointStore::State::stPeriod2) ||
      (orbit_->pointData.store->state==MandelPointStore::State::stPeriod3))
  {
    int circ_x, circ_y;
    painter.setPen(QColor(0, 0xff, 0xff)); //paint root
    reimToPixel(&circ_x, &circ_y, &orbit_->pointData.root, &orbit_->tmp);
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
    orbit_->bulb.foundMult=0;
    orbit_->bulb.is_card=false;

    orbit_->bulb.valid=orbit_->evaluator.bulb.findBulbBase(orbit_->pointData.store->period,
        &orbit_->evaluator.currentParams.c, &orbit_->bulb.cb, &orbit_->bulb.rb,
        &orbit_->bulb.xc, &orbit_->bulb.baseZC, &orbit_->bulb.baseCC,
        &orbit_->bulb.is_card, &orbit_->bulb.foundMult);
    if (orbit_->bulb.valid)
    {
      painter.setBrush(QBrush(QColor(0, 0xff, 0xff)));
      reimToPixel(&circ_x, &circ_y, &orbit_->evaluator.bulb.dbg_first_cb, &orbit_->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0x80, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &orbit_->evaluator.bulb.dbg_first_rb, &orbit_->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0x80, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &orbit_->bulb.xc, &orbit_->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        if (orbit_->bulb.is_card)
          painter.setPen(QColor(0, 0xff, 0xff)); //card center
        else
          painter.setPen(QColor(0x80, 0xc0, 0xc0)); //bulb center
        //painter.setBrush(Qt::BrushStyle::SolidPattern);
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &orbit_->bulb.cb, &orbit_->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0xff, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &orbit_->bulb.rb, &orbit_->tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0xff, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
    };
  }

  if (!orbit_->currentWorker->is0(orbit_->lagu_c.re) ||
      !orbit_->currentWorker->is0(orbit_->lagu_c.im))
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0, 0xff, 0)); //paint c
    reimToPixel(&circ_x, &circ_y, &orbit_->lagu_c, &orbit_->tmp);
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
    reimToPixel(&circ_x, &circ_y, &orbit_->lagu_r, &orbit_->tmp);
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
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      bool knownenum=false;
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
      //MandelPoint data_(&pointStore_[y*imageWidth+x], &allo);
      wtiPoint->store=&pointStore_[y*imageWidth+x];
      {
        int first, last;
        wtiPointAllocator->_getRange(first, last);
        currentWorker->assign_block(first, storeWorker, (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN);
      }
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleKind:
        {
          switch (wtiPoint->store->state)
          {
            case MandelPointStore::State::stUnknown:
              image.image->setPixel(x, y, 0xffffffff);
              //image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
              break;
            case MandelPointStore::State::stOutside:
            case MandelPointStore::State::stOutAngle:
            {
              int b=0x9f+floor(0x60*cos((wtiPoint->store->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b<<0));
              knownenum=true;
            } break;
            case MandelPointStore::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
              knownenum=true;
            } break;
            case MandelPointStore::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
              knownenum=true;
            } break;
            case MandelPointStore::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
              knownenum=true;
            } break;
            case MandelPointStore::State::stPeriod2:
            case MandelPointStore::State::stPeriod3:
            {
              /*int r;
              switch (data->period)
              {
                case 1: r=0xc0; break;
                case 2: r=0xff; break;
                case 3: r=0x80; break;
                default: r=0xe0;
              }*/
              if (wtiPoint->store->period>wtiPoint->store->near0iter)
                image.image->setPixel(x, y, 0xffff00ff); //seems to only happen by mistake, not in reality
              else
              {
                int index=periodToIndex(wtiPoint->store->period);
                //reverse bottom 7 bits:
                int rh=0x73516240>>((index&7)<<2);
                int rl=0x73516240>>((index&0x70)>>2);
                rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
                image.image->setPixel(x, y, 0xff000000+(rh<<16));
              }
              knownenum=true;
            } break;
            case MandelPointStore::State::stMaxIter:
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
          switch (wtiPoint->store->state)
          {
            case MandelPointStore::State::stUnknown:
              image.image->setPixel(x, y, 0x00906090);
              knownenum=true;
              break;
            case MandelPointStore::State::stOutside:
            case MandelPointStore::State::stOutAngle:
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
              double re=position_->worker->toDouble(wtiPoint->f.re);
              double im=position_->worker->toDouble(wtiPoint->f.im);
              double iter=wtiPoint->store->iter+6-log2(log2(re*re+im*im)); //+6 to match integer coloring
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
            case MandelPointStore::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPointStore::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPointStore::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPointStore::State::stPeriod2:
            case MandelPointStore::State::stPeriod3:
            {
              int index=periodToIndex(wtiPoint->store->period);
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
            case MandelPointStore::State::stMaxIter:
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
          switch (wtiPoint->store->state)
          {
            case MandelPointStore::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::State::stOutside:
            case MandelPointStore::State::stOutAngle:
            {
              double tf;
              if ((wtiPoint->store->exterior_avoids>10000) || (wtiPoint->store->exterior_avoids<=0))
                tf=0;
              else if (wtiPoint->store->exterior_avoids>=1)
                tf=(1-wtiPoint->store->exterior_avoids)*1;
              else
                tf=sqrt(1-log(wtiPoint->store->exterior_avoids))*2-2;
              int r=0x9f+qRound(0x60*sin(tf*2.828)); //red middle
              int g=0x9f+qRound(0x60*sin(tf*6.928)); //green fastest
              int b=0x9f+qRound(0x60*sin(tf)); //blue slowest
              if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              knownenum=true;
            } break;
            case MandelPointStore::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stPeriod2:
            case MandelPointStore::State::stPeriod3:
            {
              double re=position_->worker->toDouble(wtiPoint->fz_r.re);
              double im=position_->worker->toDouble(wtiPoint->fz_r.im);
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
            case MandelPointStore::State::stMaxIter:
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
          switch (wtiPoint->store->state)
          {
            case MandelPointStore::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::State::stOutside:
            case MandelPointStore::State::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiPoint->store->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPointStore::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stPeriod2:
            case MandelPointStore::State::stPeriod3:
            {
              int ti=30;
              if ((wtiPoint->store->interior>1) || (wtiPoint->store->interior<=0))
                ti=0;
              else
                ti=(qRound(-log(wtiPoint->store->interior/4)*300)+12*0xc0) % (6*0xc0);
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
            case MandelPointStore::State::stMaxIter:
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
          switch (wtiPoint->store->state)
          {
            case MandelPointStore::State::stUnknown:
              image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
              break;
            case MandelPointStore::State::stOutside:
            case MandelPointStore::State::stOutAngle:
            {
              int ti=wtiPoint->store->near0iter;
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
            case MandelPointStore::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
              knownenum=true;
            } break;
            case MandelPointStore::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
              knownenum=true;
            } break;
            case MandelPointStore::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
              knownenum=true;
            } break;
            case MandelPointStore::State::stPeriod2:
            case MandelPointStore::State::stPeriod3:
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

              int ti=wtiPoint->store->near0iter;
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
            case MandelPointStore::State::stMaxIter:
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
          switch (wtiPoint->store->state)
          {
            case MandelPointStore::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::State::stOutside:
            case MandelPointStore::State::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiPoint->store->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPointStore::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stPeriod2:
            case MandelPointStore::State::stPeriod3:
            {
              double re=position_->worker->toDouble(wtiPoint->fz_r.re);
              double im=position_->worker->toDouble(wtiPoint->fz_r.im);
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
            case MandelPointStore::State::stMaxIter:
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
          switch (wtiPoint->store->state)
          {
            case MandelPointStore::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPointStore::State::stOutside:
            case MandelPointStore::State::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiPoint->store->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPointStore::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPointStore::State::stPeriod2:
            case MandelPointStore::State::stPeriod3:
            {
              if (!wtiPoint->store->has_fc_r)
                image.image->setPixel(x, y, 0xffc0c0c0);
              else
              {
                double re=position_->worker->toDouble(wtiPoint->fc_c.re);
                double im=position_->worker->toDouble(wtiPoint->fc_c.im);
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
            case MandelPointStore::State::stMaxIter:
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

void MandelModel::giveWork(MandelEvaluator *evaluator)
{
  int retryEffortFrom=0;
  while (retryEffortFrom>=0)
  {
    retryEffortFrom=-1;
    int quickrun=0;
    for (int pi=0; pi<imageWidth*imageHeight; pi++)
    {
      int pointIndex=(lastGivenPointIndex_+pi)%(imageWidth*imageHeight);
      //if ((lastGivenPointIndex_!=0) && (pointIndex==0))
        //dbgPoint();
      MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), pointIndex*MandelPoint::LEN, MandelPoint::LEN, nullptr);
      MandelPoint pointData_(&pointStore_[pointIndex], &allo);
      bool needsEval=(pointData_.store->state==MandelPointStore::State::stUnknown);
      if (!needsEval)
        needsEval=(_selectedPaintStyle==paintStyleFC) &&
                  (!pointData_.store->has_fc_r) &&
                  ((pointData_.store->state==MandelPointStore::State::stPeriod2) ||
                   (pointData_.store->state==MandelPointStore::State::stPeriod3));
      if (needsEval)
      {
        bool found=false;
        for (int t=0; t<threadCount; t++)
          if (threads[t]->currentParams.pixelIndex==pointIndex)
          {
            found=true;
            break;
          };
        if (evaluator->currentParams.pixelIndex!=-1)
          dbgPoint();
        assert(evaluator->currentParams.pixelIndex==-1);
        if (!found)
        {
          int phasex=(pointIndex%imageWidth-imageWidth/2+position_->cached_center_re_mod+32768)%32768;
          int phasey=(pointIndex/imageWidth-imageHeight/2+position_->cached_center_im_mod+32768)%32768;
          //int effort=ctz16(phasex)+ctz16(phasey);
          int effort=MandelMath::ctz16(phasex | phasey);
          if (effort>8)
            effort=8;
          effort+=effortBonus_;
          if (effort>=MAX_EFFORT)
            effort=MAX_EFFORT;
          evaluator->currentParams.maxiter_=1<<effort;
          if (pointData_.store->iter>=evaluator->currentParams.maxiter_)
          {
            if (effort>=MAX_EFFORT)
              pointData_.store->state=MandelPointStore::State::stMaxIter;
            else if (retryEffortFrom<0)
              retryEffortFrom=pointIndex;
          }
          else
          {
            if (evaluator->currentWorker->ntype()!=position_->worker->ntype())
              dbgPoint();
            //evaluator->switchType(position.worker);
            position_->pixelXtoRE(pointIndex%imageWidth - imageWidth/2, evaluator->currentParams.c.re);
            position_->pixelYtoIM(imageHeight/2-pointIndex/imageWidth, evaluator->currentParams.c.im);
            evaluator->currentParams.epoch=epoch;
            evaluator->currentParams.pixelIndex=pointIndex;
            evaluator->currentParams.want_fc_r=(_selectedPaintStyle==paintStyleFC);
            if (evaluator->startCompute(&pointData_, quickrun>=100?-1:0))
            //if (worker->startCompute(pointData, true))
            {
              _threadsWorking++;
              evaluator->timeOuter.start();
              lastGivenPointIndex_=pointIndex;
              return;
            }
            else
            {
              donePixel1(evaluator);
              quickrun++;
            }
          }
        }
      }
        //MandelEvaluator::simple(cr, ci, pointStore[y*imageWidth+x]);
    }
    if ((retryEffortFrom>=0) && (effortBonus_<MAX_EFFORT))
    {
      effortBonus_++;
      lastGivenPointIndex_=retryEffortFrom;
    }
    else
      retryEffortFrom=-1;
  }
}

void MandelModel::donePixel1(MandelEvaluator *me)
{
  if ((me->currentParams.epoch==epoch) && (me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
  {
    MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), me->currentParams.pixelIndex*MandelPoint::LEN, MandelPoint::LEN, nullptr);
    MandelPoint point_(&pointStore_[me->currentParams.pixelIndex], &allo);
    /* it's OK now as fc_r is computed later
    if (point->state!=MandelPoint::State::stUnknown)
      qDebug()<<"Finished pixel finished again";
    else*/
    {
      if (position_->worker==nullptr)
        dbgPoint();
      else
        point_.assign(me->currentData);
      if ((point_.store->state==MandelPointStore::State::stUnknown) &&
          (point_.store->iter>=(1<<MAX_EFFORT)))
        point_.store->state=MandelPointStore::State::stMaxIter;
    }
  }
  else if (me->currentParams.epoch!=epoch)
  { }//qDebug()<<"Old pixel finished";
  else
    qWarning()<<"Invalid pixel finished";
  me->currentParams.pixelIndex=-1;
}

void MandelModel::donePixel(MandelEvaluator *me)
{
  me->timeOuterTotal+=me->timeOuter.nsecsElapsed();
  _threadsWorking--;
  donePixel1(me);
  giveWork(me);
}

void MandelModel::selectedPrecisionChanged()
{
  MandelMath::worker_multi *newWorker=nullptr;
  MandelMath::worker_multi *newStoreWorker=nullptr;
  MandelMath::worker_multi::Allocator *newStoreAllocator=nullptr;
  if (currentWorker!=nullptr)
  {
    for (int t=threadCount-1; t>=0; t--)
    {
      threads[t]->wantStop=true;
      threads[t]->quit();
    }
    for (int t=threadCount-1; t>=0; t--)
    {
      threads[t]->wait(1000);
      delete threads[t];
    }
    delete threads;
    threads=nullptr;
    //keep old number because reallocating worker would not reflect new length//threadCount=0;

    //bool empty=imageWidth<=0 || imageHeight<=0; //Qt begins with width=0, height=-13
    switch (_selectedPrecision)
    {
      case precisionDouble:
        newWorker=new MandelMath::worker_multi_double(currentWorker->getAllocator());
        newStoreWorker=new MandelMath::worker_multi_double(storeWorker->getAllocator());
        break;
#if !ONLY_DOUBLE_WORKER
      case precisionFloat128:
        newWorker=new MandelMath::worker_multi_float128(currentWorker_);
        newStoreWorker=new MandelMath::worker_multi_float128(storeWorker);
        break;
      case precisionDDouble:
        newWorker=new MandelMath::worker_multi_ddouble(currentWorker_);
        newStoreWorker=new MandelMath::worker_multi_ddouble(storeWorker);
        break;
      case precisionQDouble:
        newWorker=new MandelMath::worker_multi_qdouble(currentWorker_);
        newStoreWorker=new MandelMath::worker_multi_qdouble(storeWorker);
        break;
#endif
    }
    newStoreAllocator=new MandelMath::worker_multi::Allocator(newStoreWorker->getAllocator(), imageWidth*imageHeight*MandelPoint::LEN);//newStoreWorker->capacity);
  }
  else
  {
    int pointCount;
    if (imageWidth<=0 || imageHeight<=0) //Qt begins with width=0, height=-13
      pointCount=0;
    else
      pointCount=imageWidth*imageHeight;
    if (threadCount!=0)
      dbgPoint();
    threadCount=1;//QThread::idealThreadCount()-1;
    if (threadCount<1)
      threadCount=1;
    switch (_selectedPrecision)
    {
      case precisionDouble:
        newWorker=new MandelMath::worker_multi_double(MandelModel::LEN+threadCount*MandelEvaluator::LEN);
        newStoreWorker=new MandelMath::worker_multi_double(pointCount*MandelPoint::LEN);
        break;
#if !ONLY_DOUBLE_WORKER
      case precisionFloat128:
        newWorker=new MandelMath::worker_multi_float128(MandelModel::LEN+threadCount*MandelEvaluator::LEN);
        newStoreWorker=pointCount==new MandelMath::worker_multi_float128(pointCount*MandelPoint::LEN);
        break;
      case precisionDDouble:
        newWorker=new MandelMath::worker_multi_ddouble(MandelModel::LEN+threadCount*MandelEvaluator::LEN);
        newStoreWorker=pointCount==new MandelMath::worker_multi_ddouble(pointCount*MandelPoint::LEN);
        break;
      case precisionQDouble:
        newWorker=new MandelMath::worker_multi_qdouble(MandelModel::LEN+threadCount*MandelEvaluator::LEN);
        newStoreWorker=pointCount==new MandelMath::worker_multi_qdouble(pointCount*MandelPoint::LEN);
        break;
#endif
    }
    newStoreAllocator=new MandelMath::worker_multi::Allocator(newStoreWorker->getAllocator(), pointCount*MandelPoint::LEN);
  }

  {
    MandelMath::worker_multi::Allocator *newshareableViewInfoAllocator=new MandelMath::worker_multi::Allocator(newWorker->getAllocator(), ShareableViewInfo::LEN);
    MandelMath::worker_multi::Allocator *newWtiPointAllocator=new MandelMath::worker_multi::Allocator(newWorker->getAllocator(), MandelPoint::LEN);
    MandelPoint *newWtiPoint=new MandelPoint(nullptr, newWtiPointAllocator);
    Position *newPosition=new Position(newWorker->getAllocator());
    Orbit *newOrbit=new Orbit(newWorker->getAllocator());
    newPosition->assign(position_);
    delete orbit_;
    delete position_;
    delete wtiPoint;
    delete wtiPointAllocator;
    delete shareableViewInfoAllocator;
    shareableViewInfoAllocator=newshareableViewInfoAllocator;
    wtiPointAllocator=newWtiPointAllocator;
    wtiPoint=newWtiPoint;
    position_=newPosition;
    orbit_=newOrbit;
  }

  threads=new MandelEvaluator *[threadCount];
  for (int t=0; t<threadCount; t++)
  {
    threads[t]=new MandelEvaluator(newWorker->getAllocator(), nullptr);
    //threads[t].setHint(t);
    QObject::connect(threads[t], &MandelEvaluator::doneCompute,
                     this, &MandelModel::donePixel,
                     Qt::ConnectionType::QueuedConnection);
  }

  delete currentWorker;
  delete storeAllocator;
  delete storeWorker;
  currentWorker=newWorker;
  storeWorker=newStoreWorker;
  storeAllocator=newStoreAllocator;
  //pointStore stays

  startNewEpoch();
}

MandelModel::Position::Position(MandelMath::worker_multi::Allocator *allocator):
  worker(allocator->worker), center(allocator)
{
  center.zero(-0.5, 0.0);
  step_log=7;
  step_size__=1.0/128;
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

void MandelModel::Position::assign(Position *src)
{
  if (src!=nullptr)
  {
    step_log=src->step_log;
    step_size__=src->step_size__;
  }
}

void MandelModel::Position::setView(const MandelMath::complex *c, double scale)
{
  step_log=-ilogb(scale);
  step_size__=ldexp(1.0, -step_log);
  //worker->zero(center.re); //worker->swap()
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
  int ccrm=cached_center_re_mod;
  int ccim=cached_center_im_mod;
  updateCachedDepth();
  if ((ccrm!=cached_center_re_mod) || (ccim!=cached_center_im_mod))
    dbgPoint();
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
    int ccrm=cached_center_re_mod;
    int ccim=cached_center_im_mod;
    updateCachedDepth();
    if ((ccrm!=cached_center_re_mod) || (ccim!=cached_center_im_mod))
      dbgPoint();
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



MandelModel::Orbit::Orbit(MandelMath::worker_multi::Allocator *allocator): currentWorker(allocator->worker),
  /*evaluatorAllocator(worker, MandelEvaluator::LEN),*/ evaluator(allocator, nullptr),
  pointAllocator(allocator, MandelPoint::LEN), pointDataStore(), pointData(&pointDataStore, &pointAllocator),
  lagu_c(allocator), lagu_r(allocator), tmp(allocator), bulb(allocator)
{
}

MandelModel::Orbit::~Orbit()
{
}

MandelModel::Orbit::Bulb::Bulb(MandelMath::worker_multi::Allocator *allocator):
  cb(allocator), rb(allocator), xc(allocator), baseZC(allocator), baseCC(allocator)
{
}

MandelModel::Orbit::Bulb::~Bulb()
{
}

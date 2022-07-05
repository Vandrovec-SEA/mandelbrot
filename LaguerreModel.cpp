#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "LaguerreModel.hpp"
#include "MandelEvaluator.hpp"

LaguerreModel::LaguerreModel(): QObject(),
  storeAllocator(nullptr), storeWorker(nullptr), pointStore_(nullptr),
  precisionRecord(nullptr)
{
  unsigned int oldcw; //524319 = 0x8001F = mask all interrupts, 80bit precision
  MandelMath::fpu_fix_start(&oldcw);
  _selectedPaintStyle=paintStyleCls;//Kind;
  _selectedPrecision=precisionDouble;
  epoch=0;
  imageWidth=0;
  imageHeight=0;
  //pointStore=nullptr;
  lastGivenPointIndex_=0;
  effortBonus=0;
  //threadCount=4;
  //threadCount=1;//so fast that sending messages is slower QThread::idealThreadCount()-1;
  _threadsWorking=0;
  QObject::connect(this, &LaguerreModel::selectedPrecisionChange,
                   this, &LaguerreModel::selectedPrecisionChanged);
  QObject::connect(this, &LaguerreModel::doneWork,
                   this, &LaguerreModel::giveWork1,
                   Qt::ConnectionType::QueuedConnection);
  selectedPrecisionChanged();
}

LaguerreModel::~LaguerreModel()
{
  delete precisionRecord;
  delete storeAllocator;
  delete storeWorker;
  storeWorker=nullptr;
  /*for (int i=imageWidth*imageHeight-1; i>=0; i--)
  {
    pointStore[i].cleanup(position.worker);
  }*/

  delete[] pointStore_;
  pointStore_=nullptr;
  imageWidth=0;
  imageHeight=0;
}

QString LaguerreModel::pixelXtoRE_str(int x)
{
  MandelMath::number num(precisionRecord->currentWorker);
  num.assign(precisionRecord->position.center.re);
  num.add_double((x - imageWidth/2)*precisionRecord->position.step_size__);
  QString result=precisionRecord->currentWorker->toString(num.ptr);
  return result;
}

QString LaguerreModel::pixelYtoIM_str(int y)
{
//  return (y - imageHeight/2)*position.step_size+position.center_im;
  MandelMath::number num(precisionRecord->currentWorker);
  num.assign(precisionRecord->position.center.im);
  num.add_double((y - imageHeight/2)*precisionRecord->position.step_size__);
  QString result=precisionRecord->currentWorker->toString(num.ptr);
  return result;
}

QString LaguerreModel::getTimes()
{
  QString result;
  /*for (int t=0; t<threadCount; t++)
    result+=QString("%1-%2[%3,%4],").
        arg((threads[t]->timeOuterTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t]->timeInnerTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t]->timeInvokePostTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t]->timeInvokeSwitchTotal)/1000000000.0, 0, 'f', 3);*/
  result="X-X";
  return result;
}

QString LaguerreModel::getTextXY()
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
    return "-";
  return precisionRecord->currentWorker->toString(precisionRecord->orbit.evaluator.currentParams.c.re)+" +i* "+
         precisionRecord->currentWorker->toString(precisionRecord->orbit.evaluator.currentParams.c.im);
}

QString doubleToString(double x)
{
  if (x>0)
    return QString("+%1").arg(x, 10, 'f');
  else
    return QString("%1").arg(x, 10, 'f');//returns like 50 digits QString::number(x, 10, 'f');
}

QString LaguerreModel::getTextInfoGen()
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number tmp(precisionRecord->currentWorker);
    tmp.assign(precisionRecord->orbit.evaluator.currentParams.c.re);
    tmp.sub(precisionRecord->position.center.re);
    orbit_x=qRound(tmp.toDouble()/precisionRecord->position.step_size__)+imageWidth/2;

    tmp.assign(precisionRecord->orbit.evaluator.currentParams.c.im);
    tmp.sub(precisionRecord->position.center.im);
    orbit_y=imageHeight/2-qRound(tmp.toDouble()/precisionRecord->position.step_size__);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
  //LaguerrePoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  precisionRecord->wtiPoint.store=&pointStore_[orbit_x+imageWidth*orbit_y];
  int wtiIndexFirst, wtiIndexLast;
  precisionRecord->wtiPoint.self_allocator._getRange(wtiIndexFirst, wtiIndexLast);
  precisionRecord->currentWorker->assign_block(wtiIndexFirst, storeWorker, (orbit_x+imageWidth*orbit_y)*LaguerrePoint::LEN, LaguerrePoint::LEN);
  LaguerrePoint *data=&precisionRecord->wtiPoint;

  QString state;
  switch (data->store->state)
  {
    case LaguerrePointStore::State::stUnknown:
      state="Unk"; break;
    case LaguerrePointStore::State::stResolved:
      state="OK"; break;
    case LaguerrePointStore::State::stFail:
      state="Err"; break;
  }

  return "Per="+QString::number(precisionRecord->params.period)+" "+state+" iter="+QString::number(data->store->iter)+
        " mu="+QString::number(precisionRecord->orbit.first_mu_re_, 'f', 3)+","+QString::number(precisionRecord->orbit.first_mu_im, 'f', 3)+
        " mum="+QString::number(precisionRecord->orbit.first_mum_re_, 'f', 3)+","+QString::number(precisionRecord->orbit.first_mum_im_, 'f', 3);
}

QString LaguerreModel::getTextInfoSpec()
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number tmp(precisionRecord->currentWorker);
    reimToPixel(&orbit_x, &orbit_y, &precisionRecord->orbit.evaluator.currentParams.c, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
  //LaguerrePoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  precisionRecord->wtiPoint.store=&pointStore_[orbit_x+imageWidth*orbit_y];
  int wtiIndexFirst, wtiIndexLast;
  precisionRecord->wtiPoint.self_allocator._getRange(wtiIndexFirst, wtiIndexLast);
  precisionRecord->currentWorker->assign_block(wtiIndexFirst, storeWorker, (orbit_x+imageWidth*orbit_y)*LaguerrePoint::LEN, LaguerrePoint::LEN);
  LaguerrePoint *data=&precisionRecord->wtiPoint;

  switch (data->store->state)
  {
    case LaguerrePointStore::State::stUnknown:
      return " ";
      break;
    case LaguerrePointStore::State::stResolved:
    {
      //MandelMath::complex fz(orbit.worker, &orbit.pointData.fz_r_re, &orbit.pointData.fz_r_im, true);
      QString naiveChoice;
      switch (precisionRecord->orbit.pointData.store->naiveChoice)
      {
        case NewtonNaiveChoice::nc03: naiveChoice="0.3"; break;
        case NewtonNaiveChoice::nc05: naiveChoice="0.5"; break;
        case NewtonNaiveChoice::nc08: naiveChoice="0.8"; break;
        case NewtonNaiveChoice::ncWide: naiveChoice="Wide"; break;
        case NewtonNaiveChoice::nc100: naiveChoice="100"; break;
        case NewtonNaiveChoice::nc90_: naiveChoice="90"; break;
        case NewtonNaiveChoice::nc80: naiveChoice="80"; break;
        case NewtonNaiveChoice::nc60: naiveChoice="60"; break;
        case NewtonNaiveChoice::ncClose: naiveChoice="Close"; break;
      }
      return QString("attr=")+QString::number(precisionRecord->orbit.pointData.fz_r.getMag_double())+
             QString(" firstM=")+QString::number(precisionRecord->orbit.pointData.store->firstM)+
             QString(" NChoice=")+naiveChoice;
    } break;
    case LaguerrePointStore::State::stFail:
      return " ";
  }
  return "-?-?-";
}

void LaguerreModel::setParams(ShareableViewInfo viewInfo)
{
  //position.setNumberType(viewInfo.worker->ntype());
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
    dbgPoint();
  precisionRecord->params.period=viewInfo.period;

  if (viewInfo.originalAllocator->worker->ntype()==precisionRecord->currentWorker->ntype())
  {
    precisionRecord->params.base.assign(&viewInfo.c);
    precisionRecord->params.root.assign(&viewInfo.root);
  }
  else
  {
    //watch this
    dbgPoint(); //probably happens quite rarely, let's debug it when it does
    MandelMath::worker_multi::Type ntype=precisionRecord->currentWorker->ntype();
    MandelMath::worker_multi *tmpworker=MandelMath::worker_multi::create(ntype, viewInfo.originalAllocator);
    {
      MandelMath::worker_multi::Allocator useTmp(tmpworker->getAllocator(), ShareableViewInfo::LEN);
      ShareableViewInfo tmpinfo(&useTmp);
      precisionRecord->params.base.assign(&tmpinfo.c);
      precisionRecord->params.root.assign(&tmpinfo.root);
    }
    delete tmpworker;
  }
  MandelMath::complex old_c(precisionRecord->currentWorker);
  old_c.assign(&precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  precisionRecord->position.setView(&precisionRecord->params.base, viewInfo.scale);

  transformStore(precisionRecord->currentWorker, storeWorker, pointStore_, 0, 0, &old_c,
                 precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &precisionRecord->position.center,
                 precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

  startNewEpoch();
}

void LaguerreModel::transformStore(MandelMath::worker_multi *old_worker, MandelMath::worker_multi *old_sworker, LaguerrePointStore *old_store, int old_width, int old_height, const MandelMath::complex *old_c,
                                   MandelMath::worker_multi *new_worker, MandelMath::worker_multi *new_sworker, LaguerrePointStore *new_store, int new_width, int new_height, const MandelMath::complex *new_c,
                                   int inlog, int new_step_log)
{
  if (precisionRecord==nullptr || precisionRecord->currentWorker==nullptr)
  {
    dbgPoint();
    return;
  };
  int indexOfWtiPoint, wtiIndexLast;
  precisionRecord->wtiPoint.self_allocator._getRange(indexOfWtiPoint, wtiIndexLast);
  assert(wtiIndexLast-indexOfWtiPoint==LaguerrePoint::LEN);
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
        {
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint::LEN, old_sworker, (oldy*old_width+oldx)*LaguerrePoint::LEN, LaguerrePoint::LEN);
        }
        else
        {
          new_worker->assign(c.re, precisionRecord->position.center.re);
          new_worker->add_double(c.re, (newx - new_width/2)*precisionRecord->position.step_size__);
          //MandelMath::worker_multi::Allocator allo(new_sworker->getAllocator(), (newy*new_width+newx)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
          //LaguerrePoint pt(&new_store[newy*new_width+newx], &allo);
          //pt.zero(&c);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(&c);
          new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint::LEN, precisionRecord->currentWorker, indexOfWtiPoint, LaguerrePoint::LEN);
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
        {
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint::LEN, old_sworker, (oldy*old_width+oldx)*LaguerrePoint::LEN, LaguerrePoint::LEN);
        }
        else
        {
          new_worker->assign(c.re, precisionRecord->position.center.re);
          new_worker->add_double(c.re, (newx - new_width/2)*precisionRecord->position.step_size__);
          //MandelMath::worker_multi::Allocator allo(new_sworker->getAllocator(), (newy*new_width+newx)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
          //LaguerrePoint pt(&new_store[newy*new_width+newx], &allo);
          //pt.zero(&c);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(&c);
          new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint::LEN, precisionRecord->currentWorker, indexOfWtiPoint, LaguerrePoint::LEN);
        }
      }
    }
  }
}

void LaguerreModel::setView_(const MandelMath::complex *c, double scale)
{
  MandelMath::complex old_c(precisionRecord->currentWorker);
  old_c.assign(&precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  precisionRecord->position.setView(c, scale);

  transformStore(precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &precisionRecord->position.center,
                 precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

  startNewEpoch();
}

void LaguerreModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::complex old_c(precisionRecord->currentWorker);
  old_c.assign(&precisionRecord->position.center);

  int dx=qRound(delta_x);
  int dy=qRound(delta_y);
  precisionRecord->position.move(dx, dy);
  //qDebug()<<"new c: re="<<position.worker->toString(&position.center_re_s)<<",im="<<position.worker->toString(&position.center_im_s);

  transformStore(precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &precisionRecord->position.center,
                 0, precisionRecord->position.step_log);

  startNewEpoch();
}

void LaguerreModel::zoom(double x, double y, int inlog)
{
  MandelMath::complex old_c(precisionRecord->currentWorker);
  old_c.assign(&precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  precisionRecord->position.scale(inlog, qRound(x)-imageWidth/2, imageHeight/2-qRound(y));

  transformStore(precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &precisionRecord->position.center,
                 precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

  startNewEpoch();
}

void LaguerreModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  int newLength=width*height;
  LaguerrePointStore *newStore_=new LaguerrePointStore[newLength];
  {
    QString size_as_text=QString::number(newLength*sizeof(LaguerrePointStore));
    for (int pos=size_as_text.length()-3; pos>0; pos-=3)
      size_as_text.insert(pos, '\'');
    qDebug()<<"laguerreStore uses"<<size_as_text.toLocal8Bit().constData()<<"B";
  }
  MandelMath::worker_multi *newStoreWorker;
  newStoreWorker=MandelMath::worker_multi::create(storeWorker->ntype(), newLength*LaguerrePoint::LEN);
  MandelMath::worker_multi::Allocator *newStoreAllocator;
  newStoreAllocator=new MandelMath::worker_multi::Allocator(newStoreWorker->getAllocator(), newLength*LaguerrePoint::LEN);


  MandelMath::complex old_c(precisionRecord->currentWorker);
  old_c.assign(&precisionRecord->position.center);

  transformStore(precisionRecord->currentWorker, storeWorker, pointStore_, imageWidth, imageHeight, &old_c,
                 precisionRecord->currentWorker, newStoreWorker, newStore_, width, height, &precisionRecord->position.center,
                 0, precisionRecord->position.step_log);

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

void LaguerreModel::startNewEpoch()
{
  epoch=(epoch+1)%2000000000;
  lastGivenPointIndex_=0;
  effortBonus=0;
  giveWorkAll();
}

void LaguerreModel::giveWorkAll()
{
  for (int t=0; t<precisionRecord->threadCount; t++)
    if (precisionRecord->threads[t]->currentParams.pixelIndex<0)
      giveWork(precisionRecord->threads[t]);
}

void LaguerreModel::reimToPixel(int *circ_x, int *circ_y, const MandelMath::complex *c, MandelMath::number *tmp)
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

void LaguerreModel::paintOrbit(ShareableImageWrapper image, int x, int y)
{
  if ((x<0) || (x>=imageWidth) || (y<0) || (y>=imageHeight))
    return;
  if (precisionRecord==nullptr)
    dbgPoint();
  else if (precisionRecord->currentWorker==nullptr)
    dbgPoint();
  else if (precisionRecord->orbit.currentWorker->ntype()!=precisionRecord->currentWorker->ntype())
    dbgPoint();
  if (precisionRecord->params.period<=0)
    return;

  QImage newOverlay(imageWidth, imageHeight, QImage::Format::Format_ARGB32);
  QPainter painter(&newOverlay);
  //QRgb what=newOverlay.pixel(0, 0);
  //if (what!=0) //0xcdcdcdcd in MSVC compiler
  {
    painter.setCompositionMode(QPainter::CompositionMode::CompositionMode_Source);
    painter.fillRect(0, 0, imageWidth, imageHeight, Qt::GlobalColor::transparent);
  };

  //LaguerrePoint *data=&pointStore[y*imageWidth+x];
  /*switch (data_->state)
  {
    case LaguerrePoint::State::stUnknown:
    case LaguerrePoint::State::stResolved:
    case LaguerrePoint::State::stFail:
    default: ;
  }*/
  //orbit.evaluator.switchType(position.worker);
  MandelMath::number tmp(precisionRecord->currentWorker);//=&orbit.evaluator.currentData.f_re;
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0, 0)); //paint c
    reimToPixel(&circ_x, &circ_y, &precisionRecord->params.base, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l3[3]={{circ_x+1, circ_y-1, circ_x-1, circ_y-1},
                   {circ_x-1, circ_y-1, circ_x-1, circ_y+1},
                   {circ_x-1, circ_y+1, circ_x+1, circ_y+1}};
      painter.drawLines(l3, 3);
    };

    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0, 0)); //paint root as +
    reimToPixel(&circ_x, &circ_y, &precisionRecord->params.root, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };
  }


  //evaluator.params.c .. saved mouse coordinates
  //evaluator.data.root .. found root (temporary)
  //pointData.f .. not used
  //pointData.fz .. result f'
  precisionRecord->position.pixelXtoRE(x-imageWidth/2, precisionRecord->orbit.evaluator.currentParams.c.re); //remember for text labels
  precisionRecord->position.pixelYtoIM(imageHeight/2-y, precisionRecord->orbit.evaluator.currentParams.c.im);
  precisionRecord->orbit.evaluator.currentParams.epoch=epoch;
  precisionRecord->orbit.evaluator.currentParams.pixelIndex=0;
  precisionRecord->orbit.pointData.zero(&precisionRecord->orbit.evaluator.currentParams.c);
  /*MandelMath::complex base(orbit.worker, &params.base_re_s_, &params.base_im_s, true);
  currentWorker->assign(&orbit.evaluator.currentData.root_re, &orbit.evaluator.currentParams.c_re);
  orbit.worker->assign(&orbit.evaluator.currentData.root_im, &orbit.evaluator.currentParams.c_im);
  MandelMath::complex root(orbit.worker, &orbit.evaluator.currentData.root_re, &orbit.evaluator.currentData.root_im, true);*/
  precisionRecord->orbit.evaluator.currentData.root.assign(&precisionRecord->orbit.evaluator.currentParams.c);
  precisionRecord->orbit.evaluator.newton(precisionRecord->params.period, &precisionRecord->params.base, &precisionRecord->orbit.evaluator.currentData.root, true, 12);
  //orbit.pointData.assign(orbit.worker, orbit.evaluator.currentData);
  precisionRecord->orbit.pointData.store->iter=precisionRecord->orbit.evaluator.newtres_.cyclesNeeded;
  precisionRecord->orbit.pointData.store->firstM=precisionRecord->orbit.evaluator.newtres_.firstMum_re_;
  precisionRecord->orbit.pointData.fz_r.assign(&precisionRecord->orbit.evaluator.newtres_.fz_r);
  precisionRecord->currentWorker->add_double(precisionRecord->orbit.pointData.fz_r.re, 1);
  precisionRecord->orbit.pointData.store->naiveChoice=precisionRecord->orbit.evaluator.newtres_.naiveChoice;
  precisionRecord->orbit.first_mu_re_=precisionRecord->orbit.evaluator.newtres_.firstMu_re_;
  precisionRecord->orbit.first_mu_im=precisionRecord->orbit.evaluator.newtres_.firstMu_im;
  precisionRecord->orbit.first_mum_re_=precisionRecord->orbit.evaluator.newtres_.firstMum_re_;
  precisionRecord->orbit.first_mum_im_=precisionRecord->orbit.evaluator.newtres_.firstMum_im_;

  int circ_x, circ_y;
  double tmp_re, tmp_im;
  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0xff, 0)); //final=yellow /
  reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.currentData.root, &tmp);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-2, circ_y+2, circ_x+2, circ_y-2);
    //size of r that maps to 1 epsilon: about the width of transition from >0 to <0
    int circ_r=ldexp(sqrt(precisionRecord->orbit.evaluator.newtres_.accy_tostop*precisionRecord->currentWorker->eps2()), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
      painter.drawEllipse(circ_x-circ_r, circ_y-circ_r, 2*circ_r, 2*circ_r);
    //size of r that maps to a few epsilon: the accuracy of root and also size of the dead pool
    circ_r=ldexp(sqrt((precisionRecord->orbit.evaluator.newtres_.accy_tostop*precisionRecord->orbit.evaluator.newtres_.accy_multiplier)*precisionRecord->currentWorker->eps2()), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
      painter.drawEllipse(circ_x-circ_r, circ_y-circ_r, 2*circ_r, 2*circ_r);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0, 0xff)); //Newton=purple \      .
  reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.newtres_.first_guess_newt, &tmp);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-2, circ_y-2, circ_x+2, circ_y+2);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush); //white = is better visible than cyan | so paint it under
  painter.setPen(QColor(0xff, 0xff, 0xff)); //Fejer=white =
  tmp_re=precisionRecord->orbit.evaluator.newtres_.first_fejer_re-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.re);
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres_.first_fejer_im-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.im);
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-3, circ_y-1, circ_x+3, circ_y-1);
    painter.drawLine(circ_x-3, circ_y+1, circ_x+3, circ_y+1);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush); //
  painter.setPen(QColor(0xff, 0x00, 0x00)); //Naive=red o
  tmp_re=precisionRecord->orbit.evaluator.newtres_.first_naive1_re_-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.re);
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres_.first_naive1_im-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.im);
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
  };
  tmp_re=precisionRecord->orbit.evaluator.newtres_.first_naive2_re-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.re);
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres_.first_naive2_im-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.im);
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };
  tmp_re=precisionRecord->orbit.evaluator.newtres_.first_naive_re-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.re);
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres_.first_naive_im-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.im);
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre=cyan |
  reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.newtres_.first_guess_lagu, &tmp);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x, circ_y-3, circ_x, circ_y+3);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre1=cyan ||
  tmp_re=precisionRecord->orbit.evaluator.newtres_.first_lagu1_re-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.re);
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres_.first_lagu1_im-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.im);
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-1, circ_y-3, circ_x-1, circ_y+3);
    painter.drawLine(circ_x+1, circ_y-3, circ_x+1, circ_y+3);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre1other=cyan o
  tmp_re=precisionRecord->orbit.evaluator.newtres_.first_lagu1o_re-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.re);
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres_.first_lagu1o_im-precisionRecord->currentWorker->toDouble(precisionRecord->position.center.im);
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };


  painter.setBrush(Qt::BrushStyle::NoBrush);
  if (precisionRecord->orbit.evaluator.newtres_.first_neumaier1_im_!=0)
  {
    painter.setPen(QColor(0xff, 0x00, 0x00)); //Neumaier1_im
    int circ_r=ldexp(abs(precisionRecord->orbit.evaluator.newtres_.first_neumaier1_im_), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  };
  if ((precisionRecord->orbit.evaluator.newtres_.first_neumaier1_im_!=0) || (precisionRecord->orbit.evaluator.newtres_.first_neumaier1_re_<0))
  {
    painter.setPen(QColor(0xff, 0x00, 0xff)); //Neumaier1_mag
    double mag=sqrt(precisionRecord->orbit.evaluator.newtres_.first_neumaier1_re_*precisionRecord->orbit.evaluator.newtres_.first_neumaier1_re_+
                    precisionRecord->orbit.evaluator.newtres_.first_neumaier1_im_*precisionRecord->orbit.evaluator.newtres_.first_neumaier1_im_);
    int circ_r=ldexp(mag, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  }
  else
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Neumaier1
    int circ_r=ldexp(precisionRecord->orbit.evaluator.newtres_.first_neumaier1_re_, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  }

  painter.setBrush(Qt::BrushStyle::NoBrush);
  if (precisionRecord->orbit.evaluator.newtres_.first_neumaier2_im!=0)
  {
    painter.setPen(QColor(0xff, 0x00, 0x00)); //Neumaier2_im
    int circ_r=ldexp(abs(precisionRecord->orbit.evaluator.newtres_.first_neumaier2_im), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  };
  if ((precisionRecord->orbit.evaluator.newtres_.first_neumaier2_im!=0) || (precisionRecord->orbit.evaluator.newtres_.first_neumaier2_re<0))
  {
    painter.setPen(QColor(0xff, 0x00, 0xff)); //Neumaier2_mag
    double mag=sqrt(precisionRecord->orbit.evaluator.newtres_.first_neumaier2_re*precisionRecord->orbit.evaluator.newtres_.first_neumaier2_re+
                    precisionRecord->orbit.evaluator.newtres_.first_neumaier2_im*precisionRecord->orbit.evaluator.newtres_.first_neumaier2_im);
    int circ_r=ldexp(mag, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }
  else
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Neumaier2
    int circ_r=ldexp(precisionRecord->orbit.evaluator.newtres_.first_neumaier2_re, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }


  image.image->swap(newOverlay);
}

int LaguerreModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return -1;
  timerWriteToImage.start();
  int indexOfWtiPoint, _discard;
  precisionRecord->wtiPoint.self_allocator._getRange(indexOfWtiPoint, _discard);
  LaguerrePointStore *wtiStore;
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      bool knownenum=false;
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
      //LaguerrePoint data_(&pointStore_[y*imageWidth+x], &allo);
      wtiStore=&pointStore_[y*imageWidth+x];
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleCls:
        {
          switch (wtiStore->state)
          {
            case LaguerrePointStore::State::stUnknown:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
              break;
            case LaguerrePointStore::State::stResolved:
            {
              if (wtiStore->iter==0)
                image.image->setPixel(x, y, 0xff606060); //the dead pool
              else
              {
                int r;
                {
                  precisionRecord->currentWorker->assign_block(indexOfWtiPoint, storeWorker, (y*imageWidth+x)*LaguerrePoint::LEN, LaguerrePoint::LEN);
                  double tr=precisionRecord->currentWorker->toDouble(precisionRecord->wtiPoint.fz_r.re);
                  double ti=precisionRecord->currentWorker->toDouble(precisionRecord->wtiPoint.fz_r.im);
                  if (tr*tr+ti*ti>1)
                    r=0xc0;
                  else
                    r=0x00;
                }
                int g;
                /*if (data->firstM>=3)
                  g=0xff;
                else if (data->firstM>=2)
                  g=0xbf+0x40*(data->firstM-2);
                else if (data->firstM>=1)
                  g=0x7f+0x40*(data->firstM-1);
                else if (data->firstM>=0)
                  g=0x3f+0x40*(data->firstM-0);
                else if (data->firstM>=-1)
                  g=0x00+0x40*(data->firstM+1);
                else
                  g=0xdf;*/
                if (wtiStore->firstM>=3)
                  g=0x7f;
                else if (wtiStore->firstM>=2)
                  g=0xbf-0x40*(wtiStore->firstM-2);
                else if (wtiStore->firstM>=1)
                  g=0xff-0x40*(wtiStore->firstM-1);
                else if (wtiStore->firstM>=0)
                  g=0x00+0x40*(wtiStore->firstM-0);
                else if (wtiStore->firstM>=-1)
                  g=0x40+0x40*(wtiStore->firstM+1);
                else
                  g=0x80;
                int b;
                switch (wtiStore->iter % 5)
                {
                  case  0: b=0x00; break;
                  case  1: b=0xff; break;
                  case  2: b=0xc0; break;
                  case  3: b=0x80; break;
                  case  4: b=0x40; break;
                  default: b=0xff;
                }
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
              knownenum=true;
            } break;
            case LaguerrePointStore::State::stFail:
            {
              image.image->setPixel(x, y, 0xffffffff);
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

void LaguerreModel::giveWork(MandelEvaluator *evaluator)
{
  int retryEffortFrom=0;
  while (retryEffortFrom>=0)
  {
    retryEffortFrom=-1;
    int quickrun=0;
    MandelMath::complex &tmpc=precisionRecord->params.base;//will assign to currentParams.c (position.worker, &params.base_re_s_, &params.base_im_s, true);
    MandelMath::complex &root=evaluator->currentData.f;//(position.worker, &evaluator->currentData.f_re, &evaluator->currentData.f_im, true);
    for (int pi=0; pi<imageWidth*imageHeight; pi++)
    {
      int pointIndex=(lastGivenPointIndex_+pi)%(imageWidth*imageHeight);
      //if ((lastGivenPointIndex_!=0) && (pointIndex==0))
        //dbgPoint();
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), pointIndex*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
      //LaguerrePoint pointData_(&pointStore_[pointIndex], &allo);
      LaguerrePointStore *storeAtIndex=&pointStore_[pointIndex];
      bool needsEval=(storeAtIndex->state==LaguerrePointStore::State::stUnknown);
      if (needsEval)
      {
        bool found=false;
        for (int t=0; t<precisionRecord->threadCount; t++)
          if (precisionRecord->threads[t]->currentParams.pixelIndex==pointIndex)
          {
            found=true;
            break;
          };
        if (evaluator->currentParams.pixelIndex!=-1)
          dbgPoint();
        assert(evaluator->currentParams.pixelIndex==-1);
        if (!found)
        {
          int phasex=(pointIndex%imageWidth-imageWidth/2+precisionRecord->position.cached_center_re_mod+32768)%32768;
          int phasey=(pointIndex/imageWidth-imageHeight/2+precisionRecord->position.cached_center_im_mod+32768)%32768;
          //int effort=ctz16(phasex)+ctz16(phasey);
          int effort=MandelMath::ctz16(phasex | phasey);
          if (effort>8)
            effort=8;
          effort+=effortBonus;
          if (effort<8)
          {
            if (retryEffortFrom<0)
              retryEffortFrom=pointIndex;
          }
          else
          {
            if (evaluator->currentWorker->ntype()!=precisionRecord->position.worker->ntype())
              dbgPoint();
            //evaluator->switchType(position.worker);
            precisionRecord->position.pixelXtoRE(pointIndex%imageWidth - imageWidth/2, root.re);
            precisionRecord->position.pixelYtoIM(imageHeight/2-pointIndex/imageWidth, root.im);
            evaluator->currentParams.epoch=epoch;
            evaluator->currentParams.pixelIndex=pointIndex;
            if (precisionRecord->params.period>10)
            {
              evaluator->startNewton(precisionRecord->params.period, &tmpc); //evaluator->currentData.f is additional hidden parameter
              _threadsWorking++;
              evaluator->timeOuter.start();
              lastGivenPointIndex_=pointIndex;
              return;
            }
            else
            {
              int result=evaluator->newton(precisionRecord->params.period, &tmpc, &root, true, 12);
              donePixel1(evaluator, result);
              quickrun++;
              if (quickrun>=100)
              {
                lastGivenPointIndex_=pointIndex;
                /*
                can't figure out how to invokeMethod with arguments
                QMetaObject::invokeMethod(this,
                                          &LaguerreModel::giveWorkAll,
                                          Qt::ConnectionType::QueuedConnection);*/
                emit doneWork(evaluator);
                return;
              };
            }
          }
        }
      }
      //MandelEvaluator::simple(cr, ci, pointStore[y*imageWidth+x]);
    }
    if ((retryEffortFrom>=0))// && (effortBonus<MAX_EFFORT))
    {
      effortBonus++;
      lastGivenPointIndex_=retryEffortFrom;
    }
    else
      retryEffortFrom=-1;
  }
}

void LaguerreModel::donePixel1(MandelEvaluator *me, int result)
{
  if ((me->currentParams.epoch==epoch) && (me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
  {
    //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), me->currentParams.pixelIndex*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
    //LaguerrePoint point_(&pointStore_[me->currentParams.pixelIndex], &allo);
    /*
    if (point->state!=LaguerrePoint::State::stUnknown)
      qDebug()<<"Finished pixel finished again";
    else*/
    {
      if (me->currentWorker==nullptr)
        dbgPoint();
      else
      {
        LaguerrePointStore *dstStore=&pointStore_[me->currentParams.pixelIndex];
        if (result>0)
        {
          dstStore->state=LaguerrePointStore::State::stResolved;
        }
        else
        {
          me->newtres_.cyclesNeeded=-1;
          dstStore->state=LaguerrePointStore::State::stFail;
        };
        //point->assign(position.worker, me->currentData);
        //in theory, LaguerrePoint can be assigned from NewtRes; in practice, they are quite different
        //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), me->currentParams.pixelIndex*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
        //LaguerrePoint point_(dstStore, &allo, nullptr); mhmm
        me->tmpLaguerrePoint.f.assign(&me->currentData.f); //root
        me->tmpLaguerrePoint.fz_r.assign(&me->newtres_.fz_r);
        me->currentWorker->add_double(me->tmpLaguerrePoint.fz_r.re, 1);
        int first, last;
        me->tmpLaguerrePoint.self_allocator._getRange(first, last);
        storeWorker->assign_block(me->currentParams.pixelIndex*LaguerrePoint::LEN, me->currentWorker, first, last-first);
        dstStore->iter=me->newtres_.cyclesNeeded;
        dstStore->naiveChoice=me->newtres_.naiveChoice;
        dstStore->firstM=me->newtres_.firstMum_re_;
      }
    }
  }
  else if (me->currentParams.epoch!=epoch)
  { }//qDebug()<<"Old pixel finished";
  else
    qWarning()<<"Invalid pixel finished";
  me->currentParams.pixelIndex=-1;
}

void LaguerreModel::donePixel(MandelEvaluator *me, int result)
{
  me->timeOuterTotal+=me->timeOuter.nsecsElapsed();
  _threadsWorking--;
  donePixel1(me, result);
  giveWork(me);
}

void LaguerreModel::selectedPrecisionChanged()
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
    newWorker=MandelMath::worker_multi::create(newPrecision, LaguerreModel::LEN);
    newStoreWorker=MandelMath::worker_multi::create(newPrecision, pointCount*LaguerrePoint::LEN);
  }

  {
    PrecisionRecord *newprecisionRecord=new PrecisionRecord(newWorker, precisionRecord, this);
    delete precisionRecord;
    precisionRecord=newprecisionRecord;
  }

  MandelMath::worker_multi::Allocator *newStoreAllocator=new MandelMath::worker_multi::Allocator(newStoreWorker->getAllocator(), pointCount*LaguerrePoint::LEN);
  delete storeAllocator;
  delete storeWorker;
  storeWorker=newStoreWorker;
  storeAllocator=newStoreAllocator;

  startNewEpoch();
}

LaguerreModel::Params::Params(MandelMath::worker_multi::Allocator *allocator):
  period(1), base(allocator), root(allocator)
{

}

void LaguerreModel::Params::assign(Params *src)
{
  if (src)
  {
    period=src->period;
  }
  else
  {
    base.zero(0, 0);
    root.zero(0, 0);
  }
}

LaguerreModel::Position::Position(MandelMath::worker_multi::Allocator *allocator, Position *src):
  worker(allocator->worker), center(allocator)
{
  if (src)
  {
    //keep preinitialized values center.zero(-0.5, 0.0);
    step_log=src->step_log;
    step_size__=src->step_size__;
  }
  else
  {
    center.zero(-0.5, 0.0);
    step_log=7;
    step_size__=1.0/128;
  }
  updateCachedDepth();
}

LaguerreModel::Position::~Position()
{
  /*if (worker!=nullptr)
  {
    worker->cleanup(&center_re_s);
    worker->cleanup(&center_im_s);
  };*/
}

void LaguerreModel::Position::assign(Position *src)
{
  if (src!=nullptr)
  {
    step_log=src->step_log;
    step_size__=src->step_size__;
  }
}

void LaguerreModel::Position::setView(const MandelMath::complex *c, double scale)
{
  step_log=-ilogb(scale);
  step_size__=ldexp(1.0, -step_log);
  //worker->zero(&center_re_s, ldexp(round(ldexp(c_re, step_log)), -step_log));
  //worker->zero(&center_re_s); //worker->swap()
  center.assign(c);
  center.lshift(step_log);
  worker->round(center.re);
  worker->round(center.im);
  center.lshift(-step_log);
  //worker->zero(&center_im_s, ldexp(round(ldexp(c_im, step_log)), -step_log));
  //worker->zero(&center_im_s); //worker->swap()
  updateCachedDepth();
}

void LaguerreModel::Position::move(int delta_x, int delta_y)
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

void LaguerreModel::Position::scale(int inlog, int center_x, int center_y)
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

void LaguerreModel::Position::updateCachedDepth()
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

void LaguerreModel::Position::pixelXtoRE(int x, MandelMath::number_pointer result)
{
  //return (x - imageWidth/2)*position.step_size+position.center_re;
  //should be already result->reinit(center_re_n.ntype());
  worker->assign(result, center.re);
  worker->add_double(result, x*step_size__);
}

void LaguerreModel::Position::pixelYtoIM(int y, MandelMath::number_pointer result)
{
  //return (y - imageHeight/2)*position.step_size+position.center_im;
  //should be already result->reinit(center_im_n.ntype());
  worker->assign(result, center.im);
  worker->add_double(result, y*step_size__);
}



LaguerreModel::Orbit::Orbit(MandelMath::worker_multi::Allocator *allocator): currentWorker(allocator->worker),
  /*evaluatorAllocator(worker, MandelEvaluator::LEN),*/ evaluator(allocator->worker->ntype(), true),
  /*pointAllocator(allocator, LaguerrePoint::LEN),*/ pointDataStore(), pointData(&pointDataStore, allocator),
  first_mu_re_(0), first_mu_im(0), first_mum_re_(0), first_mum_im_(0)
{

}

LaguerreModel::Orbit::~Orbit()
{

}


LaguerreModel::PrecisionRecord::PrecisionRecord(MandelMath::worker_multi *newWorker, PrecisionRecord *source, LaguerreModel *doneReceiver):
  currentWorker(newWorker),
  paramsAllocator(newWorker->getAllocator(), Params::LEN), params(&paramsAllocator),
  wtiPoint(nullptr, newWorker->getAllocator()),
  position(newWorker->getAllocator(), source?&source->position:nullptr), orbit(newWorker->getAllocator()),
  threadCount(source?source->threadCount:0), threads(nullptr)
{
  if (source)
  {
    params.assign(&source->params);
    position.assign(&source->position);
  }
  else
  {
    params.assign(nullptr);
    position.assign(nullptr);
  }

  if (threadCount<=0)
  {
    threadCount=QThread::idealThreadCount()-1;
    if (threadCount<1)
      threadCount=1;
  };
  threads=new MandelEvaluator *[threadCount];
  for (int t=0; t<threadCount; t++)
  {
    threads[t]=new MandelEvaluator(newWorker->ntype(), false);
    //threads[t].setHint(t);
    QObject::connect(threads[t], &MandelEvaluator::doneNewton,
                     doneReceiver, &LaguerreModel::donePixel,
                     Qt::ConnectionType::QueuedConnection);
  }
}

LaguerreModel::PrecisionRecord::~PrecisionRecord()
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
  threadCount=0;
}

#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "MandelModel.hpp"
#include "MandelEvaluator.hpp"

MandelModel::MandelModel(): QObject(), position()
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
  MandelMath::number_worker_ddouble wdd;
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
  dd1.split(0x20000010); //-> 0x
  dd1.split(0x20000011); //-> 0x20000000 + 0x1
  dd1.split(-0.00043541188577694845);

  _selectedPaintStyle=paintStyleCls;//Kind;
  _selectedPrecision=precisionDouble;
  epoch=0;
  imageWidth=0;
  imageHeight=0;
  pointStore=nullptr;
  lastGivenPointIndex_=0;
  effortBonus_=0;
  orbit.worker=nullptr;
  //threadCount=4;
  _threadsWorking=0;
  threadCount=QThread::idealThreadCount()-1;
  if (threadCount<1)
    threadCount=1;
  threads=new MandelEvaluator[threadCount];
  for (int t=0; t<threadCount; t++)
  {
    //threads[t].setHint(t);
    QObject::connect(&threads[t], &MandelEvaluator::doneCompute,
                     this, &MandelModel::donePixel,
                     Qt::ConnectionType::QueuedConnection);
  }
  QObject::connect(this, &MandelModel::selectedPrecisionChange,
                   this, &MandelModel::selectedPrecisionChanged);
}

MandelModel::~MandelModel()
{
  for (int t=0; t<threadCount; t++)
  {
    threads[t].wantStop=true;
    threads[t].quit();
    threads[t].wait(1000);
  }
  delete[] threads;
  threadCount=0;
  threads=nullptr;
  //orbit.evaluator.switchType(nullptr);
  orbit.evaluator.wantStop=true;
  orbit.evaluator.quit();
  orbit.evaluator.wait(1000);
  if (orbit.worker!=nullptr)
  {
    orbit.pointData.cleanup(orbit.worker);
    orbit.worker->cleanup(&orbit.tmp);
    orbit.worker->cleanup(&orbit.lagu_r_im);
    orbit.worker->cleanup(&orbit.lagu_r_re);
    orbit.worker->cleanup(&orbit.lagu_c_im);
    orbit.worker->cleanup(&orbit.lagu_c_re_);

    orbit.worker->cleanup(&orbit.bulb.cb_re);
    orbit.worker->cleanup(&orbit.bulb.cb_im);
    orbit.worker->cleanup(&orbit.bulb.rb_re_);
    orbit.worker->cleanup(&orbit.bulb.rb_im);
    orbit.worker->cleanup(&orbit.bulb.xc_re);
    orbit.worker->cleanup(&orbit.bulb.xc_im);
    orbit.worker->cleanup(&orbit.bulb.baseZC_re);
    orbit.worker->cleanup(&orbit.bulb.baseZC_im);
    orbit.worker->cleanup(&orbit.bulb.baseCC_re);
    orbit.worker->cleanup(&orbit.bulb.baseCC_im);
  };

  for (int i=imageWidth*imageHeight-1; i>=0; i--)
  {
    pointStore[i].cleanup(position.worker);
  }

  delete[] pointStore;
  pointStore=nullptr;
  imageWidth=0;
  imageHeight=0;
}

QString MandelModel::pixelXtoRE_str(int x)
{
  MandelMath::number_store num;
  position.worker->assign(&num, &position.center_re_s);
  position.worker->add_double(&num, (x - imageWidth/2)*position.step_size__);
  QString result=position.worker->toString(&num);
  position.worker->cleanup(&num);
  return result;
}

QString MandelModel::pixelYtoIM_str(int y)
{
//  return (y - imageHeight/2)*position.step_size+position.center_im;
  MandelMath::number_store num;
  position.worker->assign(&num, &position.center_im_s);
  position.worker->add_double(&num, (y - imageHeight/2)*position.step_size__);
  QString result=position.worker->toString(&num);
  position.worker->cleanup(&num);
  return result;
}

QString MandelModel::getTimes()
{
  QString result;
  for (int t=0; t<threadCount; t++)
    result+=QString("%1-%2[%3,%4],").
        arg((threads[t].timeOuterTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t].timeInnerTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t].timeInvokePostTotal)/1000000000.0, 0, 'f', 3).
        arg((threads[t].timeInvokeSwitchTotal)/1000000000.0, 0, 'f', 3);
  return result;
}

QString MandelModel::getTextXY()
{
  if (orbit.worker==nullptr)
    return "-";
  return orbit.worker->toString(&orbit.evaluator.currentParams.c_re)+" +i* "+orbit.worker->toString(&orbit.evaluator.currentParams.c_im);
}

QString doubleToString(double x)
{
  if (x>0)
    return QString("+%1").arg(x, 10, 'f');
  else
    return QString("%1").arg(x, 10, 'f');//returns like 50 digits QString::number(x, 10, 'f');
}

QString MandelModel::getTextInfoGen()
{
  if (orbit.worker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number_store tmp;
    MandelMath::number_place tmp_p;
    orbit.worker->init_(&tmp, &tmp_p);
    orbit.worker->assign(&tmp, &orbit.evaluator.currentParams.c_re);
    orbit.worker->sub(&tmp, &position.center_re_s);
    orbit_x=qRound(orbit.worker->toDouble(&tmp)/position.step_size__)+imageWidth/2;

    orbit.worker->assign(&tmp, &orbit.evaluator.currentParams.c_im);
    orbit.worker->sub(&tmp, &position.center_im_s);
    orbit_y=imageHeight/2-qRound(orbit.worker->toDouble(&tmp)/position.step_size__);
    orbit.worker->cleanup(&tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  MandelPoint *data=&this->pointStore[orbit_x+imageWidth*orbit_y];

  QString state;
  switch (data->state)
  {
    case MandelPoint::State::stUnknown:
      state="Unk"; break;
    case MandelPoint::State::stOutside:
      state="Out"; break;
    case MandelPoint::State::stOutAngle:
      state="OutA"; break;
    case MandelPoint::State::stBoundary:
      state="Bound"; break;
    case MandelPoint::State::stDiverge:
      state="Diver"; break;
    case MandelPoint::State::stMisiur:
      state="Misiur"; break;
    case MandelPoint::State::stPeriod2:
      state="Per2"; break;
    case MandelPoint::State::stPeriod3:
      state="Per3"; break;
    case MandelPoint::State::stMaxIter:
      state="Max"; break;
  }

  return state+" iter="+QString::number(data->iter)+" near="+QString::number(data->near0iter)+
      " fc="+doubleToString(orbit.worker->toDouble(&data->fc_c_re_))+doubleToString(orbit.worker->toDouble(&data->fc_c_im_))+"i";
}

QString MandelModel::getTextInfoSpec()
{
  if (orbit.worker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number_store tmp;
    MandelMath::number_place tmp_p;
    orbit.worker->init_(&tmp, &tmp_p);
    orbit.worker->assign(&tmp, &orbit.evaluator.currentParams.c_re);
    orbit.worker->sub(&tmp, &position.center_re_s);
    orbit_x=qRound(orbit.worker->toDouble(&tmp)/position.step_size__)+imageWidth/2;

    orbit.worker->assign(&tmp, &orbit.evaluator.currentParams.c_im);
    orbit.worker->sub(&tmp, &position.center_im_s);
    orbit_y=imageHeight/2-qRound(orbit.worker->toDouble(&tmp)/position.step_size__);
    orbit.worker->cleanup(&tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  MandelPoint *data=&this->pointStore[orbit_x+imageWidth*orbit_y];

  switch (data->state)
  {
    case MandelPoint::State::stUnknown:
      return " ";
      break;
    case MandelPoint::State::stOutside:
      return QString("ext=")+QString::number(data->exterior_hits); break;
    case MandelPoint::State::stOutAngle:
      return QString("ext=")+QString::number(data->exterior_hits); break;
    case MandelPoint::State::stBoundary:
      return " ";
    case MandelPoint::State::stDiverge:
      return " ";
    case MandelPoint::State::stMisiur:
      return " ";
    case MandelPoint::State::stPeriod2:
      return QString("per=")+QString::number(data->period)+" int="+QString::number(data->interior)   +" mult="+QString::number(orbit.evaluator.bulb.dbg_guessmult); break;
    case MandelPoint::State::stPeriod3:
      return QString("per=")+QString::number(data->period)+" int="+QString::number(data->interior)   +" mult="+QString::number(orbit.evaluator.bulb.dbg_guessmult); break;
    case MandelPoint::State::stMaxIter:
      return " "; break;
  }
  return "-?-?-";
}

ShareableViewInfo MandelModel::getViewInfo()
{
  ShareableViewInfo result;
  result.worker=orbit.worker;
  result.period=orbit.pointData.near0iter;//evaluator.currentData.lookper_lastGuess;//orbit.pointData.period;
  if (result.period<1)
    result.period=1;
  result.scale=position.step_size__;
  orbit.worker->init_(&result.re_, &result.re_p);
  orbit.worker->init_(&result.im, &result.im_p);
  orbit.worker->init_(&result.root_re, &result.rre_p);
  orbit.worker->init_(&result.root_im, &result.rim_p);
  orbit.worker->assign(&result.re_, &orbit.evaluator.currentParams.c_re);
  orbit.worker->assign(&result.im, &orbit.evaluator.currentParams.c_im);
  orbit.worker->assign(&result.root_re, &orbit.evaluator.currentData.root_re);
  orbit.worker->assign(&result.root_im, &orbit.evaluator.currentData.root_im);
  orbit.worker->assign(&orbit.lagu_c_re_, &orbit.evaluator.currentParams.c_re);
  orbit.worker->assign(&orbit.lagu_c_im, &orbit.evaluator.currentParams.c_im);
  orbit.worker->assign(&orbit.lagu_r_re, &orbit.evaluator.currentData.root_re);
  orbit.worker->assign(&orbit.lagu_r_im, &orbit.evaluator.currentData.root_im);
  return result;
}

void MandelModel::transformStore(MandelPoint *old_store, int old_width, int old_height, MandelMath::number_store *old_cre, MandelMath::number_store *old_cim,
                                 MandelPoint *new_store, int new_width, int new_height, const MandelMath::number_store *new_cre, const MandelMath::number_store *new_cim,
                                 int inlog, int new_step_log)
{
  if (position.worker==nullptr)
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
    position.worker->sub(old_cim, new_cim); //and reversing y at the last minute
    position.worker->lshift(old_cim, new_step_log+step_scale_n_shift);
    delta_y_int=position.worker->toRound(old_cim);
    delta_y_int-=(new_height/2)<<step_scale_n_shift;

    //double delta_x=-new_width/2+(new_cre-old_cre)/new_step;
    //delta_x*=(1<<step_scale_n_shift);
    //if (fabs(delta_x-qRound(delta_x))>0.0001)
    //  dbgPoint();
    //int delta_x_int=qRound(delta_x);
    position.worker->rsub(old_cre, new_cre);
    position.worker->lshift(old_cre, new_step_log+step_scale_n_shift);
    delta_x_int=position.worker->toRound(old_cre);
    delta_x_int-=(new_width/2)<<step_scale_n_shift;
  }

  MandelMath::number_store c_im;
  MandelMath::number_store c_re;
  MandelMath::number_place cre_p, cim_p;
  position.worker->init_(&c_re, &cre_p);
  position.worker->init_(&c_im, &cim_p);
  for (int newy=0; newy<new_height; newy++)
  {
    int oldy;
    if ((step_scale_d_mask==0) || ((((newy<<step_scale_n_shift)+delta_y_int)&step_scale_d_mask)==0))
      oldy=(old_height/2) + (((newy<<step_scale_n_shift)+delta_y_int)>>step_scale_d_shift);
    else
      oldy=-1;//call reset() in the second loop, we may still need the points  =imageHeight;
    position.pixelYtoIM(new_height/2-newy, &c_im);
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
          new_store[newy*new_width+newx].assign(position.worker, old_store[oldy*old_width+oldx]);
        else
        {
          position.pixelXtoRE(newx - new_width/2, &c_re);
          new_store[newy*new_width+newx].zero(position.worker, &c_re, &c_im);
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
    position.pixelYtoIM(new_height/2-newy, &c_im);
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
          new_store[newy*new_width+newx].assign(position.worker, old_store[oldy*old_width+oldx]);
        else
        {
          position.pixelXtoRE(newx - new_width/2, &c_re);
          new_store[newy*new_width+newx].zero(position.worker, &c_re, &c_im);
        }
      }
    }
  }
  position.worker->cleanup(&c_re);
  position.worker->cleanup(&c_im);
}

void MandelModel::setView_double(double c_re, double c_im, double scale)
{
  MandelMath::number_store cre_s, cim_s;
  MandelMath::number_place cre_p, cim_p;
  position.worker->init_(&cre_s, &cre_p, c_re);
  position.worker->init_(&cim_s, &cim_p, c_im);
  setView_(&cre_s, &cim_s, scale);
  position.worker->cleanup(&cim_s);
  position.worker->cleanup(&cre_s);
}

void MandelModel::setView_(const MandelMath::number_store *c_re, const MandelMath::number_store *c_im, double scale)
{
  MandelMath::number_store old_cre_s, old_cim_s;
  MandelMath::number_place old_cre_p, old_cim_p;
  position.worker->init_(&old_cre_s, &old_cre_p);
  position.worker->init_(&old_cim_s, &old_cim_p);
  position.worker->assign(&old_cre_s, &position.center_re_s);
  position.worker->assign(&old_cim_s, &position.center_im_s);
  int old_step_log=position.step_log;

  position.setView(c_re, c_im, scale);

  transformStore(pointStore, imageWidth, imageHeight, &old_cre_s, &old_cim_s,
                 pointStore, imageWidth, imageHeight, &position.center_re_s, &position.center_im_s,
                 position.step_log-old_step_log, position.step_log);

  position.worker->cleanup(&old_cim_s);
  position.worker->cleanup(&old_cre_s);

  startNewEpoch();
}

void MandelModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::number_store old_cre_s, old_cim_s;
  MandelMath::number_place old_cre_p, old_cim_p;
  position.worker->init_(&old_cre_s, &old_cre_p);
  position.worker->init_(&old_cim_s, &old_cim_p);
  position.worker->assign(&old_cre_s, &position.center_re_s);
  position.worker->assign(&old_cim_s, &position.center_im_s);

  int dx=qRound(delta_x);
  int dy=qRound(delta_y);
  position.move(dx, dy);
  //qDebug()<<"new c: re="<<position.worker->toString(&position.center_re_s)<<",im="<<position.worker->toString(&position.center_im_s);

  transformStore(pointStore, imageWidth, imageHeight, &old_cre_s, &old_cim_s,
                 pointStore, imageWidth, imageHeight, &position.center_re_s, &position.center_im_s,
                 0, position.step_log);

  position.worker->cleanup(&old_cim_s);
  position.worker->cleanup(&old_cre_s);

  startNewEpoch();
}

void MandelModel::zoom(double x, double y, int inlog)
{
  MandelMath::number_store old_cre_s, old_cim_s;
  MandelMath::number_place old_cre_p, old_cim_p;
  position.worker->init_(&old_cre_s, &old_cre_p);
  position.worker->init_(&old_cim_s, &old_cim_p);
  position.worker->assign(&old_cre_s, &position.center_re_s);
  position.worker->assign(&old_cim_s, &position.center_im_s);
  int old_step_log=position.step_log;

  position.scale(inlog, qRound(x)-imageWidth/2, imageHeight/2-qRound(y));

  transformStore(pointStore, imageWidth, imageHeight, &old_cre_s, &old_cim_s,
                 pointStore, imageWidth, imageHeight, &position.center_re_s, &position.center_im_s,
                 position.step_log-old_step_log, position.step_log);

  position.worker->cleanup(&old_cim_s);
  position.worker->cleanup(&old_cre_s);

  startNewEpoch();
}

void MandelModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  int newLength=width*height;
  MandelPoint *newStore=new MandelPoint[newLength];
  QString size_as_text=QString::number(sizeof(MandelPoint)*newLength);
  for (int pos=size_as_text.length()-3; pos>0; pos-=3)
  {
    size_as_text.insert(pos, '\'');
  }
  qDebug()<<"pointStore uses"<<size_as_text.toLocal8Bit().constData()<<"B"; //lots of work to skip those quotes... can't skip spaces at all
  for (int i=0; i<newLength; i++)
    newStore[i].init(position.worker);

  MandelMath::number_store old_cre_s, old_cim_s;
  MandelMath::number_place old_cre_p, old_cim_p;
  position.worker->init_(&old_cre_s, &old_cre_p);
  position.worker->init_(&old_cim_s, &old_cim_p);
  position.worker->assign(&old_cre_s, &position.center_re_s);
  position.worker->assign(&old_cim_s, &position.center_im_s);

  transformStore(pointStore, imageWidth, imageHeight, &old_cre_s, &old_cim_s,
                 newStore, width, height, &position.center_re_s, &position.center_im_s,
                 0, position.step_log);

  position.worker->cleanup(&old_cim_s);
  position.worker->cleanup(&old_cre_s);

  for (int i=imageWidth*imageHeight-1; i>=0; i--)
    pointStore[i].cleanup(position.worker);
  delete[] pointStore;
  pointStore=newStore;
  imageWidth=width;
  imageHeight=height;

  startNewEpoch();
}

void MandelModel::setWorker(MandelMath::number_worker *newWorker)
{
  if ((imageWidth<=0) || (imageHeight<=0)) //Qt begins with width=0, height=-13
    return;
  /* cannot sutract oldc-newc because different types... so I use a shortcut
  int newLength=imageWidth*imageHeight;
  MandelPoint *newStore=new MandelPoint[newLength];
  QString size_as_text=QString::number(sizeof(MandelPoint)*newLength);
  for (int pos=size_as_text.length()-3; pos>0; pos-=3)
  {
    size_as_text.insert(pos, '\'');
  }
  qDebug()<<"pointStore uses"<<size_as_text.toLocal8Bit().constData()<<"B"; //lots of work to skip those quotes... can't skip spaces at all
  for (int i=0; i<newLength; i++)
    newStore[i].init(position.worker);

  MandelMath::number_store old_cre_s, old_cim_s;
  position.worker->init(&old_cre_s);
  position.worker->init(&old_cim_s);
  position.worker->assign(&old_cre_s, &position.center_re_s);
  position.worker->assign(&old_cim_s, &position.center_im_s);

  MandelMath::number_worker *oldWorker=position.worker;
  position.setNumberType(newWorker->ntype());

  transformStore(pointStore, imageWidth, imageHeight, &old_cre_s, &old_cim_s,
                 newStore, imageWidth, imageHeight, &position.center_re_s, &position.center_im_s,
                 oldWorker, 0, position.step_log);

  oldWorker->cleanup(&old_cim_s);
  oldWorker->cleanup(&old_cre_s);

  for (int i=imageWidth*imageHeight-1; i>=0; i--)
    pointStore[i].cleanup(oldWorker);
  delete[] pointStore;
  pointStore=newStore;
  */

  MandelMath::number_worker::Type oldType=position.worker->ntype();
  position.setNumberType(newWorker->ntype());
  for (int i=imageWidth*imageHeight-1; i>=0; i--)
    pointStore[i].promote(oldType, position.worker->ntype());

  startNewEpoch();
}

void MandelModel::startNewEpoch()
{
  epoch=(epoch+1)%2000000000;
  lastGivenPointIndex_=0;
  effortBonus_=0;
  for (int t=0; t<threadCount; t++)
    if (threads[t].currentParams.pixelIndex<0)
      giveWork(&threads[t]);
}

void MandelModel::reimToPixel(int *circ_x, int *circ_y, const MandelMath::number_store *re, const MandelMath::number_store *im, MandelMath::number_store *tmp)
{
  position.worker->assign(tmp, re);
  position.worker->sub(tmp, &position.center_re_s);
  position.worker->lshift(tmp, position.step_log);
  *circ_x=position.worker->toDouble(tmp)+imageWidth/2;
  position.worker->assign(tmp, im);
  position.worker->sub(tmp, &position.center_im_s);
  position.worker->lshift(tmp, position.step_log);
  *circ_y=imageHeight/2-position.worker->toDouble(tmp);
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

  if (orbit.worker!=position.worker)
  { //should be switched before getting here... not from null to some
    MandelMath::number_worker::Type oldType, newType=position.worker->ntype();
    if (orbit.worker==nullptr)
      oldType=MandelMath::number_worker::Type::typeEmpty;
    else
      oldType=orbit.worker->ntype();
    orbit.pointData.promote(oldType, newType);
    orbit.lagu_c_re_.promote_(oldType, newType, &orbit.lagu_c_re_p);
    orbit.lagu_c_im.promote_(oldType, newType, &orbit.lagu_c_im_p);
    orbit.lagu_r_re.promote_(oldType, newType, &orbit.lagu_r_re_p);
    orbit.lagu_r_im.promote_(oldType, newType, &orbit.lagu_r_im_p);
    orbit.tmp.promote_(oldType, newType, &orbit.tmp_p);

    orbit.bulb.cb_re.promote_(oldType, newType, &orbit.bulb.cb_re_p);
    orbit.bulb.cb_im.promote_(oldType, newType, &orbit.bulb.cb_im_p);
    orbit.bulb.rb_re_.promote_(oldType, newType, &orbit.bulb.rb_re_p);
    orbit.bulb.rb_im.promote_(oldType, newType, &orbit.bulb.rb_im_p);
    orbit.bulb.xc_re.promote_(oldType, newType, &orbit.bulb.xc_re_p);
    orbit.bulb.xc_im.promote_(oldType, newType, &orbit.bulb.xc_im_p);
    orbit.bulb.baseZC_re.promote_(oldType, newType, &orbit.bulb.baseZC_re_p);
    orbit.bulb.baseZC_im.promote_(oldType, newType, &orbit.bulb.baseZC_im_p);
    orbit.bulb.baseCC_re.promote_(oldType, newType, &orbit.bulb.baseCC_re_p);
    orbit.bulb.baseCC_im.promote_(oldType, newType, &orbit.bulb.baseCC_im_p);

    orbit.worker=position.worker;
  };
  MandelPoint *data=&pointStore[y*imageWidth+x];
  switch (data->state)
  {
    case MandelPoint::State::stOutside:
    case MandelPoint::State::stOutAngle:
    {
      int exterior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0xff, 0xff, 0));
      exterior=qRound(data->exterior_hits/position.step_size__);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
      painter.setPen(QColor(0xc0, 0xc0, 0));
      exterior=qRound(data->exterior_avoids/position.step_size__);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
    } break;
    case MandelPoint::State::stPeriod2:
    case MandelPoint::State::stPeriod3:
    {
      int interior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0, 0xff, 0xff));
      interior=qRound(data->interior/position.step_size__);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
      painter.setPen(QColor(0, 0xc0, 0xc0));
      interior=qRound(data->interior/4/position.step_size__);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
    } break;
    default: ;
  }
  orbit.evaluator.switchType(position.worker);
  position.pixelXtoRE(x-imageWidth/2, &orbit.evaluator.currentParams.c_re);
  position.pixelYtoIM(imageHeight/2-y, &orbit.evaluator.currentParams.c_im);
  orbit.evaluator.currentParams.epoch=epoch;
  orbit.evaluator.currentParams.pixelIndex=0;
  orbit.pointData.zero(position.worker, &orbit.evaluator.currentParams.c_re, &orbit.evaluator.currentParams.c_im);
  /*for (int effort=0; effort<=MAX_EFFORT; effort++)
  {
    orbit.evaluator.currentParams.maxiter=1<<effort;
    orbit.evaluator.startCompute(&orbit.pointData, +1);
    orbit.pointData.assign(orbit.worker, orbit.evaluator.currentData);
    if (orbit.pointData.state!=MandelPoint::State::stUnknown)
      break;
  }*/
  orbit.evaluator.currentParams.breakOnNewNearest=true;
  orbit.evaluator.currentParams.maxiter_=1<<MAX_EFFORT;
  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0xff, 0xff)); //paint path
  while ((orbit.pointData.state==MandelPoint::State::stUnknown) &&
         (orbit.pointData.iter<(1<<MAX_EFFORT)))
  {
    int line_sx, line_sy, line_ex, line_ey;
    reimToPixel(&line_sx, &line_sy, &orbit.pointData.f_re, &orbit.pointData.f_im, &orbit.tmp);

    if ((data->state==MandelPoint::State::stPeriod2 || data->state==MandelPoint::State::stPeriod3) &&
        orbit.pointData.iter<data->period)
    {
      orbit.evaluator.currentParams.maxiter_=orbit.pointData.iter+1;
    }
    else if (orbit.evaluator.currentData.lookper_lastGuess==0)
      orbit.evaluator.currentParams.maxiter_=1<<MAX_EFFORT;
    else
      orbit.evaluator.currentParams.maxiter_=(orbit.pointData.iter/orbit.evaluator.currentData.lookper_lastGuess+1)*orbit.evaluator.currentData.lookper_lastGuess;
    orbit.evaluator.startCompute(&orbit.pointData, +1);
    orbit.pointData.assign(orbit.worker, orbit.evaluator.currentData);

    reimToPixel(&line_ex, &line_ey, &orbit.pointData.f_re, &orbit.pointData.f_im, &orbit.tmp);
    painter.drawLine(line_sx, line_sy, line_ex, line_ey);
  }
  if ((orbit.pointData.state==MandelPoint::State::stPeriod2) ||
      (orbit.pointData.state==MandelPoint::State::stPeriod3))
  {
    int circ_x, circ_y;
    painter.setPen(QColor(0, 0xff, 0xff)); //paint root
    reimToPixel(&circ_x, &circ_y, &orbit.pointData.root_re, &orbit.pointData.root_im, &orbit.tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };

    //orbit.bulbinfo
    MandelMath::complex c(orbit.worker, &orbit.evaluator.currentParams.c_re, &orbit.evaluator.currentParams.c_im, true);
    MandelMath::complex cb(orbit.worker, &orbit.bulb.cb_re, &orbit.bulb.cb_im, true);
    MandelMath::complex rb(orbit.worker, &orbit.bulb.rb_re_, &orbit.bulb.rb_im, true);
    MandelMath::complex xc(orbit.worker, &orbit.bulb.xc_re, &orbit.bulb.xc_im, true);
    MandelMath::complex baseZC(orbit.worker, &orbit.bulb.baseZC_re, &orbit.bulb.baseZC_im, true);
    MandelMath::complex baseCC(orbit.worker, &orbit.bulb.baseCC_re, &orbit.bulb.baseCC_im, true);
    orbit.bulb.foundMult=0;
    orbit.bulb.is_card=false;

    orbit.bulb.valid=orbit.evaluator.findBulbBase(orbit.pointData.period, &c, &cb, &rb, &xc, &baseZC, &baseCC, &orbit.bulb.is_card, &orbit.bulb.foundMult);
    if (orbit.bulb.valid)
    {
      painter.setBrush(QBrush(QColor(0, 0xff, 0xff)));
      reimToPixel(&circ_x, &circ_y, &orbit.evaluator.bulb.dbg_first_cb_re, &orbit.evaluator.bulb.dbg_first_cb_im, &orbit.tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0x80, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &orbit.evaluator.bulb.dbg_first_rb_re, &orbit.evaluator.bulb.dbg_first_rb_im, &orbit.tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0x80, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &orbit.bulb.xc_re, &orbit.bulb.xc_im, &orbit.tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        if (orbit.bulb.is_card)
          painter.setPen(QColor(0, 0xff, 0xff)); //card center
        else
          painter.setPen(QColor(0x80, 0xc0, 0xc0)); //bulb center
        //painter.setBrush(Qt::BrushStyle::SolidPattern);
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &orbit.bulb.cb_re, &orbit.bulb.cb_im, &orbit.tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0xff, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &orbit.bulb.rb_re_, &orbit.bulb.rb_im, &orbit.tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0xff, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
    };
  }

  if (!orbit.worker->is0(&orbit.lagu_c_re_) ||
      !orbit.worker->is0(&orbit.lagu_c_im))
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0, 0xff, 0)); //paint c
    reimToPixel(&circ_x, &circ_y, &orbit.lagu_c_re_, &orbit.lagu_c_im, &orbit.tmp);
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
    reimToPixel(&circ_x, &circ_y, &orbit.lagu_r_re, &orbit.lagu_r_im, &orbit.tmp);
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
  wtiElapsed.start();
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      bool knownenum=false;
      MandelPoint *data=&pointStore[y*imageWidth+x];
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleKind:
        {
          switch (data->state)
          {
            case MandelPoint::State::stUnknown:
              image.image->setPixel(x, y, 0xffffffff);
              //image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
              break;
            case MandelPoint::State::stOutside:
            case MandelPoint::State::stOutAngle:
            {
              int b=0x9f+floor(0x60*cos((data->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b<<0));
              knownenum=true;
            } break;
            case MandelPoint::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
              knownenum=true;
            } break;
            case MandelPoint::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
              knownenum=true;
            } break;
            case MandelPoint::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
              knownenum=true;
            } break;
            case MandelPoint::State::stPeriod2:
            case MandelPoint::State::stPeriod3:
            {
              /*int r;
              switch (data->period)
              {
                case 1: r=0xc0; break;
                case 2: r=0xff; break;
                case 3: r=0x80; break;
                default: r=0xe0;
              }*/
              if (data->period>data->near0iter)
                image.image->setPixel(x, y, 0xffff00ff); //seems to only happen by mistake, not in reality
              else
              {
                int index=periodToIndex(data->period);
                //reverse bottom 7 bits:
                int rh=0x73516240>>((index&7)<<2);
                int rl=0x73516240>>((index&0x70)>>2);
                rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
                image.image->setPixel(x, y, 0xff000000+(rh<<16));
              }
              knownenum=true;
            } break;
            case MandelPoint::State::stMaxIter:
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
          switch (data->state)
          {
            case MandelPoint::State::stUnknown:
              image.image->setPixel(x, y, 0x00906090);
              knownenum=true;
              break;
            case MandelPoint::State::stOutside:
            case MandelPoint::State::stOutAngle:
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
              double re=position.worker->toDouble(&data->f_re);
              double im=position.worker->toDouble(&data->f_im);
              double iter=data->iter+6-log2(log2(re*re+im*im)); //+6 to match integer coloring
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
            case MandelPoint::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPoint::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPoint::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffffff);
              knownenum=true;
            } break;
            case MandelPoint::State::stPeriod2:
            case MandelPoint::State::stPeriod3:
            {
              int index=periodToIndex(data->period);
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
            case MandelPoint::State::stMaxIter:
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
          switch (data->state)
          {
            case MandelPoint::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPoint::State::stOutside:
            case MandelPoint::State::stOutAngle:
            {
              double tf;
              if ((data->exterior_avoids>10000) || (data->exterior_avoids<=0))
                tf=0;
              else if (data->exterior_avoids>=1)
                tf=(1-data->exterior_avoids)*1;
              else
                tf=sqrt(1-log(data->exterior_avoids))*2-2;
              int r=0x9f+qRound(0x60*sin(tf*2.828)); //red middle
              int g=0x9f+qRound(0x60*sin(tf*6.928)); //green fastest
              int b=0x9f+qRound(0x60*sin(tf)); //blue slowest
              if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              knownenum=true;
            } break;
            case MandelPoint::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stPeriod2:
            case MandelPoint::State::stPeriod3:
            {
              double re=position.worker->toDouble(&data->fz_r_re);
              double im=position.worker->toDouble(&data->fz_r_im);
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
            case MandelPoint::State::stMaxIter:
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
          switch (data->state)
          {
            case MandelPoint::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPoint::State::stOutside:
            case MandelPoint::State::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((data->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPoint::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stPeriod2:
            case MandelPoint::State::stPeriod3:
            {
              int ti=30;
              if ((data->interior>1) || (data->interior<=0))
                ti=0;
              else
                ti=(qRound(-log(data->interior/4)*300)+12*0xc0) % (6*0xc0);
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
            case MandelPoint::State::stMaxIter:
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
          switch (data->state)
          {
            case MandelPoint::State::stUnknown:
              image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
              break;
            case MandelPoint::State::stOutside:
            case MandelPoint::State::stOutAngle:
            {
              int ti=data->near0iter;
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
            case MandelPoint::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
              knownenum=true;
            } break;
            case MandelPoint::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
              knownenum=true;
            } break;
            case MandelPoint::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
              knownenum=true;
            } break;
            case MandelPoint::State::stPeriod2:
            case MandelPoint::State::stPeriod3:
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

              int ti=data->near0iter;
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
            case MandelPoint::State::stMaxIter:
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
          switch (data->state)
          {
            case MandelPoint::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPoint::State::stOutside:
            case MandelPoint::State::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((data->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPoint::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stPeriod2:
            case MandelPoint::State::stPeriod3:
            {
              double re=position.worker->toDouble(&data->fz_r_re);
              double im=position.worker->toDouble(&data->fz_r_im);
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
            case MandelPoint::State::stMaxIter:
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
          switch (data->state)
          {
            case MandelPoint::State::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              knownenum=true;
              break;
            case MandelPoint::State::stOutside:
            case MandelPoint::State::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((data->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
              knownenum=true;
            } break;
            case MandelPoint::State::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
              knownenum=true;
            } break;
            case MandelPoint::State::stPeriod2:
            case MandelPoint::State::stPeriod3:
            {
              if (!data->has_fc_r)
                image.image->setPixel(x, y, 0xffc0c0c0);
              else
              {
                double mag=sqrt(MandelMath::sqr_double(position.worker->toDouble(&data->fc_c_re_))+
                           MandelMath::sqr_double(position.worker->toDouble(&data->fc_c_im_)));
                if (mag<0) mag=0;
                else if (mag>1) mag=1;
                //int magi=qRound(mag*127.49)+128;
                double mag2=mag*127.49+128;
                double phi=std::atan2(position.worker->toDouble(&data->fc_c_im_), position.worker->toDouble(&data->fc_c_re_))/(2*M_PI);
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
            case MandelPoint::State::stMaxIter:
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
  return wtiElapsed.elapsed();
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
        MandelPoint *pointData=&pointStore[pointIndex];
        bool needsEval=(pointData->state==MandelPoint::State::stUnknown);
        if (!needsEval)
          needsEval=(_selectedPaintStyle==paintStyleFC) &&
                    ((pointData->state==MandelPoint::State::stPeriod2) ||
                     (pointData->state==MandelPoint::State::stPeriod3)) &&
                    (!pointData->has_fc_r);
        if (needsEval)
        {
          bool found=false;
          for (int t=0; t<threadCount; t++)
            if (threads[t].currentParams.pixelIndex==pointIndex)
            {
              found=true;
              break;
            };
          if (evaluator->currentParams.pixelIndex!=-1)
            dbgPoint();
          assert(evaluator->currentParams.pixelIndex==-1);
          if (!found)
          {
            int phasex=(pointIndex%imageWidth-imageWidth/2+position.cached_center_re_mod+32768)%32768;
            int phasey=(pointIndex/imageWidth-imageHeight/2+position.cached_center_im_mod+32768)%32768;
            //int effort=ctz16(phasex)+ctz16(phasey);
            int effort=MandelMath::ctz16(phasex | phasey);
            if (effort>8)
              effort=8;
            effort+=effortBonus_;
            if (effort>=MAX_EFFORT)
              effort=MAX_EFFORT;
            evaluator->currentParams.maxiter_=1<<effort;
            if (pointData->iter>=evaluator->currentParams.maxiter_)
            {
              if (effort>=MAX_EFFORT)
                pointData->state=MandelPoint::State::stMaxIter;
              else if (retryEffortFrom<0)
                retryEffortFrom=pointIndex;
            }
            else
            {
              evaluator->switchType(position.worker);
              position.pixelXtoRE(pointIndex%imageWidth - imageWidth/2, &evaluator->currentParams.c_re);
              position.pixelYtoIM(imageHeight/2-pointIndex/imageWidth, &evaluator->currentParams.c_im);
              evaluator->currentParams.epoch=epoch;
              evaluator->currentParams.pixelIndex=pointIndex;
              evaluator->currentParams.want_fc_r=(_selectedPaintStyle==paintStyleFC);
              if (evaluator->startCompute(pointData, quickrun>=100?-1:0))
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
    MandelPoint *point=&pointStore[me->currentParams.pixelIndex];
    /* it's OK now as fc_r is computed later
    if (point->state!=MandelPoint::State::stUnknown)
      qDebug()<<"Finished pixel finished again";
    else*/
    {
      if (position.worker==nullptr)
        dbgPoint();
      else
        point->assign(position.worker, me->currentData);
      if ((point->state==MandelPoint::State::stUnknown) &&
          (point->iter>=(1<<MAX_EFFORT)))
        point->state=MandelPoint::State::stMaxIter;
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
  switch (_selectedPrecision)
  {
    case precisionDouble: setWorker(&position.number_worker_double_template); break;
    case precisionDDouble: setWorker(&position.number_worker_ddouble_template); break;
    case precisionQDouble: setWorker(&position.number_worker_ddouble_template); break;
    case precisionMulti: setWorker(&position.number_worker_multi_template); break;
  }
  if (orbit.worker!=nullptr)
  {
    MandelMath::number_worker::Type oldType, newType=position.worker->ntype();
    if (orbit.worker==nullptr)
      oldType=MandelMath::number_worker::Type::typeEmpty;
    else
      oldType=orbit.worker->ntype();
    orbit.pointData.promote(oldType, newType);
    orbit.lagu_c_re_.promote_(oldType, newType, &orbit.lagu_c_re_p);
    orbit.lagu_c_im.promote_(oldType, newType, &orbit.lagu_c_im_p);
    orbit.lagu_r_re.promote_(oldType, newType, &orbit.lagu_r_re_p);
    orbit.lagu_r_im.promote_(oldType, newType, &orbit.lagu_r_im_p);
    orbit.tmp.promote_(oldType, newType, &orbit.tmp_p);

    orbit.bulb.cb_re.promote_(oldType, newType, &orbit.bulb.cb_re_p);
    orbit.bulb.cb_im.promote_(oldType, newType, &orbit.bulb.cb_im_p);
    orbit.bulb.rb_re_.promote_(oldType, newType, &orbit.bulb.rb_re_p);
    orbit.bulb.rb_im.promote_(oldType, newType, &orbit.bulb.rb_im_p);
    orbit.bulb.xc_re.promote_(oldType, newType, &orbit.bulb.xc_re_p);
    orbit.bulb.xc_im.promote_(oldType, newType, &orbit.bulb.xc_im_p);
    orbit.bulb.baseZC_re.promote_(oldType, newType, &orbit.bulb.baseZC_re_p);
    orbit.bulb.baseZC_im.promote_(oldType, newType, &orbit.bulb.baseZC_im_p);
    orbit.bulb.baseCC_re.promote_(oldType, newType, &orbit.bulb.baseCC_re_p);
    orbit.bulb.baseCC_im.promote_(oldType, newType, &orbit.bulb.baseCC_im_p);

    orbit.evaluator.switchType(position.worker);
    orbit.worker=position.worker;
  };
}

MandelModel::Position::Position():
  worker(nullptr),
  center_re_s(),
  center_im_s()
{
  setNumberType(MandelMath::number_worker::Type::typeDouble);
  worker->zero(&center_re_s, -0.5);
  worker->zero(&center_im_s, 0.0);
  //center_re_n.reinit(MandelMath::number::Type::typeDDouble);
  //center_im_n.reinit(MandelMath::number::Type::typeDDouble);
  step_log=7;
  step_size__=1.0/128;
  updateCachedDepth();
}

MandelModel::Position::~Position()
{
  if (worker!=nullptr)
  {
    worker->cleanup(&center_re_s);
    worker->cleanup(&center_im_s);
  };
}

void MandelModel::Position::setNumberType(MandelMath::number_worker::Type ntype)
{
  MandelMath::number_worker::Type oldType;
  if (worker!=nullptr)
    oldType=worker->ntype();
  else
    oldType=MandelMath::number_worker::Type::typeEmpty;
  center_re_s.promote_(oldType, ntype, &center_re_p);
  center_im_s.promote_(oldType, ntype, &center_im_p);
  /*if (worker!=nullptr)
  {
    worker->cleanup(&center_re_s);
    worker->cleanup(&center_im_s);
  }*/
  switch (ntype)
  {
#if NUMBER_DOUBLE_EXISTS
    case MandelMath::number_worker::Type::typeDouble:
      worker=&number_worker_double_template;
      break;
#endif //NUMBER_DOUBLE_EXISTS
    case MandelMath::number_worker::Type::typeDDouble:
      worker=&number_worker_ddouble_template;
      break;
    case MandelMath::number_worker::Type::typeMulti:
      worker=&number_worker_multi_template;
      break;
    case MandelMath::number_worker::Type::typeEmpty:
      worker=nullptr;
  }
  /*if (worker)
  {
    worker->init(&center_re_s);
    worker->init(&center_im_s);
  };*/
}

void MandelModel::Position::setView(const MandelMath::number_store *c_re, const MandelMath::number_store *c_im, double scale)
{
  step_log=-ilogb(scale);
  step_size__=ldexp(1.0, -step_log);
  //worker->zero(&center_re_s, ldexp(round(ldexp(c_re, step_log)), -step_log));
  worker->zero(&center_re_s); //worker->swap()
  worker->assign(&center_re_s, c_re);
  worker->lshift(&center_re_s, step_log);
  worker->round(&center_re_s);
  worker->lshift(&center_re_s, -step_log);
  //worker->zero(&center_im_s, ldexp(round(ldexp(c_im, step_log)), -step_log));
  worker->zero(&center_im_s); //worker->swap()
  worker->assign(&center_im_s, c_im);
  worker->lshift(&center_im_s, step_log);
  worker->round(&center_im_s);
  worker->lshift(&center_im_s, -step_log);
  updateCachedDepth();
}

void MandelModel::Position::move(int delta_x, int delta_y)
{
  //qDebug()<<"move ("<<delta_x<<","<<delta_y<<")";
  worker->add_double(&center_re_s, -delta_x*step_size__);
  worker->add_double(&center_im_s, +delta_y*step_size__);
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
    worker->add_double(&center_re_s, center_x*(old_step_size-step_size__));
    worker->add_double(&center_im_s, center_y*(old_step_size-step_size__));
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
    worker->add_double(&center_re_s, (center_x-adjust_x)*(old_step_size-step_size__)); //(old_step_size-step_size__)=(1-(1<<-inlog))*old_step_size
    worker->add_double(&center_im_s, (center_y-adjust_y)*(old_step_size-step_size__));
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
  //number_any d_re_n(MandelMath::number_store::DbgType::typeEmpty, &d_re_s);
  MandelMath::number_store d_re_s;
  MandelMath::number_place d_re_p;
  worker->init_(&d_re_s, &d_re_p);
  worker->assign(&d_re_s, &center_re_s);
  worker->lshift(&d_re_s, step_log-15);
  worker->mod1(&d_re_s);
  worker->lshift(&d_re_s, 15);
  cached_center_re_mod=worker->toRound(&d_re_s);
  worker->cleanup(&d_re_s);

  MandelMath::number_store d_im_s;
  MandelMath::number_place d_im_p;
  worker->init_(&d_im_s, &d_im_p);
  worker->assign(&d_im_s, &center_im_s);
  worker->lshift(&d_im_s, step_log-15);
  worker->mod1(&d_im_s);
  worker->lshift(&d_im_s, 15);
  cached_center_im_mod=worker->toRound(&d_im_s);
  worker->cleanup(&d_im_s);
}

void MandelModel::Position::pixelXtoRE(int x, MandelMath::number_store *result)
{
  //return (x - imageWidth/2)*position.step_size+position.center_re;
  //should be already result->reinit(center_re_n.ntype());
  worker->assign(result, &center_re_s);
  worker->add_double(result, x*step_size__);
}

void MandelModel::Position::pixelYtoIM(int y, MandelMath::number_store *result)
{
  //return (y - imageHeight/2)*position.step_size+position.center_im;
  //should be already result->reinit(center_im_n.ntype());
  worker->assign(result, &center_im_s);
  worker->add_double(result, y*step_size__);
}



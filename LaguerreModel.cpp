#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "LaguerreModel.hpp"
#include "MandelEvaluator.hpp"

LaguerreModel::LaguerreModel(): QObject(), position()
{
  unsigned int oldcw; //524319 = 0x8001F = mask all interrupts, 80bit precision
  MandelMath::fpu_fix_start(&oldcw);
  _selectedPaintStyle=paintStyleCls;//Kind;
  epoch=0;
  imageWidth=0;
  imageHeight=0;
  pointStore=nullptr;
  lastGivenPointIndex_=0;
  //effortBonus=0;
  orbit.worker=nullptr;
  params.period=1;
  position.worker->init(&params.base_re_s_, 0);
  position.worker->init(&params.base_im_s, 0);
  position.worker->init(&params.root_re_s, 0);
  position.worker->init(&params.root_im_s, 0);
  //threadCount=4;
  threadCount=1;//so fast that sending messages is slower QThread::idealThreadCount()-1;
  if (threadCount<1)
    threadCount=1;
  threads=new MandelEvaluator[threadCount];
  for (int t=0; t<threadCount; t++)
  {
    //threads[t].setHint(t);
    QObject::connect(&threads[t], &MandelEvaluator::doneCompute,
                     this, &LaguerreModel::donePixel,
                     Qt::ConnectionType::QueuedConnection);
  }
}

LaguerreModel::~LaguerreModel()
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
  orbit.evaluator.wantStop=true;
  orbit.evaluator.quit();
  orbit.evaluator.wait(1000);
  if (orbit.worker!=nullptr)
    orbit.pointData.cleanup(orbit.worker);

  for (int i=imageWidth*imageHeight-1; i>=0; i--)
  {
    pointStore[i].cleanup(position.worker);
  }

  delete[] pointStore;
  pointStore=nullptr;
  position.worker->cleanup(&params.root_im_s);
  position.worker->cleanup(&params.root_re_s);
  position.worker->cleanup(&params.base_im_s);
  position.worker->cleanup(&params.base_re_s_);
  imageWidth=0;
  imageHeight=0;
}

QString LaguerreModel::pixelXtoRE_str(int x)
{
  MandelMath::number_store num;
  position.worker->assign(&num, &position.center_re_s);
  position.worker->add_double(&num, (x - imageWidth/2)*position.step_size__);
  QString result=position.worker->toString(&num);
  position.worker->cleanup(&num);
  return result;
}

QString LaguerreModel::pixelYtoIM_str(int y)
{
//  return (y - imageHeight/2)*position.step_size+position.center_im;
  MandelMath::number_store num;
  position.worker->assign(&num, &position.center_im_s);
  position.worker->add_double(&num, (y - imageHeight/2)*position.step_size__);
  QString result=position.worker->toString(&num);
  position.worker->cleanup(&num);
  return result;
}

QString LaguerreModel::getTimes()
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

QString LaguerreModel::getTextXY()
{
  if (orbit.worker==nullptr)
    return "-";
  return orbit.worker->toString(&orbit.evaluator.currentParams.c_re)+" +i* "+orbit.worker->toString(&orbit.evaluator.currentParams.c_im);
}

QString LaguerreModel::getTextInfoGen()
{
  if (orbit.worker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number_store tmp;
    orbit.worker->init(&tmp);
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
  LaguerrePoint *data=&this->pointStore[orbit_x+imageWidth*orbit_y];

  QString state;
  switch (data->state)
  {
    case LaguerrePoint::State::stUnknown:
      state="Unk"; break;
    case LaguerrePoint::State::stResolved:
      state="OK"; break;
    case LaguerrePoint::State::stFail:
      state="Err"; break;
  }

  return "Per="+QString::number(params.period)+" "+state+" iter="+QString::number(data->iter)+
        " mu="+QString::number(orbit.first_mu_re, 'f', 3)+","+QString::number(orbit.first_mu_im, 'f', 3);
}

QString LaguerreModel::getTextInfoSpec()
{
  if (orbit.worker==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number_store tmp;
    orbit.worker->init(&tmp);
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
  LaguerrePoint *data=&this->pointStore[orbit_x+imageWidth*orbit_y];

  switch (data->state)
  {
    case LaguerrePoint::State::stUnknown:
      return " ";
      break;
    case LaguerrePoint::State::stResolved:
    {
      MandelMath::complex fz(orbit.worker, &orbit.pointData.fz_r_re, &orbit.pointData.fz_r_im, true);
      QString naiveChoice;
      switch (orbit.pointData.naiveChoice)
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
      return QString("attr=")+QString::number(orbit.worker->toDouble(fz.getMagTmp()))+
             QString(" firstM=")+QString::number(orbit.pointData.firstM)+
             QString(" NChoice=")+naiveChoice;
    } break;
    case LaguerrePoint::State::stFail:
      return " ";
  }
  return "-?-?-";
}

void LaguerreModel::setParams(ShareableViewInfo viewInfo)
{
  position.setNumberType(viewInfo.worker->ntype());
  params.period=viewInfo.period;
  if (position.worker==nullptr)
    dbgPoint();
  else
  {
    position.worker->assign(&params.base_re_s_, &viewInfo.re_);
    position.worker->assign(&params.base_im_s, &viewInfo.im);
    position.worker->assign(&params.root_re_s, &viewInfo.root_re);
    position.worker->assign(&params.root_im_s, &viewInfo.root_im);

    MandelMath::number_store old_cre_s, old_cim_s;
    position.worker->init(&old_cre_s);
    position.worker->init(&old_cim_s);
    position.worker->assign(&old_cre_s, &position.center_re_s);
    position.worker->assign(&old_cim_s, &position.center_im_s);
    int old_step_log=position.step_log;

    position.setView(position.worker->toDouble(&viewInfo.re_), position.worker->toDouble(&viewInfo.im), viewInfo.scale);

    transformStore(pointStore, 0, 0, &old_cre_s, &old_cim_s,
                   pointStore, imageWidth, imageHeight, &position.center_re_s, &position.center_im_s,
                   position.step_log-old_step_log, position.step_log);
    position.worker->cleanup(&old_cim_s);
    position.worker->cleanup(&old_cre_s);

    viewInfo.worker->cleanup(&viewInfo.root_im);
    viewInfo.worker->cleanup(&viewInfo.root_re);
    viewInfo.worker->cleanup(&viewInfo.im);
    viewInfo.worker->cleanup(&viewInfo.re_);

    startNewEpoch();
  }
}

void LaguerreModel::transformStore(LaguerrePoint *old_store, int old_width, int old_height, MandelMath::number_store *old_cre, MandelMath::number_store *old_cim,
                                   LaguerrePoint *new_store, int new_width, int new_height, const MandelMath::number_store *new_cre, const MandelMath::number_store *new_cim,
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
  position.worker->init(&c_im);
  position.worker->init(&c_re);
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

void LaguerreModel::setView(double c_re, double c_im, double scale)
{
  MandelMath::number_store old_cre_s, old_cim_s;
  position.worker->init(&old_cre_s);
  position.worker->init(&old_cim_s);
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

void LaguerreModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::number_store old_cre_s, old_cim_s;
  position.worker->init(&old_cre_s);
  position.worker->init(&old_cim_s);
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

void LaguerreModel::zoom(double x, double y, int inlog)
{
  MandelMath::number_store old_cre_s, old_cim_s;
  position.worker->init(&old_cre_s);
  position.worker->init(&old_cim_s);
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

void LaguerreModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  int newLength=width*height;
  LaguerrePoint *newStore=new LaguerrePoint[newLength];
  QString size_as_text=QString::number(sizeof(LaguerrePoint)*newLength);
  for (int pos=size_as_text.length()-3; pos>0; pos-=3)
  {
    size_as_text.insert(pos, '\'');
  }
  qDebug()<<"laguerreStore uses"<<size_as_text.toLocal8Bit().constData()<<"B";
  for (int i=0; i<newLength; i++)
    newStore[i].init(position.worker);

  MandelMath::number_store old_cre_s, old_cim_s;
  position.worker->init(&old_cre_s);
  position.worker->init(&old_cim_s);
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

void LaguerreModel::startNewEpoch()
{
  epoch=(epoch+1)%2000000000;
  lastGivenPointIndex_=0;
  //effortBonus=0;
  for (int t=0; t<threadCount; t++)
    if (threads[t].currentParams.pixelIndex<0)
      giveWork(&threads[t]);
}

void LaguerreModel::giveWorkAll()
{
  for (int t=0; t<threadCount; t++)
    if (threads[t].currentParams.pixelIndex<0)
      giveWork(&threads[t]);
}

void LaguerreModel::paintOrbit(ShareableImageWrapper image, int x, int y)
{
  if ((x<0) || (x>=imageWidth) || (y<0) || (y>=imageHeight))
    return;
  if (params.period<=0)
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
  if (orbit.worker!=position.worker)
  {
    if (orbit.worker)
    {
      orbit.pointData.cleanup(orbit.worker);
    };
    orbit.pointData.init(position.worker);
    orbit.worker=position.worker;
  };
  /*switch (data->state)
  {
    case LaguerrePoint::State::stUnknown:
    case LaguerrePoint::State::stResolved:
    case LaguerrePoint::State::stFail:
    default: ;
  }*/
  orbit.evaluator.switchType(position.worker);
  {
    int circ_x, circ_y;
    MandelMath::number_store *tmp=&orbit.evaluator.currentData.f_re;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0, 0)); //paint c
    position.worker->assign(tmp, &params.base_re_s_);
    position.worker->sub(tmp, &position.center_re_s);
    position.worker->lshift(tmp, position.step_log);
    circ_x=position.worker->toDouble(tmp)+imageWidth/2;
    position.worker->assign(tmp, &params.base_im_s);
    position.worker->sub(tmp, &position.center_im_s);
    position.worker->lshift(tmp, position.step_log);
    circ_y=imageHeight/2-position.worker->toDouble(tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l3[3]={{circ_x+1, circ_y-1, circ_x-1, circ_y-1},
                   {circ_x-1, circ_y-1, circ_x-1, circ_y+1},
                   {circ_x-1, circ_y+1, circ_x+1, circ_y+1}};
      painter.drawLines(l3, 3);
    };

    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0, 0)); //paint root as /\    .
    position.worker->assign(tmp, &params.root_re_s);
    position.worker->sub(tmp, &position.center_re_s);
    position.worker->lshift(tmp, position.step_log);
    circ_x=position.worker->toDouble(tmp)+imageWidth/2;
    position.worker->assign(tmp, &params.root_im_s);
    position.worker->sub(tmp, &position.center_im_s);
    position.worker->lshift(tmp, position.step_log);
    circ_y=imageHeight/2-position.worker->toDouble(tmp);
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
  position.pixelXtoRE(x-imageWidth/2, &orbit.evaluator.currentParams.c_re); //remember for text labels
  position.pixelYtoIM(imageHeight/2-y, &orbit.evaluator.currentParams.c_im);
  orbit.evaluator.currentParams.epoch=epoch;
  orbit.evaluator.currentParams.pixelIndex=0;
  orbit.pointData.zero(position.worker, &orbit.evaluator.currentParams.c_re, &orbit.evaluator.currentParams.c_im);
  MandelMath::complex base(orbit.worker, &params.base_re_s_, &params.base_im_s, true);
  orbit.worker->assign(&orbit.evaluator.currentData.root_re, &orbit.evaluator.currentParams.c_re);
  orbit.worker->assign(&orbit.evaluator.currentData.root_im, &orbit.evaluator.currentParams.c_im);
  MandelMath::complex root(orbit.worker, &orbit.evaluator.currentData.root_re, &orbit.evaluator.currentData.root_im, true);
  orbit.evaluator.newton(params.period, &base, &root, true, 12);
  //orbit.pointData.assign(orbit.worker, orbit.evaluator.currentData);
  orbit.pointData.iter=orbit.evaluator.newtres_.cyclesNeeded;
  orbit.pointData.firstM=orbit.evaluator.newtres_.firstM;
  orbit.worker->assign(&orbit.pointData.fz_r_re, &orbit.evaluator.newtres_.fz_r_re);
  orbit.worker->add_double(&orbit.pointData.fz_r_re, 1);
  orbit.worker->assign(&orbit.pointData.fz_r_im, &orbit.evaluator.newtres_.fz_r_im);
  orbit.pointData.naiveChoice=orbit.evaluator.newtres_.naiveChoice;
  orbit.first_mu_re=orbit.evaluator.newtres_.firstMu_re;
  orbit.first_mu_im=orbit.evaluator.newtres_.firstMu_im;

  int circ_x, circ_y;
  double tmp_re, tmp_im;
  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0xff, 0)); //final=yellow /
  position.worker->sub(root.re_s, &position.center_re_s);
  position.worker->lshift(root.re_s, position.step_log);
  circ_x=position.worker->toDouble(root.re_s)+imageWidth/2;
  position.worker->sub(root.im_s, &position.center_im_s);
  position.worker->lshift(root.im_s, position.step_log);
  circ_y=imageHeight/2-position.worker->toDouble(root.im_s);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-2, circ_y+2, circ_x+2, circ_y-2);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0, 0xff)); //Newton=purple \      .
  position.worker->sub(&orbit.evaluator.newtres_.first_guess_newt_re, &position.center_re_s);
  position.worker->lshift(&orbit.evaluator.newtres_.first_guess_newt_re, position.step_log);
  circ_x=position.worker->toDouble(&orbit.evaluator.newtres_.first_guess_newt_re)+imageWidth/2;
  position.worker->sub(&orbit.evaluator.newtres_.first_guess_newt_im, &position.center_im_s);
  position.worker->lshift(&orbit.evaluator.newtres_.first_guess_newt_im, position.step_log);
  circ_y=imageHeight/2-position.worker->toDouble(&orbit.evaluator.newtres_.first_guess_newt_im);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-2, circ_y-2, circ_x+2, circ_y+2);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush); //white = is better visible than cyan | so paint it under
  painter.setPen(QColor(0xff, 0xff, 0xff)); //Fejer=white =
  tmp_re=orbit.evaluator.newtres_.first_fejer_re-position.worker->toDouble(&position.center_re_s);
  tmp_re=ldexp(tmp_re, position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=orbit.evaluator.newtres_.first_fejer_im-position.worker->toDouble(&position.center_im_s);
  tmp_im=ldexp(tmp_im, position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-3, circ_y-1, circ_x+3, circ_y-1);
    painter.drawLine(circ_x-3, circ_y+1, circ_x+3, circ_y+1);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush); //
  painter.setPen(QColor(0xff, 0x00, 0x00)); //Naive=red o
  tmp_re=orbit.evaluator.newtres_.first_naive1_re_-position.worker->toDouble(&position.center_re_s);
  tmp_re=ldexp(tmp_re, position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=orbit.evaluator.newtres_.first_naive1_im-position.worker->toDouble(&position.center_im_s);
  tmp_im=ldexp(tmp_im, position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
  };
  tmp_re=orbit.evaluator.newtres_.first_naive2_re-position.worker->toDouble(&position.center_re_s);
  tmp_re=ldexp(tmp_re, position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=orbit.evaluator.newtres_.first_naive2_im-position.worker->toDouble(&position.center_im_s);
  tmp_im=ldexp(tmp_im, position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };
  tmp_re=orbit.evaluator.newtres_.first_naive_re-position.worker->toDouble(&position.center_re_s);
  tmp_re=ldexp(tmp_re, position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=orbit.evaluator.newtres_.first_naive_im-position.worker->toDouble(&position.center_im_s);
  tmp_im=ldexp(tmp_im, position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre=cyan |
  position.worker->sub(&orbit.evaluator.newtres_.first_guess_lagu_re, &position.center_re_s);
  position.worker->lshift(&orbit.evaluator.newtres_.first_guess_lagu_re, position.step_log);
  circ_x=position.worker->toDouble(&orbit.evaluator.newtres_.first_guess_lagu_re)+imageWidth/2;
  position.worker->sub(&orbit.evaluator.newtres_.first_guess_lagu_im, &position.center_im_s);
  position.worker->lshift(&orbit.evaluator.newtres_.first_guess_lagu_im, position.step_log);
  circ_y=imageHeight/2-position.worker->toDouble(&orbit.evaluator.newtres_.first_guess_lagu_im);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x, circ_y-3, circ_x, circ_y+3);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre1=cyan ||
  tmp_re=orbit.evaluator.newtres_.first_lagu1_re-position.worker->toDouble(&position.center_re_s);
  tmp_re=ldexp(tmp_re, position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=orbit.evaluator.newtres_.first_lagu1_im-position.worker->toDouble(&position.center_im_s);
  tmp_im=ldexp(tmp_im, position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-1, circ_y-3, circ_x-1, circ_y+3);
    painter.drawLine(circ_x+1, circ_y-3, circ_x+1, circ_y+3);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre1other=cyan o
  tmp_re=orbit.evaluator.newtres_.first_lagu1o_re-position.worker->toDouble(&position.center_re_s);
  tmp_re=ldexp(tmp_re, position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=orbit.evaluator.newtres_.first_lagu1o_im-position.worker->toDouble(&position.center_im_s);
  tmp_im=ldexp(tmp_im, position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };


  painter.setBrush(Qt::BrushStyle::NoBrush);
  if (orbit.evaluator.newtres_.first_neumaier1_im_!=0)
  {
    painter.setPen(QColor(0xff, 0x00, 0x00)); //Batra_im
    int circ_r=ldexp(abs(orbit.evaluator.newtres_.first_neumaier1_im_), position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  };
  if ((orbit.evaluator.newtres_.first_neumaier1_im_!=0) || (orbit.evaluator.newtres_.first_neumaier1_re_<0))
  {
    painter.setPen(QColor(0xff, 0x00, 0xff)); //Batra
    double mag=sqrt(orbit.evaluator.newtres_.first_neumaier1_re_*orbit.evaluator.newtres_.first_neumaier1_re_+
                    orbit.evaluator.newtres_.first_neumaier1_im_*orbit.evaluator.newtres_.first_neumaier1_im_);
    int circ_r=ldexp(mag, position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  }
  else
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Batra
    int circ_r=ldexp(orbit.evaluator.newtres_.first_neumaier1_re_, position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  }

  painter.setBrush(Qt::BrushStyle::NoBrush);
  if (orbit.evaluator.newtres_.first_neumaier2_im!=0)
  {
    painter.setPen(QColor(0xff, 0x00, 0x00)); //Batra_im
    int circ_r=ldexp(abs(orbit.evaluator.newtres_.first_neumaier2_im), position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  };
  if ((orbit.evaluator.newtres_.first_neumaier2_im!=0) || (orbit.evaluator.newtres_.first_neumaier2_re<0))
  {
    painter.setPen(QColor(0xff, 0x00, 0xff)); //Batra
    double mag=sqrt(orbit.evaluator.newtres_.first_neumaier2_re*orbit.evaluator.newtres_.first_neumaier2_re+
                    orbit.evaluator.newtres_.first_neumaier2_im*orbit.evaluator.newtres_.first_neumaier2_im);
    int circ_r=ldexp(mag, position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }
  else
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Batra
    int circ_r=ldexp(orbit.evaluator.newtres_.first_neumaier2_re, position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }


  image.image->swap(newOverlay);
}

void LaguerreModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return;
  MandelMath::number_store point_re;
  MandelMath::number_store point_im;
  position.worker->init(&point_re);
  position.worker->init(&point_im);
  MandelMath::complex tmp_point(position.worker, &point_re, &point_im, false);
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      bool knownenum=false;
      LaguerrePoint *data=&pointStore[y*imageWidth+x];
      position.pixelXtoRE(x-imageWidth/2, tmp_point.re_s);
      position.pixelYtoIM(imageHeight/2-y, tmp_point.im_s);
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleCls:
        {
          switch (data->state)
          {
            case LaguerrePoint::State::stUnknown:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xff000000);
              knownenum=true;
              break;
            case LaguerrePoint::State::stResolved:
            {
              int r;
              {
                double tr=position.worker->toDouble(&data->fz_r_re);
                double ti=position.worker->toDouble(&data->fz_r_im);
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
              if (data->firstM>=3)
                g=0x7f;
              else if (data->firstM>=2)
                g=0xbf-0x40*(data->firstM-2);
              else if (data->firstM>=1)
                g=0xff-0x40*(data->firstM-1);
              else if (data->firstM>=0)
                g=0x00+0x40*(data->firstM-0);
              else if (data->firstM>=-1)
                g=0x40+0x40*(data->firstM+1);
              else
                g=0x80;
              int b;
              switch (data->iter % 5)
              {
                case  0: b=0x00; break;
                case  1: b=0xff; break;
                case  2: b=0xc0; break;
                case  3: b=0x80; break;
                case  4: b=0x40; break;
                default: b=0xff;
              }
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              knownenum=true;
            } break;
            case LaguerrePoint::State::stFail:
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
}

void LaguerreModel::giveWork(MandelEvaluator *evaluator)
{
  int quickrun=0;
  MandelMath::complex tmpc(position.worker, &params.base_re_s_, &params.base_im_s, true);
  MandelMath::complex root(position.worker, &evaluator->currentData.f_re, &evaluator->currentData.f_im, true);
  for (int pi=0; pi<imageWidth*imageHeight; pi++)
  {
    int pointIndex=(lastGivenPointIndex_+pi)%(imageWidth*imageHeight);
    //if ((lastGivenPointIndex_!=0) && (pointIndex==0))
      //dbgPoint();
    LaguerrePoint *pointData=&pointStore[pointIndex];
    if (pointData->state==LaguerrePoint::State::stUnknown)
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
        evaluator->switchType(position.worker);
        position.pixelXtoRE(pointIndex%imageWidth - imageWidth/2, root.re_s);
        position.pixelYtoIM(imageHeight/2-pointIndex/imageWidth, root.im_s);
        evaluator->currentParams.epoch=epoch;
        evaluator->currentParams.pixelIndex=pointIndex;
        if (evaluator->newton(params.period, &tmpc, &root, true, 12)>0)
        {
          pointData->state=LaguerrePoint::State::stResolved;
        }
        else
        {
          evaluator->newtres_.cyclesNeeded=-1;
          pointData->state=LaguerrePoint::State::stFail;
        };
        donePixel1(evaluator);
        quickrun++;
        if (quickrun>=1000)
        {
          lastGivenPointIndex_=pointIndex;
          QMetaObject::invokeMethod(this,
                                    &LaguerreModel::giveWorkAll,
                                    Qt::ConnectionType::QueuedConnection);
          return;
        };
      }
    }
    //MandelEvaluator::simple(cr, ci, pointStore[y*imageWidth+x]);
  }
}

void LaguerreModel::donePixel1(MandelEvaluator *me)
{
  if ((me->currentParams.epoch==epoch) && (me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
  {
    LaguerrePoint *point=&pointStore[me->currentParams.pixelIndex];
    /*if (point->state!=LaguerrePoint::State::stUnknown)
      qDebug()<<"Finished pixel finished again";
    else*/
    {
      if (position.worker==nullptr)
        dbgPoint();
      else
      {
        //point->assign(position.worker, me->currentData);
        position.worker->assign(&point->f_re, &me->currentData.f_re); //root
        position.worker->assign(&point->f_im, &me->currentData.f_im);
        position.worker->assign(&point->fz_r_re, &me->newtres_.fz_r_re);
        position.worker->add_double(&point->fz_r_re, 1);
        position.worker->assign(&point->fz_r_im, &me->newtres_.fz_r_im);
        point->iter=me->newtres_.cyclesNeeded;
        point->firstM=me->newtres_.firstM;
      }
    }
  }
  else if (me->currentParams.epoch!=epoch)
  { }//qDebug()<<"Old pixel finished";
  else
    qWarning()<<"Invalid pixel finished";
  me->currentParams.pixelIndex=-1;
}

void LaguerreModel::donePixel(MandelEvaluator *me)
{
  me->timeOuterTotal+=me->timeOuter.nsecsElapsed();
  donePixel1(me);
  giveWork(me);
}

LaguerreModel::Position::Position():
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

LaguerreModel::Position::~Position()
{
  if (worker!=nullptr)
  {
    worker->cleanup(&center_re_s);
    worker->cleanup(&center_im_s);
  };
}

void LaguerreModel::Position::setNumberType(MandelMath::number_worker::Type ntype)
{
  //TODO: try to convert old value to new
  if (worker!=nullptr)
  {
    worker->cleanup(&center_re_s);
    worker->cleanup(&center_im_s);
  }
  switch (ntype)
  {
    case MandelMath::number_worker::Type::typeDouble:
      worker=&number_worker_double_template;
      break;
    case MandelMath::number_worker::Type::typeDDouble:
      worker=&number_worker_ddouble_template;
      break;
    case MandelMath::number_worker::Type::typeMulti:
      worker=&number_worker_multi_template;
      break;
    case MandelMath::number_worker::Type::typeEmpty:
      worker=nullptr;
  }
  if (worker)
  {
    worker->init(&center_re_s);
    worker->init(&center_im_s);
  };
}

void LaguerreModel::Position::setView(double c_re, double c_im, double scale)
{
  step_log=-ilogb(scale);
  step_size__=ldexp(1.0, -step_log);
  worker->zero(&center_re_s, ldexp(round(ldexp(c_re, step_log)), -step_log));
  worker->zero(&center_im_s, ldexp(round(ldexp(c_im, step_log)), -step_log));
  updateCachedDepth();
}

void LaguerreModel::Position::move(int delta_x, int delta_y)
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
    if (step_log+inlog>MAX_ZOOM_IN_DOUBLE)
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

void LaguerreModel::Position::updateCachedDepth()
{
  //number_any d_re_n(MandelMath::number_store::DbgType::typeEmpty, &d_re_s);
  MandelMath::number_store d_re_s;
  worker->init(&d_re_s);
  worker->assign(&d_re_s, &center_re_s);
  worker->lshift(&d_re_s, step_log-15);
  worker->frac_pos(&d_re_s);
  worker->lshift(&d_re_s, 15);
  cached_center_re_mod=worker->toRound(&d_re_s);
  worker->cleanup(&d_re_s);

  MandelMath::number_store d_im_s;
  worker->init(&d_im_s);
  worker->assign(&d_im_s, &center_im_s);
  worker->lshift(&d_im_s, step_log-15);
  worker->frac_pos(&d_im_s);
  worker->lshift(&d_im_s, 15);
  cached_center_im_mod=worker->toRound(&d_im_s);
  worker->cleanup(&d_im_s);
}

void LaguerreModel::Position::pixelXtoRE(int x, MandelMath::number_store *result)
{
  //return (x - imageWidth/2)*position.step_size+position.center_re;
  //should be already result->reinit(center_re_n.ntype());
  worker->assign(result, &center_re_s);
  worker->add_double(result, x*step_size__);
}

void LaguerreModel::Position::pixelYtoIM(int y, MandelMath::number_store *result)
{
  //return (y - imageHeight/2)*position.step_size+position.center_im;
  //should be already result->reinit(center_im_n.ntype());
  worker->assign(result, &center_im_s);
  worker->add_double(result, y*step_size__);
}



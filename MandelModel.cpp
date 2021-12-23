#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "MandelModel.hpp"
#include "MandelEvaluator.hpp"

int ctz16(int x)
{
  int ctzidx=(((0xF65*(x&-x))&0x7800)>>10);
  int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);
  return ctz1;
  /* shift, append i, rotate to be max; is i followed by 0 only? -> negate else not negate
  ctz 8 bit
  00011101
  000-00i-i00-y   1
   001-01i-1i0-y   1
    011-11i-11i-y   1
     111-11i-11i-y   0
      110-10i-i10-n   1
       101-01i-1i0-y   0
        010-10i-i10-n   0
         100-00i-i00-y   0
          000
   0x1D * 0=0x0000 >>4&7=0
   0x1D * 1=0x001D >>4&7=1
   0x1D * 2=0x003A >>4&7=3   MAGIC[3]:=ctz(2)
   0x1D * 4=0x0074 >>4&7=7   MAGIC[7]:=ctz(4)
   0x1D * 8=0x00E8 >>4&7=6
   0x1D *16=0x01D0 >>4&7=5
   0x1D *32=0x03A0 >>4&7=2
   0x1D *64=0x0740 >>4&7=4
   0x1D*128=0x0E80 >>4&7=0   MAGIC[0]=ctz(128)
   0x1D*256=0x1D00 >>4&7=0
   (0x23461507 >> (((0x1D*(i&-i))&0x70)>>2))&0x07 ?= ctz(0..255)

   ctz 16 bit
   0000-000i-i000-y 1
    0001-001i-1i00-y 1
     0011-011i-11i0-y 1
      0111-111i-111i-y 1
       1111-111i-111i-y 0
        1110-110i-i110-n 1
         1101-101i-1i10-n 1
          1011-011i-11i0-y 0
           0110-110i-i110-n 0
            1100-100i-i100-n 1
             1001-001i-1i00-y 0
              0010-010i-10i0 or i010-y or n-1 or 0  (10i0 correct)
               0101-101i-1i10-n 0
                1010-010i-10i0 or i010-y or n-0 or 1 (10i0 correct)
                 0100-100i-i100-n 0
                  1000-000i-i000-y 0
                   0000
   0000111101100101 = 0xF65
   (0x... >> (((0xF65*(i&-i))&0x7800)>>10))&0x0F ?= ctz(0..255)
         FEDCBA9876543210
   MAGIC=34586C9E27BD1A0F
   MAGIC23=0x59EC6C8C; //00 0101 1001 1110 1100 0110 1100 1000 1100
   MAGIC01=0xC486BD63; // 1100 0100 1000 0110 1011 1101 0110 0011
    int ctzidx=(((0xF65*(i&-i))&0x7800)>>10);
    int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);

   ctz 32 bit
   00000-0000i-i0000-y 1
    00001-0001i-1i000-y 1                                        0
     00011-0011i-11i00-y 1                                       1
      00111-0111i-111i0-y 1                                      2
       01111-1111i-1111i-y 1                                     3
        11111-1111i-1111i-y 0                                    4
         11110-1110i-i1110-n 1                                   5
          11101-1101i-1i110-n 1                                  6
           11011-1011i-11i10-n 1                                 7
            10111-0111i-111i0-y 0                                8
             01110-1110i-i1110-n 0                               9
              11100-1100i-i1100-n 1                              a
               11001-1001i-1i100-n 1                             b
                10011-0011i-11i00-y 0                            c
                 00110-0110i-110i0-y 1
                  01101-1101i-1i110-n 0
                   11010-1010i-i1010-n 1
                    10101-0101i-1i010-n 1
                     01011-1011i-11i10-n 0
                      10110-0110i-110i0-y 0
                       01100-1100i-i1100-n 0
                        11000-1000i-i1000-n 1
                         10001-0001i-1i000-y 0
                          00010-0010i-10i00-y 1
                           00101-0101i-1i010-n 0
                            01010-1010i-i1010-n 0
                             10100-0100i-i0100-n 1
                              01001-1001i-1i100-n 0
                               10010-0010i-10i00-y 0
                                00100-0100i-i0100-n 0
                                 01000-1000i-i1000-n 0
                                  10000-0000i-i0000-y 0
                                   00000
   00000111110111001101011000101001 = 0000 0111 1101 1100 1101 0110 0010 1001 = 0x7DCD629
         1f 1e 1d 1c 1b 1a 19 18 17 16 15 14 13 12 11 10 0f 0e 0d 0c 0b 0a 09 08 07 06 05 04 03 02 01 00
   MAGIC= 4  5  6  a  7  f  b 14  8 12 10 19  c 1b 15 1e  3  9  e 13 11 18 1a 1d  2  d 17 1c  1 16  0 1f
   MAGIC4=00000001011101110001111100110101 = 0000 0001 0111 0111 0001 1111 0011 0101 0000 = 0x1771f350
   MAGIC3=00010110100111010110011101010001 = 000 1011 0100 1110 1011 0011 1010 1000 1000 =  0xb4eb3a88
   MAGIC2=11101101000010110010000101110101 = .. 1110 1101 0000 1011 0010 0001 0111 0101 = 0xed0b2175; 0x3b42c85d4 won't fit
   MAGIC1=00111110010001011011001010100101 = 0 0111 1100 1000 1011 0110 0101 0100 1010 = 0x7c8b654a
   MAGIC0=01001110000101101101100101101001 = 0100 1010 0001 0110 1101 1001 0110 1001 = 0x4e16d969

    int ctzidx=(((0x7DCD629*(i&-i))&0x7C000000)>>26);
    int ctz1=((0x1771f350 >> ctzidx)&0x10) |
             ((0xb4eb3a88 >> ctzidx)&0x08) |
             ((0xed0b2175 >> ctzidx)&1)<<2 |
             ((0x7c8b654a >> ctzidx)&0x02) |
             ((0x4e16d969 >> ctzidx)&0x01);
   */

  /*for (int i=0; i<256; i++) //0 returns 7
  {
    int ctz1=(0x23461507 >> (((0x1D*(i&-i))&0x70)>>2))&0x07;
    int ctz2=0;
    if (i&1) ctz2=0;
    else if (i&0x03) ctz2=1;
    else if (i&0x07) ctz2=2;
    else if (i&0x0F) ctz2=3;
    else if (i&0x1F) ctz2=4;
    else if (i&0x3F) ctz2=5;
    else if (i&0x7F) ctz2=6;
    else if (i&0xFF) ctz2=7;
    else ctz2=7;
    if (ctz2!=ctz1)
      dbgPoint();
  }

  for (int i=0; i<65536; i++) //0 returns 15
  {
    int ctzidx=(((0xF65*(i&-i))&0x7800)>>10);
    int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);
    int ctz2=0;
    if (i&1) ctz2=0;
    else if (i&0x03) ctz2=1;
    else if (i&0x07) ctz2=2;
    else if (i&0x0F) ctz2=3;
    else if (i&0x1F) ctz2=4;
    else if (i&0x3F) ctz2=5;
    else if (i&0x7F) ctz2=6;
    else if (i&0xFF) ctz2=7;
    else if (i&0x1FF) ctz2=8;
    else if (i&0x3FF) ctz2=9;
    else if (i&0x7FF) ctz2=10;
    else if (i&0xFFF) ctz2=11;
    else if (i&0x1FFF) ctz2=12;
    else if (i&0x3FFF) ctz2=13;
    else if (i&0x7FFF) ctz2=14;
    else if (i&0xFFFF) ctz2=15;
    else ctz2=15;
    if (ctz2!=ctz1)
      dbgPoint();
  }

  for (int ii=0; ii<32; ii++)
  {
    unsigned int i=(1u<<ii);
    int ctzidx=(((0x7DCD629*(i&-i))&0x7C000000)>>26);
    int ctz1=((0x1771f350 >> ctzidx)&0x10) |
             ((0xb4eb3a88 >> ctzidx)&0x08) |
             ((0xed0b2175 >> ctzidx)&1)<<2 |
             ((0x7c8b654a >> ctzidx)&0x02) |
             ((0x4e16d969 >> ctzidx)&0x01);
    if (ii==32)
    {
      if (ctz1!=31)
        dbgPoint();
    }
    else if (ctz1!=ii)
      dbgPoint();
  }*/
}

MandelModel::MandelModel(): QObject(), position()
{
  unsigned int oldcw; //524319 = 0x8001F = mask all interrupts, 80bit precision
  MandelMath::fpu_fix_start(&oldcw);
  _selectedPaintStyle=paintStyleCls;//Kind;
  epoch=0;
  imageWidth=0;
  imageHeight=0;
  pointStore=nullptr;
  lastGivenPointIndex_=0;
  effortBonus=0;
  orbit.worker=nullptr;
  //threadCount=4;
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

QString MandelModel::getTextInfoGen()
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

  return state+" iter="+QString::number(data->iter)+" near="+QString::number(data->near0iter);
}

QString MandelModel::getTextInfoSpec()
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
      return QString("per=")+QString::number(data->period)+" int="+QString::number(data->interior); break;
    case MandelPoint::State::stPeriod3:
      return QString("per=")+QString::number(data->period)+" int="+QString::number(data->interior); break;
    case MandelPoint::State::stMaxIter:
      return " "; break;
  }
  return "-?-?-";
}

ShareableViewInfo MandelModel::getViewInfo()
{
  ShareableViewInfo result;
  result.worker=orbit.worker;
  result.period=orbit.pointData.period;
  result.scale=position.step_size__;
  orbit.worker->init(&result.re_);
  orbit.worker->init(&result.im);
  orbit.worker->init(&result.root_re);
  orbit.worker->init(&result.root_im);
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

void MandelModel::setView(double c_re, double c_im, double scale)
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

void MandelModel::drag(double delta_x, double delta_y)
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

void MandelModel::zoom(double x, double y, int inlog)
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

void MandelModel::startNewEpoch()
{
  epoch=(epoch+1)%2000000000;
  lastGivenPointIndex_=0;
  effortBonus=0;
  for (int t=0; t<threadCount; t++)
    if (threads[t].currentParams.pixelIndex<0)
      giveWork(&threads[t]);
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

/*  {
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

  /*/fillRect(transparent) does nothing, does not clear the image
  painter.setPen(QColor("cyan"));
  painter.setBrush(QBrush(QColor("yellow")));
  painter.drawEllipse(150, 50, 100, 50);
  painter.fillRect(0, 0, imageWidth, imageHeight, Qt::GlobalColor::transparent);
  painter.drawEllipse(150, 150, 100, 50);
  //painter.setBackground(QBrush(Qt::BrushStyle::NoBrush));
  painter.setBackgroundMode(Qt::BGMode::TransparentMode);
  painter.eraseRect(200, 75, 20, 100);*/

  MandelPoint *data=&pointStore[y*imageWidth+x];
  if (orbit.worker!=position.worker)
  {
    if (orbit.worker)
    {
      orbit.pointData.cleanup(orbit.worker);
      orbit.worker->cleanup(&orbit.tmp);
      orbit.worker->cleanup(&orbit.lagu_r_im);
      orbit.worker->cleanup(&orbit.lagu_r_re);
      orbit.worker->cleanup(&orbit.lagu_c_im);
      orbit.worker->cleanup(&orbit.lagu_c_re_);
    };
    orbit.pointData.init(position.worker);
    orbit.worker=position.worker;
    orbit.worker->init(&orbit.lagu_c_re_);
    orbit.worker->init(&orbit.lagu_c_im);
    orbit.worker->init(&orbit.lagu_r_re);
    orbit.worker->init(&orbit.lagu_r_im);
    orbit.worker->init(&orbit.tmp);
  };
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
  for (int effort=0; effort<=MAX_EFFORT; effort++)
  {
    orbit.evaluator.currentParams.maxiter=1<<effort;
    orbit.evaluator.startCompute(&orbit.pointData, +1);
    orbit.pointData.assign(orbit.worker, orbit.evaluator.currentData);
    if (orbit.pointData.state!=MandelPoint::State::stUnknown)
      break;
  }

  if (!orbit.worker->is0(&orbit.lagu_c_re_) ||
      !orbit.worker->is0(&orbit.lagu_c_im))
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0, 0xff, 0)); //paint c
    position.worker->assign(&orbit.tmp, &orbit.lagu_c_re_);
    position.worker->sub(&orbit.tmp, &position.center_re_s);
    position.worker->lshift(&orbit.tmp, position.step_log);
    circ_x=position.worker->toDouble(&orbit.tmp)+imageWidth/2;
    position.worker->assign(&orbit.tmp, &orbit.lagu_c_im);
    position.worker->sub(&orbit.tmp, &position.center_im_s);
    position.worker->lshift(&orbit.tmp, position.step_log);
    circ_y=imageHeight/2-position.worker->toDouble(&orbit.tmp);
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
    position.worker->assign(&orbit.tmp, &orbit.lagu_r_re);
    position.worker->sub(&orbit.tmp, &position.center_re_s);
    position.worker->lshift(&orbit.tmp, position.step_log);
    circ_x=position.worker->toDouble(&orbit.tmp)+imageWidth/2;
    position.worker->assign(&orbit.tmp, &orbit.lagu_r_im);
    position.worker->sub(&orbit.tmp, &position.center_im_s);
    position.worker->lshift(&orbit.tmp, position.step_log);
    circ_y=imageHeight/2-position.worker->toDouble(&orbit.tmp);
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

void MandelModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return;
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
              image.image->setPixel(x, y, 0xff808080);
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
              double re=position.worker->toDouble(&data->fc_c_re);
              double im=position.worker->toDouble(&data->fc_c_im);
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
      }

    }
  //return result;
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
        if (pointData->state==MandelPoint::State::stUnknown)
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
            int effort=ctz16(phasex | phasey);
            if (effort>8)
              effort=8;
            effort+=effortBonus;
            if (effort>=MAX_EFFORT)
              effort=MAX_EFFORT;
            evaluator->currentParams.maxiter=1<<effort;
            if (pointData->iter>=evaluator->currentParams.maxiter)
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
              if (evaluator->startCompute(pointData, quickrun>=1000?-1:0))
              //if (worker->startCompute(pointData, true))
              {
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
    if ((retryEffortFrom>=0) && (effortBonus<MAX_EFFORT))
    {
      effortBonus++;
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
    if (point->state!=MandelPoint::State::stUnknown)
      qDebug()<<"Finished pixel finished again";
    else
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
  donePixel1(me);
  giveWork(me);
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

void MandelModel::Position::setView(double c_re, double c_im, double scale)
{
  step_log=-ilogb(scale);
  step_size__=ldexp(1.0, -step_log);
  worker->zero(&center_re_s, ldexp(round(ldexp(c_re, step_log)), -step_log));
  worker->zero(&center_im_s, ldexp(round(ldexp(c_im, step_log)), -step_log));
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

void MandelModel::Position::updateCachedDepth()
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



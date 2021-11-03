#include <math.h>
#include <QObject>
#include <QDebug>

#include "MandelModel.hpp"
#include "MandelEvaluator.hpp"
#include "double_double.hpp"

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
   MAGIC2=11101101000010110010000101110101 = 1110 1101 0000 1011 0010 0001 0111 0101 = 0xed0b2175; 0x3b42c85d4 won't fit
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
  MandelMath::fpu_fix_start(nullptr);
  epoch=0;
  imageWidth=0;
  imageHeight=0;
  pointStore=nullptr;
  lastGivenPointIndex_=0;
  effortBonus=0;
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

  MandelMath::number::Type ntype=position.center_re_n.ntype();
  for (int i=imageWidth*imageHeight-1; i>=0; i--)
  {
    pointStore[i].cleanup(ntype);
  }

  delete[] pointStore;
  pointStore=nullptr;
  imageWidth=0;
  imageHeight=0;
}

QString MandelModel::pixelXtoRE_str(int x)
{
  MandelMath::number_any num(&position.center_re_n);
  num.impl->add_double((x - imageWidth/2)*position.step_size__);
  return num.impl->toString();
}

QString MandelModel::pixelYtoIM_str(int y)
{
//  return (y - imageHeight/2)*position.step_size+position.center_im;
  MandelMath::number_any num(&position.center_im_n);
  num.impl->add_double((y - imageHeight/2)*position.step_size__);
  return num.impl->toString();
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

void MandelModel::transformStore(MandelPoint *old_store, int old_width, int old_height, MandelMath::number *old_cre, MandelMath::number *old_cim,
                                 MandelPoint *new_store, int new_width, int new_height, const MandelMath::number_store *new_cre, const MandelMath::number_store *new_cim,
                                 int inlog, int new_step_log)
{
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
    old_cim->sub(new_cim); //and reversing y at the last minute
    old_cim->lshift_(new_step_log+step_scale_n_shift);
    delta_y_int=old_cim->toRound();
    delta_y_int-=(new_height/2)<<step_scale_n_shift;

    //double delta_x=-new_width/2+(new_cre-old_cre)/new_step;
    //delta_x*=(1<<step_scale_n_shift);
    //if (fabs(delta_x-qRound(delta_x))>0.0001)
    //  dbgPoint();
    //int delta_x_int=qRound(delta_x);
    old_cre->rsub(new_cre);
    old_cre->lshift_(new_step_log+step_scale_n_shift);
    delta_x_int=old_cre->toRound();
    delta_x_int-=(new_width/2)<<step_scale_n_shift;
  }

  MandelMath::number::Type ntype=position.center_re_n.ntype();
  for (int newy=0; newy<new_height; newy++)
  {
    int oldy;
    if ((step_scale_d_mask==0) || ((((newy<<step_scale_n_shift)+delta_y_int)&step_scale_d_mask)==0))
      oldy=(old_height/2) + (((newy<<step_scale_n_shift)+delta_y_int)>>step_scale_d_shift);
    else
      oldy=-1;//call reset() in the second loop, we may still need the points  =imageHeight;
    for (int newx=0; newx<new_width; newx++)
    {
      int oldx;
      if ((step_scale_d_mask==0) || ((((newx<<step_scale_n_shift)+delta_x_int)&step_scale_d_mask)==0))
        oldx=(old_width/2) + (((newx<<step_scale_n_shift)+delta_x_int)>>step_scale_d_shift);
      else
        oldx=-1;//call reset() in the second loop, we may still need the points  =imageWidth;
      if ((oldy>newy) || ((oldy==newy) && (oldx>newx)))
      {
        switch (ntype)
        {
          case MandelMath::number::Type::typeEmpty:
          case MandelMath::number::Type::typeDouble:
          {
            if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
              new_store[newy*new_width+newx].assign_double(old_store[oldy*old_width+oldx]);
            else
              new_store[newy*new_width+newx].zero(ntype);
          } break;
          case MandelMath::number::Type::typeDDouble:
          {
            if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
              new_store[newy*new_width+newx].assign_ddouble(old_store[oldy*old_width+oldx]);
            else
              new_store[newy*new_width+newx].zero(ntype);
          } break;
          case MandelMath::number::Type::typeMulti:
          {
            if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
              new_store[newy*new_width+newx].assign_multi(old_store[oldy*old_width+oldx]);
            else
              new_store[newy*new_width+newx].zero(ntype);
          } break;
        }

      }
    }
  }
  for (int newy=(new_height-1)&0xfffffff; newy>=0; newy--) //avoid Clang warning about newy possibly ~ 2^31
  {
    int oldy;
    if ((step_scale_d_mask==0) || ((((newy<<step_scale_n_shift)+delta_y_int)&step_scale_d_mask)==0))
      oldy=(old_height/2) + (((newy<<step_scale_n_shift)+delta_y_int)>>step_scale_d_shift);
    else
      oldy=-1;
    for (int newx=new_width-1; newx>=0; newx--)
    {
      int oldx;
      if ((step_scale_d_mask==0) || ((((newx<<step_scale_n_shift)+delta_x_int)&step_scale_d_mask)==0))
        oldx=(old_width/2) + (((newx<<step_scale_n_shift)+delta_x_int)>>step_scale_d_shift);
      else
        oldx=-1;
      if ((oldy<newy) || ((oldy==newy) && (oldx<=newx)))
      {
        switch (ntype)
        {
          case MandelMath::number::Type::typeEmpty:
          case MandelMath::number::Type::typeDouble:
          {
            if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
              new_store[newy*new_width+newx].assign_double(old_store[oldy*old_width+oldx]);
            else
              new_store[newy*new_width+newx].zero(ntype);
          } break;
          case MandelMath::number::Type::typeDDouble:
          {
            if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
              new_store[newy*new_width+newx].assign_ddouble(old_store[oldy*old_width+oldx]);
            else
              new_store[newy*new_width+newx].zero(ntype);
          } break;
          case MandelMath::number::Type::typeMulti:
          {
            if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
              new_store[newy*new_width+newx].assign_multi(old_store[oldy*old_width+oldx]);
            else
              new_store[newy*new_width+newx].zero(ntype);
          } break;
        }
      }
    }
  }
}

void MandelModel::setView(double c_re, double c_im, double scale)
{
  MandelMath::number_any old_cre_n(&position.center_re_n);
  MandelMath::number_any old_cim_n(&position.center_im_n);
  int old_step_log=position.step_log;

  position.setView(c_re, c_im, scale);

  transformStore(pointStore, imageWidth, imageHeight, old_cre_n.impl, old_cim_n.impl,
                 pointStore, imageWidth, imageHeight, position.center_re_n.impl->store, position.center_im_n.impl->store,
                 position.step_log-old_step_log, position.step_log);

  startNewEpoch();
}

void MandelModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::number_any old_cre_n(&position.center_re_n);
  MandelMath::number_any old_cim_n(&position.center_im_n);

  int dx=qRound(delta_x);
  int dy=qRound(delta_y);
  position.move(dx, dy);
  qDebug()<<"new c: re="<<position.center_re_n.impl->toString()<<",im="<<position.center_im_n.impl->toString();

  transformStore(pointStore, imageWidth, imageHeight, old_cre_n.impl, old_cim_n.impl,
                 pointStore, imageWidth, imageHeight, position.center_re_n.impl->store, position.center_im_n.impl->store,
                 0, position.step_log);

  startNewEpoch();
}

void MandelModel::zoom(double x, double y, int inlog)
{
  MandelMath::number_any old_cre_n(&position.center_re_n);
  MandelMath::number_any old_cim_n(&position.center_im_n);
  int old_step_log=position.step_log;

  position.scale(inlog, qRound(x)-imageWidth/2, imageHeight/2-qRound(y));

  transformStore(pointStore, imageWidth, imageHeight, old_cre_n.impl, old_cim_n.impl,
                 pointStore, imageWidth, imageHeight, position.center_re_n.impl->store, position.center_im_n.impl->store,
                 position.step_log-old_step_log, position.step_log);

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
  MandelMath::number::Type ntype=position.center_re_n.ntype();
  for (int i=0; i<newLength; i++)
    newStore[i].init(ntype);

  MandelMath::number_any old_cre_n(&position.center_re_n);
  MandelMath::number_any old_cim_n(&position.center_im_n);

  transformStore(pointStore, imageWidth, imageHeight, old_cre_n.impl, old_cim_n.impl,
                 newStore, width, height, position.center_re_n.impl->store, position.center_im_n.impl->store,
                 0, position.step_log);

  for (int i=imageWidth*imageHeight-1; i>=0; i--)
    pointStore[i].cleanup(ntype);
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


void MandelModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return;
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      MandelPoint *data=&pointStore[y*imageWidth+x];
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (data->state)
      {
        case MandelPoint::State::stUnknown:
          image.image->setPixel(x, y, 0x00906090);
          break;
        case MandelPoint::State::stOutside:
        {
          /*int r=128+floor(127*cos(data->iter/10.0*2*3.1415926535));
          int g=128+floor(127*cos((data->iter/10.0+0.333)*2*3.1415926535));
          int b=128+floor(127*cos((data->iter/10.0+0.667)*2*3.1415926535));
          if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
            result.setPixel(x, y, 0xffffffff);
          else
            result.setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));*/
          switch (data->iter % 12)
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
        } break;
        case MandelPoint::State::stMaxIter:
        {
          image.image->setPixel(x, y, 0xff808080);
        } break;
      }
    }
  //return result;
}

void MandelModel::giveWork(MandelEvaluator *worker)
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
          if (worker->currentParams.pixelIndex!=-1)
            dbgPoint();
          assert(worker->currentParams.pixelIndex==-1);
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
            worker->currentParams.maxiter=1<<effort;
            if (pointData->iter>=worker->currentParams.maxiter)
            {
              if (effort>=MAX_EFFORT)
                pointData->state=MandelPoint::State::stMaxIter;
              else if (retryEffortFrom<0)
                retryEffortFrom=pointIndex;
            }
            else
            {
              worker->switchType(position.center_re_n.ntype());
              position.pixelXtoRE(pointIndex%imageWidth - imageWidth/2, &worker->currentParams.cr_n);
              position.pixelYtoIM(imageHeight/2-pointIndex/imageWidth, &worker->currentParams.ci_n);
              worker->currentParams.epoch=epoch;
              worker->currentParams.pixelIndex=pointIndex;
              if (worker->startCompute(pointData, quickrun>=1000))
              //if (worker->startCompute(pointData, true))
              {
                worker->timeOuter.start();
                lastGivenPointIndex_=pointIndex;
                return;
              }
              else
              {
                donePixel1(worker);
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
      switch (position.center_re_n.ntype())
      {
        case MandelMath::number::Type::typeDouble:
          point->assign_double(me->currentData);
          break;
        case MandelMath::number::Type::typeDDouble:
          point->assign_ddouble(me->currentData);
          break;
        case MandelMath::number::Type::typeMulti:
          point->assign_multi(me->currentData);
          break;
        case MandelMath::number::Type::typeEmpty:
          dbgPoint();
      }
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
  center_re_n(),
  center_im_n()
{
  center_re_n.reinit(MandelMath::number::Type::typeDouble);
  center_im_n.reinit(MandelMath::number::Type::typeDouble);
  //center_re_n.reinit(MandelMath::number::Type::typeDDouble);
  //center_im_n.reinit(MandelMath::number::Type::typeDDouble);
  center_re_n.impl->zero(-0.5);
  center_im_n.impl->zero(0);
  step_log=7;
  step_size__=1.0/128;
  updateCachedDepth();
}

void MandelModel::Position::setNumberType(MandelMath::number::Type ntype)
{
  center_re_n.reinit(ntype);
  center_im_n.reinit(ntype);
}

void MandelModel::Position::setView(double c_re, double c_im, double scale)
{
  step_log=-ilogb(scale);
  step_size__=ldexp(1.0, -step_log);
  center_re_n.impl->zero(ldexp(round(ldexp(c_re, step_log)), -step_log));
  center_im_n.impl->zero(ldexp(round(ldexp(c_im, step_log)), -step_log));
  updateCachedDepth();
}

void MandelModel::Position::move(int delta_x, int delta_y)
{
  //qDebug()<<"move ("<<delta_x<<","<<delta_y<<")";
  center_re_n.impl->add_double(-delta_x*step_size__);
  center_im_n.impl->add_double(+delta_y*step_size__);
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
    center_re_n.impl->add_double(center_x*(old_step_size-step_size__));
    center_im_n.impl->add_double(center_y*(old_step_size-step_size__));
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
    center_re_n.impl->add_double((center_x-adjust_x)*(old_step_size-step_size__)); //(old_step_size-step_size__)=(1-(1<<-inlog))*old_step_size
    center_im_n.impl->add_double((center_y-adjust_y)*(old_step_size-step_size__));
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
  MandelMath::number_any d_re_n(&center_re_n);
  d_re_n.impl->lshift_(step_log-15);
  d_re_n.impl->frac_pos();
  d_re_n.impl->lshift_(15);
  cached_center_re_mod=d_re_n.impl->toRound();

  MandelMath::number_any d_im_n(&center_im_n);
  d_im_n.impl->lshift_(step_log-15);
  d_im_n.impl->frac_pos();
  d_im_n.impl->lshift_(15);
  cached_center_im_mod=d_im_n.impl->toRound();
}

void MandelModel::Position::pixelXtoRE(int x, MandelMath::number_any *result)
{
  //return (x - imageWidth/2)*position.step_size+position.center_re;
  //should be already result->reinit(center_re_n.ntype());
  result->impl->assign(center_re_n.impl->store);
  result->impl->add_double(x*step_size__);
}

void MandelModel::Position::pixelYtoIM(int y, MandelMath::number_any *result)
{
  //return (y - imageHeight/2)*position.step_size+position.center_im;
  //should be already result->reinit(center_im_n.ntype());
  result->impl->assign(center_im_n.impl->store);
  result->impl->add_double(y*step_size__);
}



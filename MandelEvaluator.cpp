#include "MandelEvaluator.hpp"
#include <cmath>

/*
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
*/


MandelPoint::MandelPoint()
{
  //reset();
  //should be overwritten before read:
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::assign(MandelMath::number_worker *worker, const MandelPoint &src)
{
  worker->assign(&f_re, &src.f_re);
  worker->assign(&f_im, &src.f_im);
  state=src.state;
  iter=src.iter;
  worker->assign(&fc_c_re, &src.fc_c_re);
  worker->assign(&fc_c_im, &src.fc_c_im);
  worker->assign(&fz_c_mag, &src.fz_c_mag);
  lookper_startiter=src.lookper_startiter;
  lookper_prevGuess=src.lookper_prevGuess;
  worker->assign(&lookper_startf_re, &src.lookper_startf_re);
  worker->assign(&lookper_startf_im, &src.lookper_startf_im);
  worker->assign(&lookper_nearr_dist, &src.lookper_nearr_dist);
  worker->assign(&lookper_totalFzmag, &src.lookper_totalFzmag);
  near0iter=src.near0iter;
  worker->assign(&near0f_re, &src.near0f_re);
  worker->assign(&near0f_im, &src.near0f_im);
  period=src.period;
  worker->assign(&root_re, &src.root_re);
  worker->assign(&root_im, &src.root_im);
  exterior_hits=src.exterior_hits;
  exterior_avoids=src.exterior_avoids;
  interior=src.interior;
}

void MandelPoint::init(MandelMath::number_worker *worker)
{
  worker->init(&f_re, 0.0);
  worker->init(&f_im, 0.0);
  worker->init(&fc_c_re, 0.0);
  worker->init(&fc_c_im, 0.0);
  worker->init(&fz_c_mag, 1);
  worker->init(&near0f_re);
  worker->init(&near0f_im);
  worker->init(&root_re);
  worker->init(&root_im);
  worker->init(&lookper_startf_re);
  worker->init(&lookper_startf_im);
  worker->init(&lookper_nearr_dist);
  worker->init(&lookper_totalFzmag);
  //should be overwritten before read:
  state=State::stUnknown;
  iter=0;
}

void MandelPoint::zero(MandelMath::number_worker *worker, const MandelMath::number_store *c_re, const MandelMath::number_store *c_im)
{
  //worker->zero(&f_re, 0);
  //worker->zero(&f_im, 0);
  worker->assign(&f_re, c_re);
  worker->assign(&f_im, c_im);
  worker->zero(&fc_c_re, 1);
  worker->zero(&fc_c_im, 0);
  worker->zero(&fz_c_mag, 1);
  lookper_prevGuess=0;
  //lookper resets at first iter
  near0iter=1;
  worker->assign(&near0f_re, c_re);
  worker->assign(&near0f_im, c_im);
  period=0;
  worker->zero(&root_re, 0);
  worker->zero(&root_im, 0);

  state=State::stUnknown;
  iter=0;
  exterior_avoids=-1;
  exterior_hits=-1;
  interior=-1;
  /*
    real exterior:=0
    real interior:=0
    initwinding(c)
    complex interiorComplex:=0
    int period:=0
    complex root:=0
  */
}

void MandelPoint::cleanup(MandelMath::number_worker *worker)
{
  if (worker==nullptr)
    dbgPoint();
  else
  {
    worker->cleanup(&lookper_totalFzmag);
    worker->cleanup(&lookper_nearr_dist);
    worker->cleanup(&lookper_startf_im);
    worker->cleanup(&lookper_startf_re);
    worker->cleanup(&root_im);
    worker->cleanup(&root_re);
    worker->cleanup(&near0f_im);
    worker->cleanup(&near0f_re);
    worker->cleanup(&fc_c_im);
    worker->cleanup(&fc_c_re);
    worker->cleanup(&fz_c_mag);
    worker->cleanup(&f_im);
    worker->cleanup(&f_re);
  }
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
      data.f_re.as.doubl=zr;
      data.f_im.as.doubl=zi;
      return;
    };
    double tmp=zr*zr-zi*zi+cr;
    zi=2*zr*zi+ci;
    zr=tmp;
  }
  //data.state=MandelPoint::State::stMaxIter;
  data.iter=maxiter;
  data.f_re.as.doubl=zr;
  data.f_im.as.doubl=zi;
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
      data.f_re.as.ddouble_.dd->assign(zr);
      data.f_im.as.ddouble_.dd->assign(zi);
      return;
    };
    t.assign(r2); t.add(-i2.hi, -i2.lo); t.add(cr->hi, cr->lo); //double tmp=zr*zr-zi*zi+cr;
    zi.mul(2*zr.hi, 2*zr.lo); zi.add(ci->hi, ci->lo);
    zr.assign(t);
  }
  data.iter=maxiter;
  data.f_re.as.ddouble_.dd->assign(zr);
  data.f_im.as.ddouble_.dd->assign(zi);
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
    currentWorker->cleanup(&interior.fz_mag);
    currentWorker->cleanup(&interior.fz_im);
    currentWorker->cleanup(&interior.fz_re);
    currentWorker->cleanup(&interior.inte_abs);
    currentWorker->cleanup(&interior.inte_im);
    currentWorker->cleanup(&interior.inte_re);

    currentWorker->cleanup(&newt.fzzf_im);
    currentWorker->cleanup(&newt.fzzf_re);
    currentWorker->cleanup(&newt.newtX_im);
    currentWorker->cleanup(&newt.newtX_re);
    currentWorker->cleanup(&newt.laguX_im);
    currentWorker->cleanup(&newt.laguX_re);
    currentWorker->cleanup(&newt.laguG2_im);
    currentWorker->cleanup(&newt.laguG2_re);
    currentWorker->cleanup(&newt.laguG_im);
    currentWorker->cleanup(&newt.laguG_re);
    currentWorker->cleanup(&newt.laguH_im);
    currentWorker->cleanup(&newt.laguH_re);
    currentWorker->cleanup(&newt.fzfix_im);
    currentWorker->cleanup(&newt.fzfix_re);
    currentWorker->cleanup(&newt.tmp2);
    currentWorker->cleanup(&newt.tmp1_im);
    currentWorker->cleanup(&newt.tmp1_re);
    currentWorker->cleanup(&newt.fzz_r_im);
    currentWorker->cleanup(&newt.fzz_r_re);
    currentWorker->cleanup(&newt.fz_r_im);
    currentWorker->cleanup(&newt.fz_r_re);
    currentWorker->cleanup(&newt.f_r_im);
    currentWorker->cleanup(&newt.f_r_re);
    currentWorker->cleanup(&newt.bestr_im);
    currentWorker->cleanup(&newt.bestr_re);

    currentWorker->cleanup(&eval.lookper_nearr_im);
    currentWorker->cleanup(&eval.lookper_nearr_re);
    currentWorker->cleanup(&eval.near0fmag);
    currentWorker->cleanup(&eval.fz_r_im);
    currentWorker->cleanup(&eval.fz_r_re);

    currentWorker->cleanup(&currentParams.c_im);
    currentWorker->cleanup(&currentParams.c_re);
    currentData.cleanup(currentWorker);
  }
  if (worker)
  if (worker)
  {
    currentData.init(worker);
    worker->init(&currentParams.c_re);
    worker->init(&currentParams.c_im);

    worker->init(&eval.fz_r_re);
    worker->init(&eval.fz_r_im);
    worker->init(&eval.near0fmag);
    worker->init(&eval.lookper_nearr_re);
    worker->init(&eval.lookper_nearr_im);

    worker->init(&newt.bestr_re);
    worker->init(&newt.bestr_im);
    worker->init(&newt.f_r_re);
    worker->init(&newt.f_r_im);
    worker->init(&newt.fz_r_re);
    worker->init(&newt.fz_r_im);
    worker->init(&newt.fzz_r_re);
    worker->init(&newt.fzz_r_im);
    worker->init(&newt.tmp1_re);
    worker->init(&newt.tmp1_im);
    worker->init(&newt.tmp2);
    worker->init(&newt.fzfix_re);
    worker->init(&newt.fzfix_im);
    worker->init(&newt.laguH_re);
    worker->init(&newt.laguH_im);
    worker->init(&newt.laguG_re);
    worker->init(&newt.laguG_im);
    worker->init(&newt.laguG2_re);
    worker->init(&newt.laguG2_im);
    worker->init(&newt.laguX_re);
    worker->init(&newt.laguX_im);
    worker->init(&newt.newtX_re);
    worker->init(&newt.newtX_im);
    worker->init(&newt.fzzf_re);
    worker->init(&newt.fzzf_im);

    worker->init(&interior.inte_re);
    worker->init(&interior.inte_im);
    worker->init(&interior.inte_abs);
    worker->init(&interior.fz_re);
    worker->init(&interior.fz_im);
    worker->init(&interior.fz_mag);
  }
  currentWorker=worker;
}

#if !COMPLEX_IS_TEMPLATE
bool MandelEvaluator::startCompute(const MandelPoint *data, int quick_route)
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
  if ((quick_route==1) ||
      ((quick_route==0) && (currentParams.maxiter-currentData.iter<=1000)))
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

//result 0..derivatives or value too large, or other fail (divide by 0)
//result>0 .. tried to return multiplicity but really returns just 1 (1 or >=3) or 2 (mult==2)
int MandelEvaluator::newton(int period, const complex *c, complex *r, const bool fastHoming, const int suggestedMultiplicity) //returns multiplicity
{ //TODO: suggestedMulti = maximumMultip ?
  double bestfm=1e10; //TODO: actually bestgm? g(z)=f(z)-z
  currentWorker->assign(&newt.bestr_re, r->re_s); //init, cleanup
  currentWorker->assign(&newt.bestr_im, r->im_s);
  bool movedOff=false;
  double accyBound=3e-28/(period*period);
  double accyBound2=3e-39*period/log(1+period)*1.5; //1.5=magic
  double order1; // 1/highest power in the polynomial, 2^period in case of mandelbrot set
  if (period<1024)
    order1=ldexp(1, -period);
  else
    order1=0;
  //int multiplicity1=1;
  int trustedMultiplicity=1;
  bool triedZeroGzrm=false;
  struct
  {
    bool didfix;
    int mult;
  } clever; //improve accuracy around point where 2 bulbs touch
  clever.didfix=false;
  if (suggestedMultiplicity>1)
    clever.mult=suggestedMultiplicity;
  else
    clever.mult=1;
  complex f_r(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  complex fz_r(currentWorker, &newt.fz_r_re, &newt.fz_r_im, true);
  complex fzz_r(currentWorker, &newt.fzz_r_re, &newt.fzz_r_im, true);
  complex tmp1(currentWorker, &newt.tmp1_re, &newt.tmp1_im, true);
  complex fzfix(currentWorker, &newt.fzfix_re, &newt.fzfix_im, true);
  complex laguH(currentWorker, &newt.laguH_re, &newt.laguH_im, true);
  complex laguG(currentWorker, &newt.laguG_re, &newt.laguG_im, true); //fzf
  complex laguG2(currentWorker, &newt.laguG2_re, &newt.laguG2_im, true); //G^2
  complex laguX(currentWorker, &newt.laguX_re, &newt.laguX_im, true); //Laguerre step
  complex newtX(currentWorker, &newt.newtX_re, &newt.newtX_im, true); //Laguerre step
  complex fzzf(currentWorker, &newt.fzzf_re, &newt.fzzf_im, true);
  for (int newtonCycle=0; newtonCycle<50; newtonCycle++)
  {
    if ((movedOff) && (newtonCycle>10) && (order1>=0))
    {                                    //  p m -> p
      order1=-1;                         //  2 2    1
      bestfm=1e10;                       //  4 3    2
      //multiplicity1=1;                   //  4 5    1
    };
    if (trustedMultiplicity>2) //tried some things and gave up
      trustedMultiplicity=1;
    currentWorker->assign(&newt.f_r_re, r->re_s);
    currentWorker->assign(&newt.f_r_im, r->im_s);
    currentWorker->zero(&newt.fz_r_re, 1.0);
    currentWorker->zero(&newt.fz_r_im, 0);
    currentWorker->zero(&newt.fzz_r_re, 0);
    currentWorker->zero(&newt.fzz_r_im, 0);
    //TODO: can we skip computing fzz_r if order1<0? and remember last valid multiplicity or set it to 1
    for (int i=0; i<period; i++)
    {
      double fzz_r_mag=currentWorker->toDouble(fzz_r.getMagTmp());
      double fz_r_mag=currentWorker->toDouble(fz_r.getMagTmp());
      double f_r_mag=currentWorker->toDouble(f_r.getMagTmp());
      if (fzz_r_mag+fz_r_mag+f_r_mag>1e60)
        return 0;
      //fzz:=2*(fz*fz + f*fzz)
      fzz_r.mul(&f_r);
      currentWorker->assign(tmp1.re_s, fz_r.re_s);
      currentWorker->assign(tmp1.im_s, fz_r.im_s);
      tmp1.sqr();
      fzz_r.add(&tmp1);
      currentWorker->lshift(tmp1.re_s, 1);
      currentWorker->lshift(tmp1.im_s, 1);
      //fz:=2*f*fz
      fz_r.mul(&f_r);
      currentWorker->lshift(fz_r.re_s, 1);
      currentWorker->lshift(fz_r.im_s, 1);
      //f:=f^2+c
      f_r.sqr();
      f_r.add(c);
    }
    double fz_r_mag=currentWorker->toDouble(fz_r.getMagTmp()); //ff1m in original code
    //g(r)=f(r)-r, gz(r)=fz(r)-1
    currentWorker->sub(f_r.re_s, r->re_s);
    currentWorker->sub(f_r.im_s, r->im_s);
    currentWorker->add_double(fz_r.re_s, -1);
    double g_r_mag=currentWorker->toDouble(f_r.getMagTmp());
    double gz_r_mag=currentWorker->toDouble(fz_r.getMagTmp());
#if CLEVER_FIX
//c=-0.7499 p=2
//  ideally, r=-0.5+-0.01i who are repelling
//  but we have f(-0.5001)=-0.49979999, f^2(-0.5001)=-0.5000999699959999
//  due to rounding errors, it looks as if we are at a root
//  and this point is attracting, so we have verified a false double period
//  there's really no way around this using finite precision
//  so we need something CLEVER
    if (!clever.didfix &&
        (g_r_mag<1e-16) &&
        (((gz_r_mag<5e-3) && !currentWorker->isle0(fz_r.re_s)) || //in bulb close to its base and at the wrong root
         (gz_r_mag<1e-9) || //so close to the base we don't know which root we have
         ((period==2) && !currentWorker->isle0(fz_r.re_s) && (currentWorker->toDouble(fz_r.re_s)<0.14)))) //we skip check for period=1 so special check for the point of attachment of bulb 1/2
    {
      clever.didfix=true;
      cleverFix(clever.mult);
      continue;
    };
#endif
    if (((trustedMultiplicity<=2) && (g_r_mag==0)) ||  //7e-33..4e-40 does not need more; much..5e-38 needs more
        ((trustedMultiplicity>2) && (gz_r_mag<7e-38))) //6.2e-38 cannot improve due to rounding errors
    { //r is good enough already      (f_c.re*f_c.re+f_c.im*f_c.im)/(f_zc.re*f_zc.re+f_zc.im*f_zc.im)
      return trustedMultiplicity;
    };
    if (gz_r_mag==0)
    {
      if (triedZeroGzrm)
        return 0;
      triedZeroGzrm=true;
    }
    else if (trustedMultiplicity<=2)
    {
      if (gz_r_mag<1e-39)
      { //one more multi-iter?
        //zero ffm is solved with triedZeroFfm
        if (clever.mult>0)
          return clever.mult;
        else
          return 1;
      };
      if (((gz_r_mag>1e-6) && (g_r_mag/gz_r_mag<accyBound)) || //((fm<NEWTON_EPSILON) and (ffm>1e-6)) or //check just f on single roots
          ((g_r_mag<1e-20) && (g_r_mag/gz_r_mag<accyBound)) || //for high period - and high derivative - we cannot minimize f given the finite precision of root
           ((gz_r_mag<1e-6) && (g_r_mag/gz_r_mag<2.1e-30)))   //check f/ff on multiple roots
      { //r is good enough already
        return trustedMultiplicity;
      };
      if ((g_r_mag<2e-38) && (gz_r_mag<3e-9)) //ffm: >8e-10
      { //close to be but not multiple
        if (clever.mult>0)
          return clever.mult;
        else
          return 1;
      };
      if ((fz_r_mag>0.5) && (g_r_mag/fz_r_mag<accyBound2))
      { //cannot improve because of limited precision
        return trustedMultiplicity;
      };
    };
    if (currentWorker->isequal(r->re_s, &newt.bestr_re) &&
        currentWorker->isequal(r->im_s, &newt.bestr_im) &&
        (bestfm<1e10))
    { //Laguerre can cycle (at least my version), e.g. per=2, c=-0.6640625-0.015625i, r=-0.614213552325963974-0,0179149806499481201i
      return 0; //just fail and try again next time
    };
    if (g_r_mag<bestfm)
    {
      bestfm=g_r_mag;
      currentWorker->assign(&newt.bestr_re, r->re_s);
      currentWorker->assign(&newt.bestr_im, r->im_s);
    };

    /* derive Laguerre's method, multiplicity m!=1, order of poly=n
    assume roots A and B, distance a=z-A away, others b=z-B away
    f(z)=C (z-A)^m (z-B)^(n-m)
    take ln, diff twice
    ln f(z) = ln C + m*ln(z-A) + (n-m)*ln(z-B)
    d/dz ln f(z) = d/dz ln C + d/dz m*ln a + d/dz (n-m)*ln b
    G = f'(z)/f(z) = m/a + (n-m)/b
    d/dz (f'(z)/f(z)) = d/dz m/a + d/dz (n-m)/b     a=z-A, (1/a)'=-1/a^2
    -H = (f''(z)*f(z) - f'(z)*f'(z)) / f(z)^2 = - m/a^2 - (n-m)/b^2
    G = m/a + (n-m)/b
    H = m/a^2 + (n-m)/b^2
    solve for aa=1/a from G=f'/f, H=G^2-f''/f   X=f'^2/(f''*f) = G^2/H
    (G-m*aa)/(n-m)=1/b
    H*(n-m) = m*(n-m)*aa^2 + (G-m*aa)^2
    0 = m*n*aa^2 - 2*m*G*aa + G^2-H*(n-m)
    0 = aa^2 - 2*G/n*aa + (G^2-H*(n-m))/(m*n)
    aa*n = G +- sqrt( (n/m-1)*(n*H-G^2) )
    Newton's step is f/f' = 1/G
    Laguerre's step is a=1/aa=n/(n*aa)=n/(G +- sqrt( (n/m-1)*(n*H-G^2) ))
    1/G/(1/n +- sqrt( (1/m-1/n)*(H/G^2-1/n) ))
    H/G^2=1-f''*f/f'^2 = 1/M
>>  a=f/f'/(1/n +- sqrt( (1/m-1/n)*(1/M-1/n) ))
      but fails if f'=0
    a=f/ffix'
    ffix'=f'/n +- sqrt( (1/m-1/n)*(f'^2-f''*f-f'^2/n) )         M=f'^2/(f'^2-f''*f)
      only fails if f'=f''=0

    when f'=0, G=0, M=inf; ideal m=1
    G=f'/f   H=G^2-f''/f   M=1/(1-f''*f/f'^2)   Re(M) rounds to 1 except for f''*f/f'^2 in circle (c=2/3 r=1/3)
>>  a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))    M=G^2/(G^2-f''/f) (G can be 0, can't divide)
    x^2+1 @ 0: f=1 f'=0 f''=2 G=0 H=-2
    +-i=1/(+- sqrt( (1-2/m) ))  m=1

    let's not talk about the derivation but
    either f'^2-f''*f=a+bi    f'^2=c+di
    or     H=f'^2/f^2-f''/f=a+bi    G^2=f'^2/f^2=c+di
    then   m=Round( (a*c+b*d)/(a*a+b*b) ) = Re(G^2*H^T)/mag(H) = Re(G^2/H)   x^T = conj(x)



    for m=n
    G = m/a
    a=n/G
    same as above even though the derivation is invalid

    when m:=M, a=f/f'/(1-f''*f/f'^2)=f/f'*M

    when laguH=0; ideal m=1
    a=1/(G/n +- sqrt( (1/m-1/n)*(-G^2/n) ))
    a=1/G/(1/n +- i*sqrt( (1/m-1/n)*(1/n) )) = n/G/(1 +- i*sqrt(n/m-1))
    x^2+x+1/2  f(0)=1/2 f'=1 f''=2  1*1-2/2=0 G=2 f''/f=4     -1/2+-i/2
    a=1/(1 +- sqrt( (1-2/m) ))     m=1
    -1/2+-i/2=1/(1 +- i*sqrt(2/m-1))
    1=m

    (x-1)^3*(x+1)^2  : H=0 at x = -1/5 +- (2 i sqrt(6))/5
    ReplaceAll [{(x - 1)^3*(x + 1)^2, D[(x - 1)^3*(x + 1)^2, {x}], D[(x - 1)^3*(x + 1)^2, {x, 2}]}, {x->-1/5 + (2*I*Sqrt[6])/5}]
    {-17856/3125 + (2112 i sqrt(6))/3125, 864/125 + (672 i sqrt(6))/125, 144/5 - (48 i sqrt(6))/5}
    f*f''=(-17856/3125 + (2112 i sqrt(6))/3125)*(144/5 - (48 i sqrt(6))/5) = -1963008/15625 + (1161216 i sqrt(6))/15625
    f'^2 = -1963008/15625 + (1161216 i sqrt(6))/15625
    ->a=4/5+(2 I sqrt(6))/5 = 5/(864/125 + (672 i sqrt(6))/125)/(1 +- I*sqrt(5/m-1))
    solve 4/5+(2*I*Sqrt(6))/5 = 5*(-17856/3125 + (2112 i sqrt(6))/3125)/(864/125 + (672*I*Sqrt(6))/125)/(1 + I*Sqrt(5/x-1))
      m=x=2, no solutions for 1-I*Sqrt


    when f''=0, f'!=0
    a=f/f'/(1/n +- sqrt( (1/m-1/n)*(1-1/n) ))
    x+1=0 f''=0 f'=1 f=1 G=1
    1=1/1/(1/1 +- sqrt( (1/m-1/1)*(1-1/1) ))    any m


    ----- how to find m ----
    simple: f=x^m
    f=x^m f'=m x^(m-1)  f''=m(m-1) x^(m-2)
    f/x^(m-2)=x^2  f'/x^(m-2)=m x  f''/x^(m-2)=m(m-1)
    f''*f/f'^2=m(m-1) x^(m-2) x^m / m/m / x^(m-1)/ x^(m-1) = (m-1)/m
    f'/f=m/x    f''/f=m(m-1)/x^2
    f''/f/f'*f=f''/f'=(m-1)/x
    1/m=1-f''*f/f'^2

    (x-1)^3*(x+1)^2 at 0.99
    f=-3.9601×10^-6  f'=0.00118405  f''=-0.23522
    f''*f/f'^2=0.6644  1/(1-...)=2.98   1/m=1-f''*f/f'^2

    full:
    f=(x-a)^m(x-b)^(n-m)
    f'=(x-a)^m(x-b)^(n-m)= m (x-a)^(m-1) (x-b)^(n-m) + (n-m) (x-a)^m (x-b)^(n-m-1)
    f''=(m-1) m (x-a)^(m-2) (x-b)^(n-m) + 2 m (n-m) (x-a)^(m-1) (x-b)^(n-m-1) + (n-m-1) (n-m) (x-a)^m (x-b)^(n-m-2)
      limit of 1/(1-D[(x-a)^m(x-b)^(n-m),{x,2}]*(x-a)^m(x-b)^(n-m)/D[(x-a)^m(x-b)^(n-m),{x,1}]^2) as b goes to infinity
        ->m
    w.l.o.g. x=0   Z=a/b
      M=1/(1-f''*f/f'^2)=(b m + a (n-m))^2/(b^2 m + a^2 (n-m))=(m + Z (n-m))^2/(m + Z^2 (n-m)) ~ m + 2*Z (n-m)
      from afar (a~b): m=n

    find b (bb) from Lagu:
    G = m/a + (n-m)/b
    H = m/a^2 + (n-m)/b^2
    solve for bb=1/b from G=f'/f, H=G^2-f''/f
    (G-(n-m)*bb)^2/m = m/a^2
    H = m/a^2 + (n-m)*bb^2
    H = (G-(n-m)*bb)^2/m + (n-m)*bb^2
    m*H = G*G-2*G*(n-m)*bb+(n-m)*(n)*bb*bb
    (n-m)*(n)*bb*bb-2*G*(n-m)*bb+(G*G-m*H) = 0
    bb1,2 = G/n*(1+-sqrt(1-1/(n-m)*(n)*(1-m*H/G^2)))
    b1,2 = n/G/(1+-sqrt(1-1/(n-m)*(n)*(1-m*H/G^2)))   H/G^2=1/M
    b1,2 = n/G/(1+-sqrt(1+n*(m/M-1)/(n-m))
    a = 1/G/(1/n + sqrt( (1/m-1/n)*(1/M-1/n) ))
    a/b=(1/n+-sqrt(1/n^2+(m/M-1)/n/(n-m))/(1/n + sqrt( (1/m-1/n)*(1/M-1/n) ))
    a/b=(1+-sqrt(1+(m/M-1)*n/(n-m))/(1 + sqrt( (n/m-1)*(n/M-1) ))
    a/b=(1+-sqrt(1+(m/M-1)*n/(n-m)))/(1 + sqrt( 1+((n-m)/M-1)*n/m ))
    m ~ (M-2*n*a/b)/(1-2*a/b)  for small a/b
    a/b=(1-sqrt(1+(m/M-1)*n/(n-m)))/(1 + sqrt( 1+((n-m)/M-1)*n/m ))
    a/b~-(m/M-1)*n/(n-m)/2/(1 + sqrt( 1+((n-m)/M-1)*n/m )) = -(m/M-1)*n/(n-m)/2/Q
    m=M+(m/M-1)*n/(1 + sqrt( 1+((n-m)/M-1)*n/m ))
      -> m:=M
      not gonna work
    a/b @ m=1: =(1-sqrt(1+(1/M-1)*n/(n-1)))/(1 + sqrt( 1+((n-1)/M-1)*n/1 ))
      m ~ (M-2*n*a/b)/(1-2*a/b)
      -> funny function that only depends on M and n; for n=5 (x=M):
        plot (x(1 + Sqrt[-(((x - 5) (-1 + 5))/x)]) + 2 5 (-1 + Sqrt[(x - 5)/(x - x 5)]))/(-1 + Sqrt[-(((x - 5) (-1 + 5))/x)] + 2 Sqrt[(x - 5)/(x - x 5)])




*/

    //1/f should be fine, or we'd be at the root
    currentWorker->assign(tmp1.re_s, f_r.re_s);
    currentWorker->assign(tmp1.im_s, f_r.im_s);
    tmp1.recip();    //1/f
    currentWorker->assign(laguG.re_s, fz_r.re_s);
    currentWorker->assign(laguG.im_s, fz_r.im_s);
    laguG.mul(&tmp1); //laguG = f'/f
    currentWorker->assign(fzzf.re_s, fzz_r.re_s);
    currentWorker->assign(fzzf.im_s, fzz_r.im_s);
    fzzf.mul(&tmp1); //f''/f

    //old algo:
    //Laguerre computed until cycle 10, using multi
    //if mag(f')==0 then step:=lagu
    //  else if mag(f)<1e-10 && mag(f')<1e-6 then trustedm:=m; step:=f/f'*m
    //         else step:=f/f'
    //which can obviously crash if we skip lagu and mag(f')=0

    //if ((trustedMultiplicity<=2) && (order1>=0))
    //{ // fzf=f'/f   fzzf=f''/f
      // laguH=fzf^2-fzzf
      // m=Round( Re(G^2*H^T)/mag(H) )
      // a=1/(G/n +- sqrt( (1/m-1/n)*(G^2-fzzf-G^2/n) ))
      currentWorker->assign(laguG2.re_s, laguG.re_s);
      currentWorker->assign(laguG2.im_s, laguG.im_s);
      laguG2.sqr();    //G^2
      currentWorker->assign(laguH.re_s, laguG2.re_s);
      currentWorker->assign(laguH.im_s, laguG2.im_s);
      currentWorker->sub(laguH.re_s, fzzf.re_s);
      currentWorker->sub(laguH.im_s, fzzf.im_s); //G^2-fzzf = f'^2/f^2-f''*f/f^2 = laguH
      //currentWorker->assign(tmp1.re_s, laguG2.re_s);
      //currentWorker->assign(tmp1.im_s, laguG2.im_s);
      double G2HT_re=currentWorker->toDouble(laguG2.mulreT(&laguH));
      double H_mag=currentWorker->toDouble(laguH.getMagTmp());
      //1.5*mag(H)>Re(G^2*H^T) ... m=1
      //300*mag(H)<Re(G^2*H^T) ... m=300
      //mag(H)<Re(G^2*H^T)*order1 ... m=1/order1
      int m=1;
      if (1.5*H_mag>=G2HT_re) //>= so we don't divide G^2/H if H=G=0
        m=1;
      else if ((clever.mult+0.5)*H_mag<=G2HT_re)
        m=1; //best practice is to use m=1 if H=0   clever.mult;
      else if (H_mag<G2HT_re*order1)
        m=qRound(1/order1);
      else
        m=qRound(G2HT_re/H_mag);
    if ((g_r_mag<1e-10) && (gz_r_mag<1e-6) && (trustedMultiplicity<=2))
    {
      trustedMultiplicity=m;
    };
    bool lagu_valid=false;
    bool newt_valid=false;
    if ((trustedMultiplicity<=2) && (order1>=0))
    {
      // a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
      // all but last few cycles can be done just in double precision
      //   but the cost of this compared to evaluation of f,f',f'' is negligible
      currentWorker->assign(laguX.re_s, laguG2.re_s);
      currentWorker->assign(laguX.im_s, laguG2.im_s);
      currentWorker->lshift(laguX.re_s, -period);
      currentWorker->lshift(laguX.im_s, -period); //G^2/n
      currentWorker->rsub(laguX.re_s, laguH.re_s);
      currentWorker->rsub(laguX.im_s, laguH.im_s); //H-G^2/n
      currentWorker->zero(&newt.tmp2, m);
      currentWorker->recip(&newt.tmp2);
      currentWorker->add_double(&newt.tmp2, -order1); //1/m-1/n
      currentWorker->mul(laguX.re_s, &newt.tmp2);
      currentWorker->mul(laguX.im_s, &newt.tmp2); //(1/m-1/n)*(H-G^2/n)
      laguX.sqrt();
      //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
      if (currentWorker->isle0(laguX.mulreT(&laguG))) //again isl0 would be nicer
      {
        currentWorker->chs(laguX.re_s);
        currentWorker->chs(laguX.im_s);
      };
      currentWorker->lshift(laguG.re_s, -period);
      currentWorker->lshift(laguG.im_s, -period); //G/n
      laguX.add(&laguG);
      //if 1/n~0: a=1/(0 +- sqrt( (1/m)*(H) )), m can still be 1..max
      //   fine if H!=0:       a=1/( sqrt( (1/m)*(H) )), m can still be 1..max
      //   if H==0: 1/G/(1/n + sqrt( (1/300-1/n)*(-1/n) ))=1/G* -i*sqrt(300*n)
      //   if H=G=0: 1/0
      //if G=0: a=1/(+- sqrt( (1/m-1/n)*(H) ))     m=1
      //   fine if H!=0: a=+-(sqrt(n/(n-1))*sqrt(f/-f''))       x^2+9 at 0: f=9 f''=2 -> +-3i
      //   if H=0: a=1/0
      //if H=0: a=1/G*m*(1 - i*sqrt(n/m-1))  m~n -> a=n/G;  m~300 -> a=-i/G*sqrt(n*300)
      //        a=1/G*m*n*(1/n - i*sqrt(1/m/n-1/n^2))
      double X_mag=currentWorker->toDouble(laguX.getMagTmp());
      if (X_mag>=1e-60)
      {
        laguX.recip_prepared();
        lagu_valid=true;
      };
      //else
      //we should move the guess a little and try again, but
      //  we can leave this to the caller
      //return 0;
    };
    if (fz_r_mag!=0)
    {
      //newton near multiroot:
      //f=(x-a)^m   f'=m*(x-a)^(m-1)  f/f'=(x-a)/m
      //Newton corrected for multiroot = f/f'*m
      //1/M=1-f''*f/f'^2   x-a=1/(f'/f-f''/f')
      currentWorker->assign(newtX.re_s, fz_r.re_s);
      currentWorker->assign(newtX.im_s, fz_r.im_s);
      newtX.recip();
      newtX.mul(&f_r); //f/f'
      currentWorker->zero(&newt.tmp2, trustedMultiplicity);
      currentWorker->mul(newtX.re_s, &newt.tmp2);
      currentWorker->mul(newtX.im_s, &newt.tmp2);
      newt_valid=true;
    };
    if (!newt_valid)
    {
      if (!lagu_valid)
      {
        return 0;
      };
      if (!movedOff)
      {
        movedOff=true;
        currentWorker->assign(newtX.re_s, laguX.re_s);
        currentWorker->assign(newtX.im_s, laguX.im_s);
      }
      else
        return 0;
    }
    else if (!lagu_valid)
    {
    }
    else
    {
      if (fastHoming && (newtonCycle<2) && (m>1))
      {
        currentWorker->assign(newtX.re_s, laguX.re_s);
        currentWorker->assign(newtX.im_s, laguX.im_s);
      }
      else
      {//take the smaller of newton and laguerre to a) avoid Newton's jumps when moving off the real axis, b) avoid Laguerre skipping to far root
        double N_mag=currentWorker->toDouble(newtX.getMagTmp());
        double L_mag=currentWorker->toDouble(laguX.getMagTmp());
        if (N_mag>L_mag)
        {
          currentWorker->assign(newtX.re_s, laguX.re_s);
          currentWorker->assign(newtX.im_s, laguX.im_s);
        };
      }
    }
    if ((g_r_mag>bestfm) && (newtonCycle>30))
    {
      currentWorker->lshift(newtX.re_s, -2);
      currentWorker->lshift(newtX.im_s, -2);
    };
    currentWorker->sub(r->re_s, newtX.re_s);
    currentWorker->sub(r->im_s, newtX.im_s);
  } //for newtonCycle
  return trustedMultiplicity;
}

//result=0 means the period check failed; -1 means the check failed and the root returned is invalid
int MandelEvaluator::periodCheck(int period/*must =eval.lookper_lastGuess*/, const complex *c)
{
  if (period<1)
  {
    dbgPoint();
    return -1;
  };
  int aroundCount; //estimate multiplicity (mult-1)
  if ((currentData.lookper_prevGuess>0) &&
      ((eval.lookper_lastGuess % currentData.lookper_prevGuess)==0))
    aroundCount=eval.lookper_lastGuess / currentData.lookper_prevGuess;
  else
    aroundCount=0;
  //look for root nearest to C - better stability of newton/laguerre
  MandelMath::complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
  currentWorker->assign(&currentData.root_re, &eval.lookper_nearr_re);
  currentWorker->assign(&currentData.root_im, &eval.lookper_nearr_im);
  //root.sqr();
  //root.add(c);

  /*checked before call to periodCheck if (currentWorker->toDouble(&eval.lookper_totalFzmag)>=MAGIC_MIN_SHRINK) //correct totalFZmag?
  { //TODO: is this correct? we're not evaluating at root, just some point around here...
    return -1;
  };*/

  if (period>MAX_PERIOD)
  {
    return -1; //special case
  };
  if (aroundCount==0)
    aroundCount=1; //a fresh nearestIteration means this is a new atom, so mult=2
  int newtRes=newton(period, c, &root, true, 1+aroundCount);
  if (newtRes<=0)
  { //this, of course, means that Newton() should be improved, not that there's a problem with the numbers
    return -1; //e.g. evaluating the initial guess mand.root leads to overflow immediately
  };

  complex f_r(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  complex fz_r(currentWorker, &newt.fz_r_re, &newt.fz_r_im, true);
  currentWorker->assign(f_r.re_s, root.re_s);
  currentWorker->assign(f_r.im_s, root.im_s);
  currentWorker->zero(fz_r.re_s, 1);
  currentWorker->zero(fz_r.im_s, 0);
  const MandelMath::number_store *fz_mag=nullptr;
  aroundCount=0;
  bool someBelow1=false;
  //int firstBelow1=-1;
  int firstBelow1dividing=-1;
  for (int i=0; i<period; i++)
  {
    if (fz_mag && currentWorker->isl1(fz_mag)) //I think we intentionally skip last f_z_r
    {
      someBelow1=true;
      //if (firstBelow1<0)
        //firstBelow1=i;
      if (firstBelow1dividing<0)
        if ((period % i)==0)
        {
          //needs more checks than that, e.g. fz_mag^(period/i) <=~ final fz_mag
          //per-actual  per-found  root-found   |   status at short
          //  long        short         x           won't happen because |fz|<1
          //  short       short       short         |f-r|<eps
          //  short       short        long         does not solve newton
          //  short        long       short         |f-r|<eps
          //  short        long        long         |f-r| big, |fz|<1 ... I think   !!
          //  long         long       short         |fz|>1
          //  long         long        long         |fz|>1
          //if (currentWorker->isle(f_r.getMagTmp(), root.getMagTmp()))
          double dist=currentWorker->toDouble(f_r.dist2_tmp(&root));
          if (dist<3.4e-28) //related to newton's accyBound=3e-28/period^2
            firstBelow1dividing=i; //short long short
          //but how to identify short long long ?
        };
    };
    //fm could hardly be >4 since it's tested in Newton (as well as f_z_r, BTW)
    //ff:=2*f*ff
    fz_r.mul(&f_r);
    currentWorker->lshift(fz_r.re_s, 1);
    currentWorker->lshift(fz_r.im_s, 1);
    fz_mag=fz_r.getMagTmp();
    if (currentWorker->toDouble(fz_mag)>LARGE_FLOAT2)
      return -1; //so is it checked or not
    //f:=f^2+c
    //use EXACTLY same code as the main loop: if temporaries are more precise at one place, it will never work!
    f_r.sqr();
    f_r.add(c);
  };
  //what's going on here -.-  normally we would check |ff|<1 or something
  //maybe that's what I'm doing but with some extra tricks for high precision?
  if (currentWorker->isle0(fz_r.re_s)) currentWorker->chs(fz_r.re_s);
  if (currentWorker->isle0(fz_r.im_s)) currentWorker->chs(fz_r.im_s);
  bool fz_r_mag_over1;
  if (!currentWorker->isle(fz_r.re_s, fz_r.im_s))
  {
    currentWorker->assign(&newt.tmp2, fz_r.re_s);
    currentWorker->assign(fz_r.re_s, fz_r.im_s);
    currentWorker->assign(fz_r.im_s, &newt.tmp2);
  };
  //im>=re
  if ((currentWorker->toDouble(fz_r.re_s)>=0.71) ||
      (currentWorker->toDouble(fz_r.im_s)>=1.01))
  { //im>=re>=0.71
    fz_r_mag_over1=true;
  }
  else if (currentWorker->toDouble(fz_r.im_s)<0.70)
    fz_r_mag_over1=false;
  else
  { //   re*re+im*im>1
    //   gets inaccurate for small re, im~1
    //   im-1>sqrt(1-re*re)-1
    //>> x=Sqrt(1-r*r)-1
    //   could try x:=((-r2/16-1/8)*r2-1/2)*r2; //-r^6/16-r^4/8-r^2/2
    //   but one cycle of Newton should work with any precision
    //   (x+1)*(x+1)=(1-r*r)
    //   x*x+2*x+r*r=0  f'=2*x+2
    //   x2=x-(x*x+2*x+r*r)/(2*x+2)
    //>> x2=(x*x-r*r)/(2*x+2)
    //   x2=(x-1)(r*r-x*x)/(1-x*x)/2
    currentWorker->assign(&newt.tmp2, fz_r.re_s);
    currentWorker->sqr(&newt.tmp2);
    currentWorker->assign(&newt.tmp1_im, &newt.tmp2); //r^2
    currentWorker->chs(&newt.tmp2);
    currentWorker->add_double(&newt.tmp2, 1);
    currentWorker->sqrt(&newt.tmp2);
    currentWorker->assign(&newt.tmp1_re, &newt.tmp2); //x+1
    currentWorker->add_double(&newt.tmp2, -1);        //x
    currentWorker->sqr(&newt.tmp2);
    currentWorker->sub(&newt.tmp2, &newt.tmp1_im);    //x^2-r^2
    currentWorker->recip(&newt.tmp1_re);
    currentWorker->mul(&newt.tmp2, &newt.tmp1_re);
    currentWorker->lshift(&newt.tmp2, -1);           //better x
    currentWorker->add_double(fz_r.im_s, -1);
    fz_r_mag_over1=!currentWorker->isle(fz_r.im_s, &newt.tmp2); //im-1>x
  };
  if (fz_r_mag_over1)
    return -1;//inevitably result:=-1
//  if (ff.re*ff.re+ff.im*ff.im)<=1 then
//  if ((f_z_r-1)<=0{1.1e-19}) then //Newton does not solve exactly (at least for multiple roots)

  //evaluate F_c^period at r ; its abs must be below 1 for the point to attract

  if (!someBelow1)
    return period; //seems to work
  //if (firstBelow1!=firstBelow1dividing)
    //dbgPoint();
  if (firstBelow1dividing<1)
    return period;
  else
    return firstBelow1dividing;
  /* lemme try to skip this nonsense, using firstBelow1 instead
  //the optimal cycle point is in mand.root
  complexOps.assign(f, mand.root);
  ffm:=1;
  SetLength(derivatives, period);
  //but...we just did this
  for i:=1 to period do //todo: find a way without the array
    begin
      //cmul(ff, f, ff); ff.re:=2*ff.re; ff.im:=2*ff.im;
      ffm:=4*ffm*complexOps.mag(f);
      //fg:=fg^2+c
      complexOps.aaMbP(f, f, mand.c);
      derivatives[i-1]:=ffm-1;//(ff.re*ff.re+ff.im*ff.im);  // f.re-r.re  f.im-r.im
    end;

  tmpre:=ffm; //=f_z_r
  found:=period;
  for prime:=1 to (period div 2) do
    begin
      if ((period mod prime)=0) and       //period/prime
         (derivatives[prime-1]<0) then
        begin
          ok:=True;
          for i:=prime+1 to period do
            begin
              if (derivatives[i-1]>derivatives[i-prime-1]) then
                begin
                  ok:=False;
                  break;
                end;
            end;
          if (ok) then
            begin
              found:=prime;
              break;
            end;
        end;
    end;
  if not someBelow1 and (found<>period) then
    found:=found; //looks like if someBelow1 then we always succeed decreasing period
  return found;*/
}

int MandelEvaluator::estimateInterior(int period, const complex *c, const complex *root)//, InteriorInfo *interior)
{
  // Initial values:  f = r;  fc = 0;  fz = 1;  fzz = 0;  fzc = 0;
  complex f(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  currentWorker->assign(f.re_s, root->re_s); //z   z^2+c  (z^2+c)^2+c
  currentWorker->assign(f.im_s, root->im_s);
  complex fz(currentWorker, &interior.fz_re, &interior.fz_im, true);
  currentWorker->zero(fz.re_s, 1);           //1   2z     4z^3+4cz=2*2z(z^2+c)
  currentWorker->zero(fz.im_s, 0);
  complex fc(currentWorker, &newt.laguG_re, &newt.laguG_im, true);
  currentWorker->zero(fc.re_s, 0);           //0   1      2z^2+2c+1=2*1*(z^2+c)+1
  currentWorker->zero(fc.im_s, 0);
  complex fzz(currentWorker, &newt.fzz_r_re, &newt.fzz_r_im, true);
  currentWorker->zero(fzz.re_s, 0);          //0   2      12z^2+4c=2(4z^2+2z^2+2c)=2*((2z)^2+(z^2+c)*2)
  currentWorker->zero(fzz.im_s, 0);
  complex fzc(currentWorker, &newt.laguG2_re, &newt.laguG2_im, true);
  currentWorker->zero(fzc.re_s, 0);          //0   0      4z=2*(2z*1+(z^2+c)0)
  currentWorker->zero(fzc.im_s, 0);
  complex fcc(currentWorker, &newt.laguH_re, &newt.laguH_im, true);
  currentWorker->zero(fcc.re_s, 0);
  currentWorker->zero(fcc.im_s, 0);
  currentWorker->zero(&interior.fz_mag, 1);
  complex tmp1(currentWorker, &newt.tmp1_re, &newt.tmp1_im, true);
  complex inte(currentWorker, &interior.inte_re, &interior.inte_im, true);
  // Start iterating.
  for (int i=0; i<period; i++)
  {
    //f^2+c -> 2f fc+1 -> 2f fcc+2fc fc
    // fcc := 2 * (fc^2 + f * fcc);
    currentWorker->assign(tmp1.re_s, fc.re_s);
    currentWorker->assign(tmp1.im_s, fc.im_s);
    tmp1.sqr();
    fcc.mul(&f);
    fcc.add(&tmp1);
    currentWorker->lshift(fcc.re_s, 1);
    currentWorker->lshift(fcc.im_s, 1);

    //f^2+c -> 2f fz -> 2fc fz+2f fzc
    // fzc := 2 * (fz * fc + f * fzc);
    currentWorker->assign(tmp1.re_s, fc.re_s);
    currentWorker->assign(tmp1.im_s, fc.im_s);
    tmp1.mul(&fz);
    fzc.mul(&f);
    fzc.add(&tmp1);
    currentWorker->lshift(fzc.re_s, 1);
    currentWorker->lshift(fzc.im_s, 1);

    //f^2+c -> 2f fz -> 2fz fz+2f fzz
    // fzz := 2 * (fz^2 + f * fzz);
    currentWorker->assign(tmp1.re_s, fz.re_s);
    currentWorker->assign(tmp1.im_s, fz.im_s);
    tmp1.sqr();
    fzz.mul(&f);
    fzz.add(&tmp1);
    currentWorker->lshift(fzz.re_s, 1);
    currentWorker->lshift(fzz.im_s, 1);

    //f^2+c -> 2f fc+1
    // fc := 2 * f * fc + 1;
    fc.mul(&f);
    currentWorker->lshift(fc.re_s, 1);
    currentWorker->lshift(fc.im_s, 1);
    currentWorker->add_double(fc.re_s, 1);

    //f^2+c -> 2f fz
    // fz := 2 * f * fz;
    fz.mul(&f);
    currentWorker->lshift(fz.re_s, 1);
    currentWorker->lshift(fz.im_s, 1);
    currentWorker->mul(&interior.fz_mag, f.getMagTmp());
    currentWorker->lshift(&interior.fz_mag, 2); //fz_mag*=4*mag(f)

    // f = f^2 + c
    f.sqr();
    f.add(c);

    if ((currentWorker->toDouble(f.getMagTmp())>LARGE_FLOAT2) ||
        (currentWorker->toDouble(fz.getMagTmp())>LARGE_FLOAT2) ||
        (currentWorker->toDouble(fc.getMagTmp())>LARGE_FLOAT2) ||
        (currentWorker->toDouble(fzz.getMagTmp())>LARGE_FLOAT2) ||
        (currentWorker->toDouble(fzc.getMagTmp())>LARGE_FLOAT2))
    {
//          DoGlobalDebug('?? when??');
      currentWorker->zero(&interior.inte_abs, -1);
      /*Result:=0; //yes sometimes it does...
          res_derivatives.f_z.re:=0; res_derivatives.f_z.im:=0;
          res_derivatives.f_c:=res_derivatives.f_z;
          res_derivatives.f_zz:=res_derivatives.f_z;
          res_derivatives.f_zc:=res_derivatives.f_z;
          res_derivatives.f_cc:=res_derivatives.f_z;
          interior.re:=0; interior.im:=0;
          Exit;*/
      return -1;
    };
  }

  //imma gonna skippa another test here
  //  of derivatives<1 -> would refine period

  currentWorker->assign(&newt.tmp2, &interior.fz_mag);
  currentWorker->add_double(&newt.tmp2, -1);
  if (abs(currentWorker->toDouble(&newt.tmp2))<6e-18)
  { //parabolic point
    currentWorker->zero(&interior.inte_abs, 0);
    currentWorker->zero(&interior.inte_re, 0);
    currentWorker->zero(&interior.inte_im, 0);
    return -1;
  };
  //                    1-|fz|^2          .    1-|fz|^2=(1-|fz|)(1+|fz|)
  // interior=  -----------------------   .
  //            | fzc + fzz fc/(1-fz) |   .
  currentWorker->chs(&newt.tmp2);
  currentWorker->assign(&interior.inte_re, &newt.tmp2);
  currentWorker->zero(&interior.inte_im, 0); //1-|fz|^2
  currentWorker->assign(tmp1.re_s, fz.re_s);
  currentWorker->assign(tmp1.im_s, fz.im_s);
  currentWorker->chs(tmp1.re_s);
  currentWorker->chs(tmp1.im_s);
  currentWorker->add_double(tmp1.re_s, 1); //1-fz
  //skip this step because it's just not right: interior:=newt.tmp2 * tmp1/|tmp1|
  tmp1.recip();
  tmp1.mul(&fc);
  tmp1.mul(&fzz);
  tmp1.add(&fzc);
  if (currentWorker->is0(tmp1.re_s) && currentWorker->is0(tmp1.im_s))
  { //probably wrong period, should not happen
    currentWorker->zero(&interior.inte_abs, 5);
    currentWorker->zero(&interior.inte_re, 5);
    currentWorker->zero(&interior.inte_im, 0);
    return period;
  };
  tmp1.recip();
  inte.mul(&tmp1);
  if (currentWorker->isle0(&newt.tmp2)) //should be isl0 again
  {
    currentWorker->assign(&interior.inte_abs, inte.getMagTmp());
    currentWorker->sqrt(&interior.inte_abs);
    currentWorker->chs(&interior.inte_abs);
  }
  else
  {
    currentWorker->assign(&interior.inte_abs, inte.getMagTmp());
    currentWorker->sqrt(&interior.inte_abs);
  }
  return period;
}

void MandelEvaluator::eval_until_bailout(complex *c, complex *f, complex *fc_c)
{
  for (int i=0; i<100; i++) //should be enough to reach 10000^2 except around (-2, 0)
  {
    const MandelMath::number_store *f_mag=f->getMagTmp();
    if (currentWorker->toDouble(f_mag)>1e8)
      return;
    //fc_c:=2*f*fc_c+1
    fc_c->mul(f);
    currentWorker->lshift(fc_c->re_s, 1);
    currentWorker->lshift(fc_c->im_s, 1);
    currentWorker->add_double(fc_c->re_s, 1);
    const MandelMath::number_store *fc_c_mag=fc_c->getMagTmp();
    if (currentWorker->toDouble(fc_c_mag)>LARGE_FLOAT2)
    {
      currentData.state=MandelPoint::State::stBoundary;
      return;
    };
    //f:=f^2+c
    f->sqr();
    f->add(c);
    currentData.iter++;
  };
}

void MandelEvaluator::evaluate()
{
  MandelMath::complex c(currentWorker, &currentParams.c_re, &currentParams.c_im, true);
  MandelMath::complex f(currentWorker, &currentData.f_re, &currentData.f_im, true);
  //currentData.fz_c_mag
  MandelMath::complex fc_c(currentWorker, &currentData.fc_c_re, &currentData.fc_c_im, true);
  //not needed for math MandelMath::complex fz_r(currentWorker, &eval.fz_r_re, &eval.fz_r_im, true);
  //currentWorker->zero(&eval.fz_r_re, 0);
  //currentWorker->zero(&eval.fz_r_im, 0);
  //currentData.near0iter
  { //near0f not needed for math, just a store
    MandelMath::complex near0f(currentWorker, &currentData.near0f_re, &currentData.near0f_im, true);
    currentWorker->assign(&eval.near0fmag, near0f.getMagTmp());
  }

  MandelMath::complex lookper_startf(currentWorker, &currentData.lookper_startf_re, &currentData.lookper_startf_im, true);
  MandelMath::complex lookper_nearr(currentWorker, &eval.lookper_nearr_re, &eval.lookper_nearr_im, true);

  for (int iter=currentData.iter; iter<currentParams.maxiter; iter++)
  {
    if (iter%(3*currentData.near0iter) ==0)
    {
      int quot=iter/(3*currentData.near0iter);
      if ((quot&(quot-1))==0) //also at iter==0
      { // //need k*iter for f' to start at the worst moment to reduce false positives; need k*iter-1 for good near0 -> switch to nearc
        currentData.lookper_startiter=iter;
        currentWorker->assign(&currentData.lookper_startf_re, f.re_s);
        currentWorker->assign(&currentData.lookper_startf_im, f.im_s);
        //MandelMath::complex lookper_bestf(currentWorker, &eval.lookper_startf_re, &eval.lookper_startf_im, true);
        //currentWorker->zero(&eval.lookper_bestf_re, 0);
        //currentWorker->zero(&eval.lookper_bestf_im, 0);
        currentWorker->assign(&eval.lookper_nearr_re, f.re_s);
        currentWorker->assign(&eval.lookper_nearr_im, f.im_s);
        if (iter<=1)
          currentWorker->assign(&currentData.lookper_nearr_dist, f.getMagTmp());
        else
          currentWorker->assign(&currentData.lookper_nearr_dist, f.dist2_tmp(&c));
        //currentWorker->zero(&eval.lookper_dist2, 1e10); //4.0 should be enough
        //mands.period stays there
        //currentData.lookper_prevGuess=0; //TODO: used for anything?
        eval.lookper_lastGuess=0;
        currentWorker->zero(&currentData.lookper_totalFzmag, 1.0);
      };
    }
    const MandelMath::number_store *f_mag=f.getMagTmp();
    if (currentWorker->toDouble(f_mag)>4)
    {
      currentData.state=MandelPoint::State::stOutside;
      currentData.iter=iter;
      //theory says the relative error in estimate is less than 3/bailout for large bailout
      //so lets move out a bit
      eval_until_bailout(&c, &f, &fc_c); //may switch state to stBoundary
      if (currentData.state!=MandelPoint::State::stOutside)
      {
        //currentWorker->zero(&currentData.exterior_avoids, 0);
        //currentWorker->zero(&currentData.exterior_hits, 0);
        currentData.exterior_avoids=0;
        currentData.exterior_hits=0;
      }
      else
      {
        //https://www.evl.uic.edu/hypercomplex/html/book/book.pdf p17, p29
        //G=ln(sqrt(f_mag))/2^iter   G'=sqrt(fc_c_mag)/(2^iter*sqrt(f_mag))       sinh(G)=(exp(G)-exp(-G))/2
        //sinh(G)/(2*exp(G)*G') < exterior < 2*sinh(G)/G'
        //(1-exp(-2G))/(4*G') < exterior < (exp(G)-exp(-G))/G'
        //for high iter (small G) we can use exp(x)=1+x, exp(-x)=1-x, sinh(G)=(1+G-(1-G))/2=G
        //  (2*G)/(4*G') < exterior < (2*G)/G'
        //  ln(f_mag)*sqrt(f_mag/fc_c_mag)/4 < exterior < ln(f_mag)*sqrt(f_mag/fc_c_mag)
        //  E/4 < exterior < E
        //otherwise, 1/G'=sqrt(f_mag/fc_c_mag)*2^iter
        //  (1-exp(-2*ln(sqrt(f_mag))/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter/4 < exterior < (exp(ln(sqrt(f_mag))/2^iter)-exp(-ln(sqrt(f_mag))/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter
        //  (1-1/exp(ln(f_mag)/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter/4 < exterior < (exp(ln(f_mag)/2/2^iter)-1/exp(ln(f_mag)/2/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter
        //  define X=ln(f_mag)/2/2^iter
        //  (1-1/exp(X)^2)*sqrt(f_mag/fc_c_mag)*2^iter/4 < exterior < (exp(X)-1/exp(X))*sqrt(f_mag/fc_c_mag)*2^iter
        //  (1-1/exp(X)^2)/X*E/8 < exterior < (exp(X)-1/exp(X))/X*E/2
        //  exp(X)=1+A*X 1/exp(X)=1-B*X 1/exp(X)^2=1-2*C*X, A,B,C~1
        //  C*ln(f_mag)*sqrt(f_mag/fc_c_mag)/4 < exterior < (A+B)/2*ln(f_mag)*sqrt(f_mag/fc_c_mag)
        //  C*E/4 < exterior < (A+B)/2*E
        //  A=(exp(X)-1)/X   B=(1-1/exp(X))/X   (A+B)/2=(exp(X)-1/exp(X))/X/2    C=(1-1/exp(X)^2)/X/2
        //  (A+B)/2 = 1 + x^2/6 + x^4/120 + x^6/5040 + x^8/362880 + x^10/39916800 + x^12/6227020800 + O(x^13)
        //  C = 1 - x + (2 x^2)/3 - x^3/3 + (2 x^4)/15 - (2 x^5)/45 + O(x^6)
        //  assuming f_mag<10000^2, approx up to x^2 should be accurate to 1e-20 with iter>26
        //  1 should be accurate to 1e-20 with iter>71
        double fm=currentWorker->toDouble(f.getMagTmp());
        double fcm=currentWorker->toDouble(fc_c.getMagTmp());
        double x=log(fm);
        currentData.exterior_hits=x*sqrt(fm/fcm);
        currentData.exterior_avoids=currentData.exterior_hits*0.25;
        if (currentData.iter>71)
        { }
        else
        {
          x=ldexp(x, -1-currentData.iter);
          if (currentData.iter>26)
          {
            currentData.exterior_hits+=x*x/6*currentData.exterior_hits;
            currentData.exterior_avoids+=x*(x*2/3-1)*currentData.exterior_avoids;
          }
          else
          {
            double ex=exp(x);
            currentData.exterior_hits*=(ex-1/ex)/x/2;
            currentData.exterior_avoids*=(1-1/(ex*ex))/x/2;
          }
        }
      }
      //already there currentWorker->assign(currentData.fc_c_re, &fc_c_re);
      //currentWorker->assign(currentData.fc_c_im, &fc_c_im);
      currentData.period=eval.lookper_lastGuess; //preliminary
      if (currentData.period<1)
        currentData.period=1;
      return;
    };
    double fc_c_mag=currentWorker->toDouble(fc_c.getMagTmp());
    if (fc_c_mag>1e57)
    {
      currentData.state=MandelPoint::State::stBoundary;
      currentData.iter=iter;
      currentData.exterior_avoids=0;
      currentData.exterior_hits=0;
      return;
    };
    double fz_c_mag=currentWorker->toDouble(&currentData.fz_c_mag);
    if (fz_c_mag>1e60)
    {
      currentData.state=MandelPoint::State::stDiverge;
      currentData.iter=iter;
      currentData.exterior_avoids=0;
      currentData.exterior_hits=0;
      return;
    };
    //TODO: similar to eval_until_bailout
    //fc_c:=2*f*fc_c+1
    fc_c.mul(&f);
    currentWorker->lshift(fc_c.re_s, 1);
    currentWorker->lshift(fc_c.im_s, 1);
    currentWorker->add_double(fc_c.re_s, 1);
    /* TODO: copy test here from above?
    fc_c_mag=currentWorker->toDouble(fc_c.getMagTmp());
    if (fc_c_mag>LARGE_FLOAT2)
    {
      currentData.state=MandelPoint::State::stBoundary;
      return;
    };*/
    f_mag=f.getMagTmp();
    //f'=2*f'*f, f'_mag=4*f'_mag*f_mag
    currentWorker->mul(&currentData.fz_c_mag, f_mag); //TODO: can use f_mag from above?
    currentWorker->mul(&currentData.lookper_totalFzmag, f_mag);
    currentWorker->lshift(&currentData.lookper_totalFzmag, 2);
    //f:=f^2+c
    f.sqr();
    f.add(&c);
    //currentData.iter++;
    f_mag=f.getMagTmp();

    if (!currentWorker->isle(&eval.near0fmag, f_mag)) //f_mag<near0fmag
    {
      currentData.near0iter=iter+2;
      currentWorker->assign(&currentData.near0f_re, f.re_s);
      currentWorker->assign(&currentData.near0f_im, f.im_s);
      currentWorker->assign(&eval.near0fmag, f_mag);
    };

    const MandelMath::number_store *lpdiff=lookper_startf.dist2_tmp(&f);
    if (!currentWorker->isle(&currentData.lookper_nearr_dist, lpdiff)) //|f-r|<best
    {
      currentWorker->assign(&eval.lookper_nearr_re, f.re_s);
      currentWorker->assign(&eval.lookper_nearr_im, f.im_s);
      currentWorker->assign(&currentData.lookper_nearr_dist, lpdiff);
      currentData.lookper_prevGuess=eval.lookper_lastGuess;
      eval.lookper_lastGuess=(iter+1-currentData.lookper_startiter);
    };

    if ((eval.lookper_lastGuess>0) &&
        (eval.lookper_lastGuess==(iter+1-currentData.lookper_startiter)) && //just found new guess
         (//actually period<=near0iter ((eval.lookper_lastGuess % currentData.near0iter)==0) ||  //either nearest divides period, or
          ((currentData.near0iter % eval.lookper_lastGuess)==0)) && //  period divides nearest, that's a fact
         ((iter>=3*currentData.near0iter)))  //speedup - don't check period too eagerly
    {
      int foundperiod=-1;
      if (f.isequal(&lookper_startf))
      { //exact match - misiurewicz or converged after too many steps
        foundperiod=eval.lookper_lastGuess;
        currentWorker->assign(&currentData.root_re, f.re_s);
        currentWorker->assign(&currentData.root_im, f.im_s);
        //TODO: still needs period cleanup... I think. Near 0+0I
      }
      else if (currentWorker->toDouble(&currentData.lookper_totalFzmag)<MAGIC_MIN_SHRINK)
      {
        //if (periodCheckHistory<>'') then
        //  periodCheckHistory:=periodCheckHistory+'('+IntToStr(mande.iter)+':'+IntToStr(mande.lookper_lastGuess)+') ';
        //periodEntered:=getTime64();
        //if (currentWorker->toDouble(c.re_s)==-0.015625 && currentWorker->toDouble(c.im_s)==0.640625)
          //dbgPoint();
        foundperiod=periodCheck(eval.lookper_lastGuess, &c); //updates iter, f, f_c, root
        //profiler.timeInPeriod:=profiler.timeInPeriod+(getTime64()-periodEntered);
        if ((foundperiod>0) && (foundperiod<eval.lookper_lastGuess))
        {
          complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
          //complex interiorComplex(currentWorker, &interior.inte_re, &interior.inte_im, true);
          foundperiod=estimateInterior(foundperiod, &c, &root); //-4.7e-22
            //foundperiod=-1; //the cycle can be exact but |f_z|>1 due to mistaken period or misplaced (rounding err) root
        }
      };
      if (foundperiod>0)
      {
        currentData.state=MandelPoint::State::stPeriod2;
        currentData.iter=iter;
        currentData.period=foundperiod;
        complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
        currentData.period=estimateInterior(foundperiod, &c, &root);
        if (currentWorker->isle0(&interior.inte_abs)) //wanted <0 here
          currentData.state=MandelPoint::State::stMisiur;
        else
        {
          currentData.interior=currentWorker->toDouble(&interior.inte_abs);
          currentWorker->assign(&currentData.fc_c_re, &interior.fz_re); //d/dz F_c(r)
          currentWorker->assign(&currentData.fc_c_im, &interior.fz_im);
          if (foundperiod!=currentData.period)
            currentData.state=MandelPoint::State::stPeriod3;
        }
        //currentWorker->assign(&currentData.fc_c_re, &eval.fz_r_re);
        //currentWorker->assign(&currentData.fc_c_im, &eval.fz_r_im);
        return;
      };
    };

  }
  //data.state=MandelPoint::State::stMaxIter;
  currentData.iter=currentParams.maxiter;
}

#else
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
    switch (currentWorker->ntype())
    {
      case MandelMath::number_worker::Type::typeDouble:
        evaluate<MandelMath::number_worker_double>();
        break;
      case MandelMath::number_worker::Type::typeDDouble:
        evaluate<MandelMath::number_worker_ddouble>();
        break;
      case MandelMath::number_worker::Type::typeMulti:
        evaluate<MandelMath::number_worker_multi>();
        break;
      case MandelMath::number_worker::Type::typeEmpty:
        currentData.state=MandelPoint::State::stMaxIter;
        break;
    }
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

void MandelEvaluator::doCompute()
{
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  //simple_double(currentParams.cr_n.impl->store->as.doubl, currentParams.ci_n.impl->store->as.doubl, currentData, currentParams.maxiter);
  switch (currentWorker->ntype())
  {
    case MandelMath::number_worker::Type::typeDouble:
      evaluate<MandelMath::number_worker_double>();
      break;
    case MandelMath::number_worker::Type::typeDDouble:
      evaluate<MandelMath::number_worker_ddouble>();
      break;
    case MandelMath::number_worker::Type::typeMulti:
      evaluate<MandelMath::number_worker_multi>();
      break;
    case MandelMath::number_worker::Type::typeEmpty:
      currentData.state=MandelPoint::State::stMaxIter;
      break;
  }
  pointsComputed++;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}

template <class NW>
void MandelEvaluator::evaluate()
{
  NW localWorker;
  MandelMath::complex<NW> c(&currentParams.cr_s, &currentParams.ci_s, true);
  MandelMath::complex<NW> z(&this->data_zr_s, &this->data_zi_s, true);
  localWorker.assign(z.re_s, &currentData.zr_);
  localWorker.assign(z.im_s, &currentData.zi_);
  for (int iter=currentData.iter; iter<currentParams.maxiter; iter++)
  {
    const MandelMath::number_store *magtmp=z.getMagTmp();
    if (localWorker.toDouble(magtmp)>4)
    {
      currentData.state=MandelPoint::State::stOutside;
      currentData.iter=iter;
      localWorker.assign(&currentData.zr_, z.re_s);
      localWorker.assign(&currentData.zi_, z.im_s);
      return;
    };
    z.sqr();
    z.add(&c);
  }
  //data.state=MandelPoint::State::stMaxIter;
  currentData.iter=currentParams.maxiter;
  localWorker.assign(&currentData.zr_, z.re_s);
  localWorker.assign(&currentData.zi_, z.im_s);
}

#endif

MandelEvaluator::ComputeParams::ComputeParams():
  c_re(),
  c_im()
{
  epoch=-1;
  pixelIndex=-1;
  maxiter=1;
}

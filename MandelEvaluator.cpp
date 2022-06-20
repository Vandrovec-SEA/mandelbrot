#include "MandelEvaluator.hpp"
#define _USE_MATH_DEFINES //some magic
#include <cmath>
//C++20 also don't have #include <numbers>
//has gcd since C++17 which I apparently don't have #include <numeric>

#define USE_GCD_FOR_CHECKPERIOD 0
#define CLEVER_FIX 0

#define assert(x) { if (!(x)) dbgPoint(); }

void nop()
{

}

LaguerrePointStore::LaguerrePointStore(): state(State::stUnknown), firstM(0), iter(0)
{

}

void LaguerrePointStore::assign(const LaguerrePointStore *src)
{
  assert(src!=nullptr);
  //TODO: if (src==nullptr) see LaguerrePoint::zero
  state=src->state;
  firstM=src->firstM;
  iter=src->iter;
}


LaguerrePoint::LaguerrePoint(LaguerrePointStore *store, MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN),
  store(store), f(&self_allocator), fz_r(&self_allocator)
{
}

void LaguerrePoint::assign(const LaguerrePoint &src)
{
  store->assign(src.store);
  f.assign(&src.f);
  fz_r.assign(&src.fz_r);
}

/*void LaguerrePoint::init(MandelMath::number_worker *worker)
{
  promote(MandelMath::number_worker::Type::typeEmpty, worker->ntype());
  worker->zero(&f_re, 0.0);
  worker->zero(&f_im, 0.0);
  worker->zero(&fz_r_re, 0.0);
  worker->zero(&fz_r_im, 0.0);

  //should be overwritten before read:
  state=State::stUnknown;
  iter=0;
  firstM=0;
}*/

void LaguerrePoint::zero(const MandelMath::complex *c)
{
  f.assign(c);
  fz_r.zero(1, 0);
  store->state=LaguerrePointStore::State::stUnknown;
  store->firstM=0;
  store->iter=0;
}

/*void LaguerrePoint::cleanup(MandelMath::number_worker *worker)
{
  if (worker==nullptr)
    dbgPoint();
  else
  {
    worker->cleanup(&fz_r_im);
    worker->cleanup(&fz_r_re);
    worker->cleanup(&f_im);
    worker->cleanup(&f_re);
  }
}

void LaguerrePoint::promote(MandelMath::number_worker::Type oldType, MandelMath::number_worker::Type newType)
{
  if (oldType==newType)
    return;
  void *old_place=place.dd;
  place.dd=nullptr;
  uint8_t *placement;
  int step;
  switch (newType)
  {
    case MandelMath::number_worker::Type::typeEmpty:
      placement=nullptr;
      step=0;
      break;
    case MandelMath::number_worker::Type::typeDouble:
      place.dd=nullptr;
      placement=nullptr;
      step=0;
      break;
    case MandelMath::number_worker::Type::typeDDouble:
      place.dd=reinterpret_cast<MandelMath::dd_real (*)[Place::LEN]> (new MandelMath::dd_real[Place::LEN]());
      placement=(uint8_t *)place.dd;
      step=sizeof(MandelMath::dd_real);
      break;
    case MandelMath::number_worker::Type::typeMulti:
      place.multi=reinterpret_cast<MandelMath::multiprec (*)[Place::LEN]> (new MandelMath::multiprec[Place::LEN]());
      placement=(uint8_t *)place.multi;
      step=sizeof(MandelMath::multiprec);
      break;
  }

  f_re.promote_(oldType, newType, placement); placement+=step;
  f_im.promote_(oldType, newType, placement); placement+=step;
  fz_r_re.promote_(oldType, newType, placement); placement+=step;
  fz_r_im.promote_(oldType, newType, placement); placement+=step;
  assert((size_t)(placement-(uint8_t *)place.dd)==step*Place::LEN);
  switch (oldType)
  {
    case MandelMath::number_worker::Type::typeEmpty:
      break;
    case MandelMath::number_worker::Type::typeDouble:
      break;
    case MandelMath::number_worker::Type::typeDDouble:
      delete[] (MandelMath::dd_real (*)[Place::LEN])old_place;
      break;
    case MandelMath::number_worker::Type::typeMulti:
      delete[] (MandelMath::multiprec (*)[Place::LEN])old_place;
      break;
  }
}*/

MandelPointStore::MandelPointStore(): state(State::stUnknown), iter(0)
{

}

void MandelPointStore::assign(const MandelPointStore *src)
{
  state=src->state;
  iter=src->iter;
  has_fc_r=src->has_fc_r;
  lookper_startiter=src->lookper_startiter;
  lookper_prevGuess_=src->lookper_prevGuess_;
  lookper_lastGuess=src->lookper_lastGuess;
  lookper_nearr_dist_touched=src->lookper_nearr_dist_touched;
  near0iter=src->near0iter;
  newton_iter=src->newton_iter;
  period=src->period;
  exterior_hits=src->exterior_hits;
  exterior_avoids=src->exterior_avoids;
  interior=src->interior;
}

MandelPoint::MandelPoint(MandelPointStore *store, MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN),
  store(store), f(&self_allocator), fc_c(&self_allocator), fz_r(&self_allocator), fz_c_mag(&self_allocator),
  lookper_startf(&self_allocator), lookper_nearr_dist(&self_allocator), lookper_totalFzmag(&self_allocator),
  near0f(&self_allocator), root(&self_allocator)
{
  //reset();
  //should be overwritten before read:
  assert(self_allocator.checkFill());
}

void MandelPoint::assign(const MandelPoint &src)
{
  store->assign(src.store);
  f.assign(&src.f);
  fc_c.assign(&src.fc_c);
  fz_r.assign(&src.fz_r);
  fz_c_mag.assign(src.fz_c_mag.ptr);
  lookper_startf.assign(&src.lookper_startf);
  lookper_nearr_dist.assign(src.lookper_nearr_dist.ptr);
  lookper_totalFzmag.assign(src.lookper_totalFzmag.ptr);
  near0f.assign(&src.near0f);
  root.assign(&src.root);
}

/*void MandelPoint::init(MandelMath::number_worker *worker)
{
  promote(MandelMath::number_worker::Type::typeEmpty, worker->ntype());
  worker->zero(&f_re, 0.0);
  worker->zero(&f_im, 0.0);
  worker->zero(&fc_c_re_, 0.0);
  worker->zero(&fc_c_im_, 0.0);
  worker->zero(&fz_r_re, 1.0);
  worker->zero(&fz_r_im, 0.0);
  worker->zero(&fz_c_mag_, 1);
  worker->zero(&near0f_re);
  worker->zero(&near0f_im);
  worker->zero(&root_re);
  worker->zero(&root_im);
  worker->zero(&lookper_startf_re);
  worker->zero(&lookper_startf_im);
  worker->zero(&lookper_nearr_dist_);
  worker->zero(&lookper_totalFzmag);

  //should be overwritten before read:
  state=State::stUnknown;
  iter=0;
  lookper_nearr_dist_touched=false;
  newton_iter=0;
  has_fc_r=false;
}*/

void MandelPoint::zero(const MandelMath::complex *c)
{
  //worker->zero(&f_re, 0);
  //worker->zero(&f_im, 0);
  f.assign(c);
  fc_c.zero(0, 0);
  fz_r.zero(1, 0);
  fz_c_mag.zero(1);
  store->lookper_prevGuess_=0;
  store->lookper_lastGuess=0;
  //lookper resets at first iter
  store->near0iter=1;
  near0f.assign(c);
  store->period=0;
  root.zero(0, 0);

  store->state=MandelPointStore::State::stUnknown;
  store->iter=0;
  store->newton_iter=0;
  store->exterior_avoids=-1;
  store->exterior_hits=-1;
  store->interior=-1;
  store->has_fc_r=false;
  /*
    real exterior:=0
    real interior:=0
    initwinding(c)
    complex interiorComplex:=0
    int period:=0
    complex root:=0
  */
}

/*void MandelPoint::cleanup(MandelMath::number_worker *worker)
{
  if (worker==nullptr)
    dbgPoint();
  else
  {
    worker->cleanup(&lookper_totalFzmag);
    worker->cleanup(&lookper_nearr_dist_);
    worker->cleanup(&lookper_startf_im);
    worker->cleanup(&lookper_startf_re);
    worker->cleanup(&root_im);
    worker->cleanup(&root_re);
    worker->cleanup(&near0f_im);
    worker->cleanup(&near0f_re);
    worker->cleanup(&fc_c_im_);
    worker->cleanup(&fc_c_re_);
    worker->cleanup(&fz_r_im);
    worker->cleanup(&fz_r_re);
    worker->cleanup(&fz_c_mag_);
    worker->cleanup(&f_im);
    worker->cleanup(&f_re);
  }
}

void MandelPoint::promote(MandelMath::number_worker::Type oldType, MandelMath::number_worker::Type newType)
{
  if (oldType==newType)
    return;

  void *old_place=place.dd;
  place.dd=nullptr;

  uint8_t *placement;
  int step;
  switch (newType)
  {
    case MandelMath::number_worker::Type::typeEmpty:
      placement=nullptr;
      step=0;
      break;
    case MandelMath::number_worker::Type::typeDouble:
      place.dd=nullptr;
      placement=nullptr;
      step=0;
      break;
    case MandelMath::number_worker::Type::typeDDouble:
      place.dd=reinterpret_cast<MandelMath::dd_real (*)[Place::LEN]> (new MandelMath::dd_real[Place::LEN]());
      placement=(uint8_t *)place.dd;
      step=sizeof(MandelMath::dd_real);
      break;
    case MandelMath::number_worker::Type::typeMulti:
      place.multi=reinterpret_cast<MandelMath::multiprec (*)[Place::LEN]> (new MandelMath::multiprec[Place::LEN]());
      placement=(uint8_t *)place.multi;
      step=sizeof(MandelMath::multiprec);
      break;
  }

  f_re.promote_(oldType, newType, placement); placement+=step;
  f_im.promote_(oldType, newType, placement); placement+=step;
  fc_c_re_.promote_(oldType, newType, placement); placement+=step;
  fc_c_im_.promote_(oldType, newType, placement); placement+=step;
  fz_r_re.promote_(oldType, newType, placement); placement+=step;
  fz_r_im.promote_(oldType, newType, placement); placement+=step;
  fz_c_mag_.promote_(oldType, newType, placement); placement+=step;
  lookper_startf_re.promote_(oldType, newType, placement); placement+=step;
  lookper_startf_im.promote_(oldType, newType, placement); placement+=step;
  lookper_nearr_dist_.promote_(oldType, newType, placement); placement+=step;
  lookper_totalFzmag.promote_(oldType, newType, placement); placement+=step;
  near0f_re.promote_(oldType, newType, placement); placement+=step;
  near0f_im.promote_(oldType, newType, placement); placement+=step;
  root_re.promote_(oldType, newType, placement); placement+=step;
  root_im.promote_(oldType, newType, placement); placement+=step;
  assert((size_t)(placement-(uint8_t *)place.dd)==step*Place::LEN);
  switch (oldType)
  {
    case MandelMath::number_worker::Type::typeEmpty:
      break;
    case MandelMath::number_worker::Type::typeDouble:
      break;
    case MandelMath::number_worker::Type::typeDDouble:
      delete[] (MandelMath::dd_real (*)[Place::LEN])old_place;
      break;
    case MandelMath::number_worker::Type::typeMulti:
      delete[] (MandelMath::multiprec (*)[Place::LEN])old_place;
      break;
  }
}*/

ShareableViewInfo::ShareableViewInfo(MandelMath::worker_multi::Allocator *allocator): QObject(),
  selfAllocator(allocator, 0, LEN, nullptr),
  originalAllocator(allocator),
  c(&selfAllocator), root(&selfAllocator), scale(1), period(0)
{
  assert(selfAllocator.checkFill());
}

ShareableViewInfo::ShareableViewInfo(ShareableViewInfo &src): QObject(),
  selfAllocator(src.originalAllocator, 0, LEN, nullptr),
  originalAllocator(src.originalAllocator),
  c(&selfAllocator), root(&selfAllocator),
  scale(src.scale), period(src.period)
{
}

ShareableViewInfo::ShareableViewInfo(const ShareableViewInfo &src): ShareableViewInfo((ShareableViewInfo &)src)
{
  //why do you need this?
}

ShareableViewInfo::ShareableViewInfo(ShareableViewInfo &&src): QObject(),
  selfAllocator(src.originalAllocator, 0, LEN, nullptr),
  originalAllocator(src.originalAllocator),
  c(&selfAllocator), root(&selfAllocator),
  scale(src.scale), period(src.period)
{
}

ShareableViewInfo &ShareableViewInfo::operator=(ShareableViewInfo &src)
{
  assert(originalAllocator==src.originalAllocator);
  scale=src.scale;
  period=src.period;
  assert(c.re.asf64==src.c.re.asf64);
  assert(c.im.asf64==src.c.im.asf64);
  assert(root.re.asf64==src.root.re.asf64);
  assert(root.im.asf64==src.root.im.asf64);
  /*c.assign(src.)
  originalAllocator->worker->assign()
  assert(c.re.asf64==nullptr);
  c.re=src.c.re;
  src.c.re.asf64=nullptr;
  assert(c.im.asf64==nullptr);
  c.im=src.c.im;
  src.c.im.asf64=nullptr;
  src.c.free_storage_on_destroy=false;
  assert(root.re.asf64==nullptr);
  root.re=src.root.re;
  src.root.re.asf64=nullptr;
  assert(root.im.asf64==nullptr);
  root.im=src.root.im;
  src.root.im.asf64=nullptr;
  src.root.free_storage_on_destroy=false;*/
  return *this;
}

ShareableViewInfo &ShareableViewInfo::operator=(ShareableViewInfo &&src)
{
  return operator=((ShareableViewInfo &)src);
}

LaguerreStep::LaguerreStep(MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN), currentWorker(allocator->worker),
  step(&self_allocator), s1(&self_allocator), s2(&self_allocator), tmp1(&self_allocator), tmp2(&self_allocator),
  laguG(&self_allocator), laguG2(&self_allocator), laguH(&self_allocator), laguX(&self_allocator), fzzf(&self_allocator)
{
  assert(self_allocator.checkFill());
}

/*void LaguerreStep::switchType(MandelMath::worker_multi *worker)
{
  if (worker==currentWorker)
    return;
  if (currentWorker!=nullptr)
  {
    currentWorker->cleanup(&tmp1_re);
    currentWorker->cleanup(&tmp1_im);
    currentWorker->cleanup(&tmp2);
    currentWorker->cleanup(&laguG_re);
    currentWorker->cleanup(&laguG_im);
    currentWorker->cleanup(&laguG2_re);
    currentWorker->cleanup(&laguG2_im);
    currentWorker->cleanup(&laguH_re);
    currentWorker->cleanup(&laguH_im);
    currentWorker->cleanup(&laguX_re);
    currentWorker->cleanup(&laguX_im);
    currentWorker->cleanup(&fzzf_re);
    currentWorker->cleanup(&fzzf_im);
    currentWorker->cleanup(&s1_re);
    currentWorker->cleanup(&s1_im);
    currentWorker->cleanup(&s2_re);
    currentWorker->cleanup(&s2_im);
    currentWorker->cleanup(&step_re);
    currentWorker->cleanup(&step_im);
    //currentWorker->cleanup(&z_re);
    //currentWorker->cleanup(&z_im);
    //currentWorker->cleanup(&c_re);
    //currentWorker->cleanup(&c_im);
  }

  if (worker)
  {
    uint8_t *placement;
    int step;
    switch (worker->ntype())
    {
      case MandelMath::number_worker::Type::typeEmpty:
        placement=nullptr;
        step=0;
        break;
      case MandelMath::number_worker::Type::typeDouble:
        place.dd=nullptr;
        placement=nullptr;
        step=0;
        break;
      case MandelMath::number_worker::Type::typeDDouble:
        place.dd=reinterpret_cast<MandelMath::dd_real (*)[Place::LEN]> (new MandelMath::dd_real[Place::LEN]());
        placement=(uint8_t *)place.dd;
        step=sizeof(MandelMath::dd_real);
        break;
      case MandelMath::number_worker::Type::typeMulti:
        place.multi=reinterpret_cast<MandelMath::multiprec (*)[Place::LEN]> (new MandelMath::multiprec[Place::LEN]());
        placement=(uint8_t *)place.multi;
        step=sizeof(MandelMath::multiprec);
        break;
    }

    //worker->init_(&z_re, placement); placement+=step;
    //worker->init_(&z_im, placement); placement+=step;
    //worker->init_(&c_re, placement); placement+=step;
    //worker->init_(&c_im, placement); placement+=step;
    worker->init_(&step_re, placement); placement+=step;
    worker->init_(&step_im, placement); placement+=step;
    worker->init_(&s1_re, placement); placement+=step;
    worker->init_(&s1_im, placement); placement+=step;
    worker->init_(&s2_re, placement); placement+=step;
    worker->init_(&s2_im, placement); placement+=step;
    worker->init_(&tmp1_re, placement); placement+=step;
    worker->init_(&tmp1_im, placement); placement+=step;
    worker->init_(&tmp2, placement); placement+=step;
    worker->init_(&laguG_re, placement); placement+=step;
    worker->init_(&laguG_im, placement); placement+=step;
    worker->init_(&laguG2_re, placement); placement+=step;
    worker->init_(&laguG2_im, placement); placement+=step;
    worker->init_(&laguH_re, placement); placement+=step;
    worker->init_(&laguH_im, placement); placement+=step;
    worker->init_(&laguX_re, placement); placement+=step;
    worker->init_(&laguX_im, placement); placement+=step;
    worker->init_(&fzzf_re, placement); placement+=step;
    worker->init_(&fzzf_im, placement); placement+=step;

    assert((size_t)(placement-(uint8_t *)place.dd)==step*Place::LEN);
  }
  currentWorker=worker;
}*/

bool LaguerreStep::eval(int lg2_degree, const MandelMath::complex *f, const MandelMath::complex *f_z, const MandelMath::complex *f_zz)
{
  if (currentWorker->is0(f->re) && currentWorker->is0(f->im))
  {
    step.zero(0, 0);
    return true;
  };

  double order1;
  int maxm;
  if (lg2_degree<5)
  {
    maxm=ldexp(1, lg2_degree-1); //in theory up to n-1 but for Mandelbrot that's rarely the case
    order1=ldexp(1, -lg2_degree);
  }
  else if (lg2_degree<1024)
  {
    maxm=15;
    order1=ldexp(1, -lg2_degree);
  }
  else
  {
    maxm=15;
    order1=0;
  }

  //1/f should be fine, or we'd be at the root
  tmp1.assign(f);
  tmp1.recip();    //1/f
  laguG.assign(f_z);
  laguG.mul(&tmp1); //laguG = f'/f
  fzzf.assign(f_zz);
  fzzf.mul(&tmp1); //f''/f


    // laguH=fzf^2-fzzf
    // m=Round( Re(G^2*H^T)/mag(H) )
    // a=1/(G/n +- sqrt( (1/m-1/n)*(G^2-fzzf-G^2/n) ))
    laguG2.assign(&laguG);
    laguG2.sqr();    //G^2
    laguH.assign(&laguG2);
    laguH.sub(&fzzf); //G^2-fzzf = f'^2/f^2-f''*f/f^2 = laguH
    //currentWorker->assign(tmp1.re_s, laguG2.re_s);
    //currentWorker->assign(tmp1.im_s, laguG2.im_s);
    int m=1;
    {
      /*double G2HT_re=currentWorker->toDouble(laguG2.mulreT(&laguH));
      double H_mag=currentWorker->toDouble(laguH.getMagTmp());
      //turns out that if mu=m then mu=m=G^2/H
      //1.5*mag(H)>Re(G^2*H^T) ... m=1
      //300*mag(H)<Re(G^2*H^T) ... m=300
      //mag(H)<Re(G^2*H^T)*order1 ... m=1/order1
      if (1.5*H_mag>=G2HT_re) //>= so we don't divide G^2/H if H=G=0
        m=1;
      else if ((clever.mult+0.5)*H_mag<=G2HT_re)
        m=1; //best practice is to use m=1 if H=0   clever.mult;
      else if (H_mag*(maxm-0.5)<G2HT_re)
        m=maxm;
      else
        m=qRound(G2HT_re/H_mag);*/

      //solve for m=mu:   m=G^2/H
      double mum_re=1, mum_im=0; //better than mu?
      double h_re=currentWorker->toDouble(laguH.re);
      double h_im=currentWorker->toDouble(laguH.im);
      double h_mag=h_re*h_re+h_im*h_im;
      double g2_re=currentWorker->toDouble(laguG2.re);
      double g2_im=currentWorker->toDouble(laguG2.im);
      if (h_mag>0.01)
      { //h_mag ok
        mum_re=(g2_re*h_re+g2_im*h_im)/h_mag;
        mum_im=(g2_im*h_re-g2_re*h_im)/h_mag;
      };
      dbg.mum_re=mum_re;
      dbg.mum_im=mum_im;
    }

    //m= some func of mu where mu is solution of ((1-1/n)*H/G^2-1/n) mu^2 + 2*mu/n -1=0
    //with m as input:                           ((1-m/n)*H/G^2-1/n) mu^2 + m/n 2*mu -m = 0
    double G2_mag=laguG2.getMag_double();
    if (G2_mag<0.01)
    { //G2_mag bad
      m=1;
      dbg.mu_re=1;
      dbg.mu_im=0;
    }
    else
    {
      laguX.assign(&laguG2);
      currentWorker->chs(laguX.im);
      laguX.mul(&laguH);
      double a_re=currentWorker->toDouble(laguX.re)/G2_mag*(1-order1)-order1;
      double a_im=currentWorker->toDouble(laguX.im)/G2_mag*(1-order1);
      double mu_re, mu_im;
      MandelMath::complex_double_quadratic(&mu_re, &mu_im, a_re, a_im, order1, 0, -1, 0);
      dbg.mu_re=mu_re;
      dbg.mu_im=mu_im;
      if (!(mu_re>=1.3)) //also m=1 if mu_re is NaN    (mu_re<1.3)
        m=1;
      else {/*if (abs(mu_im)>mu_re/2)
        m=1;
      else
      {
        double mu_mag=mu_re*mu_re+mu_im*mu_im;
        m=qRound(sqrt(mu_mag)); //or just round mu_re?
        */
        m=qRound(mu_re);
        if (m>maxm)
          m=maxm;
      }
    }
#if 0
    if (newtonCycle==0)
    {
      //Fejer bound: smaller solution x of
      //fzz/(n-1) x^2+2 fz x + n f=0
      //x=y*n
      //fzz*n/(n-1) y^2+2 fz y + f=0

      double r_re=currentWorker->toDouble(r->re_s);
      double r_im=currentWorker->toDouble(r->im_s);
      //numbers are small but don't need precision so let's do it in double
      double a_re=currentWorker->toDouble(fzz_r.re_s)/(1-order1);
      double a_im=currentWorker->toDouble(fzz_r.im_s)/(1-order1);
      double fz_re=currentWorker->toDouble(fz_r.re_s);
      double fz_im=currentWorker->toDouble(fz_r.im_s);
      double f_re=currentWorker->toDouble(f_r.re_s);
      double f_im=currentWorker->toDouble(f_r.im_s);
      MandelMath::complex_double_quadratic(
            &newtres_.first_fejer_re, &newtres_.first_fejer_im,
            a_re, a_im,
            fz_re, fz_im,
            f_re, f_im);
      newtres_.first_fejer_re=r_re+ldexp(newtres_.first_fejer_re, period);
      newtres_.first_fejer_im=r_im+ldexp(newtres_.first_fejer_im, period);

      //Batra's bound https://www.tuhh.de/ti3/paper/rump/Ru03c.pdf theorem 3.8
        //but only for real coefficients
      //|fz r|-|f + fzz/2 r^2|=0, find r
      //sqrt(fz fz^T) r=sqrt((f + fzz/2 r^2)(f^T + fzz^T/2 r^2))
      //sqrt(fz fz^T) r=sqrt((|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4))
      //(a+bi)(c-di)+(a-bi)(c+di)=2ac+2bd=2 Re(f fzz^T)
      //|fz|^2 r^2=|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4
      //0=|f|^2+ (Re(f^T fzz)-|fz|^2) rr + |fzz|^2/4 rr^2    r=sqrt(rr)

      /*MandelMath::complex_double_quadratic(&newtres_.first_batra, &a_im,
          currentWorker->toDouble(fzz_r.getMagTmp())/4, 0,
          (currentWorker->toDouble(f_r.mulreT(&fzz_r))-currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
          currentWorker->toDouble(f_r.getMagTmp()), 0);
      if (newtres_.first_batra>=0)
        newtres_.first_batra=sqrt(newtres_.first_batra);*/

      //https://ur.booksc.eu/book/5736333/a5b588
      //ZAMM - Journal of Applied Mathematics and Mechanics / Zeitschrift für Angewandte Mathematik und Mechanik
      //1988 Vol. 68; Iss. 6
      //Dr. A. Neumaier: An Existence Test for Root Clusters and Multiple Roots
      //fi(c, r, alpha)=r abs(re((f(c+r e^ialpha)-f(c))/(c+r e^ialpha)))-abs(f(c))
      //  addition from https://ur.booksc.eu/book/5736333/a5b588 remark 3:
      //  f needs to be divided (or rotated) by f' first to make f' real
      //for all alpha, which r makes fi==0 ?
      //abs(re(f'*r+f''/2 r^2 e^ialpha))=abs(f)
      //for max re(f'*r+f''/2 r^2 e^ialpha), we need max re(f'+f''/2 r e^ialpha) because r is real
      //f'' e^ialpha=real
      //e^ialpha=f''^T/sqrt(f'' f''^T)=sqrt(f''^T/f'')
      //abs(re(f'*r+ r^2/2 sqrt(f'' f''^T)))-abs(f)=0
      //r*abs(re(f'))+ r^2/2 sqrt(f'' f''^T)-abs(f)=0
      /*if (currentWorker->isle0(fz_r.re_s))
        MandelMath::complex_double_quadratic(&newtres_.first_batra, &a_im,
            +sqrt(currentWorker->toDouble(fzz_r.getMagTmp()))/2, 0,
            //currentWorker->toDouble(fz_r.re_s)/2, 0,
            sqrt(currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
            +sqrt(currentWorker->toDouble(f_r.getMagTmp())), 0);
      else*/
      MandelMath::complex_double_quadratic(&newtres_.first_neumaier1_re_, &newtres_.first_neumaier1_im_,
          -sqrt(currentWorker->toDouble(fzz_r.getMagTmp()))/2, 0,
          sqrt(currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
          -sqrt(currentWorker->toDouble(f_r.getMagTmp())), 0);

      /* Neumaier for k=2:
      Re(f''(z)/2) > |f| r^-2 + |f'| r^-1      r real, r>0
      f''(z)=f''+(z-z0)f'''+...
      |f''+r|f'''|+...|/2 r^2 > |f| + |f'| r
      ...+|f'''|r^3/2+|f''| r^2/2 > |f| + |f'|r
                      |f''| r^2/2 > |f| + |f'|r
      gives always r<0 but that's the wrong root
      |f''| r^2/2 - |f'|r - |f| =0
      2*|f'|/|f''|+-sqrt(4*|f'|^2/|f''|^2-8*|f|/|f''|)
      |f'|/|f''|+-sqrt(|f'|^2/|f''|^2+2*|f|/|f''|)
      works if |f'|^2/|f''|+2*|f|>0 i.e. always
      but r1 always<0, r2>2*|f'|/|f''|
      r2=(|f'|+sqrt(|f'|^2+2*|f|*|f''|))/|f''|

      test x(x-1) at 2+i
      f=1+3i f'=2x-1=3+2i f''=2
      r2=(|3+2*I|+sqrt(|3+2*I|^2+2*2*|1+3*I|))/2
      4.33502318885498454
      correct is 2.236
      */
      double fm=sqrt(currentWorker->toDouble(f_r.getMagTmp()));
      double fzm=sqrt(currentWorker->toDouble(fz_r.getMagTmp()));
      double fzzm=sqrt(currentWorker->toDouble(fzz_r.getMagTmp()));
      newtres_.first_neumaier2_re=(fzm + sqrt(fzm*fzm+2*fm*fzzm))/fzzm;
      newtres_.first_neumaier2_im=0;

      /* naive: approximate f with c(x-a)^m
      m=f'^2/(f'^2-f f'') = f'^2/f^2/(f'^2/f^2-f''/f)=G^2/H
      x-a=m/(f'/f)=m/G=G/H    looks good if |m_im|<|m_re|
      m*(x-a)=G^3/H^2

      trouble: singularities when f f''=f'^2 -> m=infinity, iteration jumps too far
                                  f'=0 -> m=0, m/(f/f') jumps too little
      */
      /*double g_re=currentWorker->toDouble(laguG.re_s);
      double g_im=currentWorker->toDouble(laguG.im_s);
      double g_mag=g_re*g_re+g_im*g_im;
      if (1e6*H_mag<=g_mag*g_mag)
      {
        newtres_.first_naive_re=currentWorker->toDouble(r->re_s);
        newtres_.first_naive_im=currentWorker->toDouble(r->im_s);
      }
      else
      {
        double g2_re=currentWorker->toDouble(laguG2.re_s);
        double g2_im=currentWorker->toDouble(laguG2.im_s);
        double h_re=currentWorker->toDouble(laguH.re_s);
        double h_im=currentWorker->toDouble(laguH.im_s);

        double m_re=(g2_re*h_re+g2_im*h_im)/H_mag;
        double m_im=(g2_im*h_re-g2_re*h_im)/H_mag;
        //couldn't find smooth function that:
        //1->1 2->2 3->3... 0->1 -1->1 i->1 -i->1
        //esp. since we need to have 1->1 exact and in neigborhood too
        if ((m_re<abs(m_im)*2))
        //if ((m_re<0.9) || (m_re<abs(m_im)*2)) //for m~0, we need something like sqrt(m): m is too small, 1 is too large
        {
          m_re=1;
          m_im=0;
        };
        newtres_.first_naive_re=currentWorker->toDouble(r->re_s)-(m_re*g_re+m_im*g_im)/g_mag;
        newtres_.first_naive_im=currentWorker->toDouble(r->im_s)-(m_im*g_re-m_re*g_im)/g_mag;
      }*/

      /* even naiver: show the 2 roots of c(x-a)(x-b) that have the same f, f', f''
      w.l.o.g. x=0
      c(x^2-(a+b)x+ab)=f''x^2/2+f'x+f
      cx^2-c(a+b)x+cab=f''x^2/2+f'x+f
      -f'/f''+-sqrt(f'^2/f''^2-2*f/f'')

      if x1 close to x2 (relative to x), use (x1+x2)/2 else use x1
      at |x1|=|x2|, 90 degrees..mult~2, use (x1+x2)/2
      at |x1|=|x2|, 60 degrees..mult~1, use x1
      at |x1|=0.8|x2|, 80% weight from x1
      at |x1|=0.5|x2|, 90% weight from x1
      at |x1|=0.3|x2|, use x1
      when x1~x2, correct guess is actually around 0.7 x1
      */
      a_re=currentWorker->toDouble(fzz_r.re_s)/2;
      a_im=currentWorker->toDouble(fzz_r.im_s)/2;
      MandelMath::complex_double_quadratic2(&newtres_.first_naive1_re_, &newtres_.first_naive1_im,
                                            &newtres_.first_naive2_re, &newtres_.first_naive2_im,
                                            a_re, a_im, fz_re/2, fz_im/2, f_re, f_im);
      double n2_rmag=1/(newtres_.first_naive2_re*newtres_.first_naive2_re+newtres_.first_naive2_im*newtres_.first_naive2_im);
      //d=naive1/naive2
      double d_re=(newtres_.first_naive1_re_*newtres_.first_naive2_re+newtres_.first_naive1_im*newtres_.first_naive2_im)*n2_rmag;
      double d_im=(newtres_.first_naive1_im*newtres_.first_naive2_re-newtres_.first_naive1_re_*newtres_.first_naive2_im)*n2_rmag;
      double d_mag=(d_re*d_re+d_im*d_im);
      double w1=1, w2=0;
      if (d_re<-0.5) //angle>120deg, even if close in magnitude
      { w1=1; w2=0; newtres_.naiveChoice=NewtonNaiveChoice::ncWide; }
      else if (d_mag<0.3*0.3)
      { w1=1; w2=0; newtres_.naiveChoice=NewtonNaiveChoice::nc03; }
      else if (d_mag<0.5*0.5)
      { w1=0.9; w2=0.1; newtres_.naiveChoice=NewtonNaiveChoice::nc05; }
      else if (d_mag<0.8*0.8)
      { w1=0.8; w2=0.2; newtres_.naiveChoice=NewtonNaiveChoice::nc08; } //or just 1;0
      else if (d_re<-0.1)
      {
        //most problematic case, some times best is [1,0] other times [1,1]
        w1=1; w2=0.0;
        newtres_.naiveChoice=NewtonNaiveChoice::nc100;
      }
      else if (d_re<0)
      {
        //most problematic case, some times best is [1,0] other times [1,1]
        //don't trust M here
        w1=1; w2=0.0;
        newtres_.naiveChoice=NewtonNaiveChoice::nc90_;
      }
      else if (d_re<0.1)
      {
        //most problematic case, some times best is [1,0] other times [1,1]
        //don't trust M here
        w1=1; w2=0.0;
        newtres_.naiveChoice=NewtonNaiveChoice::nc80;
      }
      else if (d_re<0.5)
      {
        //can (try) use M here
        w1=1; w2=0;
        newtres_.naiveChoice=NewtonNaiveChoice::nc60;
      }
      else
      {
        //can (try) use M here
        w1=newtres_.firstMum_re_-1; w2=0;
        newtres_.naiveChoice=NewtonNaiveChoice::ncClose;
      }
      newtres_.first_naive_re=w1*newtres_.first_naive1_re_+w2*newtres_.first_naive2_re;
      newtres_.first_naive_im=w1*newtres_.first_naive1_im+w2*newtres_.first_naive2_im;
      newtres_.first_naive1_re_=r_re+newtres_.first_naive1_re_;
      newtres_.first_naive1_im=r_im+newtres_.first_naive1_im;
      newtres_.first_naive2_re=r_re+newtres_.first_naive2_re;
      newtres_.first_naive2_im=r_im+newtres_.first_naive2_im;
      newtres_.first_naive_re=r_re+newtres_.first_naive_re;
      newtres_.first_naive_im=r_im+newtres_.first_naive_im;

      //Laguerre is the solution of
      //   c=-n  b=f'/f  a=f'^2/f^2*(1-n/m+1/m)-f''/f*(1-n/m)=H*(1-n/m)+G^2/m
      //   G=f'/f   H=G^2-f''/f
      //>> a=H*(1-m/n)-G^2/n  b=m*G/n  c=-m    ok
      /*
      a_re=currentWorker->toDouble(laguH.re_s)*(1-m*order1)-currentWorker->toDouble(laguG2.re_s)*order1;
      a_im=currentWorker->toDouble(laguH.im_s)*(1-m*order1)-currentWorker->toDouble(laguG2.im_s)*order1;
      double b_re=currentWorker->toDouble(laguG.re_s)*m*order1;
      double b_im=currentWorker->toDouble(laguG.im_s)*m*order1;
      MandelMath::complex_double_quadratic(
            &newtres_.first_lagum_re, &newtres_.first_lagum_im,
            a_re, a_im,
            b_re, b_im,
            -m, 0);
      newtres_.first_lagum_re=currentWorker->toDouble(r->re_s)-newtres_.first_lagum_re;
      newtres_.first_lagum_im=currentWorker->toDouble(r->im_s)-newtres_.first_lagum_im;
      */
      a_re=currentWorker->toDouble(laguH.re_s)*(1-order1)-currentWorker->toDouble(laguG2.re_s)*order1;
      a_im=currentWorker->toDouble(laguH.im_s)*(1-order1)-currentWorker->toDouble(laguG2.im_s)*order1;
      double b_re=currentWorker->toDouble(laguG.re_s)*order1;
      double b_im=currentWorker->toDouble(laguG.im_s)*order1;
      MandelMath::complex_double_quadratic2(
            &newtres_.first_lagu1_re, &newtres_.first_lagu1_im,
            &newtres_.first_lagu1o_re, &newtres_.first_lagu1o_im,
            a_re, a_im,
            b_re, b_im,
            -1, 0);
      newtres_.first_lagu1_re=r_re-newtres_.first_lagu1_re;
      newtres_.first_lagu1_im=r_im-newtres_.first_lagu1_im;
      newtres_.first_lagu1o_re=r_re-newtres_.first_lagu1o_re;
      newtres_.first_lagu1o_im=r_im-newtres_.first_lagu1o_im;
    };
#endif
  dbg.lastm=m;
  bool lagu_valid=false;
  bool newt_valid=false;
  if (order1>=0)
  {
    // a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
    // all but last few cycles can be done just in double precision
    //   but the cost of this compared to evaluation of f,f',f'' is negligible
    laguX.assign(&laguG2);
    laguX.lshift(-lg2_degree); //G^2/n
    laguX.rsub(&laguH); //H-G^2/n
    tmp2.zero(m);
    tmp2.recip();
    tmp2.add_double(-order1); //1/m-1/n
    currentWorker->mul(laguX.re, tmp2.ptr);
    currentWorker->mul(laguX.im, tmp2.ptr); //(1/m-1/n)*(H-G^2/n)
    laguX.sqrt();
    //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
    if (currentWorker->isle0(laguX.mulreT_tmp(&laguG))) //again isl0 would be nicer
    {
      currentWorker->chs(laguX.re);
      currentWorker->chs(laguX.im);
    };
    laguG.lshift(-lg2_degree); //G/n
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
    double X_mag=laguX.getMag_double();
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
  if (!currentWorker->is0(f_z->re) || !currentWorker->is0(f_z->im)) //gz_r_mag!=0)
  {
    //newton near multiroot:
    //f=(x-a)^m   f'=m*(x-a)^(m-1)  f/f'=(x-a)/m
    //Newton corrected for multiroot = f/f'*m
    //1/M=1-f''*f/f'^2   x-a=1/(f'/f-f''/f')
    step.assign(f_z);
    step.recip();
    step.mul(f); //f/f'
    if (m!=1)
    {
      tmp2.zero(m);
      currentWorker->mul(step.re, tmp2.ptr);
      currentWorker->mul(step.im, tmp2.ptr);
    };
    newt_valid=true;
  };
#if 0
  if (newtonCycle==0)
  {
    currentWorker->assign(&newtres_.first_guess_newt_re, r->re_s);
    currentWorker->assign(&newtres_.first_guess_newt_im, r->im_s);
    if (newt_valid)
    {
      currentWorker->sub(&newtres_.first_guess_newt_re, newtX.re_s);
      currentWorker->sub(&newtres_.first_guess_newt_im, newtX.im_s);
    };

    currentWorker->assign(&newtres_.first_guess_lagu_re, r->re_s);
    currentWorker->assign(&newtres_.first_guess_lagu_im, r->im_s);
    if (lagu_valid)
    {
      currentWorker->sub(&newtres_.first_guess_lagu_re, laguX.re_s);
      currentWorker->sub(&newtres_.first_guess_lagu_im, laguX.im_s);
    };
  };
#endif
  if (!newt_valid)
  {
    if (!lagu_valid)
    {
      return false;
    };
    step.assign(&laguX);
  }
  else if (!lagu_valid)
  {
    //keep newtX
  }
  else
  {
    if (m>1)//(fastHoming && (newtonCycle<2) && (m>1))
    {
      step.assign(&laguX);
    }
    else
    {//take the smaller of newton and laguerre to a) avoid Newton's jumps when moving off the real axis, b) avoid Laguerre skipping to far root
      double N_mag=step.getMag_double();
      double L_mag=laguX.getMag_double();
      if (N_mag*1.05>L_mag) //5% will do no harm, and switch to Lagu can speed up convergence
      {
        step.assign(&laguX);
      };
    }
  }
#if 0
  if ((g_r_mag>bestfm) && (newtonCycle>30))
  {
    currentWorker->lshift(newtX.re_s, -2);
    currentWorker->lshift(newtX.im_s, -2);
  };
#endif
  //currentWorker->sub(r->re_s, newtX.re_s);
  //currentWorker->sub(r->im_s, newtX.im_s);
  return true;
}

MandelLoopEvaluator::MandelLoopEvaluator(MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN), currentWorker(allocator->worker),
  f(&self_allocator), f_z(&self_allocator), f_c(&self_allocator),
  f_zz(&self_allocator), f_zc(&self_allocator), f_cc(&self_allocator), f_zzc(&self_allocator),
  multi(0), first_multi(&self_allocator), sumA(&self_allocator), s1(&self_allocator), s2(&self_allocator)
{
  assert(self_allocator.checkFill());
}

/*void MandelLoopEvaluator::switchType(MandelMath::number_worker *worker)
{
  if (worker==currentWorker)
    return;
  if (currentWorker!=nullptr)
  {
    currentWorker->cleanup(&first_multi_re);
    currentWorker->cleanup(&first_multi_im);
    currentWorker->cleanup(&s1_re);
    currentWorker->cleanup(&s1_im);
    currentWorker->cleanup(&s2_re);
    currentWorker->cleanup(&s2_im);
    currentWorker->cleanup(&sumA_re);
    currentWorker->cleanup(&sumA_im);
    currentWorker->cleanup(&f_zz_re);
    currentWorker->cleanup(&f_zz_im);
    currentWorker->cleanup(&f_zc_re);
    currentWorker->cleanup(&f_zc_im);
    currentWorker->cleanup(&f_cc_re);
    currentWorker->cleanup(&f_cc_im);
    currentWorker->cleanup(&f_zzc_re);
    currentWorker->cleanup(&f_zzc_im);
    currentWorker->cleanup(&f_z_re);
    currentWorker->cleanup(&f_z_im);
    currentWorker->cleanup(&f_c_re);
    currentWorker->cleanup(&f_c_im);
    currentWorker->cleanup(&f_re);
    currentWorker->cleanup(&f_im);
  }

  if (worker)
  {
    uint8_t *placement;
    int step;
    switch (worker->ntype())
    {
      case MandelMath::number_worker::Type::typeEmpty:
        placement=nullptr;
        step=0;
        break;
      case MandelMath::number_worker::Type::typeDouble:
        place.dd=nullptr;
        placement=nullptr;
        step=0;
        break;
      case MandelMath::number_worker::Type::typeDDouble:
        place.dd=reinterpret_cast<MandelMath::dd_real (*)[Place::LEN]> (new MandelMath::dd_real[Place::LEN]());
        placement=(uint8_t *)place.dd;
        step=sizeof(MandelMath::dd_real);
        break;
      case MandelMath::number_worker::Type::typeMulti:
        place.multi=reinterpret_cast<MandelMath::multiprec (*)[Place::LEN]> (new MandelMath::multiprec[Place::LEN]());
        placement=(uint8_t *)place.multi;
        step=sizeof(MandelMath::multiprec);
        break;
    }

    //worker->init_(&z_re, placement); placement+=step;
    //worker->init_(&z_im, placement); placement+=step;
    //worker->init_(&c_re, placement); placement+=step;
    //worker->init_(&c_im, placement); placement+=step;
    worker->init_(&s1_re, placement); placement+=step;
    worker->init_(&s1_im, placement); placement+=step;
    worker->init_(&s2_re, placement); placement+=step;
    worker->init_(&s2_im, placement); placement+=step;
    worker->init_(&sumA_re, placement); placement+=step;
    worker->init_(&sumA_im, placement); placement+=step;
    worker->init_(&f_re, placement); placement+=step;
    worker->init_(&f_im, placement); placement+=step;
    worker->init_(&f_z_re, placement); placement+=step;
    worker->init_(&f_z_im, placement); placement+=step;
    worker->init_(&f_c_re, placement); placement+=step;
    worker->init_(&f_c_im, placement); placement+=step;
    worker->init_(&f_zz_re, placement); placement+=step;
    worker->init_(&f_zz_im, placement); placement+=step;
    worker->init_(&f_zc_re, placement); placement+=step;
    worker->init_(&f_zc_im, placement); placement+=step;
    worker->init_(&f_cc_re, placement); placement+=step;
    worker->init_(&f_cc_im, placement); placement+=step;
    worker->init_(&f_zzc_re, placement); placement+=step;
    worker->init_(&f_zzc_im, placement); placement+=step;
    worker->init_(&first_multi_re, placement); placement+=step;
    worker->init_(&first_multi_im, placement); placement+=step;

    assert((size_t)(placement-(uint8_t *)place.dd)==step*Place::LEN);
  }
  currentWorker=worker;
}*/

bool MandelLoopEvaluator::evalg(int period, const MandelMath::complex *c)
{
  f.assign(c);
  f_c.zero(1, 0);
  f_cc.zero(0, 0);
  //also using s2
  for (int i=0; i<period; i++)
  {
    //g_cc=2*(g_cc*g + g_c*g_c)
    f_cc.mul(&f);
    s2.assign(&f_c);
    s2.sqr();
    f_cc.add(&s2);
    f_cc.lshift(1);
    //g_c=2*g_c*g+1
    f_c.mul(&f);
    f_c.lshift(1);
    currentWorker->add_double(f_c.re, 1);
    //g=g^2+c
    f.sqr();
    f.add(c);
    double f_mag=f.getMag_double();
    double allmag=f_mag+
                  f_c.getMag_double()+
                  f_cc.getMag_double();
    if (allmag>1e60)
      return false;
  }
  return true;
}

bool MandelLoopEvaluator::eval2(int period, const MandelMath::complex *c, const MandelMath::complex *z)
{
  f.assign(z);
  f_z.zero(1, 0);
  f_c.zero(0, 0);
  f_zz.zero(0, 0);
  f_zc.zero(0, 0);
  f_cc.zero(0, 0);
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double()+
                 f_c.getMag_double()+
                 f_zc.getMag_double()+
                 f_cc.getMag_double();
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_cc=2 * (f_c^2 + f * f_cc)
    f_cc.mul(&f);
    s2.assign(&f_c);
    s2.sqr();
    f_cc.add(&s2);
    f_cc.lshift(1);
    // f_zc = 2 * (f_z * f_c + f * f_zc);
    f_zc.mul(&f);
    s2.assign(&f_c);
    s2.mul(&f_z);
    f_zc.add(&s2);
    f_zc.lshift(1);
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    s2.assign(&f_z);
    s2.sqr();
    f_zz.add(&s2);
    f_zz.lshift(1);
    // f_c = 2 * f * f_c + 1;
    f_c.mul(&f);
    f_c.lshift(1);
    currentWorker->add_double(f_c.re, 1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    f_z.lshift(1);
    //f:=f^2+c
    f.sqr();
    f.add(c);
  }
  return true;
}

bool MandelLoopEvaluator::eval_zz(int period, const MandelMath::complex *c, const MandelMath::complex *z)
{
  f.assign(z);
  f_z.zero(1, 0);
  f_zz.zero(0, 0);
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double();
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    s2.assign(&f_z);
    s2.sqr();
    f_zz.add(&s2);
    f_zz.lshift(1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    f_z.lshift(1);
    //bigger than 3 is common
    //f:=f^2+c
    f.sqr();
    f.add(c);
  }
  return true;
}

bool MandelLoopEvaluator::eval_multi(int period, const MandelMath::complex *c, const MandelMath::complex *z, const MandelMath::complex *f_z_target)
{
  f.assign(z);
  f_z.zero(1, 0);
  f_zz.zero(0, 0);
  int near1=0;
  int sumnear1=0;
  //bool dangerzone=false;
  int reducedbymag=period;
  double closest_accepted=10000, closest_rejected=10000;
  int sumA_cnt=0;
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double();
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    s2.assign(&f_z);
    s2.sqr();
    f_zz.add(&s2);
    f_zz.lshift(1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    f_z.lshift(1);
    //f:=f^2+c
    f.sqr();
    f.add(c);

    /*
    assume all nearby roots are on a circle and we visit them all during period (plus other points that we want to filter out)
    z is on the circle at angle f_z_target
    f is on the circle at angle f_z (normalized to |f_z|=1 ?)
    center of circle=A, circle at angle 0 =D
    (z-A)/(D-A)=f_z_target
    (f-A)/(D-A)=f_z
    A=? (D=?)
    z-A=D*f_z_target-A*f_z_target
    f-A=D*f_z-A*f_z
    z*f_z=D*f_z_target*f_z+A*f_z*(1-f_z_target)
    f*f_z_target=D*f_z_target*f_z+A*f_z_target*(1-f_z)
    (z*f_z-f*f_z_target)/(f_z-f_z_target)=A
    */
    if (f_z.dist2_double(f_z_target)*period>0.5) //assuming all points are spaced regularly at 1/period around the circle,
    {                                                                  //skip those close to target to avoid div by 0
      s1.assign(z);
      s1.mul(&f_z);
      s2.assign(&f);
      s2.mul(f_z_target);
      s2.rsub(&s1);
      s1.assign(&f_z);
      s1.sub(f_z_target);
      s1.recip();
      s1.mul(&s2); //s1=A

      double dist_to_A=s1.dist2_double(z);
      double f_z_err=f_z.dist2_double(f_z_target);
      double f_zz_mag=f_zz.getMag_double();
      double expected=f_z_err/f_zz_mag;
      (void)expected;
      if (dist_to_A*3<closest_accepted)
      {
        closest_rejected=closest_accepted;
        closest_accepted=dist_to_A;
        sumA_cnt=1;
        sumA.assign(&s1);
        reducedbymag=MandelMath::gcd(period, i+1);
        first_multi.assign(&f_z);
      }
      else if (dist_to_A<3*closest_accepted)
      {
        if (dist_to_A<closest_accepted)
          closest_accepted=dist_to_A;
        reducedbymag=MandelMath::gcd(reducedbymag, i+1);
        sumA.add(&s1);
        sumA_cnt++;
      }
      else if (dist_to_A<closest_rejected)
        closest_rejected=dist_to_A;
    }
    /*double f_z_mag=currentWorker->toDouble(f_z.getMagTmp());
    if (f_z_mag<0.99)
    {
      //dangerzone=true;
      reducedbymag=MandelMath::gcd(reducedbymag, i+1);
    }
    else if (f_z_mag<1.01)
    {
      near1++;
      sumnear1+=i;
      if (near1==1)
      {
        currentWorker->assign(&first_multi_re, f_z.re_s);
        currentWorker->assign(&first_multi_im, f_z.im_s);
      };
    }
    else if (f_z_mag<2.59) //period 15/5 iter 0+3k: 2.682..2.687
      dangerzone=true;     //period 9/3 iter 0: 2.5905
                           //
                           //
    //bigger than 3 is common
    */
  }
  s1.assign(f_z_target);
  double f_z_err=f_z.getMag_double()-s1.getMag_double();
  double f_zz_mag=f_zz.getMag_double();
  double expected=std::abs(f_z_err)/f_zz_mag;
  if (reducedbymag>=period)
    multi=1;
  else if (closest_accepted>expected*1000)
    multi=1;
  else if (closest_accepted>expected*3)
  { dbgPoint(); multi=1; }
  else if (closest_rejected<13*closest_accepted) //12.27 at second 76/38/19
    multi=1;
  else if (closest_rejected>100000*closest_accepted)
  {
    multi=period/reducedbymag;
    currentWorker->zero(s1.re, sumA_cnt);
    currentWorker->recip(s1.re);
    currentWorker->mul(sumA.re, s1.re);
    currentWorker->mul(sumA.im, s1.re);
  }
  else if (closest_rejected<100*closest_accepted)
    multi=1;
  else if (closest_rejected>1000*closest_accepted)
  {
    multi=period/reducedbymag;
    currentWorker->zero(s1.re, sumA_cnt);
    currentWorker->recip(s1.re);
    currentWorker->mul(sumA.re, s1.re);
    currentWorker->mul(sumA.im, s1.re);
  }
  else if (near1<1)// || dangerzone)
  {
    multi=0;
  }
  else if (near1==1)
  {
    multi=1;
  }
  else
  {
    //expect near1 numbers, each a multiple of period/near1, and that's 1,2,3..period/near1 multiple
    //(period/near1)*near1*(near1+1)/2=period*(near1+1)/2
    //ex: p=15 near1=3 5+10+15=30=15*4/2
    //except we add 2,5,8... not 3,6,9 so period*(near1+1)/2-near1
    if (period%near1!=0)
    { dbgPoint(); multi=0; }
    else if (sumnear1==period*(near1+1)/2-near1)
    {
      multi=near1;
    }
    else
    { dbgPoint(); multi=0; }
  }
  return true;
}

bool MandelLoopEvaluator::eval2zzc(int period, const MandelMath::complex *c, const MandelMath::complex *z)
{
  f.assign(z);
  f_z.zero(1, 0);
  f_c.zero(0, 0);
  f_zz.zero(0, 0);
  f_zc.zero(0, 0);
  f_cc.zero(0, 0);
  f_zzc.zero(0, 0);
  f_zzc.zero(0, 0);
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double()+
                 f_c.getMag_double()+
                 f_zc.getMag_double()+
                 f_cc.getMag_double()+
                 f_zzc.getMag_double();
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_zzc=d/dc 2*(f_z*f_z + f*f_zz)=2*(2*f_z*f_zc + f_c*f_zz+f*f_zzc)
    f_zzc.mul(&f);
    s2.assign(&f_c);
    s2.mul(&f_zz);
    f_zzc.add(&s2);
    s2.assign(&f_z);
    s2.mul(&f_zc);
    s2.lshift(1);
    f_zzc.add(&s2);
    f_zzc.lshift(1);
    //f_cc=2 * (f_c^2 + f * f_cc)
    f_cc.mul(&f);
    s2.assign(&f_c);
    s2.sqr();
    f_cc.add(&s2);
    f_cc.lshift(1);
    // f_zc = 2 * (f_z * f_c + f * f_zc);
    f_zc.mul(&f);
    s2.assign(&f_c);
    s2.mul(&f_z);
    f_zc.add(&s2);
    f_zc.lshift(1);
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    s2.assign(&f_z);
    s2.sqr();
    f_zz.add(&s2);
    f_zz.lshift(1);
    // f_c = 2 * f * f_c + 1;
    f_c.mul(&f);
    f_c.lshift(1);
    currentWorker->add_double(f_c.re, 1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    f_z.lshift(1);
    //f:=f^2+c
    f.sqr();
    f.add(c);
  }
  return true;
}




MandelEvaluator::MandelEvaluator(MandelMath::worker_multi::Type ntype): QThread(nullptr),
  currentWorker(MandelMath::worker_multi::create(ntype, LEN)),
  params_allocator(currentWorker->getAllocator(), ComputeParams::LEN), currentParams(&params_allocator),
  /*currentDataAllocator(currentWorker->getAllocator(), MandelPoint::LEN),*/ currentData(&currentDataStore, currentWorker->getAllocator(), nullptr),
  tmpLaguerrePoint(nullptr, currentWorker->getAllocator(), nullptr),
  /*newtres_allocator(&self_allocator, NewtRes::LEN),*/ newtres_(currentWorker->getAllocator(), nullptr),
  /*eval_allocator(&self_allocator, Eval::LEN),*/ eval(currentWorker->getAllocator(), nullptr),
  /*newt_allocator(&self_allocator, Newt::LEN),*/ newt(currentWorker->getAllocator(), nullptr),
  /*interior_allocator(&self_allocator, InteriorInfo::LEN),*/ interior(currentWorker->getAllocator(), nullptr),
  /*bulb_allocator(&self_allocator, Bulb::LEN),*/ bulb(currentWorker->getAllocator(), nullptr)
{
  QThread::start(QThread::Priority::LowestPriority);
  assert(currentWorker->getAllocator()->checkFill());
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
}

#if NUMBER_DOUBLE_EXISTS
void MandelEvaluator::simple_double(double cr, double ci, MandelPoint &data, int maxiter)
{
  double zr=0;
  double zi=0;
  for (int iter=0; iter<maxiter; iter++)
  {
    if (zr*zr+zi*zi>4)
    {
      data.store->state=MandelPointStore::State::stOutside;
      data.store->iter=iter;
      data.f.zero(zr, zi);
      return;
    };
    double tmp=zr*zr-zi*zi+cr;
    zi=2*zr*zi+ci;
    zr=tmp;
  }
  //data.state=MandelPoint::State::stMaxIter;
  data.store->iter=maxiter;
  data.f.zero(zr, zi);
}
#endif //NUMBER_DOUBLE_EXISTS

#if 0
maybe with local worker
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
    t.assign(r2); t.add(i2.hi, i2.lo_);
    if (t.hi>4)
    {
      data.state=MandelPoint::State::stOutside;
      data.iter=iter;
      data.f_re.as.ddouble_.dd->assign(zr);
      data.f_im.as.ddouble_.dd->assign(zi);
      return;
    };
    t.assign(r2); t.add(-i2.hi, -i2.lo_); t.add(cr->hi, cr->lo_); //double tmp=zr*zr-zi*zi+cr;
    zi.mul(2*zr.hi, 2*zr.lo_); zi.add(ci->hi, ci->lo_);
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
#endif

/*void MandelEvaluator::switchType(MandelMath::number_worker *worker)
{
  if (worker==currentWorker)
    return;
  //TODO: use .promote()
  bulb.bulbe.switchType(worker);
  bulb.lagu.switchType(worker);
  bulb.dbg_guessmult=0;
  if (currentWorker!=nullptr)
  {
    currentWorker->cleanup(&bulb.rb_re_);
    currentWorker->cleanup(&bulb.rb_im);
    currentWorker->cleanup(&bulb.cb_re);
    currentWorker->cleanup(&bulb.cb_im);
    currentWorker->cleanup(&bulb.xc_re);
    currentWorker->cleanup(&bulb.xc_im);
    currentWorker->cleanup(&bulb.baseZC_re);
    currentWorker->cleanup(&bulb.baseZC_im);
    currentWorker->cleanup(&bulb.baseCC_re);
    currentWorker->cleanup(&bulb.baseCC_im);
    currentWorker->cleanup(&bulb.s1_re);
    currentWorker->cleanup(&bulb.s1_im);
    currentWorker->cleanup(&bulb.s2_re_);
    currentWorker->cleanup(&bulb.s2_im_);
    currentWorker->cleanup(&bulb.s3_re);
    currentWorker->cleanup(&bulb.s3_im_);
    currentWorker->cleanup(&bulb.deltac_re);
    currentWorker->cleanup(&bulb.deltac_im);
    currentWorker->cleanup(&bulb.deltar_re);
    currentWorker->cleanup(&bulb.deltar_im);
    currentWorker->cleanup(&bulb.target_f_z_re);
    currentWorker->cleanup(&bulb.target_f_z_im);
    //currentWorker->cleanup(&bulb.test_x0_re);
    //currentWorker->cleanup(&bulb.test_x0_im);
    //currentWorker->cleanup(&bulb.test_xn_re);
    //currentWorker->cleanup(&bulb.test_xn_im);
    currentWorker->cleanup(&bulb.cbx_re);
    currentWorker->cleanup(&bulb.cbx_im);
    currentWorker->cleanup(&bulb.rbx_re);
    currentWorker->cleanup(&bulb.rbx_im);
    currentWorker->cleanup(&bulb.B_re);
    currentWorker->cleanup(&bulb.B_im);
    currentWorker->cleanup(&bulb.C_re);
    currentWorker->cleanup(&bulb.C_im);
    currentWorker->cleanup(&bulb.dbg_first_rb_re);
    currentWorker->cleanup(&bulb.dbg_first_rb_im);
    currentWorker->cleanup(&bulb.dbg_first_cb_re);
    currentWorker->cleanup(&bulb.dbg_first_cb_im);

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
    //currentWorker->cleanup(&newt.fzfix_im);
    //currentWorker->cleanup(&newt.fzfix_re);
    currentWorker->cleanup(&newt.tmp2);
    currentWorker->cleanup(&newt.tmp1_im);
    currentWorker->cleanup(&newt.tmp1_re);
    currentWorker->cleanup(&newt.fzz_r_im);
    currentWorker->cleanup(&newt.fzz_r_re);
    currentWorker->cleanup(&newtres_.first_guess_lagu_re);
    currentWorker->cleanup(&newtres_.first_guess_lagu_im);
    currentWorker->cleanup(&newtres_.first_guess_newt_re);
    currentWorker->cleanup(&newtres_.first_guess_newt_im);
    currentWorker->cleanup(&newtres_.fz_r_im_);
    currentWorker->cleanup(&newtres_.fz_r_re_);
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
  {
    currentData.init(worker);

    uint8_t *placement;
    int step;
    switch (worker->ntype())
    {
      case MandelMath::number_worker::Type::typeEmpty:
        placement=nullptr;
        step=0;
        break;
      case MandelMath::number_worker::Type::typeDouble:
        place.dd=nullptr;
        placement=nullptr;
        step=0;
        break;
      case MandelMath::number_worker::Type::typeDDouble:
        place.dd=reinterpret_cast<MandelMath::dd_real (*)[Place::LEN]> (new MandelMath::dd_real[Place::LEN]());
        placement=(uint8_t *)place.dd;
        step=sizeof(MandelMath::dd_real);
        break;
      case MandelMath::number_worker::Type::typeMulti:
        place.multi=reinterpret_cast<MandelMath::multiprec (*)[Place::LEN]> (new MandelMath::multiprec[Place::LEN]());
        placement=(uint8_t *)place.multi;
        step=sizeof(MandelMath::multiprec);
        break;
    }

    worker->init_(&currentParams.c_re, placement); placement+=step;
    worker->init_(&currentParams.c_im, placement); placement+=step;

    worker->init_(&eval.fz_r_re, placement); placement+=step;
    worker->init_(&eval.fz_r_im, placement); placement+=step;
    worker->init_(&eval.near0fmag, placement); placement+=step;
    worker->init_(&eval.lookper_nearr_re, placement); placement+=step;
    worker->init_(&eval.lookper_nearr_im, placement); placement+=step;

    worker->init_(&newt.bestr_re, placement); placement+=step;
    worker->init_(&newt.bestr_im, placement); placement+=step;
    worker->init_(&newt.f_r_re, placement); placement+=step;
    worker->init_(&newt.f_r_im, placement); placement+=step;
    worker->init_(&newtres_.fz_r_im_, placement); placement+=step;
    worker->init_(&newtres_.fz_r_re_, placement); placement+=step;
    worker->init_(&newtres_.first_guess_newt_re, placement); placement+=step;
    worker->init_(&newtres_.first_guess_newt_im, placement); placement+=step;
    worker->init_(&newtres_.first_guess_lagu_re, placement); placement+=step;
    worker->init_(&newtres_.first_guess_lagu_im, placement); placement+=step;
    worker->init_(&newt.fzz_r_re, placement); placement+=step;
    worker->init_(&newt.fzz_r_im, placement); placement+=step;
    worker->init_(&newt.tmp1_re, placement); placement+=step;
    worker->init_(&newt.tmp1_im, placement); placement+=step;
    worker->init_(&newt.tmp2, placement); placement+=step;
    //worker->init_(&newt.fzfix_re, placement); placement+=step;
    //worker->init_(&newt.fzfix_im, placement); placement+=step;
    worker->init_(&newt.laguH_re, placement); placement+=step;
    worker->init_(&newt.laguH_im, placement); placement+=step;
    worker->init_(&newt.laguG_re, placement); placement+=step;
    worker->init_(&newt.laguG_im, placement); placement+=step;
    worker->init_(&newt.laguG2_re, placement); placement+=step;
    worker->init_(&newt.laguG2_im, placement); placement+=step;
    worker->init_(&newt.laguX_re, placement); placement+=step;
    worker->init_(&newt.laguX_im, placement); placement+=step;
    worker->init_(&newt.newtX_re, placement); placement+=step;
    worker->init_(&newt.newtX_im, placement); placement+=step;
    worker->init_(&newt.fzzf_re, placement); placement+=step;
    worker->init_(&newt.fzzf_im, placement); placement+=step;

    worker->init_(&interior.inte_re, placement); placement+=step;
    worker->init_(&interior.inte_im, placement); placement+=step;
    worker->init_(&interior.inte_abs, placement); placement+=step;
    worker->init_(&interior.fz_re, placement); placement+=step;
    worker->init_(&interior.fz_im, placement); placement+=step;
    worker->init_(&interior.fz_mag, placement); placement+=step;

    worker->init_(&bulb.rb_re_, placement); placement+=step;
    worker->init_(&bulb.rb_im, placement); placement+=step;
    worker->init_(&bulb.cb_re, placement); placement+=step;
    worker->init_(&bulb.cb_im, placement); placement+=step;
    worker->init_(&bulb.xc_re, placement); placement+=step;
    worker->init_(&bulb.xc_im, placement); placement+=step;
    worker->init_(&bulb.baseZC_re, placement); placement+=step;
    worker->init_(&bulb.baseZC_im, placement); placement+=step;
    worker->init_(&bulb.baseCC_re, placement); placement+=step;
    worker->init_(&bulb.baseCC_im, placement); placement+=step;
    worker->init_(&bulb.s1_re, placement); placement+=step;
    worker->init_(&bulb.s1_im, placement); placement+=step;
    worker->init_(&bulb.s2_re_, placement); placement+=step;
    worker->init_(&bulb.s2_im_, placement); placement+=step;
    worker->init_(&bulb.s3_re, placement); placement+=step;
    worker->init_(&bulb.s3_im_, placement); placement+=step;
    worker->init_(&bulb.deltac_re, placement); placement+=step;
    worker->init_(&bulb.deltac_im, placement); placement+=step;
    worker->init_(&bulb.deltar_re, placement); placement+=step;
    worker->init_(&bulb.deltar_im, placement); placement+=step;
    worker->init_(&bulb.target_f_z_re, placement); placement+=step;
    worker->init_(&bulb.target_f_z_im, placement); placement+=step;
    //worker->init_(&bulb.test_x0_re, placement); placement+=step;
    //worker->init_(&bulb.test_x0_im, placement); placement+=step;
    //worker->init_(&bulb.test_xn_re, placement); placement+=step;
    //worker->init_(&bulb.test_xn_im, placement); placement+=step;
    worker->init_(&bulb.cbx_re, placement); placement+=step;
    worker->init_(&bulb.cbx_im, placement); placement+=step;
    worker->init_(&bulb.rbx_re, placement); placement+=step;
    worker->init_(&bulb.rbx_im, placement); placement+=step;
    worker->init_(&bulb.B_re, placement); placement+=step;
    worker->init_(&bulb.B_im, placement); placement+=step;
    worker->init_(&bulb.C_re, placement); placement+=step;
    worker->init_(&bulb.C_im, placement); placement+=step;
    worker->init_(&bulb.dbg_first_rb_re, placement); placement+=step;
    worker->init_(&bulb.dbg_first_rb_im, placement); placement+=step;
    worker->init_(&bulb.dbg_first_cb_re, placement); placement+=step;
    worker->init_(&bulb.dbg_first_cb_im, placement); placement+=step;

    assert((size_t)(placement-(uint8_t *)place.dd)==step*Place::LEN);
  }
  currentWorker=worker;
}*/

bool MandelEvaluator::startCompute(/*const MandelPoint *data,*/ int quick_route)
{
  //currentParams=params;
  //currentData.assign(*data);
  if ((quick_route==1) ||
      ((quick_route==0) && (currentParams.maxiter_-currentData.store->iter<=1000)))
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

void MandelEvaluator::startNewton(int period, const MandelMath::complex *c /*, currentData.f const *root, */)
{
  currentData.store->period=period;
  currentParams.c.assign(c);
  QMetaObject::invokeMethod(this, &MandelEvaluator::doNewton, Qt::ConnectionType::QueuedConnection);
}

void MandelEvaluator::doNewton()
{
  //MandelMath::complex tmpc(currentWorker, &currentParams.c_re, &currentParams.c_im, true);
  //MandelMath::complex root(currentWorker, &currentData.f_re, &currentData.f_im, true);
  int result=newton(currentData.store->period, &currentParams.c, &currentData.f, true, 12);
  emit doneNewton(this, result);
}

bool MandelEvaluator::Bulb::findBulbBase(int period2, const MandelMath::complex *c, MandelMath::complex *cb, MandelMath::complex *rb, MandelMath::complex *xc, MandelMath::complex *baseZC, MandelMath::complex *baseCC, bool *is_card, int *foundMult)
//on input foundMult=0 -> guess rb here; =1 -> rb already set
//xc: bulb center, both z and c
//cb: bulb base c, rb: bulb base root (final point)
{ //"findBulbBaseOri"
  //note: cb, rb, xc are not this->cb, rb, xc when called from Orbit
  if (period2==1)
  {
    cb->zero(0.25, 0);
    rb->zero(0.5, 0);
    xc->zero(0, 0);
    baseZC->zero(0, 0);
    baseCC->zero(0, 0);
    *is_card=true;
    *foundMult=2;
    return true;
  };
  int period=period2;
  xc->assign(c);
  *is_card=false;
  bool did_reduce_period=false;

  /*MandelMath::complex f(currentWorker, &bulb.bulbe.f_re, &bulb.bulbe.f_im, true);
  MandelMath::complex f_z(currentWorker, &bulb.bulbe.f_z_re, &bulb.bulbe.f_z_im, true);
  MandelMath::complex f_c(currentWorker, &bulb.bulbe.f_c_re, &bulb.bulbe.f_c_im, true);
  MandelMath::complex f_zz(currentWorker, &bulb.bulbe.f_zz_re, &bulb.bulbe.f_zz_im, true);
  MandelMath::complex f_zc(currentWorker, &bulb.bulbe.f_zc_re, &bulb.bulbe.f_zc_im, true);
  MandelMath::complex f_cc_(currentWorker, &bulb.bulbe.f_cc_re, &bulb.bulbe.f_cc_im, true);
  MandelMath::complex f_zzc(currentWorker, &bulb.bulbe.f_zzc_re, &bulb.bulbe.f_zzc_im, true);
  MandelMath::complex s1(currentWorker, &bulb.s1_re, &bulb.s1_im, true);
  MandelMath::complex s2_(currentWorker, &bulb.s2_re_, &bulb.s2_im_, true);
  MandelMath::complex s3(currentWorker, &bulb.s3_re, &bulb.s3_im_, true);
  MandelMath::complex deltac(currentWorker, &bulb.deltac_re, &bulb.deltac_im, true);
  MandelMath::complex deltar(currentWorker, &bulb.deltar_re, &bulb.deltar_im, true);
  MandelMath::complex target_f_z(currentWorker, &bulb.target_f_z_re, &bulb.target_f_z_im, true);
  currentWorker->zero(target_f_z.re_s, 1);
  currentWorker->zero(target_f_z.im_s, 0);*/
  target_f_z.zero(1, 0);
  *foundMult=1;
  /*MandelMath::complex g(currentWorker, &bulb.bulbe.g_re, &bulb.g_im, true);
  MandelMath::complex g_c(currentWorker, &bulb.g_c_re, &bulb.g_c_im, true);
  MandelMath::complex g_c2(currentWorker, &bulb.g_c2_re, &bulb.g_c2_im, true);
  MandelMath::complex g_cc(currentWorker, &bulb.g_cc_re, &bulb.g_cc_im, true);*/
  //1) find bulb center
  for (;;)
  {
    for (int cyc=0; cyc<10; cyc++)
    {
      if (!bulbe.evalg(period, xc))
        return false;
      bulbe.f.rsub(xc);
      currentWorker->add_double(bulbe.f_c.re, -1);
      double g_mag=bulbe.f.getMag_double();
      //ideally copy the Lagu code from newton() but newton should be enough
      double g_c_mag=bulbe.f_c.getMag_double();
      double g_cc_mag;
      if (g_mag>currentWorker->eps2()*1000)
      { //xc=xc-2*g/g_c
        g_cc_mag=1;
        bulbe.f_c.recip_prepared();
        bulbe.f.mul(&bulbe.f_c);
        bulbe.f.lshift(1);
        xc->add(&bulbe.f);
      }
      else
      { //xc=xc-g_c/g_cc
        g_cc_mag=bulbe.f_cc.getMag_double();
        bulbe.f_cc.recip_prepared();
        bulbe.f_c.mul(&bulbe.f_cc);
        xc->sub(&bulbe.f_c);
      }
      if (g_c_mag<g_cc_mag*currentWorker->eps2()*2)
        break;
      if (cyc==9)
        cyc=cyc;
    }
    //cleanup period?
    break;
  }
  //f=0 f_c~0 f_cc=-0.00860858-i0.03615690
  //deus ex machina   bulb 1/4: cent~0.2822713907669139+i0.5300606175785253 base c~0.249+i0.500 r~-0.01+i0.499
  /*currentWorker->assign(rb->re_s, xc->re_s);
  currentWorker->assign(rb->im_s, xc->im_s);
  currentWorker->assign(cb->re_s, xc->re_s);
  currentWorker->assign(cb->im_s, xc->im_s);
  bulb.bulbe.eval2(period, cb, rb);*/ //currentWorker->sub(f.re_s, rb->re_s); currentWorker->sub(f.im_s, rb->im_s); currentWorker->add_double(f_z.re_s, -1);
  /*currentWorker->assign(s1.re_s, &bulb.bulbe.f_zc_re);
  currentWorker->assign(s1.im_s, &bulb.bulbe.f_zc_im);
  currentWorker->add(s1.re_s, &bulb.bulbe.f_zz_re);
  currentWorker->add(s1.im_s, &bulb.bulbe.f_zz_im);
  s1.recip();
  cb->add(&s1);
  currentWorker->lshift(s1.re_s, 1);
  currentWorker->lshift(s1.im_s, 1);
  rb->add(&s1);*/

  /*bulb.bulbe.eval2(period, cb, rb);
  currentWorker->assign(s1.re_s, cb->re_s);
  currentWorker->assign(s1.im_s, cb->im_s);
  currentWorker->sub(s1.re_s, xc->re_s); //s1=cb-xc
  currentWorker->sub(s1.im_s, xc->im_s);
  currentWorker->assign(s3.re_s, &bulb.bulbe.f_zz_re);
  currentWorker->assign(s3.im_s, &bulb.bulbe.f_zz_im);
  s3.mul(&f_z);
  currentWorker->rsub(s3.re_s, &bulb.bulbe.f_zc_re); //s3=f_zc-f_z*f_zz
  currentWorker->rsub(s3.im_s, &bulb.bulbe.f_zc_im);
  currentWorker->assign(s2_.re_s, f_z.re_s); //s2=f_z+1
  currentWorker->assign(s2_.im_s, f_z.im_s);
  s2_.recip();
  s2_.mul(&s1);
  s2_.mul(&s3);
  bulb.dbg_guessmult=currentWorker->toDouble(s2_.re_s); //1+1/x*/


  /*
  find r, c where f=0 fz=1
  we already have f=0 so just keep that: fz*(r-xc)+fc*(c-xc)=0
      //but fz(xc,xc)==0 fc(xc,xc)==1 so we need 2nd derivatives
      //fz*(r-xc)+fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+fc*(c-xc)+fcc/2*(c-xc)^2=0
      //fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+(c-xc)+fcc/2*(c-xc)^2=0
      actually fz(xc,xc)=-1, we need fz=0
      fzz*(r-xc)+fzc*(c-xc)=1 //move from fz=-1 to fz=0
      fz*(r-xc)+fc*(c-xc)=0   //keep f=0
      solve
      fzz*(r-xc)+fzc*(c-xc)=1
      (r-xc)=(c-xc)
      (fzz+fzc)*(c-xc)=1 correct but uhh delta r=delta c ?
      2nd derivatives
      fzz*(r-xc)+fzc*(c-xc)=1 //move from fz=-1 to fz=0
      fz*(r-xc)+fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+fc*(c-xc)+fcc/2*(c-xc)^2=0   //keep f=0
      solve, S=c-xc
      fzz*(r-xc)=(1-fzc*S)    8*0.5=1+4*0.25  4=2
      -fzz*(r-xc)/fzz+fzz/2*fzz*(r-xc)^2/fzz+fzc*fzz*(r-xc)/fzz*S+S+fcc/2*S^2=0
      -(1-fzc*S)/fzz+1/2*(1-fzc*S)^2/fzz+fzc*(1-fzc*S)/fzz*S+S+fcc/2*S^2=0   |*2*fzz
      -2*(1-fzc*S)+(1-fzc*S)^2+2*fzc*(1-fzc*S)*S+2*fzz*S+fzz*fcc*S^2=0
      2*(fzz+fzc)*S+(fzz*fcc-fzc*fzc)*S^2-1=0
      fzz*fcc-fzc*fzc=0 at r=c=xc
      2*(fzz+fzc)*S=1

  being at xr, xc assume fzz=fz=0 at target r0, c0: replace fz*(r-xr)+fzz/2*(r-xr)^2 with fzzz*(r-r0)^3 =fzzz*(r-r0)*(r-r0)^2=fzz*(r-r0)^2
  fz~3*fzzz*(xr-r0)^2  fzz~6*fzzz*(xr-r0)

  xc=-0.76 -> xr=1/10 (5 - sqrt(101))
  fz=0.02004975155164389191229403490 fc=-0.00997512422417805404385298255 fzz=0.02014925465493167573688210469  fzc=-2.0199502484483561080877059651038304 fcc=2
  at xc=-0.76 r0=-0.5: f=0.0001  r0-xr=0.004987562112089027021926491275957618694502347002637729057282829...

  S=c0-xc
  fzz/2*(r0-xr)+fzc*S=-fz   //move from fz to fz=0       4*0.5-4*0.25=1  2-1=1
  fzz*(r0-xr)=2*(-fz-fzc*S)
  fz*(r0-xr)+fzz/2*(r0-xr)^2+fzc*(r0-xr)*S+fc*S+fcc/2*S^2=0
  -fzz*(r0-xr)^2+fzc*(r0-xr)*S+fc*S+fcc/2*S^2=0
  -4*(-fz-fzc*S)*(-fz-fzc*S)+fzc*2*(-fz-fzc*S)*S+fc*S*fzz+fzz*fcc/2*S^2=0
  -4*fz^2-4*2*fz*fzc*S-4*fzc*fzc*S*S-2*fzc*fz*S-2*fzc*fzc*S*S+fc*S*fzz+fzz*fcc/2*S^2=0
  -4*fz^2+(fc*fzz-10*fz*fzc)*S+(fzz*fcc/2-6*fzc*fzc)*S^2=0
  4+(8-10*4)*S+(8*2-6*16)*S^2 = 4-32S-80S^2 = 1-8S-10S^2



                        1-|fz+1|^2
  // c0-c=  -----------------------
  //            | fzc + fzz fc/fz |
  if?
  fzz/2*(r-r0)+fzc*(c-c0)=-fz


  change in fz = 1 = fzz*(r-xc)+fzc*(c-xc)
  solve
    fzz*(r-xc)+fzc*(c-xc)=1
    fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+(c-xc)+fcc/2*(c-xc)^2=0
    |  S=c-xc
    fzz*(r-xc)=(1-fzc*S)
    1/fzz/2*fzz^2*(r-xc)^2+fzc*fzz*(r-xc)/fzz*S+S+fcc/2*S^2=0
    1/fzz/2*(1-fzc*S)^2+fzc*(1-fzc*S)/fzz*S+S+fcc/2*S^2=0    |*2*fzz
    1+2*fzz*S+(fzz*fcc-fzc*fzc)*S^2=0
    bulb 1/2: fz=-1  fc=1  fzz=8  fzc=-4  fcc=2  fzcc=0  fzzc=4  fzzz=-24
    correct guess 1/(fzz+fzc)=1/(8-4)=1/4   -1+1/4=-0.75 correct
            8*0.5=(1+4*0.25)   4=1+1

    base of bulb 1/2: c=-0.75 r=-0.5 fz=0 fc=0 fzz=0 fzc=-2 fcc=2  fzcc=0  fzzc=4  fzzz=-12   evaluate d^2/dz^2 (z^2+c)^2+c-z at c=-0.75 and z=-0.5



    (r-xc)=(1-fzc*S)/fzz
    fzz/2*(1-fzc*S)/fzz*(1-fzc*S)/fzz+fzc*(1-fzc*S)/fzz*S+S+fcc/2*S^2=0  |*2*fzz
    (1-fzc*S)*(1-fzc*S)+2*fzc*(1-fzc*S)*S+2*fzz*S+fzz*fcc*S^2=0
    1+2*fzz*S+(fzz*fcc-fzc*fzc)*S^2=0
    1+2*fzz*S+(fzz*fcc-fzc*fzc)*S^2=0

    fzz^2*(1-fzc*S)/fzz*(1-fzc*S)/fzz+2*fzz*fzc*(1-fzc*S)/fzz*S+S*fzz*2+fcc/2*S^2=0
    1/2*(1-fzc*S)*(1-fzc*S)/fzz+2*fzc*(1-fzc*S)*S+S+fcc/2*S^2=0
    (fcc*fzz-fzc*fzc)*S^2+2*(fzc+fzz)*S+1=0
      proof? that fcc*fzz-fzc^2=0
      w.l.o.g f=a*x^2+b*x*y+c*y^2
      fxx=2a fyy=2c fxy=b
      not really



    fzz/2*R+fzc*C+fz=0    f=fzzz/6*R^3  fz=3*fzzz/6*R^2  fzz=6*fzzz/6*R  fzzz=fzzz  2*fz/fzz=R  2*-1/8=-1/4!=R
    fc=1+fcc*C+fzc*R  0=1+2*0.25-(4+2)/2*0.5=1+0.5-1.5
    lower case at xr, xc; upper case at r0, c0  C=c-c0 R=r-r0
    f=0+FC*C+FZC*R*C+FCC*C^2/2+FZZZ*R^3/6+FZZC*R^2*C/2  (FZ=FZZ=0  (FZZZ),FZZC,FZCC,FCC don't change much)
    fz=FZC*C+FZZZ*R^2/2+FZZC*R*C    -1=-2*-0.25+(-12..-24)*0.5^2/2+4*-0.5*-0.25=1+(-1.5..-3)=-0.5..-2~-1
    fzz=FZZZ*R+FZZC*C               8=(-12..-24)*-0.5+4*-0.25=6..12-1=5..11~8
    fc=FC+FZC*R+FCC*C+FZZC*R^2/2    1=0+-2*-0.5+2*-0.25+4*0.5^2/2=1-0.5+1/2=1
    fzc=FZC+FZZC*R                  -4=-2+4*-0.5=-2-2=-4
    fcc=FCC            -> FC FZC FCC FZZZ FZZC R C
    fzzc=FZZC
    now find R, C from fz..fzzc and f=0
    0=FC*C+(fzc-fzzc*R)*R*C+fcc*C^2/2+(fzz-fzzc*C)*R^2/6+fzzc*R^2*C/2
    fz=(fzc-fzzc*R)*C+(fzz-fzzc*C)*R/2+fzzc*R*C
    //(fzz-fzzc*C)=FZZZ*R
    fc=FC+(fzc-fzzc*R)*R+fcc*C+fzzc*R^2/2  |*C-(1)
    //(fzc-fzzc*R)=FZC
    -
    //fc*C= FC*C+(fzc-fzzc*R)*R*C+fcc*C*C+fzzc*R^2*C/2
    //0   =-FC*C-(fzc-fzzc*R)*R*C-fcc*C*C/2-(fzz-fzzc*C)*R^2/6-fzzc*R^2*C/2
    fc*C=fcc*C*C/2-fzz*R^2/6+fzzc*R^2*C/6
    fz=fzc*C+fzz*R/2-fzzc*C*R/2
    -
    2*fzc*C-2*fz=(-fzz+fzzc*C)*R   (2*fzc*C-2*fz)/(-fzz+fzzc*C)=R
    fc*C=fcc*C*C/2+(-fzz+fzzc*C)*R*R/6
    -
    fc*C=fcc*C*C/2+(2*fzc*C-2*fz)*(2*fzc*C-2*fz)/(-fzz+fzzc*C)/6
    -
    0=4*fz*fz+(6*fc*fzz-8*fz*fzc)*C+(4*fzc*fzc-3*fcc*fzz-6*fc*fzzc)*C*C+3*fzzc*fcc*C*C*C
    at 1/2 bulb center: ->C=-0.212873 R=-0.418346
    first newton step: 4*fz*fz/(6*fc*fzz-8*fz*fzc), at bulb center 2/(3*fzz+4*fzc)

    FC=d FZC=e FCC=f FZZZ=g FZZC=h
    f=0+d*C+e*R*C+f*C^2/2+g*R^3/6+h*R^2*C/2
    fz=e*C+g*R^2/2+h*R*C
    fzz=g*R+h*C
    fc=d+e*R+f*C+h*R^2/2
    fzc=e+h*R
    fcc=f
    fzzc=h
    -> C≈-0.212873 R≈-0.418346 fc≈0.102389 fzc≈-2.32662 fcc=2 fzzz≈-21.1583 fzzc=4
    correct C=-0.25 R=-0.5     fc=0        fzc=-2       fcc=2 fzzz=-12      fzzc=4
  */

  //so let's find r0, c0 aka rb, cb
  rb->assign(xc);
  cb->assign(xc);
  //double test_x0_mag=0;
  for (int cycle=0; cycle<5; cycle++)
  {
#if 0 //attempt at 2-nd order approximation... quite a fail, and not needed any more since we reduce period later
    bulb.bulbe.eval2zzc(period, cb, rb);
    currentWorker->sub(f.re_s, rb->re_s);
    currentWorker->sub(f.im_s, rb->im_s);
    currentWorker->add_double(f_z.re_s, -1);
    //0=4*fz*fz+(6*fc*fzz-8*fz*fzc)*C+(4*fzc*fzc-3*fcc*fzz-6*fc*fzzc)*C*C+3*fzzc*fcc*C*C*C
    currentWorker->assign(s2_.re_s, f_zzc.re_s);
    currentWorker->assign(s2_.im_s, f_zzc.im_s);
    s2_.mul(&f_cc_);
    double A3_re=3*currentWorker->toDouble(s2_.re_s);
    double A3_im=3*currentWorker->toDouble(s2_.im_s);

    currentWorker->assign(s2_.re_s, f_zzc.re_s);
    currentWorker->assign(s2_.im_s, f_zzc.im_s);
    s2_.mul(&f_c);
    currentWorker->lshift(s2_.re_s, 1);
    currentWorker->lshift(s2_.im_s, 1);
    currentWorker->assign(s1.re_s, f_zz.re_s);
    currentWorker->assign(s1.im_s, f_zz.im_s);
    s1.mul(&f_cc_);
    s2_.add(&s1);
    currentWorker->mul_double(s2_.re_s, -3);
    currentWorker->mul_double(s2_.im_s, -3);
    currentWorker->assign(s1.re_s, f_zc.re_s);
    currentWorker->assign(s1.im_s, f_zc.im_s);
    currentWorker->lshift(s1.re_s, 1);
    currentWorker->lshift(s1.im_s, 1);
    s1.sqr();
    s2_.add(&s1);
    double B3_re=currentWorker->toDouble(s2_.re_s);
    double B3_im=currentWorker->toDouble(s2_.im_s);

    currentWorker->assign(s2_.re_s, f_zc.re_s);
    currentWorker->assign(s2_.im_s, f_zc.im_s);
    s2_.mul(&f_z);
    currentWorker->lshift(s2_.re_s, 3);
    currentWorker->lshift(s2_.im_s, 3);
    currentWorker->assign(s1.re_s, f_zz.re_s);
    currentWorker->assign(s1.im_s, f_zz.im_s);
    s1.mul(&f_c);
    currentWorker->mul_double(s1.re_s, 6);
    currentWorker->mul_double(s1.im_s, 6);
    currentWorker->sub(s1.re_s, s2_.re_s);
    currentWorker->sub(s1.im_s, s2_.im_s);
    double C3_re=currentWorker->toDouble(s1.re_s);
    double C3_im=currentWorker->toDouble(s1.im_s);

    currentWorker->assign(s2_.re_s, f_z.re_s);
    currentWorker->assign(s2_.im_s, f_z.im_s);
    currentWorker->lshift(s2_.re_s, 1);
    currentWorker->lshift(s2_.im_s, 1);
    s2_.sqr();
    double D3_re=currentWorker->toDouble(s2_.re_s);
    double D3_im=currentWorker->toDouble(s2_.im_s);

    //initial guess between (1-(fz+1)^2)/(fzc-fzz*fc/fz) and 1/2x that -> mul by (1+|fz|)/2
    //-4*fz-4*fz*fz-fz^3    (2-1)*(1-(-1+1)^2)/(-4+8*1/-1)=1/(-4+8*1/-1)
    //(1+fzmag)/2*fz^2*(2+fz))/(fzc*fz-fzz*fc)  (1+1)/2*1*(2-1)/(-4*-1-8*1)=1/-4
    //deltac:=fz^2*(4+4*fz+fz^2)/(fz*fzc-fc*fzz)   -(4-4+1)/(-1*-4-1*8)=-1/(-4)=0.25
    currentWorker->assign(s1.re_s, f_zc.re_s);
    currentWorker->assign(s1.im_s, f_zc.im_s);
    s1.mul(&f_z);
    currentWorker->assign(s2_.re_s, f_zz.re_s);
    currentWorker->assign(s2_.im_s, f_zz.im_s);
    s2_.mul(&f_c);
    currentWorker->sub(s1.re_s, s2_.re_s);
    currentWorker->sub(s1.im_s, s2_.im_s);
    s1.recip();
    currentWorker->assign(s2_.re_s, f_z.re_s);
    currentWorker->assign(s2_.im_s, f_z.im_s);
    s2_.sqr();
    s1.mul(&s2_); //fz^2/(fz*fzc-fc*fzz)
    /*currentWorker->assign(s3.re_s, f_z.re_s);
    currentWorker->assign(s3.im_s, f_z.im_s);
    currentWorker->add_double(s3.re_s, 1);
    currentWorker->lshift(s3.re_s, 2);
    currentWorker->lshift(s3.im_s, 2);
    s2_.add(&s3);*/
    currentWorker->assign(s2_.re_s, f_z.re_s);
    currentWorker->assign(s2_.im_s, f_z.im_s);
    currentWorker->add_double(s2_.re_s, 2);
    s1.mul(&s2_); // *(2+fz)
    currentWorker->assign(s2_.re_s, f_z.getMagTmp());
    currentWorker->add_double(s2_.re_s, 1);
    currentWorker->lshift(s2_.re_s, -1);
    currentWorker->zero(s2_.im_s);
    s1.mul(&s2_); // *(1+fz_mag)/2   not sure what function is best but we need fzmag=1 -> *1 fzmag=0 -> *0.5
    double Z_re=currentWorker->toDouble(s1.re_s); //-0.158439 -> -0.03827124; 0.00449088903772943->/4; -4.8888e-8->-1.2222e-8
    double Z_im=currentWorker->toDouble(s1.im_s);

    //now let's solve A3*Z^3+B3*Z^2+C3*Z+D3=0
    for (int cycle3=0; cycle3<9; cycle3++)
    {
      double F_re=A3_re, F_im=A3_im;
      double t;
      t=F_re*Z_re-F_im*Z_im+B3_re;
      F_im=F_im*Z_re+F_re*Z_im+B3_im;
      F_re=t; //F=A3*F+B3
      t=F_re*Z_re-F_im*Z_im+C3_re;
      F_im=F_im*Z_re+F_re*Z_im+C3_im;
      F_re=t; //F=(A3*Z+B3)*Z+C3
      t=F_re*Z_re-F_im*Z_im+D3_re;
      F_im=F_im*Z_re+F_re*Z_im+D3_im;
      F_re=t; //F=((A3*Z+B3)*Z+C3)*Z+D3
      double Fmag=(F_re*F_re+F_im*F_im);
      if (Fmag==0)
        break;
      t=1.0/Fmag;
      double f1_re=F_re*t;
      double f1_im=-F_im*t; //f1=1/F

      double laguG_re=3*A3_re;
      double laguG_im=3*A3_im; //laguG=3*A3
      t=laguG_re*Z_re-laguG_im*Z_im+2*B3_re;
      laguG_im=laguG_im*Z_re+laguG_re*Z_im+2*B3_im;
      laguG_re=t; //laguG=3*A3*Z+2*B3
      t=laguG_re*Z_re-laguG_im*Z_im+C3_re;
      laguG_im=laguG_im*Z_re+laguG_re*Z_im+C3_im;
      laguG_re=t; //laguG=(3*A3*Z+2*B3)*Z+C3=f' = 3*A3*Z^2+2*B3*Z+C3
      t=laguG_re*f1_re-laguG_im*f1_im;
      laguG_im=laguG_im*f1_re+laguG_re*f1_im;
      laguG_re=t; //laguG=f'/f

      double fzzf_re=6*A3_re;
      double fzzf_im=6*A3_im;
      t=fzzf_re*Z_re-fzzf_im*Z_im+2*B3_re;
      fzzf_im=fzzf_im*Z_re+fzzf_re*Z_im+2*B3_im;
      fzzf_re=t; //fzzf=6*A3*Z+2*B3=f'' = 6*A3*Z+2*B3
      t=fzzf_re*f1_re-fzzf_im*f1_im;
      fzzf_im=fzzf_im*f1_re+fzzf_re*f1_im;
      fzzf_re=t; //fzzf=f''/f

      t=laguG_re*laguG_re-laguG_im*laguG_im;
      double laguG2_im=2*laguG_re*laguG_im;
      double laguG2_re=t; //laguG2=laguG^2
      // laguH=fzf^2-fzzf
      // m=Round( Re(G^2*H^T)/mag(H) )
      // a=1/(G/n +- sqrt( (1/m-1/n)*(G^2-fzzf-G^2/n) ))
      double laguH_re=laguG2_re-fzzf_re;
      double laguH_im=laguG2_im-fzzf_im;

      // a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
      double laguX_re=(laguH_re-laguG2_re/3.0)*(2/3.0);
      double laguX_im=(laguH_im-laguG2_im/3.0)*(2/3.0); //laguX=(1/m-1/n)*(H-G^2/n)  m=1 n=3
      MandelMath::complex_double_sqrt(&laguX_re, &laguX_im, laguX_re, laguX_im);
      //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
      t=laguX_re*laguG_re+laguX_im*laguG_im; //laguX.mulreT(&laguG)
      if (t<0)
      {
        laguX_re=-laguX_re;
        laguX_im=-laguX_im;
      };
      laguX_re+=laguG_re/3.0;
      laguX_im+=laguG_im/3.0; //(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
      t=1/(laguX_re*laguX_re+laguX_im*laguX_im);
      laguX_re=laguX_re*t;
      laguX_im=-laguX_im*t; //1/(...)

      double oldzre=Z_re; Z_re-=laguX_re;
      double oldzim=Z_im; Z_im-=laguX_im;
      if ((Fmag<1e-37) || (Z_re==oldzre && Z_im==oldzim))
        break;
    }
    //now Z is the step in C (in cb)
    //(2*fzc*C-2*fz)/(fzzc*C-fzz)=R
    currentWorker->zero(deltac.re_s, Z_re);
    currentWorker->zero(deltac.im_s, Z_im);
    currentWorker->assign(deltar.re_s, deltac.re_s);
    currentWorker->assign(deltar.im_s, deltac.im_s);
    deltar.mul(&f_zc);
    currentWorker->sub(deltar.re_s, f_z.re_s);
    currentWorker->sub(deltar.im_s, f_z.im_s);
    currentWorker->lshift(deltar.re_s, 1);
    currentWorker->lshift(deltar.im_s, 1);
    currentWorker->assign(s1.re_s, deltac.re_s);
    currentWorker->assign(s1.im_s, deltac.im_s);
    s1.mul(&f_zzc);
    currentWorker->sub(s1.re_s, f_zz.re_s);
    currentWorker->sub(s1.im_s, f_zz.im_s);
    s1.recip();
    deltar.mul(&s1); //R
#else //4-card c~-0.154723+I*1.031046 r~-0.153526+I*1.029663
    bulbe.eval2(period, cb, rb); //TODO: dont' need f_cc here
    bulbe.f.sub(rb);
    currentWorker->add_double(bulbe.f_z.re, -1);
#endif
    //for cardioid, we need to use f_zz as well in the equation for f because it does not go to 0
    //0=f=0+C*fc+R*fz+R*R*fzz/2    C=-R*(fz+R*fzz/2)/fc
    //target_fz-1=fz+C*fzc+R*fzz   target_fz-1=fz-R*(fz+R*fzz/2)/fc*fzc+R*fzz    target_fz-1-fz=-R*R*fzz/fc*fzc/2+R*(fzz-fz/fc*fzc)
    //0=R*R*fzz/2/fc*fzc+R*(fz/fc*fzc-fzz)+target_fz-fz-1
    //x1,2= -(2*c)/(b+-sqrt(b^2-2*a*2*c))         good except both b,c small e.g. 0
    s1.assign(&bulbe.f_c);
    s1.recip();
    s1.mul(&bulbe.f_zc);
    s2.assign(&s1); //f_zc/f_c
    s1.mul(&bulbe.f_z);
    s1.sub(&bulbe.f_zz); //s1=b
    s2.mul(&bulbe.f_zz); //s2=2*a
    s3.assign(&bulbe.f_z);
    s3.sub(&target_f_z);
    currentWorker->add_double(s3.re, 1); //(fz-target_fz+1)
    s3.lshift(1); //s3=-2*c
    s2.mul(&s3);
    deltar.assign(&s1);
    deltar.sqr();
    deltar.add(&s2); //b^2-4ac
    deltar.sqrt();
    if (currentWorker->toDouble(deltar.mulreT_tmp(&s1))<0)
    {
      currentWorker->chs(deltar.re);
      currentWorker->chs(deltar.im);
    };
    deltar.add(&s1);
    deltar.recip();
    deltar.mul(&s3);
    //C=-R*(fz+R*fzz/2)/fc
    deltac.assign(&bulbe.f_zz);
    deltac.mul(&deltar);
    deltac.lshift(-1);
    deltac.add(&bulbe.f_z);
    deltac.mul(&deltar);
    s1.assign(&bulbe.f_c);
    s1.recip();
    deltac.mul(&s1);
    currentWorker->chs(deltar.re);
    currentWorker->chs(deltar.im);
    s1.assign(&deltac);
    s2.assign(&deltar);

    //check the easy way:
    //        0=f=0+C*fc+R*fz -> R=-C*fc/fz
    //target_fz-1=fz+C*fzc+R*fzz   fz/(fc/fz*fzz-fzc)=C=-deltaC    (fz-target_fz+1)/(fzc-fc/fz*fzz)=deltaC
    deltac.assign(&bulbe.f_z);
    deltac.recip();
    deltac.mul(&bulbe.f_c);
    deltar.assign(&deltac); // fc/fz
    deltac.mul(&bulbe.f_zz);
    deltac.rsub(&bulbe.f_zc);
    {
      double mag=deltac.getMag_double();
      if (mag<1e-8)
      {
        nop();
        return false;
      }
      else if (mag<0.1)
        nop();
    }
    deltac.recip(); // 1/(fzc-fc/fz*fzz)
    s3.assign(&bulbe.f_z);
    s3.sub(&target_f_z);
    currentWorker->add_double(s3.re, 1); //(fz-target_fz+1)
    deltac.mul(&s3);  //deltac
    deltar.mul(&deltac);
    currentWorker->chs(deltar.re); // deltar
    currentWorker->chs(deltar.im);

#if 0
    //at card, we are not at result yet so f_zz/f_z does not blow clearly enough
    if (!did_reduce_period && !*is_card) //after reduction, we can arrive at cardioid but must not report it any more
    { //preferably we need to decide at cycle 1 because otherwise it never converges
      currentWorker->sub(&bulb.test_xn_re, &bulb.test_x0_re);
      currentWorker->sub(&bulb.test_xn_im, &bulb.test_x0_im);
      double test_xn=MandelMath::sqr_double(currentWorker->toDouble(&bulb.test_xn_re))+MandelMath::sqr_double(currentWorker->toDouble(&bulb.test_xn_im));
      if (test_xn>test_x0_mag*100)
        *is_card=true;
      else if (test_xn>test_x0_mag*4.1)
        nop(); //?
      else if (test_xn*4.1>test_x0_mag)
        nop(); //xn~x0
      else
        nop(); //xn<x0 should not be
      /*double test_x1=currentWorker->toDouble(f_zz.getMagTmp());
      double test_x2=currentWorker->toDouble(f_zz.getMagTmp())/currentWorker->toDouble(f_z.getMagTmp())*currentWorker->toDouble(f_c.getMagTmp());
        //does not work, maybe f_zz^2
      double test_zc=currentWorker->toDouble(f_zc.getMagTmp());
      if (test_x1*1000<test_zc)
        *is_card=true;
      else if (test_x1<test_zc) // *0.25 is bulb
        nop();
      else
        nop();
      (void)test_x1;
      (void)test_x2;*/
    }
#endif
    if (*is_card)
    {
      //no idea why but the quadratic approximation gives deltac=easy_deltac/2, deltar=easy_deltar
      //easy/2 actually seems better than quadratic approx
      /*currentWorker->assign(deltac.re_s, s1.re_s);
      currentWorker->assign(deltac.im_s, s1.im_s);
      currentWorker->assign(deltar.re_s, s2_.re_s);
      currentWorker->assign(deltar.im_s, s2_.im_s);*/
      deltac.lshift(-1);
    }



    //s2:=f_zc^2-f_zz*f_cc
    /* currentWorker->assign(s2_.re_s, f_zz.re_s);
    currentWorker->assign(s2_.im_s, f_zz.im_s);
    s2_.mul(&f_cc_);
    currentWorker->assign(s1.re_s, f_zc.re_s);
    currentWorker->assign(s1.im_s, f_zc.im_s);
    s1.sqr();
    currentWorker->sub(s2_.re_s, s1.re_s);
    currentWorker->sub(s2_.im_s, s1.im_s);
    //g always=0+0i at the bulb center (or it should) and 1+0i at cardioid
    double cardioid_discrim=currentWorker->toDouble(s2_.getMagTmp());
    double precis=currentWorker->toDouble(f_zc.getMagTmp());
    if (cardioid_discrim<precis*1e-10)
      *is_card=false;
    else if (cardioid_discrim<1-precis*1e-10)
      *is_card=false; //?
    else if (cardioid_discrim<1+precis*1e-10)
      *is_card=true;
    else
      *is_card=true; //? */


    cb->sub(&deltac);
    rb->sub(&deltar);
    if (cycle==0)
    {
      dbg_first_cb.assign(cb);
      dbg_first_rb.assign(rb);
    };


/*
    //verification
    would have to compute all uppercase, which I don't really need for anything else
    f=0+FC*C+FZC*R*C+FCC*C^2/2+FZZZ*R^3/6+FZZC*R^2*C/2  (FZ=FZZ=0  (FZZZ),FZZC,FZCC,FCC don't change much)
    fz=FZC*C+FZZZ*R^2/2+FZZC*R*C    -1=-2*-0.25+(-12..-24)*0.5^2/2+4*-0.5*-0.25=1+(-1.5..-3)=-0.5..-2~-1
    fzz=FZZZ*R+FZZC*C               8=(-12..-24)*-0.5+4*-0.25=6..12-1=5..11~8
    fc=FC+FZC*R+FCC*C+FZZC*R^2/2    1=0+-2*-0.5+2*-0.25+4*0.5^2/2=1-0.5+1/2=1
    fzc=FZC+FZZC*R                  -4=-2+4*-0.5=-2-2=-4

           f                                        fz                    fzz        fc                   fzc       fcc  fzzc
    solve [0=0+d*C+e*R*C+f*C^2/2+g*R^3/6+h*R^2*C/2, -1=e*C+g*R^2/2+h*R*C, 8=g*R+h*C, 1=d+e*R+f*C+h*R^2/2, -4=e+h*R, 2=f, 4=h]
solve [0=0+d*C+e*R*C+f*C^2/2+g*R^3/6+h*R^2*C/2, -1=e*C+g*R^2/2+h*R*C, -124.796127+I*6.194576=g*R+h*C, 1=d+e*R+f*C+h*R^2/2, 48.31925+I*48.506833=e+h*R, 2.0051-I*37.4627558=f, -1020.0523+I*3920.8781217=h]
*/

    s1.assign(rb);

    //now we guessed step in cb and rb, need to tune rb so that f(cb, rb)=0 again
    for (int lcycle=0; lcycle<10; lcycle++)
    {
      //x^5=1 -> x:=x-(x^5-1)/(5*x^4)   NestList[#-(#^5-1)/(5*#^4) &, 0.6+0.8I, 10]
      //x=a^(1/5) |y=x/a| a*y=a^(1/5) y=a^(-4/5)  y^5=a^-4  y:=y-(y^5-a^-4)/(5y^4)  y:=(4/5)*y+a^-4/(5y^4)
      //                                          y^-5=a^4  y:=y-(y^-5-a^4)/(-5*y^-6)=y*(6-a^4*y^5)/5   x=a*y
      if (!bulbe.eval_zz(period, cb, rb))
      {
        *foundMult=1;
        return false;
      };
      bulbe.f.sub(rb);
      currentWorker->add_double(bulbe.f_z.re, -1);
#if 0 //worse than without; leave fixing of f_z to the 2D newton above
      if (currentWorker->toDouble(f.getMagTmp())<3*currentWorker->eps2())
      {
        //half-step to improve f'==target: rb:=rb-(f'-target)/(f'')
        currentWorker->assign(s2_.re_s, f_z.re_s);
        currentWorker->assign(s2_.im_s, f_z.im_s);
        currentWorker->rsub(s2_.re_s, target_f_z.re_s);
        currentWorker->rsub(s2_.im_s, target_f_z.im_s);
        currentWorker->add_double(s2_.re_s, -1);
        currentWorker->assign(s3.re_s, f_zz.re_s);
        currentWorker->assign(s3.im_s, f_zz.im_s);
        s3.recip();
        s2_.mul(&s3);
        rb->add(&s2_);
        break;
      }
#endif
      if (currentWorker->is0(bulbe.f.re) && currentWorker->is0(bulbe.f.im))
        break;
      if (!lagu.eval(period, &bulbe.f, &bulbe.f_z, &bulbe.f_zz))
      {
        *foundMult=1;
        return false;
      };
      rb->sub(&lagu.step);
      if (bulbe.f.getMag_double()<3*currentWorker->eps2())
        break;
    }

    if (cycle==0)
    { //no idea why but first guess, after fixing rb, has f_z either 0+0i at bulb, or 0+-i at cardioid
      s2.assign(&bulbe.f_z);
      double dist_to_0=s2.getMag_double();
      currentWorker->add_double(s2.im, 1);
      double dist_to_ni=s2.getMag_double();
      currentWorker->assign(s2.im, bulbe.f_z.im);
      currentWorker->add_double(s2.im, -1);
      double dist_to_pi=s2.getMag_double();
      if (dist_to_0<0.01)
        nop();
      else if (dist_to_ni<0.01 || dist_to_pi<0.01)
        *is_card=true;
      else if (dist_to_0<0.5)
        nop();
      else if (dist_to_ni<0.5 || dist_to_pi<0.5)
        *is_card=true;
      else
        nop();
    };

    double f_error, fz_error;
    s2.assign(&bulbe.f_z);
    s2.sub(&target_f_z);
    currentWorker->add_double(s2.re, 1);
    f_error=bulbe.f.getMag_double(); //supposed to be 0 but just in case
    fz_error=s2.getMag_double();
    if (!*is_card && f_error<1e-15 && fz_error<1e-4)
    {
      if (!bulbe.eval_multi(period, cb, rb, &target_f_z))
      {
        *foundMult=1;
        return false;
      };
      if (bulbe.multi>1)
      {
        rb->assign(&bulbe.sumA);
        /* don't need them anyway
        currentWorker->sub(f.re_s, rb->re_s);
        currentWorker->sub(f.im_s, rb->im_s);
        currentWorker->add_double(f_z.re_s, -1);
        */
        //reduce period
        *foundMult *= bulbe.multi;
        period /= bulbe.multi;
        did_reduce_period=true;

        //fix rb after moving to guessed root position, we rely on f(cb, rb)==0
        for (int lcycle=0; lcycle<7; lcycle++)
        {
          //x^5=1 -> x:=x-(x^5-1)/(5*x^4)   NestList[#-(#^5-1)/(5*#^4) &, 0.6+0.8I, 10]
          //x=a^(1/5) |y=x/a| a*y=a^(1/5) y=a^(-4/5)  y^5=a^-4  y:=y-(y^5-a^-4)/(5y^4)  y:=(4/5)*y+a^-4/(5y^4)
          //                                          y^-5=a^4  y:=y-(y^-5-a^4)/(-5*y^-6)=y*(6-a^4*y^5)/5   x=a*y
          if (!bulbe.eval_zz(period, cb, rb))
          {
            *foundMult=1;
            return false;
          };
          bulbe.f.sub(rb);
          currentWorker->add_double(bulbe.f_z.re, -1);
          if (currentWorker->is0(bulbe.f.re) && currentWorker->is0(bulbe.f.im))
            break;
          if (!lagu.eval(period, &bulbe.f, &bulbe.f_z, &bulbe.f_zz))
          {
            *foundMult=1;
            return false;
          };
          rb->sub(&lagu.step);
          if (bulbe.f.getMag_double()<3*currentWorker->eps2())
            break;
        }

        //change target f_z to multi-th root of 1, using root near first_multi
        //could do a lot of tricks but let's just use newton to find the root of newf_z^multi=1
        //or rather new_target^multi=old_target
        //starting from bulb.bulbe.first_multi
        //"next_target" = "bulbe.first_multi" MandelMath::complex next_target(currentWorker, &bulbe.first_multi);
        for (int multicyc=0; multicyc<10; multicyc++)
        {
          //a=target_f_z  y=next_target  y:=y+y*(1-a^4*y^5)/5   x=a*y
          s3.assign(&bulbe.first_multi);
          s3.mul(&target_f_z);
          s2.assign(&s3);
          for (int i=2; i<bulbe.multi; i++)
            s3.mul(&s2);
          s3.mul(&bulbe.first_multi); //a^4*y^5
          currentWorker->chs(s3.re);
          currentWorker->chs(s3.im);
          currentWorker->add_double(s3.re, 1); //(1-a^4*y^5)
          currentWorker->zero(s2.re, bulbe.multi);
          currentWorker->recip(s2.re);
          currentWorker->mul(s3.re, s2.re);
          currentWorker->mul(s3.im, s2.re); //(1-a^4*y^5)/5
          double dist1=s3.getMag_double();
          s3.mul(&bulbe.first_multi); //y*(1-a^4*y^5)/5
          bulbe.first_multi.add(&s3); //y+=y*(1-a^4*y^5)/5
          if (dist1<3*currentWorker->eps2())
            break;
        }
        target_f_z.mul(&bulbe.first_multi);
      }
    }
    /*s1 does not hold old deltar any more
    currentWorker->sub(s1.re_s, rb->re_s);
    currentWorker->sub(s1.im_s, rb->im_s);
    currentWorker->sub(s1.re_s, deltar.re_s);
    currentWorker->sub(s1.im_s, deltar.im_s); //should be around 0*/
    if (f_error<3*currentWorker->eps2() &&
        fz_error<3*(1+bulbe.f_zz.getMag_double())*currentWorker->eps2() &&
        (did_reduce_period || *is_card))
    {
      nop(); //ok
      break;
    };
    nop();
  }
  baseZC->assign(&bulbe.f_zc);
  baseCC->assign(&bulbe.f_cc);
  return (*foundMult > 1) || *is_card;







/*

  //suppose the center is exact enough
  //2) find bulb base guess c = xc+1/(f_zc+f_zz)   , derivatives at (z=xc,c=xc) (from estimateInterior)
  bulb.bulbe.eval2(period, xc, xc);
  currentWorker->sub(&bulb.bulbe.f_re, xc->re_s);//is 0
  currentWorker->sub(&bulb.bulbe.f_im, xc->im_s);
  currentWorker->add_double(&bulb.bulbe.f_z_re, -1); //is -1

  //s2:=f_zc^2-f_zz*f_cc   always 0 at bulb and card center
  currentWorker->assign(s2_.re_s, f_zz.re_s);
  currentWorker->assign(s2_.im_s, f_zz.im_s);
  s2_.mul(&f_cc_);
  currentWorker->assign(s1.re_s, f_zc.re_s);
  currentWorker->assign(s1.im_s, f_zc.im_s);
  s1.sqr();
  currentWorker->sub(s2_.re_s, s1.re_s);
  currentWorker->sub(s2_.im_s, s1.im_s);
  //g always=0+0i at the bulb center (or it should) and 1+0i at cardioid
  double cardioid_discrim=currentWorker->toDouble(s2_.getMagTmp());
  double precis=currentWorker->toDouble(f_zc.getMagTmp());
  if (cardioid_discrim<precis*1e-10)
    *is_card=false;
  else if (cardioid_discrim<1-precis*1e-10)
    *is_card=false; //?
  else if (cardioid_discrim<1+precis*1e-10)
    *is_card=true;
  else
    *is_card=true; //?



  //taken from estimateInterior: center-base=1/(f_zc+f_zz) ... unless it's a cardioid...
  currentWorker->assign(s1.re_s, f_zz.re_s);
  currentWorker->assign(s1.im_s, f_zz.im_s);
  s1.add(&f_zc);
  double recip_mag=currentWorker->toDouble(s1.getMagTmp());
  //if (recip_mag<1) and (recip_mag>=1e-30) then
  //  recip_mag:=recip_mag; //sure//really?
  if (recip_mag<1e-30) //test for about 1.0 would be enough; it gets larger for small mandels
  {
    *foundMult=1;
    return false;
  };
  s1.recip_prepared();
  currentWorker->assign(cb->re_s, xc->re_s);
  currentWorker->assign(cb->im_s, xc->im_s);
  cb->add(&s1); //cb=xc-1/(f_zc+f_zz), sign is a bit unclear
  if (*foundMult==0)
  {
    currentWorker->assign(rb->re_s, cb->re_s);
    currentWorker->assign(rb->im_s, cb->im_s);
  }
  *foundMult=1;

  //fix r a little
  //move c to estimate of bulb base
  //  and adjust r accordingly
  //repeat
  MandelMath::complex B(currentWorker, &bulb.B_re, &bulb.B_im, true);
  MandelMath::complex C(currentWorker, &bulb.C_re, &bulb.C_im, true);
  MandelMath::complex inte(currentWorker, &interior.inte_re, &interior.inte_im, true);
  for (int cycle=0; cycle<10; cycle++)
  {
    for (int lcycle=0; lcycle<5; lcycle++) //TODO: stop at convergence or just improve everything
    {
      if (!bulb.bulbe.eval2(period, cb, rb))
      {
        *foundMult=1;
        return false;
      };
      currentWorker->sub(f.re_s, rb->re_s);
      currentWorker->sub(f.im_s, rb->im_s);
      currentWorker->add_double(f_z.re_s, -1);
      if (!bulb.lagu.eval(period, &f, &f_z, &f_zz))
      {
        *foundMult=1;
        return false;
      };
      currentWorker->sub(rb->re_s, &bulb.lagu.step_re);
      currentWorker->sub(rb->im_s, &bulb.lagu.step_im);
    }
    int ei=estimateInterior(period, cb, rb);
    if (ei==0)
    {
      *foundMult=1;
      return true;
    }
    else if (ei<=0)
    {
      *foundMult=1;
      return false;
    };
    if (currentWorker->toDouble(&interior.inte_abs)>=4.5)
    {
      *foundMult=1;
      return false;
    };
    //correct step seems to be inte/2 (correct step is guaranteed between inte and inte/4)
    currentWorker->lshift(&interior.inte_re, -1);
    currentWorker->lshift(&interior.inte_im, -1);
    currentWorker->add(cb->re_s, &interior.inte_re);
    currentWorker->add(cb->im_s, &interior.inte_im);

    //f-xc=(cb-xc)*f_c+(rb-xc)*f_z
    //rb-xc=(cb-xc)*f_c+(rb-xc)*f_z
    //(rb-xc)=(cb-xc)*f_c/(1-f_z)   will crash when we reach bulb base because f_z==0

    //rb-xc=(cb-xc)*f_c+(rb-xc)*f_z+(cb-xc)^2*f_cc/2+(rb-xc)^2*f_zz/2+f_zc*(rb-xc)*(cb-xc)
    //-(rb-xc)^2*f_zz/2+(rb-xc)*(1-f_z-f_zc*(cb-xc))=(cb-xc)*f_c+(cb-xc)^2*f_cc/2
    //-(rb-xc)^2*f_zz/2+(rb-xc)*(1-f_z-f_zc*(cb-xc))=(cb-xc)*f_c+(cb-xc)^2*f_cc/2
    //A=f_zz/2  B=f_z+f_zc*(cb-xc)-1  C=(cb-xc)*f_c+(cb-xc)^2*f_cc/2
    currentWorker->assign(B.re_s, f_zc.re_s);
    currentWorker->assign(B.im_s, f_zc.im_s);
    B.mul(&inte);
    B.add(&f_z);
    currentWorker->add_double(B.re_s, -1);

    currentWorker->assign(C.re_s, f_cc_.re_s);
    currentWorker->assign(C.im_s, f_cc_.im_s);
    currentWorker->lshift(C.re_s, -1);
    currentWorker->lshift(C.im_s, -1);
    C.mul(&inte);
    currentWorker->add(C.re_s, &bulb.bulbe.f_c_re);
    currentWorker->add(C.im_s, &bulb.bulbe.f_c_im);
    C.mul(&inte);

    //let's try in doubles first
    //TODO: pretty much completely wrong
    double r_step_re, r_step_im;
    MandelMath::complex_double_quadratic(&r_step_re, &r_step_im,
                                         currentWorker->toDouble(f_zz.re_s)/2, currentWorker->toDouble(f_zz.im_s)/2,
                                         currentWorker->toDouble(B.re_s)/2, currentWorker->toDouble(B.im_s)/2,
                                         currentWorker->toDouble(C.re_s), currentWorker->toDouble(C.im_s));
    currentWorker->add_double(rb->re_s, -r_step_re);
    currentWorker->add_double(rb->im_s, -r_step_im);
  }
  currentWorker->assign(baseZC->re_s, &bulb.bulbe.f_zc_re);
  currentWorker->assign(baseZC->im_s, &bulb.bulbe.f_zc_im);
  currentWorker->assign(baseCC->re_s, &bulb.bulbe.f_cc_re);
  currentWorker->assign(baseCC->im_s, &bulb.bulbe.f_cc_im);
  return true;
  */
}

MandelEvaluator::Bulb::Bulb(MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN), currentWorker(allocator->worker),
  bulbe(&self_allocator, nullptr), rb(&self_allocator), cb(&self_allocator), xc(&self_allocator),
  cbx(&self_allocator), rbx(&self_allocator), baseZC(&self_allocator), baseCC(&self_allocator),
  s1(&self_allocator), s2(&self_allocator), s3(&self_allocator), deltac(&self_allocator), deltar(&self_allocator),
  B(&self_allocator), C(&self_allocator), dbg_first_cb(&self_allocator), dbg_first_rb(&self_allocator),
  target_f_z(&self_allocator), dbg_guessmult(0), lagu(&self_allocator, nullptr)
{
  assert(self_allocator.checkFill());
}

void MandelEvaluator::Bulb::fixRnearBase(MandelMath::complex *r, const MandelMath::complex *c, int period, int *mult)
{ //"cleverFix" in old code
  //TODO: cycles unused?
  /*MandelMath::complex rb(currentWorker, &bulb.rb_re_, &bulb.rb_im, false);
  MandelMath::complex cb(currentWorker, &bulb.cb_re, &bulb.cb_im, false);
  MandelMath::complex xc(currentWorker, &bulb.xc_re, &bulb.xc_im, false);
  MandelMath::complex baseZC_(currentWorker, &bulb.baseZC_re, &bulb.baseZC_im, false);
  MandelMath::complex baseCC_(currentWorker, &bulb.baseCC_re, &bulb.baseCC_im, false);
  MandelMath::complex s1(currentWorker, &bulb.s1_re, &bulb.s1_im, false);
  MandelMath::complex s2(currentWorker, &bulb.s2_re_, &bulb.s2_im_, false);
  MandelMath::complex cbx_(currentWorker, &bulb.cbx_re, &bulb.cbx_im, false);
  MandelMath::complex rbx_(currentWorker, &bulb.rbx_re, &bulb.rbx_im, false);*/
  if (*mult<=1)
    return;
  rb.assign(r);
  bool is_card=false;
  int foundMult=1;
  bool baseFound=findBulbBase(period, c, &cb, &rb, &xc, &baseZC, &baseCC, &is_card, &foundMult);
  if (!baseFound)
    return;
  if (foundMult<=1)
    return;
  if (*mult!=foundMult)
    *mult=foundMult;
  //s1=c-cb
  s1.assign(c);
  s1.sub(&cb);
  /*currentWorker->assign(rbx.re_s, s1.re_s);
  currentWorker->assign(rbx.im_s, s1.im_s);
  rbx.add(&xc);*/
  //double ratio=currentWorker->toDouble(s1.getMagTmp())/currentWorker->toDouble(cbx.getMagTmp());
  //if (ratio<0.05)
  //s2=s1*baseZC ~= root_z
  s2.assign(&s1);
  s2.mul(&baseZC);
  r->assign(&rb);
  if (!currentWorker->isl0(s2.re) || (*mult==2)) //=0 for baseZC=0 at period=1
  { //above base
    cbx.assign(&cb);
    cbx.sub(&xc);
    rbx.assign(&rb);
    rbx.sub(&xc);
    s2.assign(&cbx);
    s2.recip();
    s2.mul(&rbx);
    double tmpre=s2.getMag_double();
    if (tmpre==0)
    {
      //r:=rb
    }
    else
    {
      double tmpim=std::atan2(currentWorker->toDouble(s2.im), currentWorker->toDouble(s2.re));
      if (tmpim<-M_PI)
        tmpim+=2*M_PI;
      else if (tmpim>=M_PI)
        tmpim-=2*M_PI;
      if (*mult==2)
      {
        tmpre=exp(log(tmpre)/4);
        tmpim/=2;
      }
      else
      {
        tmpre=exp(log(tmpre)/(2*(*mult-1)));
        tmpim/=(*mult-1);
      }
      s2.zero(-tmpre*cos(tmpim), -tmpre*sin(tmpim));
      s2.mul(&rbx);
      //r:=rb+s2
      r->add(&s2);
    }
  }
  else
  { //we're below the bulb base, move close to the base's root
    //dz f_zc+dc/2 f_cc=0
    //dz=-dc/2 f_cc/f_zc
    s1.mul(&baseCC);
    s2.assign(&baseZC);
    s2.recip();
    s2.mul(&s1);
    s2.lshift(-1);
    r->sub(&s2);
  }
}

//result 0..derivatives or value too large, or other fail (divide by 0)
//result>0 .. tried to return multiplicity but really returns just 1 (1 or >=3) or 2 (mult==2)
int MandelEvaluator::newton(int period, const MandelMath::complex *c, MandelMath::complex *r, const bool fastHoming, const int suggestedMultiplicity) //returns multiplicity
{ //TODO: suggestedMulti = maximumMultip ?
  double bestfm=1e10; //TODO: actually bestgm? g(z)=f(z)-z
  newt.bestr.assign(r);
  bool movedOff=false;
  //double accyBound=3e-28/(period*period);
  //was for 80b floats double accyBound2=3e-39*period/log(1+period)*1.5; //1.5=magic
  //double accyBound2=1.23e-32*period/log(1+period)*1.5; //1.5=magic
  double order1; // 1/highest power in the polynomial, 2^period in case of mandelbrot set
  double r_mag_rough;
  int maxm;
  {
    newtres_.first_guess_lagu.assign(r);
    newtres_.first_guess_newt.assign(r);
    double r_re=currentWorker->toDouble(r->re);
    double r_im=currentWorker->toDouble(r->im);
    newtres_.first_fejer_re=r_re; newtres_.first_fejer_im=r_im;
    newtres_.first_naive1_re_=r_re; newtres_.first_naive1_im=r_im;
    newtres_.first_naive2_re=r_re; newtres_.first_naive2_im=r_im;
    newtres_.first_naive_re=r_re; newtres_.first_naive_im=r_im;
    newtres_.naiveChoice=NewtonNaiveChoice::ncClose;
    newtres_.first_neumaier1_re_=r_re; newtres_.first_neumaier1_im_=r_im;
    newtres_.first_neumaier2_re=r_re; newtres_.first_neumaier2_im=r_im;
    newtres_.first_lagu1_re=r_re; newtres_.first_lagu1_im=r_im;
    newtres_.first_lagu1o_re=r_re; newtres_.first_lagu1o_im=r_im;
    newtres_.firstMu_re_=1; newtres_.firstMu_im=0; //newtres_.firstM=1;
    newtres_.firstMum_re_=1; newtres_.firstMum_im_=0;
    newtres_.accy_tostop=1;
    newtres_.accy_multiplier=1;

    //ideally, I want max(abs(r_re), abs(c_re)) then round up to next power of 2
    //but neither ilogb or frexp can do that so I round 1 to 2, 2 to 4
    int lor=std::ilogb(r_re); //3->1 2->1 1->0  0.75->-1  0->-max
    /*does not seem to affect precision of root? int loc=std::ilogb(currentWorker->toDouble(c->re_s));
    if (lor<loc)
      lor=loc;*/
    if (lor<-2)
      lor=-2;
    r_re=ldexp(2, lor); //1->4  0->2  -2->0.5
    lor=std::ilogb(r_im);
    /*loc=std::ilogb(currentWorker->toDouble(c->im_s));
    if (lor<loc)
      lor=loc;*/
    if (lor<-2)
      lor=-2;
    r_im=ldexp(2, lor);
    r_mag_rough=r_re*r_re+r_im*r_im;
  }
  //double accyBound=r_mag_rough*currentWorker->eps2()*period; //eps*sqrt(period) as eps bleeds out// 3e-28/(period*period);
  if (period<5)
  {
    maxm=ldexp(1, period-1); //in theory up to n-1 but for Mandelbrot that's rarely the case
    order1=ldexp(1, -period);
  }
  else if (period<1024)
  {
    maxm=15;
    order1=ldexp(1, -period);
  }
  else
  {
    maxm=15;
    order1=0;
  }
  //maxm=1;
  //int multiplicity1=1;
  int lastm=1;
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
  /*complex f_r(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  complex fz_r(currentWorker, &newtres_.fz_r_re_, &newtres_.fz_r_im_, true);
  complex fzz_r(currentWorker, &newt.fzz_r_re, &newt.fzz_r_im, true);
  complex tmp1(currentWorker, &newt.tmp1_re, &newt.tmp1_im, true);
  //complex fzfix(currentWorker, &newt.fzfix_re, &newt.fzfix_im, true);
  complex laguH(currentWorker, &newt.laguH_re, &newt.laguH_im, true);
  complex laguG(currentWorker, &newt.laguG_re, &newt.laguG_im, true); //fzf
  complex laguG2(currentWorker, &newt.laguG2_re, &newt.laguG2_im, true); //G^2
  complex laguX(currentWorker, &newt.laguX_re, &newt.laguX_im, true); //Laguerre step
  complex newtX(currentWorker, &newt.newtX_re, &newt.newtX_im, true); //Laguerre step
  complex fzzf(currentWorker, &newt.fzzf_re, &newt.fzzf_im, true);*/
  for (int newtonCycle=0; newtonCycle<50; newtonCycle++)
  {
    newtres_.cyclesNeeded=newtonCycle;
    if ((movedOff) && (newtonCycle>10) && (order1>=0))
    {                                    //  p m -> p
      order1=-1;                         //  2 2    1
      bestfm=1e10;                       //  4 3    2
      //multiplicity1=1;                   //  4 5    1
    };
    newt.f_r.assign(r);
    newtres_.fz_r.zero(1, 0);
    newt.fzz_r.zero(0, 0);
    //TODO: can we skip computing fzz_r if order1<0? and remember last valid multiplicity or set it to 1
    //always half of eps_cumul10   double eps_cumul05=0.5;
    double eps_cumul10=1;
    double eps_cumul=1; //calc itself looks good, lemme start with 1 instead of 0
    /* { //I need to increase eps_cumul somewhere but don't know where, all looks good
      double f_rough=currentWorker->radixfloor(f_r.re_s, c->re_s);
      eps_cumul+=f_rough*f_rough;
      f_rough=currentWorker->radixfloor(f_r.im_s, c->im_s);
      eps_cumul+=f_rough*f_rough;
    }*/
    //very close to eps_cumul10   double err_cumul=0; //maybe sum(ln fz)/ln(final_fz) but similar to sum(1/fz)*final_fz
    int err_simple=0;
    double fc_re=0, fc_im=0;
    for (int i=0; i<period; i++)
    {
      double fzz_r_mag=newt.fzz_r.getMag_double();
      double fz_r_mag=newtres_.fz_r.getMag_double();
      double f_r_mag=newt.f_r.getMag_double();
      if (fzz_r_mag+fz_r_mag+f_r_mag>1e60)
        return 0;
      //fzz:=2*(fz*fz + f*fzz)
      newt.fzz_r.mul(&newt.f_r);
      newt.tmp1.assign(&newtres_.fz_r);
      newt.tmp1.sqr();
      newt.fzz_r.add(&newt.tmp1);
      newt.fzz_r.lshift(1);
      //fz:=2*f*fz
      newtres_.fz_r.mul(&newt.f_r);
      newtres_.fz_r.lshift(1);

      //eps_cumul05=4*eps_cumul05*f_r_mag+0.5;
      eps_cumul10=4*eps_cumul10*f_r_mag+1;
      eps_cumul=4*eps_cumul*f_r_mag;
      //err_cumul+=1/currentWorker->toDouble(fz_r.getMagTmp());
      if (fz_r_mag<=1)
        err_simple++;;
      //fc=2*fc*f+1
      double f_re=currentWorker->toDouble(newt.f_r.re);
      double f_im=currentWorker->toDouble(newt.f_r.im);
      double tmp=2*(f_re*fc_re-f_im*fc_im)+1;
      fc_im=2*(f_im*fc_re+f_re*fc_im);
      fc_re=tmp;

      //f:=f^2+c
      newt.f_r.sqr();
      newt.f_r.add(c);
      double f_rough=currentWorker->radixfloor(newt.f_r.re, c->re);
      eps_cumul+=f_rough*f_rough;
      f_rough=currentWorker->radixfloor(newt.f_r.im, c->im);
      eps_cumul+=f_rough*f_rough;
      /*
      (a+-c+bi+-di)^2=
        (a+bi)^2=aa+abi+abi+bbii=a^2-b^2+2abi
        (a+c+bi+di)^2=a^2-b^2+c^2-d^2+2ac-2bd +2i(ab+ad+bc+cd)
        (a+c+bi-di)^2=a^2-b^2+c^2-d^2+2ac+2bd +2i(ab-ad+bc-cd)
        (a-c+bi+di)^2=a^2-b^2+c^2-d^2-2ac-2bd +2i(ab+ad-bc-cd)
        (a-c+bi-di)^2=a^2-b^2+c^2-d^2-2ac+2bd +2i(ab-ad-bc+cd)
      for c=d:
        (a+c+bi+di)^2=a^2-b^2+2ac-2bd +2i(ab+ad+bc+cd)  new c=max(+-2ac+-2bd)=maxabs(2ac+-2bd)=abs(2ac)+abs(2bd)
        (a+c+bi-di)^2=a^2-b^2+2ac+2bd +2i(ab-ad+bc-cd)  new d=maxabs(ad+bc+cd; ad-bc+cd; ad-bc-cd; ad+bc-cd)
        (a-c+bi+di)^2=a^2-b^2-2ac-2bd +2i(ab+ad-bc-cd)  new d=abs(ad)+maxabs(bc+-cd)=abs(ad)+abs(bc)+abs(cd)  cd is small
        (a-c+bi-di)^2=a^2-b^2-2ac+2bd +2i(ab-ad-bc+cd)
        new c=2c(abs(a)+abs(b))
        new d=c(abs(a)+abs(b))
      */
    }
    //double fz_r_mag=currentWorker->toDouble(fz_r.getMagTmp()); //ff1m in original code
    //g(r)=f(r)-r, gz(r)=fz(r)-1
    newt.f_r.sub(r);
    currentWorker->add_double(newtres_.fz_r.re, -1);
    double g_r_mag=newt.f_r.getMag_double();
    double gz_r_mag=newtres_.fz_r.getMag_double();
    //err_cumul*=fz_r_mag; //err_cumul->25774, eps_cumul10=25775, eps_cumul05=12888, simple=84
    if (eps_cumul<0.25) //only happens near c=0 r=0 where the bound keeps falling faster than g_r_mag
      eps_cumul=0.25;
    newtres_.accy_tostop=eps_cumul;//r_mag_rough*eps_cumul10;
    newtres_.accy_multiplier=1/gz_r_mag;
#if CLEVER_FIX
//c=-0.7499 p=2
//  ideally, r=-0.5+-0.01i who are repelling  (and +0.5+-sqrt(0.9999) who are repel and attr)
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
      fixRnearBase(r, c, 0, period, &clever.mult);
      continue;
    };
#endif
    if (g_r_mag==0)  //7e-33..4e-40 does not need more; much..5e-38 needs more
    { //r is good enough already      (f_c.re*f_c.re+f_c.im*f_c.im)/(f_zc.re*f_zc.re+f_zc.im*f_zc.im)
      return lastm;
    };
    if (gz_r_mag==0)
    {
      if (triedZeroGzrm)
        return 0;
      triedZeroGzrm=true;
    }
    else
    {
      if (gz_r_mag<1e-39)
      { //one more multi-iter?
        //zero ffm is solved with triedZeroFfm
        if (clever.mult>0)
          return clever.mult;
        else
          return 1;
      };

      //new conditions
      //the one legit reason to end: step<2^-53/|f'| (for |f'|<1) exactly because step*|f'|=2^-53
      //    |f|/|f'|<2^-53/|f'|
      //    |f|^2<2^-106=1.23e-32
      //if (g_r_mag<1.0*r_mag_rough*currentWorker->eps2()) //maybe up to (1+gz_r_mag)*r_mag*eps2*log2(period)
      //if (g_r_mag<eps_cumul10*r_mag_rough*currentWorker->eps2()) //maybe up to (1+gz_r_mag)*r_mag*eps2*log2(period)
      if (g_r_mag<newtres_.accy_tostop*currentWorker->eps2())
        return lastm;
        //for per=11, rough=1, gz_m=0.13: around 6.3<x/eps<3.2 (biggest that won't improve < x <smallest that will improve)
        //for per=20, rough=0.5, gz=0.0216: 0.063<x<0.5
        //for per=94, rough=1, gz_m=0.56: 9.3<x<~61; gz_m=6e-5: 9.1<x<7.3
        //for per=63(7), rough=4, gz_m=1.5e-6: 9.9<x<5.6

#if 0
      if (//maybe eps2*period but 1e-29 is too much ((gz_r_mag>1e-6) && (g_r_mag/gz_r_mag<accyBound)) || //((fm<NEWTON_EPSILON) and (ffm>1e-6)) or //check just f on single roots
          //(should be handled with newer eps computation   (g_r_mag<1e-20) && (g_r_mag/gz_r_mag<accyBound)) || //for high period - and high derivative - we cannot minimize f given the finite precision of root
           ((gz_r_mag<1e-6) && (g_r_mag/gz_r_mag<2.1e-30)))   //check f/ff on multiple roots
      { //r is good enough already
        return lastm;
      };
      if ((g_r_mag<2e-38) && (gz_r_mag<3e-9)) //ffm: >8e-10
      { //close to be but not multiple
        if (clever.mult>0)
          return clever.mult;
        else
          return 1;
      };
#endif
      /*if ((fz_r_mag>0.5) && (g_r_mag/fz_r_mag<accyBound2))
      { //cannot improve because of limited precision
        return lastm;
      };*/
    };
    if (currentWorker->isequal(r->re, newt.bestr.re) &&
        currentWorker->isequal(r->im, newt.bestr.im) &&
        (bestfm<1e10) && (newtonCycle>2)) //Lagu can cycle in first 2 cycles
    { //Laguerre can cycle (at least my version), e.g. per=2, c=-0.6640625-0.015625i, r=-0.614213552325963974-0,0179149806499481201i
      if (g_r_mag<6*eps_cumul10*r_mag_rough*currentWorker->eps2()) //should be tested above but maybe use different margin here?
        return lastm;
      return 0; //just fail and try again next time
    };
    if (g_r_mag<bestfm)
    {
      bestfm=g_r_mag;
      newt.bestr.assign(r);
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
    newt.tmp1.assign(&newt.f_r);
    newt.tmp1.recip();    //1/f
    newt.laguG.assign(&newtres_.fz_r);
    newt.laguG.mul(&newt.tmp1); //laguG = f'/f
    newt.fzzf.assign(&newt.fzz_r);
    newt.fzzf.mul(&newt.tmp1); //f''/f

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
      newt.laguG2.assign(&newt.laguG);
      newt.laguG2.sqr();    //G^2
      newt.laguH.assign(&newt.laguG2);
      newt.laguH.sub(&newt.fzzf); //G^2-fzzf = f'^2/f^2-f''*f/f^2 = laguH
      //currentWorker->assign(tmp1.re_s, laguG2.re_s);
      //currentWorker->assign(tmp1.im_s, laguG2.im_s);
      int m=1;
      {
        /*double G2HT_re=currentWorker->toDouble(laguG2.mulreT(&laguH));
        double H_mag=currentWorker->toDouble(laguH.getMagTmp());
        //turns out that if mu=m then mu=m=G^2/H
        //1.5*mag(H)>Re(G^2*H^T) ... m=1
        //300*mag(H)<Re(G^2*H^T) ... m=300
        //mag(H)<Re(G^2*H^T)*order1 ... m=1/order1
        if (1.5*H_mag>=G2HT_re) //>= so we don't divide G^2/H if H=G=0
          m=1;
        else if ((clever.mult+0.5)*H_mag<=G2HT_re)
          m=1; //best practice is to use m=1 if H=0   clever.mult;
        else if (H_mag*(maxm-0.5)<G2HT_re)
          m=maxm;
        else
          m=qRound(G2HT_re/H_mag);*/

        //solve for m=mu:   m=G^2/H
        double mum_re=1, mum_im=0; //better than mu?
        double h_re=currentWorker->toDouble(newt.laguH.re);
        double h_im=currentWorker->toDouble(newt.laguH.im);
        double h_mag=h_re*h_re+h_im*h_im;
        double g2_re=currentWorker->toDouble(newt.laguG2.re);
        double g2_im=currentWorker->toDouble(newt.laguG2.im);
        if (h_mag>0.01)
        { //h_mag ok
          mum_re=(g2_re*h_re+g2_im*h_im)/h_mag;
          mum_im=(g2_im*h_re-g2_re*h_im)/h_mag;
        };
        if (newtonCycle==0)
        {
          newtres_.firstMum_re_=mum_re;
          newtres_.firstMum_im_=mum_im;

          /*if ((2*maxm+1)*H_mag<=G2HT_re)
            newtres_.firstM=2*maxm+1; //need to see a bit above order//2*maxm ~ order but without overflow
          else
            newtres_.firstM=G2HT_re/H_mag;*/
        };
      }

      //m= some func of mu where mu is solution of ((1-1/n)*H/G^2-1/n) mu^2 + 2*mu/n -1=0
      //with m as input:                           ((1-m/n)*H/G^2-1/n) mu^2 + m/n 2*mu -m = 0
      double G2_mag=newt.laguG2.getMag_double();
      if (G2_mag<0.01)
      { //G2_mag bad
        m=1;
        if (newtonCycle==0)
        {
          newtres_.firstMu_re_=1;
          newtres_.firstMu_im=0;
        };
      }
      else
      {
        newt.laguX.assign(&newt.laguG2);
        currentWorker->chs(newt.laguX.im);
        newt.laguX.mul(&newt.laguH);
        double a_re=currentWorker->toDouble(newt.laguX.re)/G2_mag*(1-order1)-order1;
        double a_im=currentWorker->toDouble(newt.laguX.im)/G2_mag*(1-order1);
        double mu_re, mu_im;
        MandelMath::complex_double_quadratic(&mu_re, &mu_im, a_re, a_im, order1, 0, -1, 0);
        if (newtonCycle==0)
        {
          newtres_.firstMu_re_=mu_re;
          newtres_.firstMu_im=mu_im;
        };
        if (!(mu_re>=1.3)) //also m=1 if mu_re is NaN    (mu_re<1.3)
          m=1;
        else {/*if (abs(mu_im)>mu_re/2)
          m=1;
        else
        {
          double mu_mag=mu_re*mu_re+mu_im*mu_im;
          m=qRound(sqrt(mu_mag)); //or just round mu_re?
          */
          m=qRound(mu_re);
          if (m>maxm)
            m=maxm;
        }
      }
      if (newtonCycle==0)
      {
        //Fejer bound: smaller solution x of
        //fzz/(n-1) x^2+2 fz x + n f=0
        //x=y*n
        //fzz*n/(n-1) y^2+2 fz y + f=0

        double r_re=currentWorker->toDouble(r->re);
        double r_im=currentWorker->toDouble(r->im);
        //numbers are small but don't need precision so let's do it in double
        double a_re=currentWorker->toDouble(newt.fzz_r.re)/(1-order1);
        double a_im=currentWorker->toDouble(newt.fzz_r.im)/(1-order1);
        double fz_re=currentWorker->toDouble(newtres_.fz_r.re);
        double fz_im=currentWorker->toDouble(newtres_.fz_r.im);
        double f_re=currentWorker->toDouble(newt.f_r.re);
        double f_im=currentWorker->toDouble(newt.f_r.im);
        MandelMath::complex_double_quadratic(
              &newtres_.first_fejer_re, &newtres_.first_fejer_im,
              a_re, a_im,
              fz_re, fz_im,
              f_re, f_im);
        newtres_.first_fejer_re=r_re+ldexp(newtres_.first_fejer_re, period);
        newtres_.first_fejer_im=r_im+ldexp(newtres_.first_fejer_im, period);

        //Batra's bound https://www.tuhh.de/ti3/paper/rump/Ru03c.pdf theorem 3.8
          //but only for real coefficients
        //|fz r|-|f + fzz/2 r^2|=0, find r
        //sqrt(fz fz^T) r=sqrt((f + fzz/2 r^2)(f^T + fzz^T/2 r^2))
        //sqrt(fz fz^T) r=sqrt((|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4))
        //(a+bi)(c-di)+(a-bi)(c+di)=2ac+2bd=2 Re(f fzz^T)
        //|fz|^2 r^2=|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4
        //0=|f|^2+ (Re(f^T fzz)-|fz|^2) rr + |fzz|^2/4 rr^2    r=sqrt(rr)

        /*MandelMath::complex_double_quadratic(&newtres_.first_batra, &a_im,
            currentWorker->toDouble(fzz_r.getMagTmp())/4, 0,
            (currentWorker->toDouble(f_r.mulreT(&fzz_r))-currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
            currentWorker->toDouble(f_r.getMagTmp()), 0);
        if (newtres_.first_batra>=0)
          newtres_.first_batra=sqrt(newtres_.first_batra);*/

        //https://ur.booksc.eu/book/5736333/a5b588
        //ZAMM - Journal of Applied Mathematics and Mechanics / Zeitschrift für Angewandte Mathematik und Mechanik
        //1988 Vol. 68; Iss. 6
        //Dr. A. Neumaier: An Existence Test for Root Clusters and Multiple Roots
        //fi(c, r, alpha)=r abs(re((f(c+r e^ialpha)-f(c))/(c+r e^ialpha)))-abs(f(c))
        //  addition from https://ur.booksc.eu/book/5736333/a5b588 remark 3:
        //  f needs to be divided (or rotated) by f' first to make f' real
        //for all alpha, which r makes fi==0 ?
        //abs(re(f'*r+f''/2 r^2 e^ialpha))=abs(f)
        //for max re(f'*r+f''/2 r^2 e^ialpha), we need max re(f'+f''/2 r e^ialpha) because r is real
        //f'' e^ialpha=real
        //e^ialpha=f''^T/sqrt(f'' f''^T)=sqrt(f''^T/f'')
        //abs(re(f'*r+ r^2/2 sqrt(f'' f''^T)))-abs(f)=0
        //r*abs(re(f'))+ r^2/2 sqrt(f'' f''^T)-abs(f)=0
        /*if (currentWorker->isle0(fz_r.re_s))
          MandelMath::complex_double_quadratic(&newtres_.first_batra, &a_im,
              +sqrt(currentWorker->toDouble(fzz_r.getMagTmp()))/2, 0,
              //currentWorker->toDouble(fz_r.re_s)/2, 0,
              sqrt(currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
              +sqrt(currentWorker->toDouble(f_r.getMagTmp())), 0);
        else*/
        MandelMath::complex_double_quadratic(&newtres_.first_neumaier1_re_, &newtres_.first_neumaier1_im_,
            -sqrt(newt.fzz_r.getMag_double())/2, 0,
            sqrt(newtres_.fz_r.getMag_double())/2, 0,
            -sqrt(newt.f_r.getMag_double()), 0);

        /* Neumaier for k=2:
        Re(f''(z)/2) > |f| r^-2 + |f'| r^-1      r real, r>0
        f''(z)=f''+(z-z0)f'''+...
        |f''+r|f'''|+...|/2 r^2 > |f| + |f'| r
        ...+|f'''|r^3/2+|f''| r^2/2 > |f| + |f'|r
                        |f''| r^2/2 > |f| + |f'|r
        gives always r<0 but that's the wrong root
        |f''| r^2/2 - |f'|r - |f| =0
        2*|f'|/|f''|+-sqrt(4*|f'|^2/|f''|^2-8*|f|/|f''|)
        |f'|/|f''|+-sqrt(|f'|^2/|f''|^2+2*|f|/|f''|)
        works if |f'|^2/|f''|+2*|f|>0 i.e. always
        but r1 always<0, r2>2*|f'|/|f''|
        r2=(|f'|+sqrt(|f'|^2+2*|f|*|f''|))/|f''|

        test x(x-1) at 2+i
        f=1+3i f'=2x-1=3+2i f''=2
        r2=(|3+2*I|+sqrt(|3+2*I|^2+2*2*|1+3*I|))/2
        4.33502318885498454
        correct is 2.236
        */
        double fm=sqrt(newt.f_r.getMag_double());
        double fzm=sqrt(newtres_.fz_r.getMag_double());
        double fzzm=sqrt(newt.fzz_r.getMag_double());
        newtres_.first_neumaier2_re=(fzm + sqrt(fzm*fzm+2*fm*fzzm))/fzzm;
        newtres_.first_neumaier2_im=0;

        /* naive: approximate f with c(x-a)^m
        m=f'^2/(f'^2-f f'') = f'^2/f^2/(f'^2/f^2-f''/f)=G^2/H
        x-a=m/(f'/f)=m/G=G/H    looks good if |m_im|<|m_re|
        m*(x-a)=G^3/H^2

        trouble: singularities when f f''=f'^2 -> m=infinity, iteration jumps too far
                                    f'=0 -> m=0, m/(f/f') jumps too little
        */
        /*double g_re=currentWorker->toDouble(laguG.re_s);
        double g_im=currentWorker->toDouble(laguG.im_s);
        double g_mag=g_re*g_re+g_im*g_im;
        if (1e6*H_mag<=g_mag*g_mag)
        {
          newtres_.first_naive_re=currentWorker->toDouble(r->re_s);
          newtres_.first_naive_im=currentWorker->toDouble(r->im_s);
        }
        else
        {
          double g2_re=currentWorker->toDouble(laguG2.re_s);
          double g2_im=currentWorker->toDouble(laguG2.im_s);
          double h_re=currentWorker->toDouble(laguH.re_s);
          double h_im=currentWorker->toDouble(laguH.im_s);

          double m_re=(g2_re*h_re+g2_im*h_im)/H_mag;
          double m_im=(g2_im*h_re-g2_re*h_im)/H_mag;
          //couldn't find smooth function that:
          //1->1 2->2 3->3... 0->1 -1->1 i->1 -i->1
          //esp. since we need to have 1->1 exact and in neigborhood too
          if ((m_re<abs(m_im)*2))
          //if ((m_re<0.9) || (m_re<abs(m_im)*2)) //for m~0, we need something like sqrt(m): m is too small, 1 is too large
          {
            m_re=1;
            m_im=0;
          };
          newtres_.first_naive_re=currentWorker->toDouble(r->re_s)-(m_re*g_re+m_im*g_im)/g_mag;
          newtres_.first_naive_im=currentWorker->toDouble(r->im_s)-(m_im*g_re-m_re*g_im)/g_mag;
        }*/

        /* even naiver: show the 2 roots of c(x-a)(x-b) that have the same f, f', f''
        w.l.o.g. x=0
        c(x^2-(a+b)x+ab)=f''x^2/2+f'x+f
        cx^2-c(a+b)x+cab=f''x^2/2+f'x+f
        -f'/f''+-sqrt(f'^2/f''^2-2*f/f'')

        if x1 close to x2 (relative to x), use (x1+x2)/2 else use x1
        at |x1|=|x2|, 90 degrees..mult~2, use (x1+x2)/2
        at |x1|=|x2|, 60 degrees..mult~1, use x1
        at |x1|=0.8|x2|, 80% weight from x1
        at |x1|=0.5|x2|, 90% weight from x1
        at |x1|=0.3|x2|, use x1
        when x1~x2, correct guess is actually around 0.7 x1
        */
        a_re=currentWorker->toDouble(newt.fzz_r.re)/2;
        a_im=currentWorker->toDouble(newt.fzz_r.im)/2;
        MandelMath::complex_double_quadratic2(&newtres_.first_naive1_re_, &newtres_.first_naive1_im,
                                              &newtres_.first_naive2_re, &newtres_.first_naive2_im,
                                              a_re, a_im, fz_re/2, fz_im/2, f_re, f_im);
        double n2_rmag=1/(newtres_.first_naive2_re*newtres_.first_naive2_re+newtres_.first_naive2_im*newtres_.first_naive2_im);
        //d=naive1/naive2
        double d_re=(newtres_.first_naive1_re_*newtres_.first_naive2_re+newtres_.first_naive1_im*newtres_.first_naive2_im)*n2_rmag;
        double d_im=(newtres_.first_naive1_im*newtres_.first_naive2_re-newtres_.first_naive1_re_*newtres_.first_naive2_im)*n2_rmag;
        double d_mag=(d_re*d_re+d_im*d_im);
        double w1=1, w2=0;
        if (d_re<-0.5) //angle>120deg, even if close in magnitude
        { w1=1; w2=0; newtres_.naiveChoice=NewtonNaiveChoice::ncWide; }
        else if (d_mag<0.3*0.3)
        { w1=1; w2=0; newtres_.naiveChoice=NewtonNaiveChoice::nc03; }
        else if (d_mag<0.5*0.5)
        { w1=0.9; w2=0.1; newtres_.naiveChoice=NewtonNaiveChoice::nc05; }
        else if (d_mag<0.8*0.8)
        { w1=0.8; w2=0.2; newtres_.naiveChoice=NewtonNaiveChoice::nc08; } //or just 1;0
        else if (d_re<-0.1)
        {
          //most problematic case, some times best is [1,0] other times [1,1]
          w1=1; w2=0.0;
          newtres_.naiveChoice=NewtonNaiveChoice::nc100;
        }
        else if (d_re<0)
        {
          //most problematic case, some times best is [1,0] other times [1,1]
          //don't trust M here
          w1=1; w2=0.0;
          newtres_.naiveChoice=NewtonNaiveChoice::nc90_;
        }
        else if (d_re<0.1)
        {
          //most problematic case, some times best is [1,0] other times [1,1]
          //don't trust M here
          w1=1; w2=0.0;
          newtres_.naiveChoice=NewtonNaiveChoice::nc80;
        }
        else if (d_re<0.5)
        {
          //can (try) use M here
          w1=1; w2=0;
          newtres_.naiveChoice=NewtonNaiveChoice::nc60;
        }
        else
        {
          //can (try) use M here
          w1=newtres_.firstMum_re_-1; w2=0;
          newtres_.naiveChoice=NewtonNaiveChoice::ncClose;
        }
        newtres_.first_naive_re=w1*newtres_.first_naive1_re_+w2*newtres_.first_naive2_re;
        newtres_.first_naive_im=w1*newtres_.first_naive1_im+w2*newtres_.first_naive2_im;
        newtres_.first_naive1_re_=r_re+newtres_.first_naive1_re_;
        newtres_.first_naive1_im=r_im+newtres_.first_naive1_im;
        newtres_.first_naive2_re=r_re+newtres_.first_naive2_re;
        newtres_.first_naive2_im=r_im+newtres_.first_naive2_im;
        newtres_.first_naive_re=r_re+newtres_.first_naive_re;
        newtres_.first_naive_im=r_im+newtres_.first_naive_im;

        //Laguerre is the solution of
        //   c=-n  b=f'/f  a=f'^2/f^2*(1-n/m+1/m)-f''/f*(1-n/m)=H*(1-n/m)+G^2/m
        //   G=f'/f   H=G^2-f''/f
        //>> a=H*(1-m/n)-G^2/n  b=m*G/n  c=-m    ok
        /*
        a_re=currentWorker->toDouble(laguH.re_s)*(1-m*order1)-currentWorker->toDouble(laguG2.re_s)*order1;
        a_im=currentWorker->toDouble(laguH.im_s)*(1-m*order1)-currentWorker->toDouble(laguG2.im_s)*order1;
        double b_re=currentWorker->toDouble(laguG.re_s)*m*order1;
        double b_im=currentWorker->toDouble(laguG.im_s)*m*order1;
        MandelMath::complex_double_quadratic(
              &newtres_.first_lagum_re, &newtres_.first_lagum_im,
              a_re, a_im,
              b_re, b_im,
              -m, 0);
        newtres_.first_lagum_re=currentWorker->toDouble(r->re_s)-newtres_.first_lagum_re;
        newtres_.first_lagum_im=currentWorker->toDouble(r->im_s)-newtres_.first_lagum_im;
        */
        a_re=currentWorker->toDouble(newt.laguH.re)*(1-order1)-currentWorker->toDouble(newt.laguG2.re)*order1;
        a_im=currentWorker->toDouble(newt.laguH.im)*(1-order1)-currentWorker->toDouble(newt.laguG2.im)*order1;
        double b_re=currentWorker->toDouble(newt.laguG.re)*order1;
        double b_im=currentWorker->toDouble(newt.laguG.im)*order1;
        MandelMath::complex_double_quadratic2(
              &newtres_.first_lagu1_re, &newtres_.first_lagu1_im,
              &newtres_.first_lagu1o_re, &newtres_.first_lagu1o_im,
              a_re, a_im,
              b_re, b_im,
              -1, 0);
        newtres_.first_lagu1_re=r_re-newtres_.first_lagu1_re;
        newtres_.first_lagu1_im=r_im-newtres_.first_lagu1_im;
        newtres_.first_lagu1o_re=r_re-newtres_.first_lagu1o_re;
        newtres_.first_lagu1o_im=r_im-newtres_.first_lagu1o_im;
      };
    lastm=m;
    bool lagu_valid=false;
    bool newt_valid=false;
    if (order1>=0)
    {
      // a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
      // all but last few cycles can be done just in double precision
      //   but the cost of this compared to evaluation of f,f',f'' is negligible
      newt.laguX.assign(&newt.laguG2);
      newt.laguX.lshift(-period); //G^2/n
      newt.laguX.rsub(&newt.laguH); //H-G^2/n
      newt.tmp2.zero(m);
      newt.tmp2.recip();
      newt.tmp2.add_double(-order1); //1/m-1/n
      currentWorker->mul(newt.laguX.re, newt.tmp2.ptr);
      currentWorker->mul(newt.laguX.im, newt.tmp2.ptr); //(1/m-1/n)*(H-G^2/n)
      newt.laguX.sqrt();
      //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
      if (currentWorker->isle0(newt.laguX.mulreT_tmp(&newt.laguG))) //again isl0 would be nicer
      {
        currentWorker->chs(newt.laguX.re);
        currentWorker->chs(newt.laguX.im);
      };
      newt.laguG.lshift(-period); //G/n
      newt.laguX.add(&newt.laguG);
      //if 1/n~0: a=1/(0 +- sqrt( (1/m)*(H) )), m can still be 1..max
      //   fine if H!=0:       a=1/( sqrt( (1/m)*(H) )), m can still be 1..max
      //   if H==0: 1/G/(1/n + sqrt( (1/300-1/n)*(-1/n) ))=1/G* -i*sqrt(300*n)
      //   if H=G=0: 1/0
      //if G=0: a=1/(+- sqrt( (1/m-1/n)*(H) ))     m=1
      //   fine if H!=0: a=+-(sqrt(n/(n-1))*sqrt(f/-f''))       x^2+9 at 0: f=9 f''=2 -> +-3i
      //   if H=0: a=1/0
      //if H=0: a=1/G*m*(1 - i*sqrt(n/m-1))  m~n -> a=n/G;  m~300 -> a=-i/G*sqrt(n*300)
      //        a=1/G*m*n*(1/n - i*sqrt(1/m/n-1/n^2))
      double X_mag=newt.laguX.getMag_double();
      if (X_mag>=1e-60)
      {
        newt.laguX.recip_prepared();
        lagu_valid=true;
      };
      //else
      //we should move the guess a little and try again, but
      //  we can leave this to the caller
      //return 0;
    };
    if (gz_r_mag!=0)
    {
      //newton near multiroot:
      //f=(x-a)^m   f'=m*(x-a)^(m-1)  f/f'=(x-a)/m
      //Newton corrected for multiroot = f/f'*m
      //1/M=1-f''*f/f'^2   x-a=1/(f'/f-f''/f')
      newt.newtX.assign(&newtres_.fz_r);
      newt.newtX.recip();
      newt.newtX.mul(&newt.f_r); //f/f'
      if (m!=1)
      {
        newt.tmp2.zero(m);
        currentWorker->mul(newt.newtX.re, newt.tmp2.ptr);
        currentWorker->mul(newt.newtX.im, newt.tmp2.ptr);
      };
      newt_valid=true;
    };
    if (newtonCycle==0)
    {
      newtres_.first_guess_newt.assign(r);
      if (newt_valid)
      {
        newtres_.first_guess_newt.sub(&newt.newtX);
      };

      newtres_.first_guess_lagu.assign(r);
      if (lagu_valid)
      {
        newtres_.first_guess_lagu.sub(&newt.laguX);
      };
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
        newt.newtX.assign(&newt.laguX);
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
        newt.newtX.assign(&newt.laguX);
      }
      else
      {//take the smaller of newton and laguerre to a) avoid Newton's jumps when moving off the real axis, b) avoid Laguerre skipping to far root
        double N_mag=newt.newtX.getMag_double();
        double L_mag=newt.laguX.getMag_double();
        if (N_mag*1.05>L_mag) //5% will do no harm, and switch to Lagu can speed up convergence
        {
          newt.newtX.assign(&newt.laguX);
        };
      }
    }
    if ((g_r_mag>bestfm) && (newtonCycle>30))
    {
      newt.newtX.lshift(-2);
    };
    r->sub(&newt.newtX);
  } //for newtonCycle
  return lastm;
}

//result=0 means the period check failed; -1 means the check failed and the root returned is invalid
int MandelEvaluator::periodCheck(int period/*must =eval.lookper_lastGuess*/, const MandelMath::complex *c)
{
  if (period<1)
  {
    dbgPoint();
    return -1;
  };
  int aroundCount; //estimate multiplicity (mult-1)
  if ((currentData.store->lookper_prevGuess_>0) &&
      ((currentData.store->lookper_lastGuess % currentData.store->lookper_prevGuess_)==0))
    aroundCount=currentData.store->lookper_lastGuess / currentData.store->lookper_prevGuess_;
  else
    aroundCount=0;
  //look for root nearest to C - better stability of newton/laguerre
  //MandelMath::complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
  currentData.root.assign(&eval.lookper_nearr);
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
  int newtRes=newton(period, c, &currentData.root, true, 1+aroundCount);
  if (newtRes<=0)
  { //this, of course, means that Newton() should be improved, not that there's a problem with the numbers
    return -1; //e.g. evaluating the initial guess mand.root leads to overflow immediately
  };

  //complex f_r(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  //complex fz_r(currentWorker, &newtres_.fz_r_re_, &newtres_.fz_r_im_, true);
  //double dist_around=5*currentWorker->eps2()/currentWorker->toDouble(fz_r.getMagTmp());
  //need at least safety factor of 14.2
  //correct root maps to .to_stop, then also 2*error can map to .to_stop
  //and we oscillate +-error so that's 4*to_stop, or 16 in dist squared
  double dist_around=16*(newtres_.accy_multiplier*newtres_.accy_tostop)*currentWorker->eps2();
  newt.f_r.assign(&currentData.root);
  double f_e_re=sqrt(dist_around/6), f_e_im=0;
  newtres_.fz_r.zero(1, 0);
  double fz_e_m=0;
  eval.fz_mag1.zero(1);
  aroundCount=0;
  bool someBelow1=false;

  //TODO: try to guess the radius of "around" roots   double dist_around=newt.fz*2/fzz  f(r+x)=r+x=r+f'x+f''x^2/2

  //int firstBelow1=-1;
  int firstBelow1dividing=-1;
  for (int i=0; i<period; i++)
  {
    //if (fz_mag1 && currentWorker->isl0(fz_mag1)) //I think we intentionally skip last fz_mag
    if (//replaced with fz_mag1.zero(1) fz_mag1.asf64 &&  //I think we intentionally skip last fz_mag
        (currentWorker->isle0(eval.fz_mag1.ptr) ||
         (MandelMath::sqr_double(currentWorker->toDouble(eval.fz_mag1.ptr))<=fz_e_m)))
    {
      someBelow1=true;
      //if (firstBelow1<0)
        //firstBelow1=i;
      if (firstBelow1dividing<0)
        if ((period % i)==0)
        {
          //needs more checks than that, e.g. fz_mag^(period/i) <=~ final fz_mag
          //per-actual  per-found  root-found   |   status at short
          //  short       short       short         no long to loop over
          //  short       short        long         does not solve newton
          //  short        long       short         |f-r|<eps
          //  short        long        long         |fz|>1
          //  long        short       short         |fz|>1
          //  long        short        long         does not solve newton
          //  long         long       short         |fz|>1
          //  long         long        long         no short to test
          //if (currentWorker->isle(f_r.getMagTmp(), root.getMagTmp()))
          double dist=newt.f_r.dist2_double(&currentData.root);
          if (dist<dist_around)//3.4e-28) //related to newton's accyBound=3e-28/period^2
            firstBelow1dividing=i; //short long short
          else if (dist>1e-7)
          { }
          else
          {
            nop();
          }
          //does "else" even exist? should be always true? 0.00804 vs 1e-31 nope
        };
    };
    //fm could hardly be >4 since it's tested in Newton (as well as f_z_r, BTW)
    //fz:=2*f*fz
    {
      double fzre=currentWorker->toDouble(newtres_.fz_r.re);
      double fzim=currentWorker->toDouble(newtres_.fz_r.im);
      double re=2*(fzre*f_e_re-fzim*f_e_im);
      double im=2*(fzim*f_e_re+fzre*f_e_im);
      fz_e_m=re*re+im*im;
    }
    newtres_.fz_r.mul(&newt.f_r);
    newtres_.fz_r.lshift(1);
    eval.fz_mag1.assign(newtres_.fz_r.getMag1_tmp());
    if (currentWorker->toDouble(eval.fz_mag1.ptr)>LARGE_FLOAT2)
      return -1; //so is it checked or not
    //f:=f^2+c
    //use EXACTLY same code as the main loop: if temporaries are more precise at one place, it will never work!
    {
      double fre=currentWorker->toDouble(newt.f_r.re);
      double fim=currentWorker->toDouble(newt.f_r.im);
      double re=fre*f_e_re-fim*f_e_im;
      f_e_im=fim*f_e_re+fre*f_e_im;
      f_e_re=re;
    }
    newt.f_r.sqr();
    newt.f_r.add(c);
  };
  //what's going on here -.-  normally we would check |ff|<1 or something
  //maybe that's what I'm doing but with some extra tricks for high precision?
  //reduce fz_re, fz_im to 2nd octant
  if (currentWorker->isle0(newtres_.fz_r.re)) currentWorker->chs(newtres_.fz_r.re);
  if (currentWorker->isle0(newtres_.fz_r.im)) currentWorker->chs(newtres_.fz_r.im);
  bool fz_r_mag_over1;
  if (!currentWorker->isle(newtres_.fz_r.re, newtres_.fz_r.im))
  {
    currentWorker->assign(newt.tmp2.ptr, newtres_.fz_r.re);
    currentWorker->assign(newtres_.fz_r.re, newtres_.fz_r.im);
    currentWorker->assign(newtres_.fz_r.im, newt.tmp2.ptr);
  };
  //im>=re, now check re^2+im^2>1
  if ((currentWorker->toDouble(newtres_.fz_r.re)>=0.71) ||
      (currentWorker->toDouble(newtres_.fz_r.im)>=1.01))
  { //im>=re>=0.71
    fz_r_mag_over1=true;
  }
  else if (currentWorker->toDouble(newtres_.fz_r.im)<0.70)
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
    currentWorker->assign(newt.tmp2.ptr, newtres_.fz_r.re);
    newt.tmp2.sqr();
    currentWorker->assign(newt.tmp1.im, newt.tmp2.ptr); //r^2
    currentWorker->chs(newt.tmp2.ptr);
    newt.tmp2.add_double(1);
    newt.tmp2.sqrt();
    currentWorker->assign(newt.tmp1.re, newt.tmp2.ptr); //x+1
    newt.tmp2.add_double(-1);        //x
    newt.tmp2.sqr();
    currentWorker->sub(newt.tmp2.ptr, newt.tmp1.im);    //x^2-r^2
    currentWorker->recip(newt.tmp1.re);
    currentWorker->mul(newt.tmp2.ptr, newt.tmp1.re);
    newt.tmp2.lshift(-1);           //better x
    currentWorker->add_double(newtres_.fz_r.im, -1);
    fz_r_mag_over1=!currentWorker->isle(newtres_.fz_r.im, newt.tmp2.ptr); //im-1>x
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

int MandelEvaluator::estimateInterior(int period, const MandelMath::complex *c, const MandelMath::complex *root)//, InteriorInfo *interior)
{
  MandelMath::complex &f=newt.f_r;
  MandelMath::complex &fz=interior.fz;
  MandelMath::complex &fc=newt.laguG;
  MandelMath::complex &fzz=newt.fzz_r;
  MandelMath::complex &fzc=newt.laguG2;
  MandelMath::complex &fcc=newt.laguH;
  // Initial values:  f = r;  fc = 0;  fz = 1;  fzz = 0;  fzc = 0;
  f.assign(root); //z   z^2+c  (z^2+c)^2+c
  interior.fz.zero(1, 0);           //1   2z     4z^3+4cz=2*2z(z^2+c)
  fc.zero(0, 0);           //0   1      2z^2+2c+1=2*1*(z^2+c)+1
  newt.fzz_r.zero(0, 0);          //0   2      12z^2+4c=2(4z^2+2z^2+2c)=2*((2z)^2+(z^2+c)*2)
  fzc.zero(0, 0);          //0   0      4z=2*(2z*1+(z^2+c)0)
  fcc.zero(0, 0);
  interior.fz_mag.zero(1);
  /*complex f(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  currentWorker->assign(f.re_s, root->re_s);
  currentWorker->assign(f.im_s, root->im_s);
  complex fz(currentWorker, &interior.fz_re, &interior.fz_im, true);
  currentWorker->zero(fz.re_s, 1);
  currentWorker->zero(fz.im_s, 0);
  complex fc(currentWorker, &newt.laguG_re, &newt.laguG_im, true);
  currentWorker->zero(fc.re_s, 0);
  currentWorker->zero(fc.im_s, 0);
  complex fzz(currentWorker, &newt.fzz_r_re, &newt.fzz_r_im, true);
  currentWorker->zero(fzz.re_s, 0);
  currentWorker->zero(fzz.im_s, 0);
  complex fzc(currentWorker, &newt.laguG2_re, &newt.laguG2_im, true);
  currentWorker->zero(fzc.re_s, 0);
  currentWorker->zero(fzc.im_s, 0);
  complex fcc(currentWorker, &newt.laguH_re, &newt.laguH_im, true);
  currentWorker->zero(fcc.re_s, 0);
  currentWorker->zero(fcc.im_s, 0);
  currentWorker->zero(&interior.fz_mag, 1);
  complex tmp1(currentWorker, &newt.tmp1_re, &newt.tmp1_im, true);
  complex inte(currentWorker, &interior.inte_re, &interior.inte_im, true);*/
  // Start iterating.
  for (int i=0; i<period; i++)
  {
    //f^2+c -> 2f fc+1 -> 2f fcc+2fc fc
    // fcc := 2 * (fc^2 + f * fcc);
    newt.tmp1.assign(&fc);
    newt.tmp1.sqr();
    fcc.mul(&f);
    fcc.add(&newt.tmp1);
    fcc.lshift(1);

    //f^2+c -> 2f fz -> 2fc fz+2f fzc
    // fzc := 2 * (fz * fc + f * fzc);
    newt.tmp1.assign(&fc);
    newt.tmp1.mul(&fz);
    fzc.mul(&f);
    fzc.add(&newt.tmp1);
    fzc.lshift(1);

    //f^2+c -> 2f fz -> 2fz fz+2f fzz
    // fzz := 2 * (fz^2 + f * fzz);
    newt.tmp1.assign(&fz);
    newt.tmp1.sqr();
    fzz.mul(&f);
    fzz.add(&newt.tmp1);
    fzz.lshift(1);

    //f^2+c -> 2f fc+1
    // fc := 2 * f * fc + 1;
    fc.mul(&f);
    fc.lshift(1);
    currentWorker->add_double(fc.re, 1);

    //f^2+c -> 2f fz
    // fz := 2 * f * fz;
    fz.mul(&f);
    fz.lshift(1);
    interior.fz_mag.mul(f.getMag_tmp_());
    interior.fz_mag.lshift(2); //fz_mag*=4*mag(f)

    // f = f^2 + c
    f.sqr();
    f.add(c);

    if ((f.getMag_double()>LARGE_FLOAT2) ||
        //(currentWorker->toDouble(fz.getMagTmp())>LARGE_FLOAT2) ||
        currentWorker->toDouble(interior.fz_mag.ptr)>LARGE_FLOAT2 ||
        (fc.getMag_double()>LARGE_FLOAT2) ||
        (fzz.getMag_double()>LARGE_FLOAT2) ||
        (fzc.getMag_double()>LARGE_FLOAT2))
    {
//          DoGlobalDebug('?? when??');
      interior.inte_abs.zero(-1);
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

  newt.tmp2.assign(interior.fz_mag.ptr);
  newt.tmp2.add_double(-1);
  if (abs(currentWorker->toDouble(newt.tmp2.ptr))<6e-18)
  { //parabolic point
    interior.inte_abs.zero(0);
    interior.inte.zero(0, 0);
    return 0;
  };
  //                    1-|fz|^2          .    1-|fz|^2=(1-|fz|)(1+|fz|)
  // interior=  -----------------------   .
  //            | fzc + fzz fc/(1-fz) |   .
  currentWorker->chs(newt.tmp2.ptr);
  currentWorker->assign(interior.inte.re, newt.tmp2.ptr);
  currentWorker->zero(interior.inte.im, 0); //1-|fz|^2
  newt.tmp1.assign(&fz);
  currentWorker->chs(newt.tmp1.re);
  currentWorker->chs(newt.tmp1.im);
  currentWorker->add_double(newt.tmp1.re, 1); //1-fz
  //skip this step because it's just not right: interior:=newt.tmp2 * tmp1/|tmp1|
  newt.tmp1.recip();
  newt.tmp1.mul(&fc);
  newt.tmp1.mul(&fzz);
  newt.tmp1.add(&fzc);
  if (currentWorker->is0(newt.tmp1.re) && currentWorker->is0(newt.tmp1.im))
  { //probably wrong period, should not happen
    interior.inte_abs.zero(5);
    interior.inte.zero(5, 0);
    return period;
  };
  newt.tmp1.recip();
  interior.inte.mul(&newt.tmp1);
  if (currentWorker->isle0(newt.tmp2.ptr)) //should be isl0 again
  {
    interior.inte_abs.assign(interior.inte.getMag_tmp_());
    interior.inte_abs.sqrt();
    currentWorker->chs(interior.inte_abs.ptr);
  }
  else
  {
    interior.inte_abs.assign(interior.inte.getMag_tmp_());
    interior.inte_abs.sqrt();
  }
  return period;
}

void MandelEvaluator::eval_until_bailout(const MandelMath::complex *c, MandelMath::complex *f, MandelMath::complex *fc_c)
{
  for (int i=0; i<100; i++) //should be enough to reach 10000^2 except around (-2, 0)
  {
    double f_mag=f->getMag_double();
    if (f_mag>1e8)
      return;
    //fc_c:=2*f*fc_c+1
    fc_c->mul(f);
    fc_c->lshift(1);
    currentWorker->add_double(fc_c->re, 1);
    double fc_c_mag=fc_c->getMag_double();
    if (fc_c_mag>LARGE_FLOAT2)
    {
      currentData.store->state=MandelPointStore::State::stBoundary;
      return;
    };
    //f:=f^2+c
    f->sqr();
    f->add(c);
    currentData.store->iter++;
  };
}

void MandelEvaluator::evaluate()
{
  /*MandelMath::complex c(currentWorker, &currentParams.c_re, &currentParams.c_im, true);
  MandelMath::complex f(currentWorker, &currentData.f_re, &currentData.f_im, true);
  //currentData.fz_c_mag
  MandelMath::complex fc_c(currentWorker, &currentData.fc_c_re_, &currentData.fc_c_im_, true);*/
  //not needed for math MandelMath::complex fz_r(currentWorker, &eval.fz_r_re, &eval.fz_r_im, true);
  //currentWorker->zero(&eval.fz_r_re, 0);
  //currentWorker->zero(&eval.fz_r_im, 0);
  //currentData.near0iter
  { //near0f not needed for math, just a store
    //MandelMath::complex near0f(currentWorker, &currentData.near0f_re, &currentData.near0f_im, true);
    eval.near0fmag.assign(currentData.near0f.getMag_tmp_()); //TODO: update on changing near0f
  }

  /*MandelMath::complex lookper_startf(currentWorker, &currentData.lookper_startf_re, &currentData.lookper_startf_im, true);
  MandelMath::complex lookper_nearr(currentWorker, &eval.lookper_nearr_re, &eval.lookper_nearr_im, true);*/

  for (; (currentData.store->iter<currentParams.maxiter_) && (currentData.store->state==MandelPointStore::State::stUnknown) && !wantStop; currentData.store->iter++)
  {
    if (currentData.store->iter%(3*currentData.store->near0iter) ==0)
    {
      int quot=currentData.store->iter/(3*currentData.store->near0iter);
      if ((quot&(quot-1))==0) //also at iter==0  //TODO: maybe better 3*(2^k-1)*near not 3*(2^k)*near
      { // //need k*iter for f' to start at the worst moment to reduce false positives; need k*iter-1 for good near0 -> switch to nearc
        currentData.store->lookper_startiter=currentData.store->iter;
        currentData.lookper_startf.assign(&currentData.f);
        //MandelMath::complex lookper_bestf(currentWorker, &eval.lookper_startf_re, &eval.lookper_startf_im, true);
        //currentWorker->zero(&eval.lookper_bestf_re, 0);
        //currentWorker->zero(&eval.lookper_bestf_im, 0);
        eval.lookper_nearr.assign(&currentData.f);
        if (currentData.store->iter<=1)
          currentData.lookper_nearr_dist.assign(currentData.f.getMag_tmp_());
        else
          currentData.store->lookper_nearr_dist_touched=false;//currentWorker->assign(&currentData.lookper_nearr_dist, f.dist2_tmp(&c));
        //currentWorker->zero(&eval.lookper_dist2, 1e10); //4.0 should be enough
        //mands.period stays there
        //currentData.lookper_prevGuess=0; //TODO: used for anything?
        currentData.store->lookper_lastGuess=0;
        currentData.lookper_totalFzmag.zero(1.0);
      };
    }
    MandelMath::number_pointer_c f_mag=currentData.f.getMag_tmp_();
    if (currentWorker->toDouble(f_mag)>4)
    {
      currentData.store->state=MandelPointStore::State::stOutside;
      //theory says the relative error in estimate is less than 3/bailout for large bailout
      //so lets move out a bit
      eval_until_bailout(&currentParams.c, &currentData.f, &currentData.fc_c); //may switch state to stBoundary
      if (currentData.store->state!=MandelPointStore::State::stOutside)
      {
        //currentWorker->zero(&currentData.exterior_avoids, 0);
        //currentWorker->zero(&currentData.exterior_hits, 0);
        currentData.store->exterior_avoids=0;
        currentData.store->exterior_hits=0;
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
        double fm=currentData.f.getMag_double();
        double fcm=currentData.fc_c.getMag_double();
        double x=log(fm);
        currentData.store->exterior_hits=x*sqrt(fm/fcm);
        currentData.store->exterior_avoids=currentData.store->exterior_hits*0.25;
        if (currentData.store->iter>71)
        { }
        else
        {
          x=ldexp(x, -1-currentData.store->iter);
          if (currentData.store->iter>26)
          {
            currentData.store->exterior_hits+=x*x/6*currentData.store->exterior_hits;
            currentData.store->exterior_avoids+=x*(x*2/3-1)*currentData.store->exterior_avoids;
          }
          else
          {
            double ex=exp(x);
            currentData.store->exterior_hits*=(ex-1/ex)/x/2;
            currentData.store->exterior_avoids*=(1-1/(ex*ex))/x/2;
          }
        }
      }
      //already there currentWorker->assign(currentData.fc_c_re, &fc_c_re);
      //currentWorker->assign(currentData.fc_c_im, &fc_c_im);
      currentData.store->period=currentData.store->lookper_lastGuess; //preliminary
      if (currentData.store->period<1)
        currentData.store->period=1;
      break;
    };
    double fc_c_mag=currentData.fc_c.getMag_double();
    if (fc_c_mag>1e57)
    {
      currentData.store->state=MandelPointStore::State::stBoundary;
      currentData.store->exterior_avoids=0;
      currentData.store->exterior_hits=0;
      break;
    };
    double fz_c_mag=currentWorker->toDouble(currentData.fz_c_mag.ptr);
    if (fz_c_mag>1e60)
    {
      currentData.store->state=MandelPointStore::State::stDiverge;
      currentData.store->exterior_avoids=0;
      currentData.store->exterior_hits=0;
      break;
    };
    //TODO: similar to eval_until_bailout
    //fc_c:=2*f*fc_c+1
    currentData.fc_c.mul(&currentData.f);
    currentData.fc_c.lshift(1);
    currentWorker->add_double(currentData.fc_c.re, 1);
    /* TODO: copy test here from above?
    fc_c_mag=currentWorker->toDouble(fc_c.getMagTmp());
    if (fc_c_mag>LARGE_FLOAT2)
    {
      currentData.state=MandelPoint::State::stBoundary;
      break;
    };*/
    f_mag=currentData.f.getMag_tmp_();
    //f'=2*f'*f, f'_mag=4*f'_mag*f_mag
    currentData.fz_c_mag.mul(f_mag); //TODO: can use f_mag from above? would need storage, tmp can't survive this long
    currentData.lookper_totalFzmag.mul(f_mag);
    currentData.lookper_totalFzmag.lshift(2);
    //f:=f^2+c
    currentData.f.sqr();
    currentData.f.add(&currentParams.c);
    //currentData.iter++;
    f_mag=currentData.f.getMag_tmp_();

    if (!currentWorker->isle(eval.near0fmag.ptr, f_mag)) //f_mag<near0fmag
    {
      currentData.store->near0iter=currentData.store->iter+2;
      currentData.near0f.assign(&currentData.f);
      eval.near0fmag.assign(f_mag);
    };

    const MandelMath::number_pointer_c lpdiff=currentData.lookper_startf.dist2_tmp_(&currentData.f);
    switch (currentWorker->compare(lpdiff, currentData.lookper_nearr_dist.ptr)) //|f-r|<best
    {
      case -1:
      {
        eval.lookper_nearr.assign(&currentData.f);
        currentData.lookper_nearr_dist.assign(lpdiff);
        currentData.store->lookper_nearr_dist_touched=false;
        currentData.store->lookper_prevGuess_=currentData.store->lookper_lastGuess;
        currentData.store->lookper_lastGuess=(currentData.store->iter+1-currentData.store->lookper_startiter);
      } break;
      case 0:
      {
        if (!currentData.store->lookper_nearr_dist_touched)
        {
          eval.lookper_nearr.assign(&currentData.f);
          currentData.lookper_nearr_dist.assign(lpdiff);
          currentData.store->lookper_nearr_dist_touched=true;
          currentData.store->lookper_prevGuess_=currentData.store->lookper_lastGuess;
          currentData.store->lookper_lastGuess=(currentData.store->iter+1-currentData.store->lookper_startiter);
        };
      } break;
    };

    if ((currentData.store->lookper_lastGuess>0) &&
        (currentData.store->lookper_lastGuess==(currentData.store->iter+1-currentData.store->lookper_startiter)) && //just found new guess
#if USE_GCD_FOR_CHECKPERIOD
#else
        ((currentData.store->near0iter % currentData.store->lookper_lastGuess)==0) && //  period divides nearest, that's a fact
#endif
        ((currentData.store->iter>=3*currentData.store->near0iter)))  //speedup - don't check period too eagerly
    {
#if USE_GCD_FOR_CHECKPERIOD
      int testperiod=MandelMath::gcd(currentData.near0iter, currentData.lookper_lastGuess);//currentData.lookper_lastGuess
#else
      int testperiod=currentData.store->lookper_lastGuess;
#endif
      int foundperiod=-1;
      if (currentData.f.isequal(&currentData.lookper_startf))
      { //exact match - misiurewicz or converged after too many steps
        foundperiod=currentData.store->lookper_lastGuess;
        currentData.root.assign(&currentData.f);
        //TODO: still needs period cleanup... I think. Near 0+0I
      }
      else if (currentWorker->toDouble(currentData.lookper_totalFzmag.ptr)<MAGIC_MIN_SHRINK)
      {
        //if (periodCheckHistory<>'') then
        //  periodCheckHistory:=periodCheckHistory+'('+IntToStr(mande.iter)+':'+IntToStr(mande.lookper_lastGuess)+') ';
        //periodEntered:=getTime64();
        //if (currentWorker->toDouble(c.re_s)==-0.015625 && currentWorker->toDouble(c.im_s)==0.640625)
          //dbgPoint();
        foundperiod=periodCheck(testperiod, &currentParams.c); //updates iter, f, f_c, root
        //profiler.timeInPeriod:=profiler.timeInPeriod+(getTime64()-periodEntered);
        if ((foundperiod>0) && (foundperiod<testperiod))
        {
          //complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
          foundperiod=estimateInterior(foundperiod, &currentParams.c, &currentData.root); //-4.7e-22
            //foundperiod=-1; //the cycle can be exact but |f_z|>1 due to mistaken period or misplaced (rounding err) root
        }
      };
      if (foundperiod>0)
      {
        currentData.store->state=MandelPointStore::State::stPeriod2;
        currentData.store->period=foundperiod;
        currentData.store->newton_iter=newtres_.cyclesNeeded;
        //complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
        currentData.store->period=estimateInterior(foundperiod, &currentParams.c, &currentData.root);
        if (currentWorker->isl0(interior.inte_abs.ptr))
          currentData.store->state=MandelPointStore::State::stMisiur;
        else
        {
          currentData.store->interior=currentWorker->toDouble(interior.inte_abs.ptr);
          currentData.fz_r.assign(&interior.fz); //d/dz F_c(r)
          if (foundperiod!=currentData.store->period)
            currentData.store->state=MandelPointStore::State::stPeriod3;
        }
        //currentWorker->assign(&currentData.fc_c_re, &eval.fz_r_re);
        //currentWorker->assign(&currentData.fc_c_im, &eval.fz_r_im);
        break;
      };
      if (currentParams.breakOnNewNearest)
      {
        currentData.store->iter++;
        break;
      }
    };

  }
  //data.state=MandelPoint::State::stMaxIter;
  if (!currentData.store->has_fc_r && currentParams.want_fc_r &&
      ((currentData.store->state==MandelPointStore::State::stPeriod2) ||
       (currentData.store->state==MandelPointStore::State::stPeriod3)))
  {
    //MandelMath::complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
    this->bulb.bulbe.eval2(currentData.store->period, &currentParams.c, &currentData.root);
    currentData.fc_c.assign(&bulb.bulbe.f_c);
    currentData.store->has_fc_r=true;
  };

}


MandelEvaluator::ComputeParams::ComputeParams(MandelMath::worker_multi::Allocator *allocator):
  c(allocator), epoch(-1), pixelIndex(-1), maxiter_(1), breakOnNewNearest(false), want_fc_r(false)
{
}

MandelEvaluator::NewtRes::NewtRes(MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN), cyclesNeeded(-1),
  fz_r(&self_allocator), first_guess_lagu(&self_allocator), first_guess_newt(&self_allocator)
{
  assert(self_allocator.checkFill());
}

MandelEvaluator::Eval::Eval(MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN),
  fz_r(&self_allocator), fz_mag1(&self_allocator),
  near0fmag(&self_allocator), lookper_nearr(&self_allocator)
{
  assert(self_allocator.checkFill());
}

MandelEvaluator::Newt::Newt(MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN),
  bestr(&self_allocator), f_r(&self_allocator), fzz_r(&self_allocator), tmp1(&self_allocator),
  laguH(&self_allocator), laguG(&self_allocator), laguG2(&self_allocator),
  laguX(&self_allocator), newtX(&self_allocator), fzzf(&self_allocator), tmp2(&self_allocator)
{
  assert(self_allocator.checkFill());
}

MandelEvaluator::InteriorInfo::InteriorInfo(MandelMath::worker_multi::Allocator *allocator, MandelMath::upgrademe *):
  self_allocator(allocator, LEN),
  inte(&self_allocator), inte_abs(&self_allocator), fz(&self_allocator), fz_mag(&self_allocator)
{
  assert(self_allocator.checkFill());
}

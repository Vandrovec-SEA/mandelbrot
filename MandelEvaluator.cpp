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

LaguerrePoint::LaguerrePoint()
{
  //reset();
  //should be overwritten before read:
  state=State::stUnknown;
  iter=0;
  firstM=0;
}

void LaguerrePoint::assign(MandelMath::number_worker *worker, const LaguerrePoint &src)
{
  worker->assign(&f_re, &src.f_re);
  worker->assign(&f_im, &src.f_im);
  worker->assign(&fz_r_re, &src.fz_r_re);
  worker->assign(&fz_r_im, &src.fz_r_im);
  state=src.state;
  iter=src.iter;
  firstM=src.firstM;
}

void LaguerrePoint::init(MandelMath::number_worker *worker)
{
  /*worker->init(&f_re, 0.0);
  worker->init(&f_im, 0.0);
  worker->init(&fz_r_re, 0.0);
  worker->init(&fz_r_im, 0.0);*/
  promote(MandelMath::number_worker::Type::typeEmpty, worker->ntype());
  worker->zero(&f_re, 0.0);
  worker->zero(&f_im, 0.0);
  worker->zero(&fz_r_re, 0.0);
  worker->zero(&fz_r_im, 0.0);

  //should be overwritten before read:
  state=State::stUnknown;
  iter=0;
  firstM=0;
}

void LaguerrePoint::zero(MandelMath::number_worker *worker, const MandelMath::number_store *c_re, const MandelMath::number_store *c_im)
{
  //worker->zero(&f_re, 0);
  //worker->zero(&f_im, 0);
  worker->assign(&f_re, c_re);
  worker->assign(&f_im, c_im);
  worker->zero(&fz_r_re, 1);
  worker->zero(&fz_r_im, 0);

  state=State::stUnknown;
  iter=0;
  firstM=0;
}

void LaguerrePoint::cleanup(MandelMath::number_worker *worker)
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
  /*if (place.dd!=nullptr)
  {
    delete[] place.dd;
    place.dd=nullptr;
  }*/
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
}



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
  worker->assign(&fc_c_re_, &src.fc_c_re_);
  worker->assign(&fc_c_im_, &src.fc_c_im_);
  has_fc_r=src.has_fc_r;
  worker->assign(&fz_r_re, &src.fz_r_re);
  worker->assign(&fz_r_im, &src.fz_r_im);
  worker->assign(&fz_c_mag_, &src.fz_c_mag_);
  lookper_startiter=src.lookper_startiter;
  lookper_prevGuess_=src.lookper_prevGuess_;
  lookper_lastGuess=src.lookper_lastGuess;
  worker->assign(&lookper_startf_re, &src.lookper_startf_re);
  worker->assign(&lookper_startf_im, &src.lookper_startf_im);
  worker->assign(&lookper_nearr_dist_, &src.lookper_nearr_dist_);
  lookper_nearr_dist_touched=src.lookper_nearr_dist_touched;
  worker->assign(&lookper_totalFzmag, &src.lookper_totalFzmag);
  near0iter=src.near0iter;
  newton_iter=src.newton_iter;
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
  /*worker->init(&f_re, 0.0);
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
  worker->init(&lookper_nearr_dist_);
  worker->init(&lookper_totalFzmag);*/
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
}

void MandelPoint::zero(MandelMath::number_worker *worker, const MandelMath::number_store *c_re, const MandelMath::number_store *c_im)
{
  //worker->zero(&f_re, 0);
  //worker->zero(&f_im, 0);
  worker->assign(&f_re, c_re);
  worker->assign(&f_im, c_im);
  worker->zero(&fc_c_re_, 0);
  worker->zero(&fc_c_im_, 0);
  worker->zero(&fz_r_re, 1);
  worker->zero(&fz_r_im, 0);
  worker->zero(&fz_c_mag_, 1);
  lookper_prevGuess_=0;
  lookper_lastGuess=0;
  //lookper resets at first iter
  near0iter=1;
  worker->assign(&near0f_re, c_re);
  worker->assign(&near0f_im, c_im);
  period=0;
  worker->zero(&root_re, 0);
  worker->zero(&root_im, 0);

  state=State::stUnknown;
  iter=0;
  newton_iter=0;
  exterior_avoids=-1;
  exterior_hits=-1;
  interior=-1;
  has_fc_r=false;
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
  /*if (place.dd!=nullptr)
  {
    delete[] place.dd;
    place.dd=nullptr;
  }*/

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
}

ShareableViewInfo::ShareableViewInfo(): worker(nullptr), re_(), im(), root_re(), root_im(), scale(1), period(0)
{

}

ShareableViewInfo::ShareableViewInfo(ShareableViewInfo &src): QObject()
{
  worker=src.worker;
  scale=src.scale;
  period=src.period;
  if (worker!=nullptr)
  {
    worker->swap(&re_, &src.re_);
    worker->swap(&im, &src.im);
    worker->swap(&root_re, &src.root_re);
    worker->swap(&root_im, &src.root_im);
  };
}

ShareableViewInfo::ShareableViewInfo(const ShareableViewInfo &src): ShareableViewInfo((ShareableViewInfo &)src)
{
  //why do you need this?
}

ShareableViewInfo::ShareableViewInfo(ShareableViewInfo &&src): QObject()
{
  worker=src.worker;
  scale=src.scale;
  period=src.period;
  if (worker!=nullptr)
  {
    worker->swap(&re_, &src.re_);
    worker->swap(&im, &src.im);
    worker->swap(&root_re, &src.root_re);
    worker->swap(&root_im, &src.root_im);
  };
}

ShareableViewInfo &ShareableViewInfo::operator=(ShareableViewInfo &src)
{
  worker=src.worker;
  scale=src.scale;
  period=src.period;
  if (worker!=nullptr)
  {
    worker->swap(&re_, &src.re_);
    worker->swap(&im, &src.im);
    worker->swap(&root_re, &src.root_re);
    worker->swap(&root_im, &src.root_im);
  };
  return *this;
}

ShareableViewInfo &ShareableViewInfo::operator=(ShareableViewInfo &&src)
{
  worker=src.worker;
  scale=src.scale;
  period=src.period;
  if (worker!=nullptr)
  {
    worker->swap(&re_, &src.re_);
    worker->swap(&im, &src.im);
    worker->swap(&root_re, &src.root_re);
    worker->swap(&root_im, &src.root_im);
  };
  return *this;
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

#if NUMBER_DOUBLE_EXISTS
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
#endif //NUMBER_DOUBLE_EXISTS

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

void MandelEvaluator::switchType(MandelMath::number_worker *worker)
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
    /*currentWorker->cleanup(&bulb.g_re);
    currentWorker->cleanup(&bulb.g_im);
    currentWorker->cleanup(&bulb.g_c_re);
    currentWorker->cleanup(&bulb.g_c_im);
    currentWorker->cleanup(&bulb.g_cc_re);
    currentWorker->cleanup(&bulb.g_cc_im);
    currentWorker->cleanup(&bulb.g_c2_re);
    currentWorker->cleanup(&bulb.g_c2_im);
    currentWorker->cleanup(&bulb.f_re);
    currentWorker->cleanup(&bulb.f_im);
    currentWorker->cleanup(&bulb.f_z_re);
    currentWorker->cleanup(&bulb.f_z_im);
    currentWorker->cleanup(&bulb.f_c_re);
    currentWorker->cleanup(&bulb.f_c_im);
    currentWorker->cleanup(&bulb.f_zz_re);
    currentWorker->cleanup(&bulb.f_zz_im);
    currentWorker->cleanup(&bulb.f_zc_re);
    currentWorker->cleanup(&bulb.f_zc_im);
    currentWorker->cleanup(&bulb.f_cc_re);
    currentWorker->cleanup(&bulb.f_cc_im);*/

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
    /*worker->init_(&bulb.g_re, placement); placement+=step;
    worker->init_(&bulb.g_im, placement); placement+=step;
    worker->init_(&bulb.g_c_re, placement); placement+=step;
    worker->init_(&bulb.g_c_im, placement); placement+=step;
    worker->init_(&bulb.g_cc_re, placement); placement+=step;
    worker->init_(&bulb.g_cc_im, placement); placement+=step;
    worker->init_(&bulb.g_c2_re, placement); placement+=step;
    worker->init_(&bulb.g_c2_im, placement); placement+=step;
    worker->init_(&bulb.f_re, placement); placement+=step;
    worker->init_(&bulb.f_im, placement); placement+=step;
    worker->init_(&bulb.f_z_re, placement); placement+=step;
    worker->init_(&bulb.f_z_im, placement); placement+=step;
    worker->init_(&bulb.f_c_re, placement); placement+=step;
    worker->init_(&bulb.f_c_im, placement); placement+=step;
    worker->init_(&bulb.f_zz_re, placement); placement+=step;
    worker->init_(&bulb.f_zz_im, placement); placement+=step;
    worker->init_(&bulb.f_zc_re, placement); placement+=step;
    worker->init_(&bulb.f_zc_im, placement); placement+=step;
    worker->init_(&bulb.f_cc_re, placement); placement+=step;
    worker->init_(&bulb.f_cc_im, placement); placement+=step;*/

    assert((size_t)(placement-(uint8_t *)place.dd)==step*Place::LEN);
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
      ((quick_route==0) && (currentParams.maxiter_-currentData.iter<=1000)))
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
  currentData.period=period;
  currentWorker->assign(&currentParams.c_re, c->re_s);
  currentWorker->assign(&currentParams.c_im, c->im_s);
  QMetaObject::invokeMethod(this, &MandelEvaluator::doNewton, Qt::ConnectionType::QueuedConnection);
}

void MandelEvaluator::doNewton()
{
  MandelMath::complex tmpc(currentWorker, &currentParams.c_re, &currentParams.c_im, true);
  MandelMath::complex root(currentWorker, &currentData.f_re, &currentData.f_im, true);
  int result=newton(currentData.period, &tmpc, &root, true, 12);
  emit doneNewton(this, result);
}

LaguerreStep::LaguerreStep(): currentWorker(nullptr)
{

}

void LaguerreStep::switchType(MandelMath::number_worker *worker)
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
}

bool LaguerreStep::eval(int lg2_degree, const MandelMath::complex *f, const MandelMath::complex *f_z, const MandelMath::complex *f_zz)
{
  MandelMath::complex tmp1(currentWorker, &tmp1_re, &tmp1_im, true);
  MandelMath::complex laguG(currentWorker, &laguG_re, &laguG_im, true);
  MandelMath::complex laguG2(currentWorker, &laguG2_re, &laguG2_im, true);
  MandelMath::complex laguH(currentWorker, &laguH_re, &laguH_im, true);
  MandelMath::complex laguX(currentWorker, &laguX_re, &laguX_im, true);
  MandelMath::complex newtX(currentWorker, &step_re, &step_im, true);
  MandelMath::complex fzzf(currentWorker, &fzzf_re, &fzzf_im, true);
  if (currentWorker->is0(f->re_s) && currentWorker->is0(f->im_s))
  {
    currentWorker->zero(&step_re);
    currentWorker->zero(&step_im);
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
  currentWorker->assign(tmp1.re_s, f->re_s);
  currentWorker->assign(tmp1.im_s, f->im_s);
  tmp1.recip();    //1/f
  currentWorker->assign(laguG.re_s, f_z->re_s);
  currentWorker->assign(laguG.im_s, f_z->im_s);
  laguG.mul(&tmp1); //laguG = f'/f
  currentWorker->assign(fzzf.re_s, f_zz->re_s);
  currentWorker->assign(fzzf.im_s, f_zz->im_s);
  fzzf.mul(&tmp1); //f''/f


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
      double h_re=currentWorker->toDouble(laguH.re_s);
      double h_im=currentWorker->toDouble(laguH.im_s);
      double h_mag=h_re*h_re+h_im*h_im;
      double g2_re=currentWorker->toDouble(laguG2.re_s);
      double g2_im=currentWorker->toDouble(laguG2.im_s);
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
    double G2_mag=currentWorker->toDouble(laguG2.getMagTmp());
    if (G2_mag<0.01)
    { //G2_mag bad
      m=1;
      dbg.mu_re=1;
      dbg.mu_im=0;
    }
    else
    {
      currentWorker->assign(laguX.re_s, laguG2.re_s);
      currentWorker->assign(laguX.im_s, laguG2.im_s);
      currentWorker->chs(laguX.im_s);
      laguX.mul(&laguH);
      double a_re=currentWorker->toDouble(laguX.re_s)/G2_mag*(1-order1)-order1;
      double a_im=currentWorker->toDouble(laguX.im_s)/G2_mag*(1-order1);
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
    currentWorker->assign(laguX.re_s, laguG2.re_s);
    currentWorker->assign(laguX.im_s, laguG2.im_s);
    currentWorker->lshift(laguX.re_s, -lg2_degree);
    currentWorker->lshift(laguX.im_s, -lg2_degree); //G^2/n
    currentWorker->rsub(laguX.re_s, laguH.re_s);
    currentWorker->rsub(laguX.im_s, laguH.im_s); //H-G^2/n
    currentWorker->zero(&tmp2, m);
    currentWorker->recip(&tmp2);
    currentWorker->add_double(&tmp2, -order1); //1/m-1/n
    currentWorker->mul(laguX.re_s, &tmp2);
    currentWorker->mul(laguX.im_s, &tmp2); //(1/m-1/n)*(H-G^2/n)
    laguX.sqrt();
    //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
    if (currentWorker->isle0(laguX.mulreT(&laguG))) //again isl0 would be nicer
    {
      currentWorker->chs(laguX.re_s);
      currentWorker->chs(laguX.im_s);
    };
    currentWorker->lshift(laguG.re_s, -lg2_degree);
    currentWorker->lshift(laguG.im_s, -lg2_degree); //G/n
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
  if (!currentWorker->is0(f_z->re_s) || !currentWorker->is0(f_z->im_s)) //gz_r_mag!=0)
  {
    //newton near multiroot:
    //f=(x-a)^m   f'=m*(x-a)^(m-1)  f/f'=(x-a)/m
    //Newton corrected for multiroot = f/f'*m
    //1/M=1-f''*f/f'^2   x-a=1/(f'/f-f''/f')
    currentWorker->assign(newtX.re_s, f_z->re_s);
    currentWorker->assign(newtX.im_s, f_z->im_s);
    newtX.recip();
    newtX.mul(f); //f/f'
    if (m!=1)
    {
      currentWorker->zero(&tmp2, m);
      currentWorker->mul(newtX.re_s, &tmp2);
      currentWorker->mul(newtX.im_s, &tmp2);
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
    currentWorker->assign(newtX.re_s, laguX.re_s);
    currentWorker->assign(newtX.im_s, laguX.im_s);
  }
  else if (!lagu_valid)
  {
    //keep newtX
  }
  else
  {
    if (m>1)//(fastHoming && (newtonCycle<2) && (m>1))
    {
      currentWorker->assign(newtX.re_s, laguX.re_s);
      currentWorker->assign(newtX.im_s, laguX.im_s);
    }
    else
    {//take the smaller of newton and laguerre to a) avoid Newton's jumps when moving off the real axis, b) avoid Laguerre skipping to far root
      double N_mag=currentWorker->toDouble(newtX.getMagTmp());
      double L_mag=currentWorker->toDouble(laguX.getMagTmp());
      if (N_mag*1.05>L_mag) //5% will do no harm, and switch to Lagu can speed up convergence
      {
        currentWorker->assign(newtX.re_s, laguX.re_s);
        currentWorker->assign(newtX.im_s, laguX.im_s);
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

MandelLoopEvaluator::MandelLoopEvaluator(): currentWorker(nullptr)
{

}

void MandelLoopEvaluator::switchType(MandelMath::number_worker *worker)
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
}

bool MandelLoopEvaluator::evalg(int period, const MandelMath::complex *c)
{
  currentWorker->assign(&f_re, c->re_s);
  currentWorker->assign(&f_im, c->im_s);
  currentWorker->zero(&f_c_re, 1);
  currentWorker->zero(&f_c_im, 0);
  currentWorker->zero(&f_cc_re, 0);
  currentWorker->zero(&f_cc_im, 0);
  MandelMath::complex f(currentWorker, &f_re, &f_im, true);
  MandelMath::complex f_c(currentWorker, &f_c_re, &f_c_im, true);
  MandelMath::complex f_cc(currentWorker, &f_cc_re, &f_cc_im, true);
  MandelMath::complex s2(currentWorker, &s2_re, &s2_im, true);
  for (int i=0; i<period; i++)
  {
    //g_cc=2*(g_cc*g + g_c*g_c)
    f_cc.mul(&f);
    currentWorker->assign(&s2_re, &f_c_re);
    currentWorker->assign(&s2_im, &f_c_im);
    s2.sqr();
    f_cc.add(&s2);
    currentWorker->lshift(f_cc.re_s, 1);
    currentWorker->lshift(f_cc.im_s, 1);
    //g_c=2*g_c*g+1
    f_c.mul(&f);
    currentWorker->lshift(f_c.re_s, 1);
    currentWorker->lshift(f_c.im_s, 1);
    currentWorker->add_double(f_c.re_s, 1);
    //g=g^2+c
    f.sqr();
    f.add(c);
    double f_mag=currentWorker->toDouble(f.getMagTmp());
    double allmag=f_mag+
                  currentWorker->toDouble(f_c.getMagTmp())+
                  currentWorker->toDouble(f_cc.getMagTmp());
    if (allmag>1e60)
      return false;
  }
  return true;
}

bool MandelLoopEvaluator::eval2(int period, const MandelMath::complex *c, const MandelMath::complex *z)
{
  currentWorker->assign(&f_re, z->re_s);
  currentWorker->assign(&f_im, z->im_s);
  currentWorker->zero(&f_z_re, 1);
  currentWorker->zero(&f_z_im);
  currentWorker->zero(&f_c_re);
  currentWorker->zero(&f_c_im);
  currentWorker->zero(&f_zz_re);
  currentWorker->zero(&f_zz_im);
  currentWorker->zero(&f_zc_re);
  currentWorker->zero(&f_zc_im);
  currentWorker->zero(&f_cc_re);
  currentWorker->zero(&f_cc_im);
  MandelMath::complex f(currentWorker, &f_re, &f_im, true);
  MandelMath::complex f_z(currentWorker, &f_z_re, &f_z_im, true);
  MandelMath::complex f_c(currentWorker, &f_c_re, &f_c_im, true);
  MandelMath::complex f_zz(currentWorker, &f_zz_re, &f_zz_im, true);
  MandelMath::complex f_zc(currentWorker, &f_zc_re, &f_zc_im, true);
  MandelMath::complex f_cc(currentWorker, &f_cc_re, &f_cc_im, true);
  //MandelMath::complex s1(currentWorker, &s1_re, &s1_im, true);
  MandelMath::complex s2(currentWorker, &s2_re, &s2_im, true);
  for (int i=0; i<period; i++)
  {
    double m_sum=currentWorker->toDouble(f.getMagTmp())+
                 currentWorker->toDouble(f_z.getMagTmp())+
                 currentWorker->toDouble(f_zz.getMagTmp())+
                 currentWorker->toDouble(f_c.getMagTmp())+
                 currentWorker->toDouble(f_zc.getMagTmp())+
                 currentWorker->toDouble(f_cc.getMagTmp());
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_cc=2 * (f_c^2 + f * f_cc)
    f_cc.mul(&f);
    currentWorker->assign(s2.re_s, f_c.re_s);
    currentWorker->assign(s2.im_s, f_c.im_s);
    s2.sqr();
    f_cc.add(&s2);
    currentWorker->lshift(f_cc.re_s, 1);
    currentWorker->lshift(f_cc.im_s, 1);
    // f_zc = 2 * (f_z * f_c + f * f_zc);
    f_zc.mul(&f);
    currentWorker->assign(s2.re_s, f_c.re_s);
    currentWorker->assign(s2.im_s, f_c.im_s);
    s2.mul(&f_z);
    f_zc.add(&s2);
    currentWorker->lshift(f_zc.re_s, 1);
    currentWorker->lshift(f_zc.im_s, 1);
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    currentWorker->assign(s2.re_s, f_z.re_s);
    currentWorker->assign(s2.im_s, f_z.im_s);
    s2.sqr();
    f_zz.add(&s2);
    currentWorker->lshift(f_zz.re_s, 1);
    currentWorker->lshift(f_zz.im_s, 1);
    // f_c = 2 * f * f_c + 1;
    f_c.mul(&f);
    currentWorker->lshift(f_c.re_s, 1);
    currentWorker->lshift(f_c.im_s, 1);
    currentWorker->add_double(f_c.re_s, 1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    currentWorker->lshift(f_z.re_s, 1);
    currentWorker->lshift(f_z.im_s, 1);
    //f:=f^2+c
    f.sqr();
    f.add(c);
  }
  return true;
}

bool MandelLoopEvaluator::eval_zz(int period, const MandelMath::complex *c, const MandelMath::complex *z)
{
  currentWorker->assign(&f_re, z->re_s);
  currentWorker->assign(&f_im, z->im_s);
  currentWorker->zero(&f_z_re, 1);
  currentWorker->zero(&f_z_im);
  currentWorker->zero(&f_zz_re);
  currentWorker->zero(&f_zz_im);
  MandelMath::complex f(currentWorker, &f_re, &f_im, true);
  MandelMath::complex f_z(currentWorker, &f_z_re, &f_z_im, true);
  MandelMath::complex f_zz(currentWorker, &f_zz_re, &f_zz_im, true);
  //MandelMath::complex s1(currentWorker, &s1_re, &s1_im, true);
  MandelMath::complex s2(currentWorker, &s2_re, &s2_im, true);
  for (int i=0; i<period; i++)
  {
    double m_sum=currentWorker->toDouble(f.getMagTmp())+
                 currentWorker->toDouble(f_z.getMagTmp())+
                 currentWorker->toDouble(f_zz.getMagTmp());
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    currentWorker->assign(s2.re_s, f_z.re_s);
    currentWorker->assign(s2.im_s, f_z.im_s);
    s2.sqr();
    f_zz.add(&s2);
    currentWorker->lshift(f_zz.re_s, 1);
    currentWorker->lshift(f_zz.im_s, 1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    currentWorker->lshift(f_z.re_s, 1);
    currentWorker->lshift(f_z.im_s, 1);
    //bigger than 3 is common
    //f:=f^2+c
    f.sqr();
    f.add(c);
  }
  return true;
}

bool MandelLoopEvaluator::eval_multi(int period, const MandelMath::complex *c, const MandelMath::complex *z, const MandelMath::complex *f_z_target)
{
  currentWorker->assign(&f_re, z->re_s);
  currentWorker->assign(&f_im, z->im_s);
  currentWorker->zero(&f_z_re, 1);
  currentWorker->zero(&f_z_im);
  currentWorker->zero(&f_zz_re);
  currentWorker->zero(&f_zz_im);
  MandelMath::complex f(currentWorker, &f_re, &f_im, true);
  MandelMath::complex f_z(currentWorker, &f_z_re, &f_z_im, true);
  MandelMath::complex f_zz(currentWorker, &f_zz_re, &f_zz_im, true);
  MandelMath::complex s1(currentWorker, &s1_re, &s1_im, true);
  MandelMath::complex s2(currentWorker, &s2_re, &s2_im, true);
  MandelMath::complex sumA(currentWorker, &sumA_re, &sumA_im, true);
  int near1=0;
  int sumnear1=0;
  //bool dangerzone=false;
  int reducedbymag=period;
  double closest_accepted=10000, closest_rejected=10000;
  int sumA_cnt=0;
  for (int i=0; i<period; i++)
  {
    double m_sum=currentWorker->toDouble(f.getMagTmp())+
                 currentWorker->toDouble(f_z.getMagTmp())+
                 currentWorker->toDouble(f_zz.getMagTmp());
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    currentWorker->assign(s2.re_s, f_z.re_s);
    currentWorker->assign(s2.im_s, f_z.im_s);
    s2.sqr();
    f_zz.add(&s2);
    currentWorker->lshift(f_zz.re_s, 1);
    currentWorker->lshift(f_zz.im_s, 1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    currentWorker->lshift(f_z.re_s, 1);
    currentWorker->lshift(f_z.im_s, 1);
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
    if (currentWorker->toDouble(f_z.dist2_tmp(f_z_target))*period>0.5) //assuming all points are spaced regularly at 1/period around the circle,
    {                                                                  //skip those close to target to avoid div by 0
      currentWorker->assign(s1.re_s, z->re_s);
      currentWorker->assign(s1.im_s, z->im_s);
      s1.mul(&f_z);
      currentWorker->assign(s2.re_s, f.re_s);
      currentWorker->assign(s2.im_s, f.im_s);
      s2.mul(f_z_target);
      currentWorker->rsub(s2.re_s, s1.re_s);
      currentWorker->rsub(s2.im_s, s1.im_s);
      currentWorker->assign(s1.re_s, f_z.re_s);
      currentWorker->assign(s1.im_s, f_z.im_s);
      currentWorker->sub(s1.re_s, f_z_target->re_s);
      currentWorker->sub(s1.im_s, f_z_target->im_s);
      s1.recip();
      s1.mul(&s2); //s1=A

      double dist_to_A=currentWorker->toDouble(s1.dist2_tmp(z));
      double f_z_err=currentWorker->toDouble(f_z.dist2_tmp(f_z_target));
      double f_zz_mag=currentWorker->toDouble(f_zz.getMagTmp());
      double expected=f_z_err/f_zz_mag;
      (void)expected;
      if (dist_to_A*3<closest_accepted)
      {
        closest_rejected=closest_accepted;
        closest_accepted=dist_to_A;
        sumA_cnt=1;
        currentWorker->assign(sumA.re_s, s1.re_s);
        currentWorker->assign(sumA.im_s, s1.im_s);
        reducedbymag=MandelMath::gcd(period, i+1);
        currentWorker->assign(&first_multi_re, f_z.re_s);
        currentWorker->assign(&first_multi_im, f_z.im_s);
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
  currentWorker->assign(s1.re_s, f_z_target->re_s);
  currentWorker->assign(s1.im_s, f_z_target->im_s);
  double f_z_err=currentWorker->toDouble(f_z.getMagTmp())-currentWorker->toDouble(s1.getMagTmp());
  double f_zz_mag=currentWorker->toDouble(f_zz.getMagTmp());
  double expected=std::abs(f_z_err)/f_zz_mag;
  if (reducedbymag>=period)
    multi=1;
  else if (closest_accepted>expected*1000)
    multi=1;
  else if (closest_accepted>expected*3)
    multi=1;
  else if (closest_rejected<13*closest_accepted) //12.27 at second 76/38/19
    multi=1;
  else if (closest_rejected>100000*closest_accepted)
  {
    multi=period/reducedbymag;
    currentWorker->zero(s1.re_s, sumA_cnt);
    currentWorker->recip(s1.re_s);
    currentWorker->mul(&sumA_re, s1.re_s);
    currentWorker->mul(&sumA_im, s1.re_s);
  }
  else if (closest_rejected<100*closest_accepted)
    multi=1;
  else if (closest_rejected>1000*closest_accepted)
  {
    multi=period/reducedbymag;
    currentWorker->zero(s1.re_s, sumA_cnt);
    currentWorker->recip(s1.re_s);
    currentWorker->mul(&sumA_re, s1.re_s);
    currentWorker->mul(&sumA_im, s1.re_s);
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
      multi=0;
    else if (sumnear1==period*(near1+1)/2-near1)
    {
      multi=near1;
    }
    else
      multi=0;
  }
  return true;
}

bool MandelLoopEvaluator::eval2zzc(int period, const MandelMath::complex *c, const MandelMath::complex *z)
{
  currentWorker->assign(&f_re, z->re_s);
  currentWorker->assign(&f_im, z->im_s);
  currentWorker->zero(&f_z_re, 1);
  currentWorker->zero(&f_z_im);
  currentWorker->zero(&f_c_re);
  currentWorker->zero(&f_c_im);
  currentWorker->zero(&f_zz_re);
  currentWorker->zero(&f_zz_im);
  currentWorker->zero(&f_zc_re);
  currentWorker->zero(&f_zc_im);
  currentWorker->zero(&f_cc_re);
  currentWorker->zero(&f_cc_im);
  currentWorker->zero(&f_zzc_re);
  currentWorker->zero(&f_zzc_im);
  MandelMath::complex f(currentWorker, &f_re, &f_im, true);
  MandelMath::complex f_z(currentWorker, &f_z_re, &f_z_im, true);
  MandelMath::complex f_c(currentWorker, &f_c_re, &f_c_im, true);
  MandelMath::complex f_zz(currentWorker, &f_zz_re, &f_zz_im, true);
  MandelMath::complex f_zc(currentWorker, &f_zc_re, &f_zc_im, true);
  MandelMath::complex f_cc_(currentWorker, &f_cc_re, &f_cc_im, true);
  MandelMath::complex f_zzc(currentWorker, &f_zzc_re, &f_zzc_im, true);
  //MandelMath::complex s1(currentWorker, &s1_re, &s1_im, true);
  MandelMath::complex s2(currentWorker, &s2_re, &s2_im, true);
  for (int i=0; i<period; i++)
  {
    double m_sum=currentWorker->toDouble(f.getMagTmp())+
                 currentWorker->toDouble(f_z.getMagTmp())+
                 currentWorker->toDouble(f_zz.getMagTmp())+
                 currentWorker->toDouble(f_c.getMagTmp())+
                 currentWorker->toDouble(f_zc.getMagTmp())+
                 currentWorker->toDouble(f_cc_.getMagTmp())+
                 currentWorker->toDouble(f_zzc.getMagTmp());
    if (m_sum>MandelEvaluator::LARGE_FLOAT2)
      return false;
    //f_zzc=d/dc 2*(f_z*f_z + f*f_zz)=2*(2*f_z*f_zc + f_c*f_zz+f*f_zzc)
    f_zzc.mul(&f);
    currentWorker->assign(s2.re_s, f_c.re_s);
    currentWorker->assign(s2.im_s, f_c.im_s);
    s2.mul(&f_zz);
    f_zzc.add(&s2);
    currentWorker->assign(s2.re_s, f_z.re_s);
    currentWorker->assign(s2.im_s, f_z.im_s);
    s2.mul(&f_zc);
    currentWorker->lshift(s2.re_s, 1);
    currentWorker->lshift(s2.im_s, 1);
    f_zzc.add(&s2);
    currentWorker->lshift(f_zzc.re_s, 1);
    currentWorker->lshift(f_zzc.im_s, 1);
    //f_cc=2 * (f_c^2 + f * f_cc)
    f_cc_.mul(&f);
    currentWorker->assign(s2.re_s, f_c.re_s);
    currentWorker->assign(s2.im_s, f_c.im_s);
    s2.sqr();
    f_cc_.add(&s2);
    currentWorker->lshift(f_cc_.re_s, 1);
    currentWorker->lshift(f_cc_.im_s, 1);
    // f_zc = 2 * (f_z * f_c + f * f_zc);
    f_zc.mul(&f);
    currentWorker->assign(s2.re_s, f_c.re_s);
    currentWorker->assign(s2.im_s, f_c.im_s);
    s2.mul(&f_z);
    f_zc.add(&s2);
    currentWorker->lshift(f_zc.re_s, 1);
    currentWorker->lshift(f_zc.im_s, 1);
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(&f);
    currentWorker->assign(s2.re_s, f_z.re_s);
    currentWorker->assign(s2.im_s, f_z.im_s);
    s2.sqr();
    f_zz.add(&s2);
    currentWorker->lshift(f_zz.re_s, 1);
    currentWorker->lshift(f_zz.im_s, 1);
    // f_c = 2 * f * f_c + 1;
    f_c.mul(&f);
    currentWorker->lshift(f_c.re_s, 1);
    currentWorker->lshift(f_c.im_s, 1);
    currentWorker->add_double(f_c.re_s, 1);
    //f_z:=2*f*f_z
    f_z.mul(&f);
    currentWorker->lshift(f_z.re_s, 1);
    currentWorker->lshift(f_z.im_s, 1);
    //f:=f^2+c
    f.sqr();
    f.add(c);
  }
  return true;
}


bool MandelEvaluator::findBulbBase(int period2, const MandelMath::complex *c, MandelMath::complex *cb, MandelMath::complex *rb, MandelMath::complex *xc, MandelMath::complex *baseZC, MandelMath::complex *baseCC, bool *is_card, int *foundMult)
//on input foundMult=0 -> guess rb here; =1 -> rb already set
//xc: bulb center, both z and c
//cb: bulb base c, rb: bulb base root (final point)
{ //"findBulbBaseOri"
  if (period2==1)
  {
    currentWorker->zero(cb->re_s, 0.25);
    currentWorker->zero(cb->im_s, 0);
    currentWorker->zero(rb->re_s, 0.5);
    currentWorker->zero(rb->im_s, 0);
    currentWorker->zero(xc->re_s, 0);
    currentWorker->zero(xc->im_s, 0);
    currentWorker->zero(baseZC->re_s, 0);
    currentWorker->zero(baseZC->im_s, 0);
    currentWorker->zero(baseCC->re_s, 0);
    currentWorker->zero(baseCC->im_s, 0);
    *is_card=true;
    *foundMult=2;
    return true;
  };
  int period=period2;
  currentWorker->assign(xc->re_s, c->re_s);
  currentWorker->assign(xc->im_s, c->im_s);
  /*enum {stLockF, stLockCard, stNewtonFFZ, stNewton2FFZ, stRepeatFFZ,
              stDecide,            stNewtonFFC, stNewton2FFC, stRepeatFFC,
                                   stNewtonFCFZ, stNewton2FCFZ, stRepeatFCFZ, //works if locked at the central root, preferably close to the base
              stLockF2,
              stNewtonX, //works inside the bulb if |ff|<1
                //inCardioid, do a 2D step to minimize f_z then switch to stNewton2
                //otherwise if |ff| small goto stNewtonFCFZ else perform some serious magic on dc and dr and go to stNewtonX2
              stNewtonX2} state; //...perform laguerre on rb until fm is small, then switch to stNewtonX*/
  /*double order1;
  if (period<1024)
    order1=std::ldexp(1, -period);
  else
    order1=0;*/
  //state=stLockF;
  //double prevFF=1e10;
  //bool Result=false;
  //bool isCardioid=false;
  *is_card=false;
  bool did_reduce_period=false;

  MandelMath::complex f(currentWorker, &bulb.bulbe.f_re, &bulb.bulbe.f_im, true);
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
  currentWorker->zero(target_f_z.im_s, 0);
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
      if (!bulb.bulbe.evalg(period, xc))
        return false;
      currentWorker->rsub(f.re_s, xc->re_s);
      currentWorker->rsub(f.im_s, xc->im_s);
      //g.add(xc);
      currentWorker->add_double(f_c.re_s, -1);
      double g_mag=currentWorker->toDouble(f.getMagTmp());
      //ideally copy the Lagu code from newton() but newton should be enough
      double g_c_mag=currentWorker->toDouble(f_c.getMagTmp());
      double g_cc_mag;
      if (g_mag>currentWorker->eps2()*1000)
      { //xc=xc-2*g/g_c
        g_cc_mag=1;
        f_c.recip_prepared();
        f.mul(&f_c);
        currentWorker->lshift(f.re_s, 1);
        currentWorker->lshift(f.im_s, 1);
        xc->add(&f);
      }
      else
      { //xc=xc-g_c/g_cc
        g_cc_mag=currentWorker->toDouble(f_cc_.getMagTmp());
        f_cc_.recip_prepared();
        f_c.mul(&f_cc_);
        currentWorker->sub(xc->re_s, f_c.re_s);
        currentWorker->sub(xc->im_s, f_c.im_s);
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
  currentWorker->assign(rb->re_s, xc->re_s);
  currentWorker->assign(rb->im_s, xc->im_s);
  currentWorker->assign(cb->re_s, xc->re_s);
  currentWorker->assign(cb->im_s, xc->im_s);
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
    bulb.bulbe.eval2(period, cb, rb); //TODO: dont' need f_cc here
    currentWorker->sub(f.re_s, rb->re_s);
    currentWorker->sub(f.im_s, rb->im_s);
    currentWorker->add_double(f_z.re_s, -1);
#endif
    //for cardioid, we need to use f_zz as well in the equation for f because it does not go to 0
    //0=f=0+C*fc+R*fz+R*R*fzz/2    C=-R*(fz+R*fzz/2)/fc
    //target_fz-1=fz+C*fzc+R*fzz   target_fz-1=fz-R*(fz+R*fzz/2)/fc*fzc+R*fzz    target_fz-1-fz=-R*R*fzz/fc*fzc/2+R*(fzz-fz/fc*fzc)
    //0=R*R*fzz/2/fc*fzc+R*(fz/fc*fzc-fzz)+target_fz-fz-1
    //x1,2= -(2*c)/(b+-sqrt(b^2-2*a*2*c))         good except both b,c small e.g. 0
    currentWorker->assign(s1.re_s, f_c.re_s);
    currentWorker->assign(s1.im_s, f_c.im_s);
    s1.recip();
    s1.mul(&f_zc);
    currentWorker->assign(s2_.re_s, s1.re_s);
    currentWorker->assign(s2_.im_s, s1.im_s); //f_zc/f_c
    s1.mul(&f_z);
    currentWorker->sub(s1.re_s, f_zz.re_s);
    currentWorker->sub(s1.im_s, f_zz.im_s);  //s1=b
    //currentWorker->lshift(s1.re_s, -1);
    //currentWorker->lshift(s1.im_s, -1);
    s2_.mul(&f_zz); //s2=2*a
    currentWorker->assign(s3.re_s, f_z.re_s);
    currentWorker->assign(s3.im_s, f_z.im_s);
    currentWorker->sub(s3.re_s, target_f_z.re_s);
    currentWorker->sub(s3.im_s, target_f_z.im_s);
    currentWorker->add_double(s3.re_s, 1); //(fz-target_fz+1)
    currentWorker->lshift(s3.re_s, 1);
    currentWorker->lshift(s3.im_s, 1); //s3=-2*c
    s2_.mul(&s3);
    currentWorker->assign(deltar.re_s, s1.re_s);
    currentWorker->assign(deltar.im_s, s1.im_s);
    deltar.sqr();
    currentWorker->add(deltar.re_s, s2_.re_s);
    currentWorker->add(deltar.im_s, s2_.im_s); //b^2-4ac
    deltar.sqrt();
    if (currentWorker->toDouble(deltar.mulreT(&s1))<0)
    {
      currentWorker->chs(deltar.re_s);
      currentWorker->chs(deltar.im_s);
    };
    deltar.add(&s1);
    deltar.recip();
    deltar.mul(&s3);
    //C=-R*(fz+R*fzz/2)/fc
    currentWorker->assign(deltac.re_s, f_zz.re_s);
    currentWorker->assign(deltac.im_s, f_zz.im_s);
    deltac.mul(&deltar);
    currentWorker->lshift(deltac.re_s, -1);
    currentWorker->lshift(deltac.im_s, -1);
    deltac.add(&f_z);
    deltac.mul(&deltar);
    currentWorker->assign(s1.re_s, f_c.re_s);
    currentWorker->assign(s1.im_s, f_c.im_s);
    s1.recip();
    deltac.mul(&s1);
    currentWorker->chs(deltar.re_s);
    currentWorker->chs(deltar.im_s);
    currentWorker->assign(s1.re_s, deltac.re_s);
    currentWorker->assign(s1.im_s, deltac.im_s);
    currentWorker->assign(s2_.re_s, deltar.re_s);
    currentWorker->assign(s2_.im_s, deltar.im_s);


    //check the easy way:
    //        0=f=0+C*fc+R*fz -> R=-C*fc/fz
    //target_fz-1=fz+C*fzc+R*fzz   fz/(fc/fz*fzz-fzc)=C=-deltaC    (fz-target_fz+1)/(fzc-fc/fz*fzz)=deltaC
    currentWorker->assign(deltac.re_s, f_z.re_s);
    currentWorker->assign(deltac.im_s, f_z.im_s);
    deltac.recip();
    deltac.mul(&f_c);
    currentWorker->assign(deltar.re_s, deltac.re_s); // fc/fz
    currentWorker->assign(deltar.im_s, deltac.im_s);
    deltac.mul(&f_zz);
    currentWorker->rsub(deltac.re_s, f_zc.re_s);
    currentWorker->rsub(deltac.im_s, f_zc.im_s);
#if 0
    if (cycle==0)
    {
      currentWorker->assign(&bulb.test_x0_re, deltac.re_s);
      currentWorker->assign(&bulb.test_x0_im, deltac.im_s);
      test_x0_mag=currentWorker->toDouble(deltac.getMagTmp());
    }
    currentWorker->assign(&bulb.test_xn_re, deltac.re_s);
    currentWorker->assign(&bulb.test_xn_im, deltac.im_s);
#endif
    deltac.recip(); // 1/(fzc-fc/fz*fzz)
    currentWorker->assign(s3.re_s, f_z.re_s);
    currentWorker->assign(s3.im_s, f_z.im_s);
    currentWorker->sub(s3.re_s, target_f_z.re_s);
    currentWorker->sub(s3.im_s, target_f_z.im_s);
    currentWorker->add_double(s3.re_s, 1); //(fz-target_fz+1)
    deltac.mul(&s3);  //deltac
    deltar.mul(&deltac);
    currentWorker->chs(deltar.re_s); // deltar
    currentWorker->chs(deltar.im_s);

#if 0
    //at card, we are not at result yet so 1/f_z does not blow clearly enough
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
      /*currentWorker->assign(deltac.re_s, s1.re_s);
      currentWorker->assign(deltac.im_s, s1.im_s);
      currentWorker->assign(deltar.re_s, s2_.re_s);
      currentWorker->assign(deltar.im_s, s2_.im_s);*/
      currentWorker->lshift(deltac.re_s, -1);
      currentWorker->lshift(deltac.im_s, -1);
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



    currentWorker->sub(cb->re_s, deltac.re_s);
    currentWorker->sub(cb->im_s, deltac.im_s);
    currentWorker->sub(rb->re_s, deltar.re_s);
    currentWorker->sub(rb->im_s, deltar.im_s);
    if (cycle==0)
    {
      currentWorker->assign(&bulb.dbg_first_cb_re, cb->re_s);
      currentWorker->assign(&bulb.dbg_first_cb_im, cb->im_s);
      currentWorker->assign(&bulb.dbg_first_rb_re, rb->re_s);
      currentWorker->assign(&bulb.dbg_first_rb_im, rb->im_s);
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

    currentWorker->assign(s1.re_s, rb->re_s);
    currentWorker->assign(s1.im_s, rb->im_s);

    //now we guessed step in cb and rb, need to tune rb so that f(cb, rb)=0 again
    for (int lcycle=0; lcycle<10; lcycle++)
    {
      //x^5=1 -> x:=x-(x^5-1)/(5*x^4)   NestList[#-(#^5-1)/(5*#^4) &, 0.6+0.8I, 10]
      //x=a^(1/5) |y=x/a| a*y=a^(1/5) y=a^(-4/5)  y^5=a^-4  y:=y-(y^5-a^-4)/(5y^4)  y:=(4/5)*y+a^-4/(5y^4)
      //                                          y^-5=a^4  y:=y-(y^-5-a^4)/(-5*y^-6)=y*(6-a^4*y^5)/5   x=a*y
      if (!bulb.bulbe.eval_zz(period, cb, rb))
      {
        *foundMult=1;
        return false;
      };
      currentWorker->sub(f.re_s, rb->re_s);
      currentWorker->sub(f.im_s, rb->im_s);
      currentWorker->add_double(f_z.re_s, -1);
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
      if (currentWorker->is0(f.re_s) && currentWorker->is0(f.im_s))
        break;
      if (!bulb.lagu.eval(period, &f, &f_z, &f_zz))
      {
        *foundMult=1;
        return false;
      };
      currentWorker->sub(rb->re_s, &bulb.lagu.step_re);
      currentWorker->sub(rb->im_s, &bulb.lagu.step_im);
      if (currentWorker->toDouble(f.getMagTmp())<3*currentWorker->eps2())
        break;
    }

    if (cycle==0)
    { //no idea why but first guess, after fixing rb, has f_z either 0+0i at bulb, or 0+-i at cardioid
      currentWorker->assign(s2_.re_s, f_z.re_s);
      currentWorker->assign(s2_.im_s, f_z.im_s);
      double dist_to_0=currentWorker->toDouble(s2_.getMagTmp());
      currentWorker->add_double(s2_.im_s, 1);
      double dist_to_ni=currentWorker->toDouble(s2_.getMagTmp());
      currentWorker->assign(s2_.im_s, f_z.im_s);
      currentWorker->add_double(s2_.im_s, -1);
      double dist_to_pi=currentWorker->toDouble(s2_.getMagTmp());
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
    currentWorker->assign(s2_.re_s, f_z.re_s);
    currentWorker->assign(s2_.im_s, f_z.im_s);
    currentWorker->sub(s2_.re_s, target_f_z.re_s);
    currentWorker->sub(s2_.im_s, target_f_z.im_s);
    currentWorker->add_double(s2_.re_s, 1);
    f_error=currentWorker->toDouble(f.getMagTmp()); //supposed to be 0 but just in case
    fz_error=currentWorker->toDouble(s2_.getMagTmp());
    if (!*is_card && f_error<1e-15 && fz_error<1e-4)
    {
      if (!bulb.bulbe.eval_multi(period, cb, rb, &target_f_z))
      {
        *foundMult=1;
        return false;
      };
      if (bulb.bulbe.multi>1)
      {
        currentWorker->assign(rb->re_s, &bulb.bulbe.sumA_re);
        currentWorker->assign(rb->im_s, &bulb.bulbe.sumA_im);
        /* don't need them anyway
        currentWorker->sub(f.re_s, rb->re_s);
        currentWorker->sub(f.im_s, rb->im_s);
        currentWorker->add_double(f_z.re_s, -1);
        */
        //reduce period
        *foundMult *= bulb.bulbe.multi;
        period/=bulb.bulbe.multi;
        did_reduce_period=true;

        //fix rb after moving to guessed root position, we rely on f(cb, rb)==0
        for (int lcycle=0; lcycle<7; lcycle++)
        {
          //x^5=1 -> x:=x-(x^5-1)/(5*x^4)   NestList[#-(#^5-1)/(5*#^4) &, 0.6+0.8I, 10]
          //x=a^(1/5) |y=x/a| a*y=a^(1/5) y=a^(-4/5)  y^5=a^-4  y:=y-(y^5-a^-4)/(5y^4)  y:=(4/5)*y+a^-4/(5y^4)
          //                                          y^-5=a^4  y:=y-(y^-5-a^4)/(-5*y^-6)=y*(6-a^4*y^5)/5   x=a*y
          if (!bulb.bulbe.eval_zz(period, cb, rb))
          {
            *foundMult=1;
            return false;
          };
          currentWorker->sub(f.re_s, rb->re_s);
          currentWorker->sub(f.im_s, rb->im_s);
          currentWorker->add_double(f_z.re_s, -1);
          if (currentWorker->is0(f.re_s) && currentWorker->is0(f.im_s))
            break;
          if (!bulb.lagu.eval(period, &f, &f_z, &f_zz))
          {
            *foundMult=1;
            return false;
          };
          currentWorker->sub(rb->re_s, &bulb.lagu.step_re);
          currentWorker->sub(rb->im_s, &bulb.lagu.step_im);
          if (currentWorker->toDouble(f.getMagTmp())<3*currentWorker->eps2())
            break;
        }

        //change target f_z to multi-th root of 1, using root near first_multi
        //could do a lot of tricks but let's just use newton to find the root of newf_z^multi=1
        //or rather new_target^multi=old_target
        //starting from bulb.bulbe.first_multi
        MandelMath::complex next_target(currentWorker, &bulb.bulbe.first_multi_re, &bulb.bulbe.first_multi_im, true);
        for (int multicyc=0; multicyc<10; multicyc++)
        {
          //a=target_f_z  y=next_target  y:=y*(6-a^4*y^5)/5   x=a*y
          currentWorker->assign(s3.re_s, next_target.re_s);
          currentWorker->assign(s3.im_s, next_target.im_s);
          s3.mul(&target_f_z);
          currentWorker->assign(s2_.re_s, s3.re_s);
          currentWorker->assign(s2_.im_s, s3.im_s);
          for (int i=2; i<bulb.bulbe.multi; i++)
            s3.mul(&s2_);
          s3.mul(&next_target); //a^4*y^5
          currentWorker->chs(s3.re_s);
          currentWorker->chs(s3.im_s);
          currentWorker->add_double(s3.re_s, 1); //(1-a^4*y^5)
          currentWorker->zero(s2_.re_s, bulb.bulbe.multi);
          currentWorker->recip(s2_.re_s);
          currentWorker->mul(s3.re_s, s2_.re_s);
          currentWorker->mul(s3.im_s, s2_.re_s); //(1-a^4*y^5)/5
          double dist1=currentWorker->toDouble(s3.getMagTmp());
          s3.mul(&next_target); //y*(1-a^4*y^5)/5
          next_target.add(&s3); //y+=y*(1-a^4*y^5)/5
          if (dist1<3*currentWorker->eps2())
            break;
        }
        target_f_z.mul(&next_target);
      }
    }
    currentWorker->sub(s1.re_s, rb->re_s);
    currentWorker->sub(s1.im_s, rb->im_s);
    currentWorker->sub(s1.re_s, deltar.re_s);
    currentWorker->sub(s1.im_s, deltar.im_s); //should be around 0
    if (f_error<3*currentWorker->eps2() &&
        fz_error<3*(1+currentWorker->toDouble(f_zz.getMagTmp()))*currentWorker->eps2() &&
        (did_reduce_period || *is_card))
    {
      nop();
      break;
    };
    nop();
  }
  currentWorker->assign(baseZC->re_s, &bulb.bulbe.f_zc_re);
  currentWorker->assign(baseZC->im_s, &bulb.bulbe.f_zc_im);
  currentWorker->assign(baseCC->re_s, &bulb.bulbe.f_cc_re);
  currentWorker->assign(baseCC->im_s, &bulb.bulbe.f_cc_im);
  return (*foundMult > 1) || *is_card;









  //suppose the center is exact enough
  //2) find bulb base guess c = xc+1/(f_zc+f_zz)   , derivatives at (z=xc,c=xc) (from estimateInterior)
  bulb.bulbe.eval2(period, xc, xc);
  currentWorker->sub(&bulb.bulbe.f_re, xc->re_s);//is 0
  currentWorker->sub(&bulb.bulbe.f_im, xc->im_s);
  currentWorker->add_double(&bulb.bulbe.f_z_re, -1); //is -1

  //s2:=f_zc^2-f_zz*f_cc
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
   /*   if complexOps.mag(xc)<1e-30 then
        begin //converged to (0,0), so the period is clearly wrong
          foundMult:=1;
          Result:=False;
          Exit;
        end
      else
        begin
          //happens if the search for xc flew to a distant root that is lower period
          //or the period cleaner didn't try hard enough
          foundMult:=1;
          Result:=False;
          Exit;
        end*/
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
}

void MandelEvaluator::fixRnearBase(MandelMath::complex *r, const MandelMath::complex *c, int period, int *mult)
{ //"cleverFix" in old code
  //TODO: cycles unused?
  MandelMath::complex rb(currentWorker, &bulb.rb_re_, &bulb.rb_im, false);
  MandelMath::complex cb(currentWorker, &bulb.cb_re, &bulb.cb_im, false);
  MandelMath::complex xc(currentWorker, &bulb.xc_re, &bulb.xc_im, false);
  MandelMath::complex baseZC_(currentWorker, &bulb.baseZC_re, &bulb.baseZC_im, false);
  MandelMath::complex baseCC_(currentWorker, &bulb.baseCC_re, &bulb.baseCC_im, false);
  MandelMath::complex s1(currentWorker, &bulb.s1_re, &bulb.s1_im, false);
  MandelMath::complex s2(currentWorker, &bulb.s2_re_, &bulb.s2_im_, false);
  MandelMath::complex cbx_(currentWorker, &bulb.cbx_re, &bulb.cbx_im, false);
  MandelMath::complex rbx_(currentWorker, &bulb.rbx_re, &bulb.rbx_im, false);
  if (*mult<=1)
    return;
  currentWorker->assign(rb.re_s, r->re_s);
  currentWorker->assign(rb.im_s, r->im_s);
  bool is_card=false;
  int foundMult=1;
  bool baseFound=findBulbBase(period, c, &cb, &rb, &xc, &baseZC_, &baseCC_, &is_card, &foundMult);
  if (!baseFound)
    return;
  if (foundMult<=1)
    return;
  if (*mult!=foundMult)
    *mult=foundMult;
  //s1=c-cb
  currentWorker->assign(s1.re_s, c->re_s);
  currentWorker->assign(s1.im_s, c->im_s);
  currentWorker->chs(s1.re_s);
  currentWorker->chs(s1.im_s);
  currentWorker->sub(s1.re_s, cb.re_s);
  currentWorker->sub(s1.im_s, cb.im_s);
  /*currentWorker->assign(rbx.re_s, s1.re_s);
  currentWorker->assign(rbx.im_s, s1.im_s);
  rbx.add(&xc);*/
  //double ratio=currentWorker->toDouble(s1.getMagTmp())/currentWorker->toDouble(cbx.getMagTmp());
  //if (ratio<0.05)
  //s2=s1*baseZC ~= root_z
  currentWorker->assign(s2.re_s, s1.re_s);
  currentWorker->assign(s2.im_s, s1.im_s);
  s2.mul(&baseZC_);
  currentWorker->assign(r->re_s, rb.re_s);
  currentWorker->assign(r->im_s, rb.im_s);
  if (!currentWorker->isl0(s2.re_s) || (*mult==2)) //=0 for baseZC=0 at period=1
  { //above base
    currentWorker->assign(cbx_.re_s, cb.re_s);
    currentWorker->assign(cbx_.im_s, cb.im_s);
    currentWorker->sub(cbx_.re_s, xc.re_s);
    currentWorker->sub(cbx_.im_s, xc.im_s);
    currentWorker->assign(rbx_.re_s, rb.re_s);
    currentWorker->assign(rbx_.im_s, rb.im_s);
    currentWorker->sub(rbx_.re_s, xc.re_s);
    currentWorker->sub(rbx_.im_s, xc.im_s);
    currentWorker->assign(s2.re_s, cbx_.re_s);
    currentWorker->assign(s2.re_s, cbx_.im_s);
    s2.recip();
    s2.mul(&rbx_);
    double tmpre=currentWorker->toDouble(s2.getMagTmp());
    if (tmpre==0)
    {
      //r:=rb
    }
    else
    {
      double tmpim=std::atan2(currentWorker->toDouble(s2.im_s), currentWorker->toDouble(s2.re_s));
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
      currentWorker->zero(s2.re_s, -tmpre*cos(tmpim));
      currentWorker->zero(s2.im_s, -tmpre*sin(tmpim));
      s2.mul(&rbx_);
      //r:=rb+s2
      r->add(&s2);
    }
  }
  else
  { //we're below the bulb base, move close to the base's root
    //dz f_zc+dc/2 f_cc=0
    //dz=-dc/2 f_cc/f_zc
    s1.mul(&baseCC_);
    currentWorker->assign(s2.re_s, baseZC_.re_s);
    currentWorker->assign(s2.im_s, baseZC_.im_s);
    s2.recip();
    s2.mul(&s1);
    currentWorker->lshift(s2.re_s, -1);
    currentWorker->lshift(s2.im_s, -1);
    currentWorker->chs(s2.re_s);
    currentWorker->chs(s2.im_s);
    r->add(&s2);
  }
}

//result 0..derivatives or value too large, or other fail (divide by 0)
//result>0 .. tried to return multiplicity but really returns just 1 (1 or >=3) or 2 (mult==2)
int MandelEvaluator::newton(int period, const complex *c, complex *r, const bool fastHoming, const int suggestedMultiplicity) //returns multiplicity
{ //TODO: suggestedMulti = maximumMultip ?
  double bestfm=1e10; //TODO: actually bestgm? g(z)=f(z)-z
  currentWorker->assign(&newt.bestr_re, r->re_s); //init, cleanup
  currentWorker->assign(&newt.bestr_im, r->im_s);
  bool movedOff=false;
  //double accyBound=3e-28/(period*period);
  //was for 80b floats double accyBound2=3e-39*period/log(1+period)*1.5; //1.5=magic
  //double accyBound2=1.23e-32*period/log(1+period)*1.5; //1.5=magic
  double order1; // 1/highest power in the polynomial, 2^period in case of mandelbrot set
  double r_mag_rough;
  int maxm;
  {
    currentWorker->assign(&newtres_.first_guess_lagu_re, r->re_s);
    currentWorker->assign(&newtres_.first_guess_lagu_im, r->im_s);
    currentWorker->assign(&newtres_.first_guess_newt_re, r->re_s);
    currentWorker->assign(&newtres_.first_guess_newt_im, r->im_s);
    double r_re=currentWorker->toDouble(r->re_s);
    double r_im=currentWorker->toDouble(r->im_s);
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
  complex f_r(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  complex fz_r(currentWorker, &newtres_.fz_r_re_, &newtres_.fz_r_im_, true);
  complex fzz_r(currentWorker, &newt.fzz_r_re, &newt.fzz_r_im, true);
  complex tmp1(currentWorker, &newt.tmp1_re, &newt.tmp1_im, true);
  //complex fzfix(currentWorker, &newt.fzfix_re, &newt.fzfix_im, true);
  complex laguH(currentWorker, &newt.laguH_re, &newt.laguH_im, true);
  complex laguG(currentWorker, &newt.laguG_re, &newt.laguG_im, true); //fzf
  complex laguG2(currentWorker, &newt.laguG2_re, &newt.laguG2_im, true); //G^2
  complex laguX(currentWorker, &newt.laguX_re, &newt.laguX_im, true); //Laguerre step
  complex newtX(currentWorker, &newt.newtX_re, &newt.newtX_im, true); //Laguerre step
  complex fzzf(currentWorker, &newt.fzzf_re, &newt.fzzf_im, true);
  for (int newtonCycle=0; newtonCycle<50; newtonCycle++)
  {
    newtres_.cyclesNeeded=newtonCycle;
    if ((movedOff) && (newtonCycle>10) && (order1>=0))
    {                                    //  p m -> p
      order1=-1;                         //  2 2    1
      bestfm=1e10;                       //  4 3    2
      //multiplicity1=1;                   //  4 5    1
    };
    currentWorker->assign(&newt.f_r_re, r->re_s);
    currentWorker->assign(&newt.f_r_im, r->im_s);
    currentWorker->zero(&newtres_.fz_r_re_, 1.0);
    currentWorker->zero(&newtres_.fz_r_im_, 0);
    currentWorker->zero(&newt.fzz_r_re, 0);
    currentWorker->zero(&newt.fzz_r_im, 0);
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
      currentWorker->lshift(fzz_r.re_s, 1);
      currentWorker->lshift(fzz_r.im_s, 1);
      //fz:=2*f*fz
      fz_r.mul(&f_r);
      currentWorker->lshift(fz_r.re_s, 1);
      currentWorker->lshift(fz_r.im_s, 1);

      //eps_cumul05=4*eps_cumul05*f_r_mag+0.5;
      eps_cumul10=4*eps_cumul10*f_r_mag+1;
      eps_cumul=4*eps_cumul*f_r_mag;
      //err_cumul+=1/currentWorker->toDouble(fz_r.getMagTmp());
      if (fz_r_mag<=1)
        err_simple++;;
      //fc=2*fc*f+1
      double f_re=currentWorker->toDouble(f_r.re_s);
      double f_im=currentWorker->toDouble(f_r.im_s);
      double tmp=2*(f_re*fc_re-f_im*fc_im)+1;
      fc_im=2*(f_im*fc_re+f_re*fc_im);
      fc_re=tmp;

      //f:=f^2+c
      f_r.sqr();
      f_r.add(c);
      double f_rough=currentWorker->radixfloor(f_r.re_s, c->re_s);
      eps_cumul+=f_rough*f_rough;
      f_rough=currentWorker->radixfloor(f_r.im_s, c->im_s);
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
    currentWorker->sub(f_r.re_s, r->re_s);
    currentWorker->sub(f_r.im_s, r->im_s);
    currentWorker->add_double(fz_r.re_s, -1);
    double g_r_mag=currentWorker->toDouble(f_r.getMagTmp());
    double gz_r_mag=currentWorker->toDouble(fz_r.getMagTmp());
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
    if (currentWorker->isequal(r->re_s, &newt.bestr_re) &&
        currentWorker->isequal(r->im_s, &newt.bestr_im) &&
        (bestfm<1e10) && (newtonCycle>2)) //Lagu can cycle in first 2 cycles
    { //Laguerre can cycle (at least my version), e.g. per=2, c=-0.6640625-0.015625i, r=-0.614213552325963974-0,0179149806499481201i
      if (g_r_mag<6*eps_cumul10*r_mag_rough*currentWorker->eps2()) //should be tested above but maybe use different margin here?
        return lastm;
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
        double h_re=currentWorker->toDouble(laguH.re_s);
        double h_im=currentWorker->toDouble(laguH.im_s);
        double h_mag=h_re*h_re+h_im*h_im;
        double g2_re=currentWorker->toDouble(laguG2.re_s);
        double g2_im=currentWorker->toDouble(laguG2.im_s);
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
      double G2_mag=currentWorker->toDouble(laguG2.getMagTmp());
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
        currentWorker->assign(laguX.re_s, laguG2.re_s);
        currentWorker->assign(laguX.im_s, laguG2.im_s);
        currentWorker->chs(laguX.im_s);
        laguX.mul(&laguH);
        double a_re=currentWorker->toDouble(laguX.re_s)/G2_mag*(1-order1)-order1;
        double a_im=currentWorker->toDouble(laguX.im_s)/G2_mag*(1-order1);
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
    lastm=m;
    bool lagu_valid=false;
    bool newt_valid=false;
    if (order1>=0)
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
    if (gz_r_mag!=0)
    {
      //newton near multiroot:
      //f=(x-a)^m   f'=m*(x-a)^(m-1)  f/f'=(x-a)/m
      //Newton corrected for multiroot = f/f'*m
      //1/M=1-f''*f/f'^2   x-a=1/(f'/f-f''/f')
      currentWorker->assign(newtX.re_s, fz_r.re_s);
      currentWorker->assign(newtX.im_s, fz_r.im_s);
      newtX.recip();
      newtX.mul(&f_r); //f/f'
      if (m!=1)
      {
        currentWorker->zero(&newt.tmp2, m);
        currentWorker->mul(newtX.re_s, &newt.tmp2);
        currentWorker->mul(newtX.im_s, &newt.tmp2);
      };
      newt_valid=true;
    };
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
        if (N_mag*1.05>L_mag) //5% will do no harm, and switch to Lagu can speed up convergence
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
  return lastm;
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
  if ((currentData.lookper_prevGuess_>0) &&
      ((currentData.lookper_lastGuess % currentData.lookper_prevGuess_)==0))
    aroundCount=currentData.lookper_lastGuess / currentData.lookper_prevGuess_;
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
  complex fz_r(currentWorker, &newtres_.fz_r_re_, &newtres_.fz_r_im_, true);
  //double dist_around=5*currentWorker->eps2()/currentWorker->toDouble(fz_r.getMagTmp());
  //need at least safety factor of 14.2
  //correct root maps to .to_stop, then also 2*error can map to .to_stop
  //and we oscillate +-error so that's 4*to_stop, or 16 in dist squared
  double dist_around=16*(newtres_.accy_multiplier*newtres_.accy_tostop)*currentWorker->eps2();
  currentWorker->assign(f_r.re_s, root.re_s);
  currentWorker->assign(f_r.im_s, root.im_s);
  double f_e_re=sqrt(dist_around/6), f_e_im=0;
  currentWorker->zero(fz_r.re_s, 1);
  currentWorker->zero(fz_r.im_s, 0);
  double fz_e_m=0;
  const MandelMath::number_store *fz_mag1=nullptr;
  aroundCount=0;
  bool someBelow1=false;

  //TODO: try to guess the radius of "around" roots   double dist_around=newt.fz*2/fzz  f(r+x)=r+x=r+f'x+f''x^2/2

  //int firstBelow1=-1;
  int firstBelow1dividing=-1;
  for (int i=0; i<period; i++)
  {
    //if (fz_mag1 && currentWorker->isl0(fz_mag1)) //I think we intentionally skip last fz_mag
    if (fz_mag1 &&  //I think we intentionally skip last fz_mag
        (currentWorker->isle0(fz_mag1) ||
         (MandelMath::sqr_double(currentWorker->toDouble(fz_mag1))<=fz_e_m)))
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
          double dist=currentWorker->toDouble(f_r.dist2_tmp(&root));
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
      double fzre=currentWorker->toDouble(fz_r.re_s);
      double fzim=currentWorker->toDouble(fz_r.im_s);
      double re=2*(fzre*f_e_re-fzim*f_e_im);
      double im=2*(fzim*f_e_re+fzre*f_e_im);
      fz_e_m=re*re+im*im;
    }
    fz_r.mul(&f_r);
    currentWorker->lshift(fz_r.re_s, 1);
    currentWorker->lshift(fz_r.im_s, 1);
    fz_mag1=fz_r.getMag1Tmp();
    if (currentWorker->toDouble(fz_mag1)>LARGE_FLOAT2)
      return -1; //so is it checked or not
    //f:=f^2+c
    //use EXACTLY same code as the main loop: if temporaries are more precise at one place, it will never work!
    {
      double fre=currentWorker->toDouble(f_r.re_s);
      double fim=currentWorker->toDouble(f_r.im_s);
      double re=fre*f_e_re-fim*f_e_im;
      f_e_im=fim*f_e_re+fre*f_e_im;
      f_e_re=re;
    }
    f_r.sqr();
    f_r.add(c);
  };
  //what's going on here -.-  normally we would check |ff|<1 or something
  //maybe that's what I'm doing but with some extra tricks for high precision?
  //reduce fz_re, fz_im to 2nd octant
  if (currentWorker->isle0(fz_r.re_s)) currentWorker->chs(fz_r.re_s);
  if (currentWorker->isle0(fz_r.im_s)) currentWorker->chs(fz_r.im_s);
  bool fz_r_mag_over1;
  if (!currentWorker->isle(fz_r.re_s, fz_r.im_s))
  {
    currentWorker->assign(&newt.tmp2, fz_r.re_s);
    currentWorker->assign(fz_r.re_s, fz_r.im_s);
    currentWorker->assign(fz_r.im_s, &newt.tmp2);
  };
  //im>=re, now check re^2+im^2>1
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
    return 0;
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
  MandelMath::complex fc_c(currentWorker, &currentData.fc_c_re_, &currentData.fc_c_im_, true);
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

  for (; (currentData.iter<currentParams.maxiter_) && (currentData.state==MandelPoint::State::stUnknown); currentData.iter++)
  {
    if (currentData.iter%(3*currentData.near0iter) ==0)
    {
      int quot=currentData.iter/(3*currentData.near0iter);
      if ((quot&(quot-1))==0) //also at iter==0  //TODO: maybe better 3*(2^k-1)*near not 3*(2^k)*near
      { // //need k*iter for f' to start at the worst moment to reduce false positives; need k*iter-1 for good near0 -> switch to nearc
        currentData.lookper_startiter=currentData.iter;
        currentWorker->assign(&currentData.lookper_startf_re, f.re_s);
        currentWorker->assign(&currentData.lookper_startf_im, f.im_s);
        //MandelMath::complex lookper_bestf(currentWorker, &eval.lookper_startf_re, &eval.lookper_startf_im, true);
        //currentWorker->zero(&eval.lookper_bestf_re, 0);
        //currentWorker->zero(&eval.lookper_bestf_im, 0);
        currentWorker->assign(&eval.lookper_nearr_re, f.re_s);
        currentWorker->assign(&eval.lookper_nearr_im, f.im_s);
        if (currentData.iter<=1)
          currentWorker->assign(&currentData.lookper_nearr_dist_, f.getMagTmp());
        else
          currentData.lookper_nearr_dist_touched=false;//currentWorker->assign(&currentData.lookper_nearr_dist, f.dist2_tmp(&c));
        //currentWorker->zero(&eval.lookper_dist2, 1e10); //4.0 should be enough
        //mands.period stays there
        //currentData.lookper_prevGuess=0; //TODO: used for anything?
        currentData.lookper_lastGuess=0;
        currentWorker->zero(&currentData.lookper_totalFzmag, 1.0);
      };
    }
    const MandelMath::number_store *f_mag=f.getMagTmp();
    if (currentWorker->toDouble(f_mag)>4)
    {
      currentData.state=MandelPoint::State::stOutside;
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
      currentData.period=currentData.lookper_lastGuess; //preliminary
      if (currentData.period<1)
        currentData.period=1;
      break;
    };
    double fc_c_mag=currentWorker->toDouble(fc_c.getMagTmp());
    if (fc_c_mag>1e57)
    {
      currentData.state=MandelPoint::State::stBoundary;
      currentData.exterior_avoids=0;
      currentData.exterior_hits=0;
      break;
    };
    double fz_c_mag=currentWorker->toDouble(&currentData.fz_c_mag_);
    if (fz_c_mag>1e60)
    {
      currentData.state=MandelPoint::State::stDiverge;
      currentData.exterior_avoids=0;
      currentData.exterior_hits=0;
      break;
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
      break;
    };*/
    f_mag=f.getMagTmp();
    //f'=2*f'*f, f'_mag=4*f'_mag*f_mag
    currentWorker->mul(&currentData.fz_c_mag_, f_mag); //TODO: can use f_mag from above?
    currentWorker->mul(&currentData.lookper_totalFzmag, f_mag);
    currentWorker->lshift(&currentData.lookper_totalFzmag, 2);
    //f:=f^2+c
    f.sqr();
    f.add(&c);
    //currentData.iter++;
    f_mag=f.getMagTmp();

    if (!currentWorker->isle(&eval.near0fmag, f_mag)) //f_mag<near0fmag
    {
      currentData.near0iter=currentData.iter+2;
      currentWorker->assign(&currentData.near0f_re, f.re_s);
      currentWorker->assign(&currentData.near0f_im, f.im_s);
      currentWorker->assign(&eval.near0fmag, f_mag);
    };

    const MandelMath::number_store *lpdiff=lookper_startf.dist2_tmp(&f);
    switch (currentWorker->compare(lpdiff, &currentData.lookper_nearr_dist_)) //|f-r|<best
    {
      case -1:
      {
        currentWorker->assign(&eval.lookper_nearr_re, f.re_s);
        currentWorker->assign(&eval.lookper_nearr_im, f.im_s);
        currentWorker->assign(&currentData.lookper_nearr_dist_, lpdiff);
        currentData.lookper_nearr_dist_touched=false;
        currentData.lookper_prevGuess_=currentData.lookper_lastGuess;
        currentData.lookper_lastGuess=(currentData.iter+1-currentData.lookper_startiter);
      } break;
      case 0:
      {
        if (!currentData.lookper_nearr_dist_touched)
        {
          currentWorker->assign(&eval.lookper_nearr_re, f.re_s);
          currentWorker->assign(&eval.lookper_nearr_im, f.im_s);
          currentWorker->assign(&currentData.lookper_nearr_dist_, lpdiff);
          currentData.lookper_nearr_dist_touched=true;
          currentData.lookper_prevGuess_=currentData.lookper_lastGuess;
          currentData.lookper_lastGuess=(currentData.iter+1-currentData.lookper_startiter);
        };
      } break;
    };

    if ((currentData.lookper_lastGuess>0) &&
        (currentData.lookper_lastGuess==(currentData.iter+1-currentData.lookper_startiter)) && //just found new guess
#if USE_GCD_FOR_CHECKPERIOD
#else
        ((currentData.near0iter % currentData.lookper_lastGuess)==0) && //  period divides nearest, that's a fact
#endif
        ((currentData.iter>=3*currentData.near0iter)))  //speedup - don't check period too eagerly
    {
#if USE_GCD_FOR_CHECKPERIOD
      int testperiod=MandelMath::gcd(currentData.near0iter, currentData.lookper_lastGuess);//currentData.lookper_lastGuess
#else
      int testperiod=currentData.lookper_lastGuess;
#endif
      int foundperiod=-1;
      if (f.isequal(&lookper_startf))
      { //exact match - misiurewicz or converged after too many steps
        foundperiod=currentData.lookper_lastGuess;
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
        foundperiod=periodCheck(testperiod, &c); //updates iter, f, f_c, root
        //profiler.timeInPeriod:=profiler.timeInPeriod+(getTime64()-periodEntered);
        if ((foundperiod>0) && (foundperiod<testperiod))
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
        currentData.period=foundperiod;
        currentData.newton_iter=newtres_.cyclesNeeded;
        complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
        currentData.period=estimateInterior(foundperiod, &c, &root);
        if (currentWorker->isle0(&interior.inte_abs)) //wanted <0 here
          currentData.state=MandelPoint::State::stMisiur;
        else
        {
          currentData.interior=currentWorker->toDouble(&interior.inte_abs);
          currentWorker->assign(&currentData.fz_r_re, &interior.fz_re); //d/dz F_c(r)
          currentWorker->assign(&currentData.fz_r_im, &interior.fz_im);
          if (foundperiod!=currentData.period)
            currentData.state=MandelPoint::State::stPeriod3;
        }
        //currentWorker->assign(&currentData.fc_c_re, &eval.fz_r_re);
        //currentWorker->assign(&currentData.fc_c_im, &eval.fz_r_im);
        break;
      };
      if (currentParams.breakOnNewNearest)
      {
        currentData.iter++;
        break;
      }
    };

  }
  //data.state=MandelPoint::State::stMaxIter;
  if (!currentData.has_fc_r && currentParams.want_fc_r &&
      ((currentData.state==MandelPoint::State::stPeriod2) ||
       (currentData.state==MandelPoint::State::stPeriod3)))
  {
    MandelMath::complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
    this->bulb.bulbe.eval2(currentData.period, &c, &root);
    currentWorker->assign(&currentData.fc_c_re_, &bulb.bulbe.f_c_re);
    currentWorker->assign(&currentData.fc_c_im_, &bulb.bulbe.f_c_im);
    currentData.has_fc_r=true;
  };

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
  maxiter_=1;
  breakOnNewNearest=false;
  want_fc_r=false;
}

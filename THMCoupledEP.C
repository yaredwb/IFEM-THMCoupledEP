// $Id$
//==============================================================================
//!
//! \file THMCoupledEP.C
//!
//! \date
//!
//! \author Yared Worku Bekele
//!
//! \brief Integrand implementations for THM coupled problems in ground freezing
//!
//==============================================================================

#include "THMCoupledEP.h"
#include "ASMbase.h"
#include "ASMmxBase.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Tensor.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Vec3Oper.h"
#include "VTF.h"
#include "StabilizationUtils.h"

//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  Uo = 0,            //!< Previous displacement
  Po = 1,            //!< Previous pore pressure
  To = 2,            //!< Previous temperature
  Uc = 3,            //!< Current displacement
  Pc = 4,            //!< Current pore pressure
  Tc = 5,            //!< Current temperature
  NSOL = 6           //!< Number of solution vectors
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Ru = 0,            //!< External RHS vector from equilibrium equation (M)
  Rp = 1,            //!< External RHS vector from mass balance equation (H)
  RT = 2,            //!< External RHS vector from energy balance equation (T)
  Rprev = 3,         //!< Internal RHS vector from previous step (THM)
  Rnow = 4,          //!< Internal RHS vector for current step (THM)
  Rres = 5,          //!< Residual RHS vector (THM)
  Rc = 6,            //!< Weak Dirichlet boundary contribution to RHS
  NVEC = 7           //!< Number of RHS vectors
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  uu = 0,            //!< Stiffness matrix
  up = 1,            //!< Mechanical-hydraulic coupling matrix
  uT = 2,            //!< Mechanical-thermal coupling matrix
  pu = 3,            //!< Hydro-mechanical coupling matrix
  pp = 4,            //!< Hydraulic matrix
  pT = 5,            //!< Hydro-thermal coupling matrix
  Tp = 6,            //!< Thermo-hydraulic coupling matrix
  TT = 7,            //!< Thermal matrix
  Kprev = 8,         //!< Fully coupled THM matrix from previous step
  Know = 9,          //!< Fully coupled THM matrix for current step
  Ktan = 10,         //!< Fully coupled THM Jacobian matrix
  NMAT = 11          //!< Number of LHS element matrices
};


THMCoupledEP::MixedElmMats::MixedElmMats()
{
  this->resize(NMAT,NVEC);
}


const Matrix& THMCoupledEP::MixedElmMats::getNewtonMatrix() const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  size_t ru = A[uu].rows();
  size_t rp = A[pp].rows();

  for (size_t i = 1; i <= ru; i++)
  {
    for (size_t j = 1; j <= ru; j++)
    {
      N(i,j) = A[uu](i,j);
    }
    for (size_t j = 1; j <= rp; j++)
    {
      size_t k = ru + 2*j - 1;
      size_t l = ru + 2*j;
      N(i,k) = A[up](i,j);
      N(k,i) = A[pu](j,i);
      N(i,l) = A[uT](i,j);
    }
  }

  for (size_t i = 1; i <= rp; i++)
  {
    for (size_t j = 1; j <= rp; j++)
    {
      size_t ki = ru + 2*i - 1;
      size_t kj = ru + 2*j - 1;
      size_t li = ru + 2*i;
      size_t lj = ru + 2*j;
      N(ki,kj) = A[pp](i,j);
      N(ki,lj) = A[pT](i,j);
      N(li,lj) = A[TT](i,j);
      N(li,kj) = A[Tp](i,j);
    }
  }

  return A[Ktan];
}


const Vector& THMCoupledEP::MixedElmMats::getRHSVector() const
{
  Vector& R = const_cast<Vector&>(b[Rres]);

  size_t ru = b[Ru].size();
  size_t rp = b[Rp].size();

  for (size_t i = 1; i <= ru; i++)
    R(i) = b[Ru](i);

  for (size_t i = 1; i <= rp; i++)
  {
    R(ru+2*i-1) = b[Rp](i);
    R(ru+2*i)   = b[RT](i);
  }

  R += b[Rprev];
  R -= b[Rnow];

  // Robin boundary contribution to the residual
  // TO DO: This should be done properly by calling finalizeElementBou after
  //        integration of the Robin boundary terms.
  R -= A[Know]*b[Rc];

  return b[Rres];
}


THMCoupledEP::THMCoupledEP(unsigned short int n, int order, bool stab) :
  nsd(n), gacc(9.81), mat(nullptr), SUPG(stab)
{
  primsol.resize(1 + order);
  tracFld = nullptr;
  fluxFld = nullptr;
}


Vec3 THMCoupledEP::getTraction(const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


LocalIntegral* THMCoupledEP::getLocalIntegral(const std::vector<size_t>& nen,
                                              size_t, bool neumann) const
{
  const size_t nedof1 = nsd*nen[0];           //!< Number of DOFs on basis 1
  const size_t nedof  = nedof1 + 2*nen[1];

  ElmMats* result = new MixedElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Ru].resize(nedof1);
  result->b[Rp].resize(nen[1]);
  result->b[RT].resize(nen[1]);
  result->b[Rprev].resize(nedof);
  result->b[Rnow].resize(nedof);
  result->b[Rres].resize(nedof);
  result->b[Rc].resize(nedof);

  if(!neumann)
  {
    result->A[uu].resize(nedof1,nedof1);
    result->A[up].resize(nedof1,nen[1]);
    result->A[uT].resize(nedof1,nen[1]);
    result->A[pu].resize(nen[1],nedof1);
    result->A[pp].resize(nen[1],nen[1]);
    result->A[pT].resize(nen[1],nen[1]);
    result->A[Tp].resize(nen[1],nen[1]);
    result->A[TT].resize(nen[1],nen[1]);
    result->A[Kprev].resize(nedof,nedof);
    result->A[Know].resize(nedof,nedof);
    result->A[Ktan].resize(nedof,nedof);
  }

  return result;
}


bool THMCoupledEP::initElement(const std::vector<int>& MNPC,
                               const std::vector<size_t>& elem_sizes,
                               const std::vector<size_t>& basis_sizes,
                               LocalIntegral& elmInt)
{
  if(primsol.front().empty())
    return true;

  // Extract element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart), nsd, primsol[0], elmInt.vec[Uc])
           + utl::gather(IntVec(fstart, MNPC.end()), 0, 2, primsol[0], elmInt.vec[Pc], nsd*basis_sizes[0], basis_sizes[0])
           + utl::gather(IntVec(fstart, MNPC.end()), 1, 2, primsol[0], elmInt.vec[Tc], nsd*basis_sizes[0], basis_sizes[0])
           + utl::gather(IntVec(MNPC.begin(), fstart), nsd, primsol[1], elmInt.vec[Uo])
           + utl::gather(IntVec(fstart, MNPC.end()), 0, 2, primsol[1], elmInt.vec[Po], nsd*basis_sizes[0], basis_sizes[0])
           + utl::gather(IntVec(fstart, MNPC.end()), 1, 2, primsol[1], elmInt.vec[To], nsd*basis_sizes[0], basis_sizes[0]);

  if (ierr == 0)
    return true;

  std::cerr << " *** THMCoupledEP::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool THMCoupledEP::initElementBou(const std::vector<int>& MNPC,
                                  const std::vector<size_t>& elem_sizes,
                                  const std::vector<size_t>& basis_sizes,
                                  LocalIntegral& elmInt)
{
  return this->IntegrandBase::initElementBou(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool THMCoupledEP::evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                             const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!mat)
  {
    std::cerr << __FUNCTION__ << ": No material data." << std::endl;
    return false;
  }

  size_t i,j,k;

  Matrix Bmat, Cmat, CB;

  // Get the strain-displacement matrix
  if(!mat->formBmatrix(Bmat,fe.grad(1),nsd))
    return false;

  // Get the updated pressure and temperature values at the current point
  double p = elMat.vec[Pc].dot(fe.basis(2));
  double T = elMat.vec[Tc].dot(fe.basis(2));
  // Rescale the pressure and temperature values
  p *= scl1;
  T *= scl2;

  // Evaluate the updated material tangent stiffness
  if(!mat->formElasticMatrix(Cmat,X,nsd,p,T))
    return false;

  // Integration of the stiffness matrix
  CB.multiply(Cmat,Bmat,false,false);
  CB *= 1.0*fe.detJxW;
  elMat.A[uu].multiply(Bmat,CB,true,false,true);

  // Get the densities of water and ice and the porosity
  double rhow = mat->getWaterDensity(X);
  double rhoi = mat->getIceDensity(X);
  double poro = mat->getPorosity(X);

  // Get the latent heat of fusion
  double Lf = mat->getLatentHeat(X);

  // Get the bulk moduli of water and ice
  double Kw = mat->getBulkWater(X);
  double Ki = mat->getBulkIce(X);

  // Evaluate the ice pressure at the current point
  double pi = mat->getIcePressure(X,p,T);

  // Evaluate the updated reference temperature
  double Tf = mat->getUpdatedRefTemp(X,pi,T);   // NB: pi is passed-in instead of p

  // Get the the reference pressure, freezing temperature and pressure melting parameter
  double pref = mat->getRefPressure(X);
  double Tf0  = mat->getFreezingTemp(X);
  double a    = mat->getPresMeltingParam(X);

  // Calculate dTf/dpi
  double dTfdpi = (Tf0 / a / pref) * pow(pi/pref + 1.0, 1.0/a - 1.0);

  // Calculate Mf
  double Mf = 1.0 + (Lf / Ki) * log(T/Tf) - (rhoi * Lf / Tf) * dTfdpi;

  // Evaluate the water capacities (Sp,ST)
  double Sp = mat->getWaterCapacity(X,p,T,true);
  double ST = mat->getWaterCapacity(X,p,T,false);
  // Evaluate the relative permeability coefficient and the unfrozen permeability
  double kr = mat->getRelPermCoeff(X,p,T);
  Vec3 perm = mat->getPermeability(X);

  // Evaluate the degrees of saturation (Sw,Si)
  double Sw = mat->getWaterSaturation(X,p,T);
  double Si = mat->getIceSaturation(X,p,T);

#if SP_DEBUG > 3
  std::cout << "Ice pressure at current point,  pi = " << pi << std::endl;
  std::cout << "Updated reference temperature,  Tf = " << Tf << std::endl;
  std::cout << "Degree of water saturation,     Sw = " << Sw << std::endl;
  std::cout << "Degree of ice saturation,       Si = " << Si << std::endl;
  std::cout << "Isothermal water capacity,      Sp = " << Sp << std::endl;
  std::cout << "Nonisothermal water capacity,   ST = " << ST << std::endl;
  std::cout << "Relative permeability coeff.,   kr = " << kr << std::endl;
#endif

  // Define the unit Voigt vector
  Vector m, Cm;
  m.resize(Cmat.rows());
  for (i = 1; i <= Cmat.rows(); i++)
    if (i <= nsd)
      m(i) = 1;

  // Evaluate Cmat*m
  Cm.resize(Cmat.rows());
  Cm = Cmat*m;

  // Integration of the mechanical-hydraulic matrix
  Matrix Cuptmp;
  const size_t nstrc = nsd*(nsd+1)/2;
  Cuptmp.resize(nstrc,fe.basis(2).size());
  // Evaluate the derivatives of total stress and phase change strain wrt pressure
  //double dsdp = Sp*p + Sw;
  double dsdp = Sp * (p - pi) + (Sw + Si * (rhoi/rhow));
  double dedp = poro * (rhoi - rhow) * Sp / (rhow*Sw + rhoi*Si) / 3.0;

#if SP_DEBUG > 3
  std::cout << "Change of stress with pore pressure,       dsdp = " << dsdp << std::endl;
  std::cout << "Change of phase strain with pore pressure, dedp = " << dedp << std::endl;
#endif

  for (i = 1; i <= nstrc; i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      Cuptmp(i,j) += -1.0 * scl1 * (dsdp*m(i) + dedp*Cm(i)) * fe.basis(2)(j) * fe.detJxW;

  elMat.A[up].multiply(Bmat,Cuptmp,true,false,true);

  // Integration of the mechanical-thermal matrix
  // Evaluate the derivatives of total stress and phase change strain wrt temperature
  //double dsdT = ST*p;
  double dsdT = ST * (p - pi) - Si * rhoi * Lf/T;
  double dedT = poro * (rhoi - rhow) * ST / (rhow*Sw + rhoi*Si) / 3.0;

#if SP_DEBUG > 3
  std::cout << "Change of stress with temperature,       dsdT = " << dsdT << std::endl;
  std::cout << "Change of phase strain with temperature, dedT = " << dedT << std::endl;
#endif

  Matrix CuTtmp;
  CuTtmp.resize(nstrc,fe.basis(2).size());
  for (i = 1; i <= nstrc; i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CuTtmp(i,j) += -1.0 * scl2 * (dsdT*m(i) + dedT*Cm(i)) * fe.basis(2)(j) * fe.detJxW;

  elMat.A[uT].multiply(Bmat,CuTtmp,true,false,true);

  // Integration of the hydraulic matrix (coefficient to constant p)
  Matrix Kpp;
  Kpp.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      for (k = 1; k <= nsd; k++)
        Kpp(i,j) += scl1 * scl1 * fe.grad(2)(i,k) * (kr/(rhow*gacc)) * perm[k-1] * fe.grad(2)(j,k) * fe.detJxW;

  // Integration of the hydro-mechanical matrix
  Matrix Cputmp;
  Cputmp.resize(fe.basis(2).size(),nstrc);
  for (i = 1; i <= fe.basis(2).size(); i++)
    for(j = 1; j <= nstrc; j++)
      Cputmp(i,j) += scl1 * (Sw + (rhoi/rhow)*Si) * fe.basis(2)(i) * m(j) * fe.detJxW;

  elMat.A[pu].multiply(Cputmp,Bmat,false,false,true);

  // Integration of the hydraulic matrix
  Matrix Cpp;
  Cpp.resize(fe.basis(2).size(),fe.basis(2).size());
  double Mp = (poro * Sw / Kw) + (poro * Si / Ki) * (rhoi/rhow) * (1.0/Mf);
  double Np = poro * (1 - (rhoi/rhow)) * Sp;
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      Cpp(i,j) += scl1 * scl1 * fe.basis(2)(i) * (Mp + Np) * fe.basis(2)(j) * fe.detJxW;

  elMat.A[pp] += Cpp;
  elMat.A[pp].add(Kpp,time.dt);

  // Integration of the hydro-thermal matrix
  Matrix CpT;
  CpT.resize(fe.basis(2).size(),fe.basis(2).size());
  double MT = -1.0 * (poro * Si / Ki) * (rhoi/rhow) * (rhoi * Lf / T) * (1.0/Mf);
  double NT = poro * (1 - (rhoi/rhow)) * ST;
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CpT(i,j) += scl1 * scl2 * fe.basis(2)(i) * (MT + NT) * fe.basis(2)(j) * fe.detJxW;

  elMat.A[pT] += CpT;

  // Integration of the thermal matrix (Coefficient to constant T)
  Matrix KTT;
  KTT.resize(fe.basis(2).size(),fe.basis(2).size());
  double cw = mat->getWaterHeatCapacity(T);
  // Evaluate the current pressure gradient
  Vector gradP;
  fe.grad(2).multiply(elMat.vec[Pc],gradP,true);
  // Evalaute the current Darcy velocity
  Vec3 grav = this->getGravity();
  Vec3 vel;
  for (k = 1; k <= nsd; k++)
    vel[k-1] = -1.0 * scl1 * (kr/(rhow*gacc)) * perm[k-1] * (gradP(k) - rhow*grav[k-1]);

#if SP_DEBUG > 3
  std::cout << "Darcy velocity at current point, vel = " << vel << std::endl;
#endif

  // Evaluate the overall thermal conductivity
  double lambda = mat->getThermalConductivity(X,p,T);
  for (i = 1; i <= fe.basis(2).size(); i++)
  {
    for (j = 1; j <= fe.basis(2).size(); j++)
    {
      double laplace = 0.0, convection = 0.0;
      for (k = 1; k <= nsd; k++)
      {
        laplace += fe.grad(2)(i,k) * fe.grad(2)(j,k);
        convection += fe.basis(2)(i) * rhow * cw * vel[k-1] * fe.grad(2)(j,k);
      }
      KTT(i,j) += scl2 * scl2 * (laplace * lambda + convection) * fe.detJxW;
    }
  }

  // Integration of the thermo-hydraulic matrix
  Matrix CTp;
  CTp.resize(fe.basis(2).size(),fe.basis(2).size());
  double xi = (poro *  rhoi) / (Sw + (rhoi/rhow) * Si);
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CTp(i,j) += scl1  * scl2 * fe.basis(2)(i) * Lf * xi * Sp * fe.basis(2)(j) * fe.detJxW;

  // Integration of the thermal matrix
  Matrix CTT;
  CTT.resize(fe.basis(2).size(),fe.basis(2).size());
  double rhoc_eff = mat->getEffHeatCapacity(X,p,T);
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CTT(i,j) += scl2 * scl2 * fe.basis(2)(i) * (rhoc_eff + Lf * xi * ST) * fe.basis(2)(j) * fe.detJxW;

  elMat.A[TT] += CTT;
  elMat.A[TT].add(KTT,time.dt);

#if SP_DEBUG > 3
  utl::zero_print_tol = 1e-30;
  std::cout << "Cup = " << elMat.A[up] << std::endl;
  std::cout << "CuT = " << elMat.A[uT] << std::endl;
  std::cout << "Cuu = " << elMat.A[uu] << std::endl;
  std::cout << "Kpp = " << Kpp << std::endl;
  std::cout << "Cpp = " << Cpp << std::endl;
  std::cout << "App = " << elMat.A[pp] << std::endl;
  std::cout << "CTp = " << CTp << std::endl;
  std::cout << "KTT = " << KTT << std::endl;
  std::cout << "CTT = " << CTT << std::endl;
  std::cout << "ATT = " << elMat.A[TT] << std::endl;
#endif

  // Evaluate elements of Kprev on basis 2
  size_t ru = elMat.A[uu].rows();
  size_t rp = elMat.A[pp].rows();

  for (i = 1; i <= rp; i++)
  {
    for (j = 1; j <= rp; j++)
    {
      size_t ki = ru + 2*i - 1;
      size_t kj = ru + 2*j - 1;
      elMat.A[Kprev](ki,kj) += Cpp(i,j);
      size_t li = ru + 2*i;
      size_t lj = ru + 2*j;
      elMat.A[Kprev](li,lj) += CTT(i,j);
    }
  }

  // If SUPG stabilization is required, evaluate the stabilizing matrices
  if (SUPG)
  {
    std::cout << "Stabilization matrices not implemented!" << std::endl;
  }

  // Add flow driving body forces to the RHS
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (k = 1; k <= nsd; k++)
      elMat.b[Rp](i) += time.dt * fe.grad(2)(i,k) * (kr/(rhow*gacc)) * perm[k-1] * rhow * grav[k-1] * fe.detJxW;

  return true;
}


bool THMCoupledEP::evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                             const TimeDomain& time, const Vec3& X,
                             const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr << "*** THMCoupledEP::evalBouMx: No fluxes/tractions." << std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the surface traction
  Vec4 Xt = static_cast<const Vec4&>(X);
  Xt.t = time.t;
  Vec3 tr2 = this->getTraction(Xt,normal);
  Xt.t -= time.dt;
  Vec3 tr1 = this->getTraction(Xt,normal);
  Vec3 dtr;
  dtr = tr2 - tr1;

  // Integration of Ru
  for (size_t i = 1; i <= fe.basis(1).size(); i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      elMat.b[Ru](nsd*(i-1)+j) += dtr[j-1] * fe.basis(1)(i) * fe.detJxW;

  // Integration of Rp and Rt
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
  {
    for (size_t k = 1; k <= nsd; k++)
    {
      elMat.b[Rp](i) += 0.0;
      elMat.b[RT](i) += 0.0;
    }
  }

  return true;
}


bool THMCoupledEP::finalizeElement(LocalIntegral& elmInt, const TimeDomain&, size_t)
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  size_t ru = elMat.A[uu].rows();
  size_t rp = elMat.A[pp].rows();

  // Evaluate the previous and updated internal matrices
  for (size_t i = 1; i <= ru; i++)
  {
    for (size_t j = 1; j <= ru; j++)
    {
      elMat.A[Kprev](i,j) = elMat.A[uu](i,j);
      elMat.A[Know](i,j)  = elMat.A[uu](i,j);
    }
    for (size_t j = 1; j <= rp; j++)
    {
      size_t k = ru + 2*j - 1;
      elMat.A[Kprev](i,k) = elMat.A[up](i,j);
      elMat.A[Know](i,k)  = elMat.A[up](i,j);
      elMat.A[Kprev](k,i) = elMat.A[pu](j,i);
      elMat.A[Know](k,i)  = elMat.A[pu](j,i);
      size_t l = ru + 2*j;
      elMat.A[Kprev](i,l) = elMat.A[uT](i,j);
      elMat.A[Know](i,l)  = elMat.A[uT](i,j);
    }
  }

  for (size_t i = 1; i <= rp; i++)
  {
    for (size_t j = 1; j <= rp; j++)
    {
      size_t ki = ru + 2*i - 1;
      size_t kj = ru + 2*j - 1;
      size_t li = ru + 2*i;
      size_t lj = ru + 2*j;
      elMat.A[Know](ki,kj)  = elMat.A[pp](i,j);
      elMat.A[Know](li,lj)  = elMat.A[TT](i,j);
      elMat.A[Kprev](ki,lj) = elMat.A[pT](i,j);
      elMat.A[Kprev](li,kj) = elMat.A[Tp](i,j);
      elMat.A[Know](ki,lj)  = elMat.A[pT](i,j);
      elMat.A[Know](li,kj)  = elMat.A[Tp](i,j);
    }
  }

  // Get the previous and current solution vectors on basis 2
  Vector PoTo, PcTc;
  PoTo.resize(2*elMat.vec[Po].size());
  PcTc.resize(2*elMat.vec[Pc].size());
  for (size_t i = 1; i <= elMat.vec[Po].size(); i++)
  {
    PoTo(2*i-1) = elMat.vec[Po](i);
    PoTo(2*i)   = elMat.vec[To](i);
    PcTc(2*i-1) = elMat.vec[Pc](i);
    PcTc(2*i)   = elMat.vec[Tc](i);
  }

  // Get the previous and current solution vectors
  Vector prevSol, currSol;
  prevSol = elMat.vec[Uo];
  prevSol.insert(prevSol.end(),PoTo.begin(),PoTo.end());

  currSol = elMat.vec[Uc];
  currSol.insert(currSol.end(),PcTc.begin(),PcTc.end());

  elMat.b[Rprev] = elMat.A[Kprev]*prevSol;
  elMat.b[Rnow]  = elMat.A[Know]*currSol;

  return true;
}


bool THMCoupledEP::evalSol(Vector& s, const MxFiniteElement& fe,
                           const Vec3& X, const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes) const
{
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  IntVec MNPC2(fstart,MNPC.end());
  for (size_t i = 0; i < MNPC2.size(); i++)
    MNPC2[i] += basis_sizes[0];

  Vector disp, pres, temp;
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart), nsd, primsol[0], disp)
           + utl::gather(MNPC2, 0, 2, primsol[0], pres, nsd*basis_sizes[0], basis_sizes[0])
           + utl::gather(MNPC2, 1, 2, primsol[0], temp, nsd*basis_sizes[0], basis_sizes[0]);

  if (ierr != 0)
    std::cerr << " *** THMCoupledEP::evalSol: Detected " << ierr/3
              << " node numbers out of range." << std::endl;

  double p = pres.dot(fe.basis(2));
  p *= scl1;
  double T = temp.dot(fe.basis(2));
  T *= scl2;

  double pi = mat->getIcePressure(X,p,T);
  double Si = mat->getIceSaturation(X,p,T);
  double Sw = mat->getWaterSaturation(X,p,T);

  s.resize(9);
  s(1) = pi;
  s(2) = Si;
  s(3) = Sw;

  Matrix Bmat, Cmat;
  if(!mat->formBmatrix(Bmat,fe.grad(1),nsd))
    return false;

  if(!mat->formElasticMatrix(Cmat,X,nsd,p,T))
    return false;

  Vector strain, stress;
  strain = Bmat*disp;
  stress = Cmat*strain;

  for (size_t i = 1; i <= strain.size(); i++) {
    s(3+i) = strain(i);
    s(6+i) = stress(i);
  }

  return true;
}


size_t THMCoupledEP::getNoFields(int fld) const
{
  size_t noSecSol = 9;
  if (fld < 2)
    return nsd+2;
  else
    return noSecSol;
}


const char* THMCoupledEP::getField1Name(size_t i, const char* prefix) const
{
  static const char* s[4] = { "u_x", "u_y", "p^w", "T" };

  if(!prefix) return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


const char* THMCoupledEP::getField2Name(size_t i, const char* prefix) const
{
  size_t noSecSol = 9;
  if (i >= noSecSol) return 0;

  static const char* s2[] = {"p^i","Si","Sw","eps_x","eps_y","eps_xy",
                             "sig_x", "sig_y","sig_xy"};
  if (!prefix) return s2[i];

  static std::string name;
  name = prefix + std::string(" ") + s2[i];

  return name.c_str();
}


THMCoupledEP::WeakDirichlet::WeakDirichlet(unsigned short int n) :
    nsd(n), flux(nullptr)
{
  primsol.resize(2);
}


LocalIntegral* THMCoupledEP::WeakDirichlet::getLocalIntegral(const std::vector<size_t>& nen,
                                                             size_t, bool neumann) const
{
  const size_t nedof1 = nsd * nen[0];
  const size_t nedof  = nedof1 + 2 * nen[1];

  ElmMats* result = new MixedElmMats();

  result->withLHS = true;
  result->b[Ru].resize(nedof1);
  result->b[Rp].resize(nen[1]);
  result->b[RT].resize(nen[1]);
  result->b[Rprev].resize(nedof);
  result->b[Rnow].resize(nedof);
  result->b[Rres].resize(nedof);
  result->b[Rc].resize(nedof);

  result->A[uu].resize(nedof1,nedof1);
  result->A[up].resize(nedof1,nen[1]);
  result->A[uT].resize(nedof1,nen[1]);
  result->A[pu].resize(nen[1],nedof1);
  result->A[pp].resize(nen[1],nen[1]);
  result->A[pT].resize(nen[1],nen[1]);
  result->A[Tp].resize(nen[1],nen[1]);
  result->A[TT].resize(nen[1],nen[1]);
  result->A[Kprev].resize(nedof,nedof);
  result->A[Know].resize(nedof,nedof);
  result->A[Ktan].resize(nedof,nedof);

  return result;
}


bool THMCoupledEP::WeakDirichlet::initElementBou(const std::vector<int>& MNPC,
                                                 const std::vector<size_t>& elem_sizes,
                                                 const std::vector<size_t>& basis_sizes,
                                                 LocalIntegral& elmInt)
{
  if(primsol.front().empty())
    return true;

  // Extract element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart), nsd, primsol[0], elmInt.vec[Uc])
           + utl::gather(IntVec(fstart,MNPC.end()), 0, 2, primsol[0], elmInt.vec[Pc], nsd*basis_sizes[0], basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()), 1, 2, primsol[0], elmInt.vec[Tc], nsd*basis_sizes[0], basis_sizes[0])
           + utl::gather(IntVec(MNPC.begin(), fstart), nsd, primsol[1], elmInt.vec[Uo])
           + utl::gather(IntVec(fstart,MNPC.end()), 0, 2, primsol[1], elmInt.vec[Po], nsd*basis_sizes[0], basis_sizes[0])
           + utl::gather(IntVec(fstart,MNPC.end()), 1, 2, primsol[1], elmInt.vec[To], nsd*basis_sizes[0], basis_sizes[0]);

  if (ierr == 0)
    return true;

  std::cerr << " *** THMCoupledEP::WeakDirichlet::initElementBou: Detected "
            << ierr/3 << " node numbers out of range." << std::endl;

  return false;
}


bool THMCoupledEP::WeakDirichlet::evalBouMx(LocalIntegral& elmInt,
                                            const MxFiniteElement& fe,
                                            const TimeDomain& time,
                                            const Vec3& X,
                                            const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the Neumann heat flux on the boundary
  double qT = 0.0;
  if (flux)
    qT = (*flux)(X);

  size_t ru = elMat.A[uu].rows();

  for (size_t i = 1; i <= fe.basis(2).size(); i++)
  {
    for (size_t j = 1; j <= fe.basis(2).size(); j++)
    {
      size_t li = ru + 2*i;
      size_t lj = ru + 2*j;
      elMat.A[TT](i,j) += time.dt * fe.basis(2)(i) * lambdae * fe.basis(2)(j) * fe.detJxW;
      elMat.A[Know](li,lj) += time.dt * fe.basis(2)(i) * lambdae * fe.basis(2)(j) * fe.detJxW;
    }
    elMat.b[RT](i) += time.dt * (qT * fe.basis(2)(i) + lambdae * Te * fe.basis(2)(i)) * fe.detJxW;
  }

  // Get the previous and current solution vectors on basis 2
  Vector PcTc;
  PcTc.resize(2*elMat.vec[Pc].size());
  for (size_t i = 1; i <= elMat.vec[Pc].size(); i++)
  {
    PcTc(2*i-1) = elMat.vec[Pc](i);
    PcTc(2*i  ) = elMat.vec[Tc](i);
  }

  elMat.b[Rc] = elMat.vec[Uc];
  elMat.b[Rc].insert(elMat.b[Rc].end(),PcTc.begin(),PcTc.end());

  return true;
}
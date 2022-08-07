// $Id$
//==============================================================================
//!
//! \file Poro3Material.C
//!
//! \date
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for a freezing three-phase porous material.
//!
//==============================================================================


#include "Poro3Material.h"

#include "Functions.h"
#include "IFEM.h"
#include "tinyxml.h"
#include "Utilities.h"
#include "Vec3.h"

template<>
RealFunc* Poro3Material::FuncConstPair<RealFunc>::parse(const char* val,
                                                        const std::string& type)
{
  return utl::parseRealFunc(val,type);
}


template<>
ScalarFunc* Poro3Material::FuncConstPair<ScalarFunc>::parse(const char* val,
                                                            const std::string& type)
{
  return utl::parseTimeFunc(val,type);
}


template<>
VecFunc* Poro3Material::FuncConstPair<VecFunc>::parse(const char* val,
                                                      const std::string& type)
{
  return utl::parseVecFunc(val,type);
}


  template<class T>
static bool propertyParse(Poro3Material::FuncConstPair<T>& data,
                          const TiXmlElement* elem,
                          const std::string& attr,
                          const std::string& tag)
{
  std::string constant;
  if (utl::getAttribute(elem,attr.c_str(),constant)) {
    std::stringstream str;
    str << constant;
    str >> data.constant;
    return true;
  }

  const TiXmlElement* child = elem->FirstChildElement(tag);
  if (child) {
    IFEM::cout <<" ";
    std::string type;
    utl::getAttribute(child,"type",type,true);
    const TiXmlNode* aval;
    if ((aval = child->FirstChild()))
      data.function = data.parse(aval->Value(),type);

    return data.function != nullptr;
  }

  return false;
}


void Poro3Material::parse(const TiXmlElement* elem)
{
  propertyParse(rhow, elem, "rhow", "waterdensity");
  propertyParse(rhos, elem, "rhos", "soliddensity");
  propertyParse(rhoi, elem, "rhoi", "icedensity");
  propertyParse(wheatcapacity, elem, "cw", "waterheatcapacity");
  propertyParse(sheatcapacity, elem, "cs", "solidheatcapacity");
  propertyParse(iheatcapacity, elem, "ci", "iceheatcapacity");
  propertyParse(wconductivity, elem, "lambdaw", "waterconductivity");
  propertyParse(sconductivity, elem, "lambdas", "solidconductivity");
  propertyParse(iconductivity, elem, "lambdai", "iceconductivity");
  propertyParse(bulkw, elem, "Kw", "waterbulk");
  propertyParse(bulks, elem, "Ks", "solidbulk");
  propertyParse(bulki, elem, "Ki", "icebulk");
  propertyParse(lheat, elem, "Lf", "latentheat");
  propertyParse(porosity, elem, "poro", "porosity");
  propertyParse(sexpansion, elem, "alphas", "solidexpansion");
  propertyParse(permeability, elem, "perm", "permeability");
  propertyParse(permparam, elem, "mk", "permeabilityparam");
  propertyParse(sstiffness, elem, "Es", "solidstiffness");
  propertyParse(istiffness, elem, "Ei", "icestiffness");
  propertyParse(spoisson, elem, "nus", "solidpoisson");
  propertyParse(ipoisson, elem, "nui", "icepoisson");
  propertyParse(strparam, elem, "eta", "strengthparam");
  propertyParse(satparam1, elem, "alpha", "saturationparam1");
  propertyParse(satparam2, elem, "beta", "saturationparam2");
  propertyParse(satparam3, elem, "gamma", "saturationparam3");
  propertyParse(frtemp, elem, "Tf0", "freezingtemp");
  propertyParse(meltingparam, elem, "a", "presmeltingparam");
  propertyParse(refpres, elem, "pref", "refpressure");
}


void Poro3Material::printLog() const
{
  IFEM::cout << "\tDensities: "
             << "\n\t\tDensity of water, rhow = " << rhow.constant
             << "\n\t\tDensity of solid, rhos = " << rhos.constant
             << "\n\t\tDensity of ice,   rhoi = " << rhoi.constant << std::endl;
  IFEM::cout << "\tHeat Capacities: "
             << "\n\t\tHeat capacity of water, cw = " << wheatcapacity.constant
             << "\n\t\tHeat capacity of solid, cs = " << sheatcapacity.constant
             << "\n\t\tHeat capacity of ice,   ci = " << iheatcapacity.constant
             << std::endl;
  IFEM::cout << "\tThermal Conductivities: "
             << "\n\t\tThermal conductivity of water, lambdaw = " << wconductivity.constant
             << "\n\t\tThermal conductivity of solid, lambdas = " << sconductivity.constant
             << "\n\t\tThermal conductivity of ice,   lambdai = " << iconductivity.constant
             << std::endl;
  IFEM::cout << "\tLaternt heat of fusion, Lf = " << lheat.constant << std::endl;
  IFEM::cout << "\tPorosity, n = " << porosity.constant << std::endl;
  IFEM::cout << "\tPermeability, k = " << permeability.constant << std::endl;
  IFEM::cout << "\tPermeability parameter, mk = " << permparam.constant << std::endl;
  IFEM::cout << "\tThermal exp. coeff. of solid, alphas = " << sexpansion.constant
             << std::endl;
  IFEM::cout << "\tYoung's moduli: "
             << "\n\t\tYoung's modulus of solid, Es = " << sstiffness.constant
             << "\n\t\tYoung's modulus of ice,   Ei = " << istiffness.constant
             << std::endl;
  IFEM::cout << "\tPoisson's ratios: "
             << "\n\t\tPoisson's ratio of solid, nus = " << spoisson.constant
             << "\n\t\tPoisson's ratio of ice,   nui = " << ipoisson.constant
             << std::endl;
  IFEM::cout << "\tStrength parameter, eta = " << strparam.constant << std::endl;
  IFEM::cout << "\tSaturation curve parameters: "
             << "\n\t\tParameter 1, alpha = " << satparam1.constant
             << "\n\t\tParameter 2,  beta = " << satparam2.constant
             << "\n\t\tParameter 3, gamma = " << satparam3.constant
             << std::endl;
  IFEM::cout << "\tFreezing temperature, Tf0 = " << frtemp.constant << std::endl;
  IFEM::cout << "\tPressure melting parameter, a = " << meltingparam.constant << std::endl;
  IFEM::cout << "\tReference Pressure, pref = " << refpres.constant << std::endl;
}


double Poro3Material::getWaterDensity(const Vec3& X) const
{
  return rhow.evaluate(X);
}


double Poro3Material::getSolidDensity(const Vec3& X) const
{
  return rhos.evaluate(X);
}


double Poro3Material::getIceDensity(const Vec3& X) const
{
  return rhoi.evaluate(X);
}


double Poro3Material::getMassDensity(const Vec3& X, double p, double T) const
{
  double n = getPorosity(X);
  double Sw = getWaterSaturation(X,p,T);
  return (1.0-n) * getSolidDensity(X) + n*Sw*getWaterDensity(X) +
          n*(1.0-Sw) * getIceDensity(X);
}


double Poro3Material::getWaterHeatCapacity(double T) const
{
  return wheatcapacity.evaluate(T);
}


double Poro3Material::getSolidHeatCapacity(double T) const
{
  return sheatcapacity.evaluate(T);
}


double Poro3Material::getIceHeatCapacity(double T) const
{
  return iheatcapacity.evaluate(T);
}


double Poro3Material::getEffHeatCapacity(const Vec3& X, double p, double T) const
{
  double n = getPorosity(X);
  double Sw = getWaterSaturation(X,p,T);
  return (1.0-n) * getSolidDensity(X) * getSolidHeatCapacity(T) +
          n*Sw*getWaterDensity(X) * getWaterHeatCapacity(T) +
          n*(1.0-Sw)*getIceDensity(X) * getIceHeatCapacity(T);
}


double Poro3Material::getWaterThermalConductivity(double T) const
{
  return wconductivity.evaluate(T);
}


double Poro3Material::getSolidThermalConductivity(double T) const
{
  return sconductivity.evaluate(T);
}


double Poro3Material::getIceThermalConductivity(double T) const
{
  return iconductivity.evaluate(T);
}


double Poro3Material::getThermalConductivity(const Vec3& X, double p,
                                             double T) const
{
  double n = getPorosity(X);
  double Sw = getWaterSaturation(X,p,T);
  return pow(getSolidThermalConductivity(T),1.0-n)*
         pow(getWaterThermalConductivity(T),n*Sw)*
         pow(getIceThermalConductivity(T),n*(1.0-Sw));
}


double Poro3Material::getBulkWater(const Vec3& X) const
{
  return bulkw.evaluate(X);
}


double Poro3Material::getBulkSolid(const Vec3& X) const
{
  return bulks.evaluate(X);
}


double Poro3Material::getBulkIce(const Vec3& X) const
{
  return bulki.evaluate(X);
}


double Poro3Material::getBulkMedium(const Vec3& X) const
{
  return 0.0;   // TO DO - compressible components is a feature to be added
}


double Poro3Material::getSolidThermalExpansion(double T) const
{
  return sexpansion.evaluate(T);
}


double Poro3Material::getLatentHeat(const Vec3& X) const
{
  return lheat.evaluate(X);
}


double Poro3Material::getPorosity(const Vec3& X) const
{
  return porosity.evaluate(X);
}


Vec3 Poro3Material::getPermeability(const Vec3& X) const
{
  return permeability.evaluate(X);
}


double Poro3Material::getPermeabilityParam(const Vec3& X) const
{
  return permparam.evaluate(X);
}


double Poro3Material::getRelPermCoeff(const Vec3& X, double p, double T) const
{
  double Sw = getWaterSaturation(X,p,T);
  double m = getPermeabilityParam(X);
  double tmp = 1.0 - pow(1.0 - pow(Sw,1.0/m),m);
  return sqrt(Sw) * pow(tmp,2.0);
}


double Poro3Material::getSolidStiffness(const Vec3& X) const
{
  return sstiffness.evaluate(X);
}


double Poro3Material::getIceStiffness(const Vec3& X) const
{
  return istiffness.evaluate(X);
}


double Poro3Material::getStrengthParam(const Vec3& X) const
{
  return strparam.evaluate(X);
}


double Poro3Material::getStiffness(const Vec3& X, double p, double T) const
{
  double Si = getIceSaturation(X,p,T);
  double eta = getStrengthParam(X);
  return pow(getIceStiffness(X)/getSolidStiffness(X), pow(Si,eta)) *
         getSolidStiffness(X);
}


double Poro3Material::getSolidPoisson(const Vec3& X) const
{
  return spoisson.evaluate(X);
}


double Poro3Material::getIcePoisson(const Vec3& X) const
{
  return ipoisson.evaluate(X);
}


double Poro3Material::getPoisson(const Vec3& X, double p, double T) const
{
  double Si = getIceSaturation(X,p,T);
  double eta = getStrengthParam(X);
  return pow(getIcePoisson(X)/getSolidPoisson(X), pow(Si,eta)) *
         getSolidPoisson(X);
}


double Poro3Material::getAlphaSWCC(const Vec3& X) const
{
  return satparam1.evaluate(X);
}


double Poro3Material::getBetaSWCC(const Vec3& X) const
{
  return satparam2.evaluate(X);
}


double Poro3Material::getGammaSWCC(const Vec3& X) const
{
  return satparam3.evaluate(X);
}


double Poro3Material::getFreezingTemp(const Vec3& X) const
{
  return frtemp.evaluate(X);
}


double Poro3Material::getPresMeltingParam(const Vec3& X) const
{
  return meltingparam.evaluate(X);
}


double Poro3Material::getRefPressure(const Vec3& X) const
{
  return refpres.evaluate(X);
}


/*double Poro3Material::getWaterSaturation(const Vec3& X, double p, double T) const
{
  double Lf = getLatentHeat(X);
  double Tf = getFreezingTemp(X);
  double ri = getIceDensity(X);
  double rw = getWaterDensity(X);
  double a = getAlphaSWCC(X);
  double b = getBetaSWCC(X);
  double g = getGammaSWCC(X);
  double tmp1 = ((ri/rw)-1.0)*p - ri*Lf*log(T/Tf);
  double tmp2;
  if (tmp1 > 0.0)
    tmp2 = 1.0 + pow(a*tmp1,b);
  else
    tmp2 = 1.0 + -1.0 * pow(-1.0*a*tmp1,b);
  double Sw = 1.0 - 0.99 * (1.0 - pow(tmp2,-1.0*g));
  return (T > Tf || Sw > 1.0) ? 1.0 : Sw;
}*/

double Poro3Material::getWaterSaturation(const Vec3& X, double p, double T) const
{
  double a = getAlphaSWCC(X);
  double b = getBetaSWCC(X);
  double g = getGammaSWCC(X);
  double ri = getIceDensity(X);
  double Lf = getLatentHeat(X);
  double pi = getIcePressure(X,p,T);

  // Evaluate the updated reference temperature
  double Tf = getUpdatedRefTemp(X,pi,T);

  //double tmp1 = pi - p;
  double tmp1 = p - ri*Lf*log(T/Tf);
  double tmp2;
  if (tmp1 > 0.0)
    tmp2 = 1.0 + pow(a*tmp1,b);
  else
    tmp2 = 1.0 + -1.0 * pow(-1.0*a*tmp1,b);

  // Calculate the degree of water saturation
  double Sw = 1.0 - 0.99 * (1.0 - pow(tmp2,-1.0*g));

  return (T > Tf || Sw > 1.0) ? 1.0 : Sw;
}


double Poro3Material::getIceSaturation(const Vec3& X, double p, double T) const
{
  double Sw = getWaterSaturation(X,p,T);
  return 1.0 - Sw;
}


/*double Poro3Material::getWaterCapacity(const Vec3& X, double p, double T, bool iso) const
{
  double Lf = getLatentHeat(X);
  double Tf = getFreezingTemp(X);
  double ri = getIceDensity(X);
  double rw = getWaterDensity(X);
  double a = getAlphaSWCC(X);
  double b = getBetaSWCC(X);
  double g = getGammaSWCC(X);
  double tmp1 = ((ri/rw)-1.0)*p - ri*Lf*log(T/Tf);
  double tmp2, tmp3;
  if (tmp1 > 0.0) {
    tmp2 = 1.0 + pow(a*tmp1,b);
    tmp3 = pow(tmp1,b-1.0);
  }
  else {
    tmp2 = 1.0 + -1.0 * pow(-1.0*a*tmp1,b);
    tmp3 = -1.0 * pow(-1.0*tmp1,b-1.0);
  }
  double Sp = -1.0 * pow(a,b) * b * g * 0.99 * ((ri/rw)-1.0) * tmp3 * pow(tmp2,-1.0-g);
  double ST = pow(a,b) * b * g * 0.99 * (ri*Lf/T) * tmp3 * pow(tmp2,-1.0-g);
  if (iso)
    return T > Tf ? 0.0 : Sp;
  else
    return T > Tf ? 0.0 : ST;
}*/


double Poro3Material::getWaterCapacity(const Vec3& X, double p, double T, bool iso) const
{
  double Lf = getLatentHeat(X);
  double ri = getIceDensity(X);
  double Ki = getBulkIce(X);
  double al = getAlphaSWCC(X);
  double b = getBetaSWCC(X);
  double g = getGammaSWCC(X);
  double a = getPresMeltingParam(X);
  double pref = getRefPressure(X);
  double Tf0 = getFreezingTemp(X);
  double pi = getIcePressure(X,p,T);
  double Tf = getUpdatedRefTemp(X,p,T);

  // Calculate dTf/dpi
  double dTfdpi = (Tf0 / a / pref) * pow(pi/pref + 1.0, 1.0/a - 1.0);

  // Calculate Mf
  double Mf = 1.0 + (Lf / Ki) * log(T/Tf) - (ri * Lf / Tf) * dTfdpi;

  // Calculate P = -alpha * rhoi * Lf * ln(T/Tf)
  double P = -al * ri * Lf * log(T/Tf);

  // Calculate P^b and P^(b-1)
  double Pb, Pb1;
  if (P > 0.0)
  {
    Pb  = pow(P,b);
    Pb1 = pow(P,b-1.0);
  }
  else
  {
    Pb  = -1.0 * pow(-1.0*P,b);
    Pb1 = -1.0 * pow(-1.0*P,b-1.0);
  }

  // Calculate dP/dpw
  double tmp1 = (Mf - 1.0) / (1.0 + Lf/Ki * log(T/Tf));
  double dPdpw = -al * tmp1 * (1.0 + (ri*Lf/Tf) * (1/Mf) * dTfdpi);

  // Calculate dP/dT
  double tmp2 = ri*Lf/Tf;
  double tmp3 = pow(tmp2,2.0) / (tmp2 * (dTfdpi - 1.0));
  double dPdT = -al * ((ri*Lf/T) - tmp3 * dTfdpi);

  // Calculate the water capacities
  double tmp4 = pow(1.0+Pb,-g-1.0) ;
  double Sp = -1.0 * b * g * tmp4 * Pb1 * dPdpw;
  double ST = -1.0 * b * g * tmp4 * Pb1 * dPdT;

  // Return isothermal/non-sothermal capacity, based on inquiry
  if (iso)
    return T > Tf ? 0.0 : Sp;
  else
    return T > Tf ? 0.0 : ST;
}


/*double Poro3Material::getIcePressure(const Vec3& X, double p, double T) const
{
  double Lf = getLatentHeat(X);
  double Tf = getFreezingTemp(X);
  double ri = getIceDensity(X);
  double rw = getWaterDensity(X);
  double pi = (ri/rw) * p - ri * Lf * log(T/Tf);
  return T > Tf ? 0.0 : pi;
}*/


double Poro3Material::getIcePressure(const Vec3& X, double p, double T) const
{
  double Lf = getLatentHeat(X);
  double Ki = getBulkIce(X);
  double ri = getIceDensity(X);
  double pref = getRefPressure(X);
  double Tf0 = getFreezingTemp(X);
  double a = getPresMeltingParam(X);
  double pi;
  double Tf, f, df;
  const double tol = 0.000001;
  const int n = 15;
  for (int i = 2; i <= n; i++)
  {
    pi = p;
    Tf = Tf0 * pow(pi/pref + 1.0, 1.0/a);
    f = pi - p + ri * Lf * log(T/Tf);
    df = 1.0 - (Lf/Ki) * log(T/Tf) - ri * Lf / Tf;
    p = pi - f/df;
    //std::cout << "pi at iteration " << i << " - " << p << std::endl;
  }

  /*if (abs(pi-p) > tol)
  {
    std::cerr << "Convergence not achieved for pi!" << std::endl;
  }*/

  return T > Tf0 ? 0.0 : pi;
}


double Poro3Material::getUpdatedRefTemp(const Vec3& X, double p, double T) const
{
  double Tf0 = getFreezingTemp(X);
  double pref = getRefPressure(X);
  double a = getPresMeltingParam(X);
  double Tf = Tf0 * pow(p/pref + 1.0, 1.0/a);
  return Tf;
}


bool Poro3Material::formBmatrix(Matrix& Bmat, const Matrix& dNdX, size_t nsd) const
{
  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
  Bmat.resize(nstrc*nsd,nenod,true);
  if (dNdX.cols() < nsd)
  {
    std::cerr << " *** Poro3Material::formBmatrix: Invalid dimension on fe.grad(1), "
              << dNdX.rows() << "x" << dNdX.cols() << "." << std::endl;
    return false;
  }

#define INDEX(i,j) i+nstrc*(j-1)

  // Strain-displacement matrix for 2D elements
  for (size_t i = 1; i <= nenod; i++)
  {
    Bmat(INDEX(1,1),i) = dNdX(i,1);
    Bmat(INDEX(2,2),i) = dNdX(i,2);
    Bmat(INDEX(3,1),i) = dNdX(i,2);
    Bmat(INDEX(3,2),i) = dNdX(i,1);
  }

#undef INDEX

  Bmat.resize(nstrc,nsd*nenod);

  return true;
}


bool Poro3Material::formElasticMatrix(Matrix& Cmat, const Vec3& X, size_t nsd,
                                        double p, double T) const
{
  double E = getStiffness(X,p,T);
  double nu = getPoisson(X,p,T);

  double C33 = E/(2.0+2.0*nu);
  double C12 = (E*nu)/((1.0+nu)*(1.0-2.0*nu));
  double C11 = C12 + 2.0*C33;

  const size_t nstrc = nsd*(nsd+1)/2;

  Cmat.resize(nstrc,nstrc,true);

  Cmat(1,1) = C11;
  Cmat(1,2) = C12;
  Cmat(2,1) = C12;
  Cmat(2,2) = C11;
  Cmat(3,3) = C33;

  Cmat.resize(nstrc,nstrc);

  return true;
}

// $Id$
//==============================================================================
//!
//! \file Poro3Material.h
//!
//! \date
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for a freezing three-phase porous material.
//!
//==============================================================================

#ifndef _PORO3_MATERIAL_H
#define _PORO3_MATERIAL_H

#include "Function.h"
#include "MatVec.h"
#include "Vec3.h"
#include "Vec3Oper.h"

class TiXmlElement;


/*!
  \brief Class representing the material properties of a freezing
         three-phase porous medium.
*/

class Poro3Material
{
public:
  //! \brief Empty constructor.
  Poro3Material() {}

  //! \brief Empty destructor
  ~Poro3Material() {}

  //! \brief Parses material parameters from an XML element.
  void parse(const TiXmlElement*);

  //! \brief Prints out material parameters to the log stream.
  void printLog() const;

  //! \brief Returns the mass density of the water at the current point.
  //! \param[in] X Cartesian coordinate of the current integration point.
  double getWaterDensity(const Vec3& X) const;
  //! \brief Returns the mass density of the solid at the current point.
  double getSolidDensity(const Vec3& X) const;
  //! \brief Returns the mass density of the ice at the current point.
  double getIceDensity(const Vec3& X) const;
  //! \brief Evaluates the mean mass density.
  //! \param[in] p Updated pore water pressure at the current point.
  //! \param[in] T Updated temperature at the current point.
  double getMassDensity(const Vec3& X, double p, double T) const;
  //! \brief Returns the water heat capacity at the current point.
  double getWaterHeatCapacity(double T) const;
  //! \brief Returns the solid heat capacity at the current point.
  double getSolidHeatCapacity(double T) const;
  //! \brief Returns the ice heat capacity at the current point.
  double getIceHeatCapacity(double T) const;
  //! \brief Evaluates the effective heat capacity at the current point.
  double getEffHeatCapacity(const Vec3& X, double p, double T) const;
  //! \brief Returns the thermal conductivity of the water at the current point.
  double getWaterThermalConductivity(double T) const;
  //! \brief Returns the thermal conductivity of the solid at the current point.
  double getSolidThermalConductivity(double T) const;
  //! \brief Returns the thermal conductivity of the ice at the current point.
  double getIceThermalConductivity(double T) const;
  //! \brief Evaluates the effective thermal conductivity at the current point.
  double getThermalConductivity(const Vec3& X, double p, double T) const;
  //! \brief Returns bulk modulus of the water at the current point.
  double getBulkWater(const Vec3& X) const;
  //! \brief Returns bulk modulus of the solid at the current point.
  double getBulkSolid(const Vec3& X) const;
  //! \brief Returns the bulk modulus of the ice at the current point.
  double getBulkIce(const Vec3& X) const;
  //! \brief Evaluates the bulk modulus of the medium at the current point.
  double getBulkMedium(const Vec3& X) const;
  //! \brief Returns the thermal expansion of the solid at the current point.
  double getSolidThermalExpansion(double T) const;
  //! \brief Returns the latent heat of fusion at the current point.
  double getLatentHeat(const Vec3& X) const;
  //! \brief Returns porosity at the current point.
  double getPorosity(const Vec3& X) const;
  //! \brief Returns the permeability at the current point.
  Vec3 getPermeability(const Vec3& X) const;
  //! \brief Returns the permeability model parameter at the current point.
  double getPermeabilityParam(const Vec3& X) const;
  //! \brief Evaluates the relative permeability coefficient at the current point.
  double getRelPermCoeff(const Vec3& X, double p, double T) const;
  //! \brief Returns the Young's modulus of the solid at the current point.
  double getSolidStiffness(const Vec3& X) const;
  //! \brief Returns the Young's modulus of the ice at the current point.
  double getIceStiffness(const Vec3& X) const;
  //! \brief Returns the strength parameter at the current point.
  double getStrengthParam(const Vec3& X) const;
  //! \brief Evaluates the overall stiffness of the porous medim.
  double getStiffness(const Vec3& X, double p, double T) const;
  //! \brief Returns the Poisson's ratio of the solid at the current point.
  double getSolidPoisson(const Vec3& X) const;
  //! \brief Returns the Poisson's ratio of the ice at the current point
  double getIcePoisson(const Vec3& X) const;
  //! \brief Evaluates the overall Poisson's ratio at the current point.
  double getPoisson(const Vec3& X, double p, double T) const;
  //! \brief Returns the first parameter of the saturation curve.
  double getAlphaSWCC(const Vec3& X) const;
  //! \brief Returns the second parameter of the saturation curve.
  double getBetaSWCC(const Vec3& X) const;
  //! \brief Returns the third parameter of the saturation curve.
  double getGammaSWCC(const Vec3& X) const;
  //! \brief Returns the freezing temperature at the current point.
  double getFreezingTemp(const Vec3& X) const;
  //! \brief Returns the temperature parameter for pressure melting
  double getPresMeltingParam(const Vec3& X) const;
  //! \brief Returns the reference pressure
  double getRefPressure(const Vec3& X) const;
  //! \brief Evaluates the water saturation at the current point.
  double getWaterSaturation(const Vec3& X, double p, double T) const;
  //! \brief Evaluates the ice saturation at the current point.
  double getIceSaturation(const Vec3& X, double p, double T) const;
  //! \brief Evaluates the isothermal/non-isothermal water capacity at the current point.
  //! \param[in] iso Boolean flag for isothermal or non-isothermal capacity
  double getWaterCapacity(const Vec3& X, double p, double T, bool iso) const;
  //! \brief Evaluates the ice pressure at the current point.
  //! \details The Modified Clausius-Clapeyron equation is used here
  double getIcePressure(const Vec3& X, double p, double T) const;
  //! \brief Evaluates the updated reference temperature.
  double getUpdatedRefTemp(const Vec3& X, double p, double T) const;

  //! \brief Calculates the strain-displacement matrix.
  //! \param[in] Bmat The strain-displacement matrix
  //! \param[in] dNdX First basis function gradients at current point
  //! \param[in] nsd Number of spatial dimensions
  bool formBmatrix(Matrix& Bmat, const Matrix& dNdX, size_t nsd) const;

  //! \brief Evalutates the constitutive matrix at an integration point
  //! \param[out] Cmat Constitutive matrix at current point
  bool formElasticMatrix(Matrix& Cmat, const Vec3& X, size_t nsd,
                         double p, double T) const;

protected:
  /*! \brief Helper template for wrapping a constant/function pair */
    template<class Function>
  struct FuncConstPair
  {
    Function* function;                   //!< Function definition
    typename Function::Output constant;   //!< Constant

    //! \brief Constructor.
    FuncConstPair() { function = nullptr; constant = 0.0; }

    //! \brief Parse an XML element. Specialized per type.
    Function* parse(const char* val, const std::string& type) { return nullptr; }

    //! \brief Evaluate function.
    //! \param[in] X The value to evaluate at.
    typename Function::Output evaluate(const typename Function::Input& X) const
    {
      return function ? (*function)(X) : constant;
    }
  };

  FuncConstPair<RealFunc> rhow;               //!< Water density
  FuncConstPair<RealFunc> rhos;               //!< Solid density
  FuncConstPair<RealFunc> rhoi;               //!< Ice density

  FuncConstPair<ScalarFunc> wheatcapacity;    //!< Specific heat capacity of water
  FuncConstPair<ScalarFunc> sheatcapacity;    //!< Specific heat capacity of solid
  FuncConstPair<ScalarFunc> iheatcapacity;    //!< Specific heat capacity of ice

  FuncConstPair<ScalarFunc> wconductivity;    //!< Thermal conductivity of water
  FuncConstPair<ScalarFunc> sconductivity;    //!< Thermal conductivity of solid
  FuncConstPair<ScalarFunc> iconductivity;    //!< Thermal conductivity of ice

  FuncConstPair<RealFunc> bulkw;              //!< Bulk modulus of water
  FuncConstPair<RealFunc> bulks;              //!< Bulk modulus of solid
  FuncConstPair<RealFunc> bulki;              //!< Bulk modulus of ice

  FuncConstPair<RealFunc> lheat;              //!< Latent heat of fusion
  FuncConstPair<RealFunc> porosity;           //!< Porosity
  FuncConstPair<ScalarFunc> sexpansion;       //!< Thermal expansion of solid
  FuncConstPair<VecFunc>  permeability;       //!< Permeability
  FuncConstPair<RealFunc>  permparam;         //!< Permeability model parameter

  FuncConstPair<RealFunc> sstiffness;         //!< Young's modulus of solid
  FuncConstPair<RealFunc> istiffness;         //!< Young's modulus of ice
  FuncConstPair<RealFunc> spoisson;           //!< Poisson's ratio of solid
  FuncConstPair<RealFunc> ipoisson;           //!< Poisson's ratio of ice
  FuncConstPair<RealFunc> strparam;           //!< Strength model parameter

  FuncConstPair<RealFunc> satparam1;          //!< Parameter 1 of saturation curve
  FuncConstPair<RealFunc> satparam2;          //!< Parameter 2 of saturation curve
  FuncConstPair<RealFunc> satparam3;          //!< Parameter 3 of saturation curve

  FuncConstPair<RealFunc> frtemp;             //!< Freezing temperature
  FuncConstPair<RealFunc> meltingparam;       //!< Pressure melting parameter
  FuncConstPair<RealFunc> refpres;            //!< Reference pressure
};

#endif

// $Id$
//==============================================================================
//!
//! \file SIMTHMCoupledEP.h
//!
//! \date
//!
//! \author Yared Worku Bekele
//!
//! \brief Simulation driver for THM coupled problems in ground freezing
//!
//==============================================================================

#ifndef _SIM_THMCOUPLED_EP_H
#define _SIM_THMCOUPLED_EP_H

#include "THMCoupledEP.h"
#include "SIMSolver.h"
#include "TimeStep.h"
#include "Poro3Material.h"
#include "NonLinSIM.h"
#include "Utilities.h"
#include "InitialConditionHandler.h"
#include "DataExporter.h"
#include "SAM.h"

/*!
  \brief Class for THM coupled ground freezing simulators
*/

template<class Dim> class SIMTHMCoupledEP : public Dim
{
public:
  //! \brief The constructor initializes the references to the integrand.
  SIMTHMCoupledEP(bool supg) : Dim({Dim::dimension,2}), thm(Dim::dimension,1,supg),
                               thmwd(Dim::dimension), nSim(*this)
  {
    Dim::myProblem = &thm;
  }

  //! \brief Destructor
  virtual ~SIMTHMCoupledEP()
  {
    Dim::myProblem = NULL;
    Dim::myInts.clear();
  }

  //! \brief Parses a subelement of the \a resultoutput XML tag
  virtual bool parseOutputTag(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),"resultpoints"))
      if (utl::getAttribute(elem,"file",pointfile) && Dim::nProc > 1)
      {
        char cPid[8];
        sprintf(cPid,"_p%04d",Dim::myPid);
        pointfile.append(cPid);
      }

    return this->Dim::parseOutputTag(elem);
  }

  //! \brief Parses a data section from an XML element
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"thmcoupledep"))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
    {
      if (!strcasecmp(child->Value(),"materialdata"))
      {
        int code = this->parseMaterialSet(child,mVec.size());
        std::cout << "\tMaterial code " << code << ":" << std::endl;
        mVec.push_back(std::unique_ptr<Poro3Material>(new Poro3Material));
        mVec.back()->parse(child);
        mVec.back()->printLog();
      }
      else if (!strcasecmp(child->Value(),"envtproperties"))
      {
        double Te = 0.0, lambdae = 0.0;
        utl::getAttribute(child,"Te",Te);
        utl::getAttribute(child,"lambdae",lambdae);
        thmwd.setEnvtTemperature(Te);
        thmwd.setEnvtConductivity(lambdae);
        std::cout << "\tEnvironment thermal properties:"
                  << "\n\t\tTe = " << Te << "\tlambdae = " << lambdae << std::endl;
      }
      else if (!strcasecmp(child->Value(),"scaling"))
      {
        double scl1 = 1.0, scl2 = 1.0;
        utl::getAttribute(child,"scl1",scl1);
        utl::getAttribute(child,"scl2",scl2);
        thm.setScalingValues(scl1,scl2);
        std::cout << "\tScaling values:"
                  << "\n\t\tscl1 (U-P) = " << scl1 << "\tscl2 (U-T) = " << scl2 << std::endl;
      }
      else if (!strcasecmp(child->Value(),"gravity"))
      {
        const char* value = utl::getValue(child,"gravity");
        Vec3 grav;
        std::stringstream str;
        str << value;
        str >> grav[0] >> grav[1] >> grav[2];
        thm.setGravity(grav);
        std::cout << "\tGravitation vector: "
                  << "\n\t\tg = " << grav << std::endl;
      }
      else if (!strcasecmp(child->Value(),"nonlinearsolver"))
        nSim.parse(child);
      else
        this->Dim::parse(child);
    }

    if (!mVec.empty())
      thm.setMaterial(mVec.front().get());

    return true;
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the generic Neumann property codes.
  virtual void preprocessA()
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    PropertyVec::iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); p++)
      if (p->pcode == Property::NEUMANN_GENERIC)
      {
        if (Dim::myInts.find(p->pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p->pindx,&thmwd));
      }
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
    typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);
    typename Dim::TracFuncMap::const_iterator tit = Dim::myTracs.find(propInd);

    if (sit != Dim::myScalars.end())
      thmwd.setFlux(sit->second);
    else if (vit != Dim::myVectors.end())
      thm.setTraction(vit->second);
    else if (tit != Dim::myTracs.end())
      thm.setTraction(tit->second);
    else
      return false;

    return true;
  }

  //! \brief Evaluates some iteration norms for convergence assessment.
  //! \param[in] u Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[out] eNorm Energy norm of solution increment
  //! \param[out] rNorm Residual norm of solution increment
  //! \param[out] dNorm Displacement norm of solution increment
  void iterationNorms(const Vector& u, const Vector& r,
          double& eNorm, double& rNorm, double& dNorm) const
  {
    eNorm = this->mySam->dot(r,u,'A');
    rNorm = this->mySam->norm2(r,'D');
    dNorm = this->mySam->norm2(u,'D');

    // Norms for field variables on basis 2
    // rNorm = this->mySam->norm2(r,'P');
    // dNorm = this->mySam->norm2(u,'P');
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const
  {
    return "THMCoupledEP";
  }

  //! \brief Obtain const reference to solution vector.
  //! \param[in] i Solution vector to get reference to
  //! \return Const reference to requested vector
  const Vector& getSolution(int i)
  {
    return nSim.getSolution(i);
  }

  //! \brief Initializes the simulator time-stepping loop
  bool init(const TimeStep& tp)
  {
    // Initialize solution vectors
    nSim.initSol();
    size_t n, nSols = this->getNoSolutions();
    std::string str = "vector1";
    for (n = 0; n < nSols; n++, str[6]++)
      this->registerField(str,nSim.getSolution(n));

    return true;
  }

  //! \brief Opens a new VTF file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF file name from
  //! \param[out] geoBlk Running geometry blcok counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    geoBlk = nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results to to a VTF file for a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (this->getNoResultPoints() > 0)
    {
      double old = utl::zero_print_tol;
      utl::zero_print_tol = 1e-16;
      this->savePoints(pointfile,nSim.getSolution(),tp.time.t,tp.step+1,3);
      utl::zero_print_tol = old;
    }

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;

    // Write solution fields
    bool result = this->writeGlvS(nSim.getSolution(), iDump, nBlock,
                                  tp.time.t, false, "vector", 89);

    return result && this->writeGlvStep(iDump, tp.time.t);
  }

  //! \brief Advances the time step one step forward
  bool advanceStep(TimeStep& tp)
  {
    return nSim.advanceStep(tp,false);
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (nSim.solveStep(tp,SIM::DYNAMIC) != SIM::CONVERGED)
      return false;

    return this->postSolve(tp);
  }

  //! \brief Post-processing of solution, if needed
  bool postSolve(const TimeStep&, bool = false)
  {
    return true;
  }

  //! \brief Register fields for data export
  void registerFields(DataExporter& exporter)
  {
    //int results = DataExporter::PRIMARY | DataExporter::SECONDARY | DataExporter::RESTART;
    int results = DataExporter::PRIMARY;
    exporter.registerField("u-p-T","primary",DataExporter::SIM,results);
    exporter.setFieldValue("u-p-T",this,&nSim.getSolution());
  }

  //! \brief Sets initial conditions.
  void setInitialConditions()
  {
    SIM::setInitialConditions(*this);
  }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    thm.setMaterial(mVec[propInd].get());
    return true;
  }

private:
  THMCoupledEP thm;                                   //!< THMCoupledEP integrand
  THMCoupledEP::WeakDirichlet thmwd;                  //!< THMCoupledEP Robin integrand
  std::string pointfile;                              //!< File name for point output data
  std::vector<std::unique_ptr<Poro3Material>> mVec;   //!< Material data
  NonLinSIM nSim;                                     //!< Nonlinear simulator
};

#endif  // _SIM_THMCOUPLED_NL_H

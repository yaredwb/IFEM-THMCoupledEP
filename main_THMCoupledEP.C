// $Id$
//==============================================================================
//!
//! \file main_THMCoupledEP.C
//!
//! \date
//!
//! \author Yared Worku Bekele
//!
//! \brief Main program for the simulator of THM coupled problems
//!
//==============================================================================

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMTHMCoupledEP.h"
#include "SIMSolver.h"
#include "ASMmxBase.h"
#include "Utilities.h"
#include "Profiler.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "VTF.h"
#include "AppCommon.h"
#include "TimeIntUtils.h"

template<class Dim>
int runSimulator(char* infile, bool supg)
{
  SIMTHMCoupledEP<Dim> model(supg);
  SIMSolver< SIMTHMCoupledEP<Dim> > solver(model);

  //! Read input file
  if (!model.read(infile) || !solver.read(infile))
    return 1;

  //! Configure finite element library
  if (!model.preprocess())
    return 2;

  //! Setup integration
  model.setQuadratureRule(model.opt.nGauss[0],true);
  model.initSystem(model.opt.solver,1,1,false);
  model.init(solver.getTimePrm());
  model.setInitialConditions();
  model.setAssociatedRHS(0,0);
  model.setMode(SIM::DYNAMIC);

  // HDF5 output
  DataExporter* exporter = nullptr;
  if (model.opt.dumpHDF5(infile))
  {
    //exporter = new DataExporter(true);
    /*exporter->registerWriter(new HDF5Writer(model->opt.hdf5,model->getProcessAdm()));
    exporter->registerWriter(new XMLWriter(model->opt.hdf5,model->getProcessAdm()));*/
    exporter = SIM::handleDataOutput(model,solver,model.opt.hdf5,
                                     false,model.getDumpInterval(),1);
  }

  if (!solver.solveProblem(infile,exporter))
    return 5;

  delete exporter;

  return 0;
}


int main(int argc, char ** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SIMoptions dummy;
  std::vector<int> ignoredPatches;

  int ndim = 3;
  char* infile = 0;
  bool supg = false;
  ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;

  IFEM::Init(argc,argv);

  for (int i = 1; i < argc; i++)
  {
    if (dummy.parseOldOptions(argc,argv,i))
      ; // Ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strcmp(argv[i],"-1D"))
      ndim = 1;
    else if (!strcasecmp(argv[i],"-SUPG"))
      supg = true;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr << "*** Unknown option ignored: " << argv[i] << std::endl;
  }

  if (!infile)
  {
    IFEM::cout << "Usage: " << argv[0]
               << " <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
               << " [-free] [-lag|-spec|-LR] [-1D|-2D] [-mixed] [-nGauss <n>]"
               << "\n       [-vtf <format> [-nviz <nviz>]"
               << " [-nu <nu> [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
               << "       [-eig <iop> [-nev <nev>] [-ncv <ncv>] [-shift <shf>]]\n"
               << "       [-ignore <p1> <p2> ...] [-fixDup]" << std::endl;
    return 0;
  }

  IFEM::cout << "\n >>> IFEM Nonlinear THM Coupled Solver <<<"
             << "\n ========================================="
             << "\n Executing command:\n";
  for (int i = 0; i < argc; i++) IFEM::cout << " " << argv[i];
  IFEM::cout << "\n\n Input file: " << infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (ndim == 3)
    return runSimulator<SIM3D>(infile,supg);
  else if (ndim == 2)
    return runSimulator<SIM2D>(infile,supg);
  else
    return runSimulator<SIM1D>(infile,supg);

  return 1;
}
// $Id: minimizer.h 1286 2013-10-28 11:54:20Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iomanip>

#include "common.h"
#include "fitpdfset.h"
using std::vector;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <NNPDF/experiments.h>
#include <NNPDF/positivity.h>
using NNPDF::Experiment;
using NNPDF::PositivitySet;

class NNPDFSettings;
/**
 *  \class Minimizer
 *  \brief Virtual minimisation base class
 */
class Minimizer
{
public:
  Minimizer(NNPDFSettings const&);
  virtual ~Minimizer();

  enum Mode {
    SetMode,
    ExpMode
  };  //!< Minimisation mode

  enum SortPDF {
    PDF_SORT,
    PDF_NOSORT
  };  //!< PDF sorting

  virtual void Init(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&);
  virtual void Iterate(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&) = 0;  //!< Perform an iteration of the minimisation

protected:
  virtual void ComputeErf(FitPDFSet*,
                  vector<Experiment*> const&,
                  vector<PositivitySet> const&,
                  Minimizer::Mode,
                  Minimizer::SortPDF); //!< Compute fChi2Mem

  real*  fChi2Mem;

  const NNPDFSettings& fSettings;
};

// *************************************************************************************

class CMAESParam
{
public:
  CMAESParam(size_t const& _n, size_t const& _lambda);
  const size_t lambda;
  const size_t mu;
  const size_t n;
  size_t eigenInterval;
  double expN;
  double mu_eff;
  double csigma;
  double dsigma;
  double cc;
  double c1;
  double cmu;
  std::vector<double> wgts;
};

/**
 *  \class CMAESMinimizer
 *  \brief CMA-ES minimiser
 */

class CMAESMinimizer : public Minimizer
{
public:
  CMAESMinimizer(NNPDFSettings const&);
  ~CMAESMinimizer();

  virtual void Init(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&);
  virtual void Iterate(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&);

private:
  std::vector<gsl_vector*> Mutation(FitPDFSet* pdf) const;
  gsl_vector* Recombination(FitPDFSet* pdf, vector<size_t> const& rank, std::vector<gsl_vector*> const& yvals) const;

  void CSA(gsl_vector const* yavg);
  void CMA(FitPDFSet*, vector<size_t> const& rank, std::vector<gsl_vector*> const& yvals, gsl_vector const* yavg);

  void GetParam(Parametrisation* const, gsl_vector*) const;
  void SetParam(gsl_vector* const, Parametrisation*) const;

  void NormVect(gsl_vector*) const; //!< Normally distributed random vector
  void ComputeEigensystem();

protected:
  size_t fNTparam;
  double fSigma;
  CMAESParam* fCMAES;
  gsl_vector *fpsigma, *fpc;
  gsl_matrix *fC, *fBD, *finvC;
  gsl_eigen_symmv_workspace *fwrkspc;
};
